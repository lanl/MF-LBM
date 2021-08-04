#include "./preprocessor.h"
!========================================================================================================================================
!---------------------- geometry related ----------------------
!========================================================================================================================================
!******************************* set walls ******************************
subroutine set_walls
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,L,M,N,k,icount,isize,out1,out2,out3,rank,status(MPI_STATUS_SIZE),nmax
    real(kind=8) :: x,y,z,r1,r2,xc,yc,zc
    LOGICAL ALIVE
    integer :: num
    character (len=200) :: flnm, dummy   !file name

    !~~~~~~~~~~~~~ temporary fix for IBM Power9 nodes. Otherwise, MPI_Bcast will be extremely slow for a large array. ~~~~~~~~~~~~~~~~~~~~
    ! cause unknown. likely have something to do with memory allocation of the MPI buffers.
    integer(kind=1), allocatable,dimension(:,:,:) :: tt
    allocate(tt(1:nx,1:ny,1:nz))
    deallocate(tt)
    !~~~~~~~~~~~~~ temporary fix for IBM Power9 nodes. Otherwise, MPI_Bcase will be extremely slow for a large array. ~~~~~~~~~~~~~~~~~~~~

    isize = nx*ny*nz

    !~~~~~~~~~~~~~~~~ initialize wall array ~~~~~~~~~~~~~~~~~~~~
    ! local
    !$OMP parallel
    !$omp do private(i,j) collapse(2)
    do k=1-overlap_walls,nz+overlap_walls
        do j=1-overlap_walls,ny+overlap_walls
            do i=1-overlap_walls,nx+overlap_walls
                walls(i,j,k)=0
            enddo
        enddo
    enddo

    ! global
    !$OMP DO private(i,j) collapse(2)
    do k=1,nzglobal
        do j=1,nyglobal
            do i=1,nxglobal
                walls_global(i,j,k)=0
            enddo
        enddo
    enddo
    !$omp end parallel

    !~~~~~~~~~~~~~~~~ read wall data ~~~~~~~~~~~~~~~~~~~~
    nx_sample=0
    ny_sample=0
    nz_sample=0
    if(external_geometry_read_cmd==1)then
        if(id==0)print*, 'This simulation uses external geometry data!'
        INQUIRE(FILE='./path_info.txt',EXIST=ALIVE)
        if(alive)then
            OPEN(UNIT=11,FILE='./path_info.txt',action='READ')
            read(11,'(A)')dummy
            read(11,'(A)')dummy
            read(11,'(A)')dummy  ! skip three lines
            read(11,'(A)')geo_file_path   ! read geometry file path
            close(11)
        else
            if(id==0)print*,'Error! path_info.txt is not found! Exiting program!'
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            call mpi_abort(MPI_COMM_WORLD,ierr) 
        endif

        INQUIRE(FILE=trim(geo_file_path),EXIST=ALIVE)
        if(alive)then
            call read_walls
        else     
            if(id==0)print*,'Error! No external geometry file found! Exiting program!'
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            call mpi_abort(MPI_COMM_WORLD,ierr)
        endif
    else
        if(id==0)print*, 'This simulation does not use external geometry data!'        ! duct flow
    endif

    !~~~~~~~~~~~~~~~~ modify geometry or create hard coded geometry ~~~~~~~~~~~~~~~~~~~~
    if(modify_geometry_cmd==1.and.id==0)then
        call modify_geometry
    endif
    
    !~~~~~~~~~~~~~~~~ specify channel walls for the global wall array ~~~~~~~~~~~~~~~~~~~~
    ! walls=1: solid;  walls=0: fluid
    if(id==0)then
        !$OMP PARALLEL DO private(i,j)
        do k=1,nzglobal
            do j=1,nyglobal
                do i=1,nxglobal
                    if(domain_wall_status_z_min==1)walls_global(i,j,1)=1

                    if(domain_wall_status_z_max==1)walls_global(i,j,nzGlobal)=1

                    if(domain_wall_status_x_min==1)walls_global(1,j,k)=1

                    if(domain_wall_status_x_max==1)walls_global(nxGlobal,j,k)=1

                    if(domain_wall_status_y_min==1)walls_global(i,1,k)=1

                    if(domain_wall_status_y_max==1)walls_global(i,nyGlobal,k)=1
                enddo
            enddo
        enddo
    endif

    !~~~~~~~~~~~~~~~~ broadcast and distribute geometry data ~~~~~~~~~~~~~~~~~~~~
    if(extreme_large_sim_cmd==0)then
        call MPI_Bcast(walls_global,nxGlobal*nyGlobal*nzGlobal,MPI_INTEGER1,0,MPI_COMM_VGRID,ierr)
    else
        ! walls_global array may be too large that (nxGlobal)*(nyGlobal)*(nzGlobal) exceeds 2^31 limit of MPI_Bcast
        nmax = 4 ! divide the array into nmax pieces (nzglobal must be completely divided by nmax)
        if(mod(nzglobal,nmax)==0)then
            do n=1,nmax
                k = 1+nzglobal/nmax*(n-1)
                out3 = nzGlobal/nmax
                call MPI_Bcast(walls_global(1,1,k),(nxGlobal)*(nyGlobal)*out3,MPI_INTEGER1,0,MPI_COMM_VGRID,ierr)
            enddo
        else
            if(id==0)print*,'Error in broadcasting walls_global array! mod(nzglobal,nmax)/=0'
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            call mpi_abort(MPI_COMM_WORLD,ierr)
        endif
    endif

    !$OMP PARALLEL DO private(i,j,l,m,n)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                L = idx*nx + i
                M = idy*ny + j
                N = idz*nz + k
                walls(i,j,k)=walls_global(l,m,n)
            enddo
        enddo
    enddo

    !~~~~~~~~~~~~~~~~ specify channel walls for local wall arrays ~~~~~~~~~~~~~~~~~~~~
    ! the reason to do this in local arrays is to specify correct solid nodes information for the 
    ! overlap layers of the local wall array. For example, walls(j>=ny)=1 means ny, ny+1,..., 
    ! ny+overlap_walls are all solid nodes.
    ! walls=1: solid;  walls=0: fluid
    !$OMP PARALLEL DO private(i,j)  collapse(2)
    do k=1-overlap_walls,nz+overlap_walls
        do j=1-overlap_walls,ny+overlap_walls
            do i=1-overlap_walls,nx+overlap_walls
                if(idz==0)then           
                    if(k<=1)then
                        !  if(domain_wall_status_z_min==0)walls(i,j,k)=0
                        if(domain_wall_status_z_min==1)walls(i,j,k)=1     !solid wall BC
                    endif
                endif
                if(idz==npz-1)then
                    if(k>=nz)then
                        !  if(domain_wall_status_z_max==0)walls(i,j,k)=0
                        if(domain_wall_status_z_max==1)walls(i,j,k)=1     !solid wall BC
                    endif
                endif
 
                if(idx==0)then
                    if(i<=1)then
                        !if(domain_wall_status_x_min==0)     !do nothing, usually accompanied by periodic BC
                        if(domain_wall_status_x_min==1)walls(i,j,k)=1     !solid wall BC
                    endif
                endif
                if(idx==npx-1)then
                    if(i>=nx)then
                        !if(domain_wall_status_x_max==0)     !do nothing, usually accompanied by periodic BC
                        if(domain_wall_status_x_max==1)walls(i,j,k)=1     !solid wall BC
                    endif
                endif
 
                if(idy==0)then                
                    if(j<=1)then
                        !if(domain_wall_status_y_min==0)     !do nothing, usually accompanied by periodic BC
                        if(domain_wall_status_y_min==1)walls(i,j,k)=1     !solid wall BC
                    endif                     
                endif
                if(idy==npy-1)then
                    if(j>=ny)then
                        !if(domain_wall_status_Y_max==0)     !do nothing, usually accompanied by periodic BC
                        if(domain_wall_status_Y_max==1)walls(i,j,k)=1     !solid wall BC                 
                    endif
                endif
            enddo
        enddo
    enddo

    !~~~~~~~~~~~~~~~~ calculate open area at inlet ~~~~~~~~~~~~~~~~~~~~
    ! A_xy_effective is not equal to A_xy
    ! A_xy_effective is used when the inlet does not cover the entire xy plane
    icount = 0
    !$OMP PARALLEL DO private(i)reduction(+:icount)
    do j=1,nyglobal
        do i=1,nxglobal
            if(walls_global(i,j,1)<=0)icount = icount + 1
        enddo
    enddo
    A_xy_effective= icount  !effecitve inlet area
    if(id==0)print*,'Inlet effective open area = ', A_xy_effective
    call pore_profile   !pore profile info along the flow direction z

    !~~~~~~~~~~~~~~~~ save fluid point information ~~~~~~~~~~~~~~~~~~~~
    ! only used when extreme_large_sim_cmd==1, to be combined with save_phi (distributed phi data)
    if(extreme_large_sim_cmd==1)then
        icount = 0
        !$OMP PARALLEL DO private(i,j) reduction(+:icount)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(walls(i,j,k)<=0)then
                        icount = icount + 1
                    endif
                enddo
            enddo
        enddo
        n_fluid_node_local = icount  ! total fluid nodes of local domain

        write(flnm,"('geometry_id',i5.5)")id
        if(id==0)print*,'Start to save fluid point information'
        open(unit=9+id, file='out3.field_data/'//trim(flnm), FORM='unformatted', status='replace',access='stream')
        write(9+id)n_fluid_node_local
        write(9+id)idx,idy,idz,nx,ny,nz
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(walls(i,j,k)==0)then
                        write(9+id)i,j,k
                    endif
                enddo
            enddo
        enddo
        close(9+id)
        if(id==0)print*,'Fluid point information saved!'
    endif

    call ztransport_walls(0,0,overlap_walls)
    call ytransport_walls(0,overlap_walls,overlap_walls)
    call xtransport_walls(overlap_walls,overlap_walls,overlap_walls)

    return
end subroutine set_walls

!************************** modify geometry *******************************
subroutine modify_geometry
    use Misc_module
    IMPLICIT NONE
    integer :: i,j,k, buffer
    double precision :: xc,yc,zc,r1,r2
    ! walls=1: solid;  walls=0: fluid
    ! below is a sample code to modify geometry

    ! domain center, used as a reference point
    xc = 0.5d0*dble(nxglobal+1)
    yc = 0.5d0*dble(nyglobal+1)
    zc = 0.5d0*dble(nzglobal+1)
    r1 = 0.15d0*nyglobal  ! used to creat simple obstacle
    r2 = nyglobal * 0.5d0  ! tube radius
    buffer = 10    ! inlet outlet reservior
    !$OMP PARALLEL DO private(i,j)
    do k=1,nzglobal
        do j=1,nyglobal
            do i=1,nxglobal
                if(( i - xc )**2+( j - yc )**2+( k - zc )**2 < r1**2)then
                    walls_global(i,j,k)=1
                endif
                ! if(( i - xc )**2+( j - yc )**2 > r2**2 .and. k>buffer .and. k< nzglobal-buffer+1)then
                !     walls_global(i,j,k)=1
                ! endif
            enddo
        enddo
    enddo
    print*, 'Internal geometry modified!'

    return
end subroutine modify_geometry

!***************************** read walls ************************************
subroutine read_walls
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k,L,M,N
    character (len=300) :: dummy
    integer :: error_signal  

    error_signal = 0
    if(id==0)then
        OPEN(UNIT=11,FILE='./path_info.txt',action='READ')
        read(11,'(A)')dummy
        read(11,'(A)')dummy
        read(11,'(A)')dummy  ! skip three lines
        read(11,'(A)')geo_file_path   ! read geometry file path
        close(11)
        OPEN(UNIT=9,FILE=trim(adjustl(geo_file_path)),STATUS='OLD',FORM='UNFORMATTED',access='stream')
        read(9)nx_sample,ny_sample,nz_sample
        write(*,"(' Porous media sample size: nx=', I5, ', ny=', I5,', nz=', I5)")nx_sample,ny_sample,nz_sample
        if(nxglobal<nx_sample.or.nyglobal<ny_sample.or.nzglobal<nz_sample)then
            print*,'Error! Domain size is smaller than porous media sample size! Exiting program!'
            error_signal = 1
        else
            read(9)(((walls_global(i,j,k),i=1,nx_sample),j=1,ny_sample),k=1,nz_sample)
        endif
        close(9)

        if(domain_wall_status_x_max==1.and.domain_wall_status_y_max==1)then
            !$OMP PARALLEL DO private(i,j)
            do k=1,nzglobal
                do j=1,nyglobal
                    do i=1,nxglobal
                        if(j>=ny_sample.or.i>=nx_sample)walls_global(i,j,k)=1   !pad with solid walls
                    enddo
                enddo
            enddo
        endif
    endif         

    call MPI_Bcast(error_signal,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)
    if(error_signal == 1)call mpi_abort(MPI_COMM_WORLD,ierr)

    call MPI_Bcast(nx_sample,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)
    call MPI_Bcast(ny_sample,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)
    call MPI_Bcast(nz_sample,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)

    return
end subroutine read_walls

!*******************************pore profile*************************************
subroutine pore_profile
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer,allocatable,dimension(:) ::pore,tt1
    integer :: i,j,k ,L,M,N,out1,out2,out3,rank,isize,status(MPI_STATUS_SIZE)

    isize =nz
    allocate(tt1(isize),pore(nz))

    !$omp parallel do private(i,j)
    do k=1,nz
        pore(k)=0
        do j=1,ny
            do i=1,nx
                if(walls(i,j,k)<=0)then
                    pore(k)=pore(k)+1
                endif
            enddo
        enddo
    enddo

    if(id.ne.0)then
        !packing the data for transfer
        !$omp parallel do
        do k=1,nz
            tt1(k)=pore(k)
        enddo
        call mpi_send(tt1,isize,MPI_INTEGER,0,600+id,MPI_COMM_vgrid,ierr)
    else
        allocate(pore_profile_z(nzglobal))
        pore_profile_z(:)=0
        !$omp parallel do
        do k=1,nz
            pore_profile_z(k)= pore(k)
        enddo

        do rank = 1, np-1
            call mpi_recv(tt1,isize,MPI_INTEGER,rank,600+rank,MPI_COMM_vgrid,status,ierr)
            CALL MPI_CART_COORDS(MPI_COMM_VGRID, rank, mpi_dim, mpi_coords, ierr)
            out1=mpi_coords(1)
            out2=mpi_coords(2)
            out3=mpi_coords(3)

            !$omp parallel do
            do k=1,nz
                pore_profile_z(out3*nz+k)= pore_profile_z(out3*nz+k) + tt1(k)
            enddo
        enddo
        pore_sum = 0
        !$omp parallel do reduction (+:pore_sum)
        do k=1,nzGlobal
            pore_sum = pore_sum + pore_profile_z(k)
        enddo

        !$omp parallel do reduction (+:pore_sum_effective)
        do k=1+n_exclude_inlet,nzGlobal-n_exclude_outlet
            pore_sum_effective = pore_sum_effective + pore_profile_z(k)
        enddo

    endif
    call MPI_Bcast(pore_sum,1,MPI_INTEGER,0,MPI_COMM_vgrid,ierr)
    call MPI_Bcast(pore_sum_effective,1,MPI_INTEGER,0,MPI_COMM_vgrid,ierr)

    deallocate(tt1,pore)
    return
end subroutine pore_profile



!======================================================================================================================================
!---------------------- compute macroscopic varaibles from PDFs ----------------------
!======================================================================================================================================
subroutine compute_macro_vars     ! u,v,w,rho  
    use Misc_module
    use Fluid_singlephase
    IMPLICIT NONE
    integer :: i,j,k
    integer*1 ::   wall_indicator
    real(kind=8) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,fx,fy,fz,tmp

    !$OMP parallel DO private(i,j,ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,wall_indicator,fx,fy,fz)
    !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
    !$acc & u,v,w,rho,walls)
    !$acc loop collapse(3) device_type(NVIDIA)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                ft0 = f0(i,j,k)
                ft1  = f1(i,j,k)
                ft2  = f2(i,j,k)
                ft3  = f3(i,j,k)
                ft4  = f4(i,j,k)
                ft5  = f5(i,j,k)
                ft6  = f6(i,j,k)
                ft7  = f7(i,j,k)
                ft8  = f8(i,j,k)
                ft9  = f9(i,j,k)
                ft10  = f10(i,j,k)
                ft11  = f11(i,j,k)
                ft12  = f12(i,j,k)
                ft13  = f13(i,j,k)
                ft14  = f14(i,j,k)
                ft15  = f15(i,j,k)
                ft16  = f16(i,j,k)
                ft17  = f17(i,j,k)
                ft18  = f18(i,j,k)

                rho(i,j,k)=(ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8+ft9+ft10+ft11+ft12+ft13+ft14+ft15+ft16+ft17+ft18)* (1-wall_indicator)

                fx = 0d0
                fy = 0d0
                fz = force_Z

                !here "- 0.5fx" is due to that PDFs are after even step, which is post collision before streaming
                !to use the post collision PDFs to calculate the velocities, one must substruct the forcing terms applied during collision step
                !thus the fomular is u = f...f + 0.5fx - fx = f...f - 0.5fx
                u(i,j,k)=  (ft1 - ft2 + ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14 - 0.5d0*fx)* (1-wall_indicator)
                v(i,j,k)=  (ft3 - ft4 + ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18 - 0.5d0*fy)* (1-wall_indicator)
                w(i,j,k)=  (ft5 - ft6 + ft11+ ft12- ft13- ft14 + ft15 + ft16 - ft17 - ft18 - 0.5d0*fz)* (1-wall_indicator)
            enddo
        enddo
    enddo
    !$acc end kernels

    return
end subroutine compute_macro_vars


!======================================================================================================================================
!---------------------- OpenACC ----------------------
!======================================================================================================================================
#ifdef _openacc
function setDevice(nprocs,myrank)
    use iso_c_binding
    use openacc
    implicit none
    include 'mpif.h'
    interface
        function gethostid() BIND(C)
            use iso_c_binding
            integer (C_INT) :: gethostid
        end function gethostid
    end interface

    integer :: nprocs, myrank
    integer, dimension(nprocs) :: hostids, localprocs
    integer :: hostid, ierr, numdev, mydev, i, numlocal
    integer :: setDevice

    ! get the hostids so we can determine what other processes are on this node
#ifdef _WIN32
    hostid = 1
#else
    hostid = gethostid()
#endif

    CALL mpi_allgather( hostid, 1, MPI_INTEGER, hostids, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

    ! determine who's on this node
    numlocal=0
    localprocs=0
    do i=1,nprocs
        if (hostid .eq. hostids(i)) then
            localprocs(i)=numlocal
            numlocal = numlocal+1
        endif
    enddo

    ! get the number of devices on this node
    numdev = acc_get_num_devices(ACC_DEVICE_NVIDIA)

    if (numdev .lt. 1) then
        print *, 'ERROR: There are no devices available on this host.  ABORTING.', myrank
        stop
    endif

    ! print a warning if the number of devices is less then the number
    ! of processes on this node.  NVIDIA does not support multiple context on a single device.
    ! Hence, having multiple processes share devices is not recommended.  It may work, but
    ! can also fail.
    if (numdev .lt. numlocal) then
        if (localprocs(myrank+1).eq.1) then
            ! print the message only once per node
            print *, 'WARNING: The number of local process is greater then the number of devices.', myrank
        endif
        mydev = mod(localprocs(myrank+1),numdev)
    else
        mydev = localprocs(myrank+1)
    endif

    call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
    call acc_init(ACC_DEVICE_NVIDIA)
    setDevice = mydev
end function setDevice
#endif



!====================================================================================================================================
!---------------------- Unsorted ----------------------
!====================================================================================================================================
!*************inlet velocity - analytical solution**********************************************************
subroutine inlet_vel_profile_rectangular(vel_avg, num_terms)    
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    real(kind=8) :: vel_avg, a, b, xx, yy, tmp1, tmp2, tmp3
    integer :: n, num_terms, i,j, x, y
    
    a = 0.5d0 * la_x
    b = 0.5d0 * la_y

    tmp1 = 0.0d0
    do n=1,num_terms,2
        tmp1 = tmp1 + (dtanh(0.5d0*dble(n)*pi*b/a)) / n**5
    enddo

    tmp2 = 1.0d0 - 192d0 / pi**5 * (a/b) * tmp1

    tmp2 = -3d0 * vel_avg  / (tmp2 * a**2)
  
    !$OMP PARALLEL DO private(i,n,x,y,xx,yy,tmp3)
    do j=1,ny
      do i=1,nx   
          x = idx*nx + i
          y = idy*ny + j
          if(x>1.and.x<nxGlobal.and.y>1.and.y<nyGlobal)then
            xx = x - 1.5d0 - a 
            yy = y - 1.5d0 - b
            tmp3 = 0d0
            do n=1,num_terms,2
              tmp3 = tmp3 + (-1.0d0)**(0.5d0*dble(n-1)) * dcos(0.5d0*n*pi*xx/a)/n**3 &
              ! * (1.0d0 - dcosh(0.5d0*n*pi*yy/a) / dcosh(0.5d0*n*pi*b/a)) =   
              * (1.0d0 - ( dexp(0.5d0*n*pi*(yy-b)/a) + dexp(0.5d0*n*pi*(-yy-b)/a)) / (1d0 + dexp(0.5d0*n*pi*(-b-b)/a)) )  
            enddo
            w_in(i,j) = tmp3 * ( -16d0 * tmp2 * a**2 * pi**(-3) ) 
          endif
      enddo
    enddo

    return
end subroutine inlet_vel_profile_rectangular


