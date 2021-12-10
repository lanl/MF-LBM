
!==============================================================================
! This code converts distributed macro data (save_macro function in MF-LBM)
! to full domain field data for easy analysis. save_macro is activated when
! extreme_large_sim_cmd = 1 in control file (99% of chance that you don't need this)
!==============================================================================

MODULE misc
    IMPLICIT NONE
    SAVE
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: u,v,w,rho
    integer :: np
    integer :: ntime0,ntime_interval,ntime_max
    integer :: nxglobal,nyglobal,nzglobal
END MODULE misc


!################################################################################################################################
!Main Start
!################################################################################################################################
program main_singlephase
    use misc
    IMPLICIT NONE
    integer :: i,j,k
    double precision :: temp
    character (len=50) :: flnm   !file name

    call system('mkdir ./output')

    open(5,file='./input_parameters.txt',status='old')
    read(5,*)np  ! total number of MPI processes in the original simulation
    read(5,*)nxglobal,nyglobal,nzglobal  !dimensions of the entire simulation domain
    read(5,*)ntime0,ntime_interval,ntime_max  !start time step, time interval, end time step
    close(5)

    allocate(u(nxglobal,nyglobal,nzglobal),v(nxglobal,nyglobal,nzglobal),w(nxglobal,nyglobal,nzglobal),rho(nxglobal,nyglobal,nzglobal))
    !$OMP parallel DO private(i,j)
    do k=1,nzglobal
        do j=1,nyglobal
            do i=1,nxglobal
                u(i,j,k)=0d0
                v(i,j,k)=0d0
                w(i,j,k)=0d0
                rho(i,j,k)=1d0
            enddo
        enddo
    enddo

    call read_save_macro  

    deallocate(u,v,w,rho)

End program main_singlephase
!################################################################################################################################
!Main End
!################################################################################################################################

!************************************* read and save macro field files *********************************
subroutine read_save_macro  
    use Misc
    implicit none
    integer :: i,j,k,l,m,n, nx,ny,nz
    integer :: n_step,ncount,num, nt,n_fluid_node_local, id,idx,idy,idz
    character (len=30) :: flnm   !file name
    double precision :: v1,v2,v3,v4,v5
    type geometry
    integer :: n_fluid_node_local
        integer, allocatable, dimension(:) :: ix,iy,iz
    end type geometry

    type(geometry), allocatable :: fluid_array(:)

    allocate(fluid_array(0:np-1))
    !$omp parallel do private(flnm,num,idx,idy,idz,nx,ny,nz,i,j,k)
    do id=0,np-1
        write(flnm,"('./data/geometry_id',i5.5)")id
        open(unit=9+id, file=trim(flnm), FORM='unformatted', status='old',access='stream')
        read(9+id)fluid_array(id)%n_fluid_node_local  !total number of fluid nodes of local MPI subdomain
        allocate(fluid_array(id)%ix(fluid_array(id)%n_fluid_node_local),fluid_array(id)%iy(fluid_array(id)%n_fluid_node_local),fluid_array(id)%iz(fluid_array(id)%n_fluid_node_local))
        read(9+id)idx,idy,idz,nx,ny,nz
        do num=1,fluid_array(id)%n_fluid_node_local
            read(9+id)i,j,k
            l = idx*nx + i
            m = idy*ny + j
            n = idz*nz + k
            fluid_array(id)%ix(num) = l
            fluid_array(id)%iy(num) = m
            fluid_array(id)%iz(num) = n
        enddo
        close(9+id)
    enddo

    ncount = (ntime_max-ntime0)/ntime_interval + 1
    print*,'number of vtk files pending processing:', ncount

    do n_step=1,ncount
        nt = ntime0 + (n_step-1)*ntime_interval
        print*, 'Start processing macro field of ntime=', nt
        !$omp parallel do private(flnm,num,v2,v3,v4,v5,i,j,k)
        do id=0,np-1
            write(flnm,"('full_nt',i9.9,'_id',i5.5)")nt,id
            open(unit=9+id, file='./data/'//trim(flnm), FORM='unformatted', status='old',access='stream')
            do num=1,fluid_array(id)%n_fluid_node_local
                read(9+id)v2,v3,v4,v5
                i = fluid_array(id)%ix(num)
                j = fluid_array(id)%iy(num)
                k = fluid_array(id)%iz(num)
                u(i,j,k) = v2
                v(i,j,k) = v3
                w(i,j,k) = v4
                rho(i,j,k) = v5
            enddo
            close(9+id)
        enddo

        !call save_macro(nt)
        !call VTK_detail_bin(nt)
        call VTK_detail_bin_half(nt)
        print*, 'End processing macro field of ntime=', nt
    enddo

    deallocate(fluid_array)

    return
end subroutine read_save_macro  
!************************************* read and save macro field files *********************************

!************************************* write data *********************************
subroutine save_macro(nt)   
    use Misc
    IMPLICIT NONE
    character*30 :: flnm
    integer :: nt

    write(flnm,'(i10.10,".dat")')nt
    open(unit=9,file='./output/rho_ntime_'//trim(flnm),FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')
    write(9)rho
    close(9)
    open(unit=9,file='./output/vel_ntime_'//trim(flnm),FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')
    write(9)u
    write(9)v
    write(9)w
    close(9)

    return
end subroutine save_macro

!************************************* VTK legacy output **************************************
subroutine VTK_detail_bin(nt)   
    use Misc
    IMPLICIT NONE
    character*30 :: flnm
    integer :: i,j,k, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9, int

    write(flnm,'(i10.10,".vtk")')nt

    open(unit=ivtk,file='./output/detail_'//flnm,FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf                                          ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nxglobal
    write(str2(1:10),'(i10)')nyglobal
    write(str3(1:10),'(i10)')nzglobal
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    num=nxglobal*nyglobal*nzglobal

    write(str1(1:10),'(i10)')num
    buffer = 'POINT_DATA '//str1//lf                                               ; write(ivtk) trim(buffer)

    !scalar - density
    buffer = 'SCALARS density double'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((rho(i,j,k),i=1,nxGlobal),j=1,nyGlobal),k=1,nzGlobal)

    buffer = 'SCALARS velocity_X double'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((u(i,j,k),i=1,nxGlobal),j=1,nyGlobal),k=1,nzGlobal)

    buffer = 'SCALARS velocity_Y double'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((v(i,j,k),i=1,nxGlobal),j=1,nyGlobal),k=1,nzGlobal)

    buffer = 'SCALARS velocity_Z double'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((w(i,j,k),i=1,nxGlobal),j=1,nyGlobal),k=1,nzGlobal)

    close(ivtk)

    return
end


subroutine VTK_detail_bin_half(nt)   
  use Misc
  IMPLICIT NONE
  character*30 :: flnm
  integer :: i,j,k, nt
  integer(kind=8) :: num
  character :: buffer*80, lf*1, str1*10, str2*10, str3*10,str4*14
  integer   :: ivtk = 9, int

  write(flnm,'(i10.10,".vtk")')nt

  open(unit=ivtk,file='./output/detail_half_'//flnm,FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

  lf = char(10) ! line feed character

  buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
  buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
  buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
  buffer = 'DATASET STRUCTURED_POINTS '//lf                                          ; write(ivtk) trim(buffer)

  write(str1(1:10),'(i10)')nxglobal/2
  write(str2(1:10),'(i10)')nyglobal/2
  write(str3(1:10),'(i10)')nzglobal/2
  buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

  write(str1(1:10),'(i10)')1
  write(str2(1:10),'(i10)')1
  write(str3(1:10),'(i10)')1
  buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

  write(str1(1:10),'(i10)')1
  write(str2(1:10),'(i10)')1
  write(str3(1:10),'(i10)')1
  buffer = 'SPACING '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

  num=int(nxGlobal,kind=8)*int(nyGlobal,kind=8)*int(nzGlobal,kind=8)/8

  write(str4(1:14),'(i14)')num
  buffer = 'POINT_DATA '//str4//lf                                               ; write(ivtk) trim(buffer)

  !scalar - density
  buffer = 'SCALARS density float'//lf
  write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk)(((real(rho(i,j,k)),i=1,nxGlobal,2),j=1,nyGlobal,2),k=1,nzGlobal,2)


  buffer = 'SCALARS velocity_X float'//lf
  write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk)(((real(u(i,j,k)),i=1,nxGlobal,2),j=1,nyGlobal,2),k=1,nzGlobal,2)

  buffer = 'SCALARS velocity_Y float'//lf
  write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk)(((real(v(i,j,k)),i=1,nxGlobal,2),j=1,nyGlobal,2),k=1,nzGlobal,2)

  buffer = 'SCALARS velocity_Z float'//lf
  write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk)(((real(w(i,j,k)),i=1,nxGlobal,2),j=1,nyGlobal,2),k=1,nzGlobal,2)

  close(ivtk)

  return
end

!************************************* VTK legacy output **************************************

