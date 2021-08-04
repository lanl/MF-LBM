#include "./preprocessor.h"
!=====================================================================================================================================
!---------------------- initialization basic ----------------------
!=====================================================================================================================================
subroutine initialization_basic
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k, lx,ly,lz,n_vin
    real(kind=8) :: omega
    LOGICAL :: ALIVE

    call read_parameter
    
    call set_MPI

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !create folders
    if(id==0)print*,''
    if(id==0)print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    if(id==0)print*,'Creating directories if not exist'
    if(id==0)call system('mkdir out2.checkpoint')
    if(id==0.and.double_bak_checkpoint_pdf_cmd==1)call system('mkdir out2.checkpoint/2rd_backup')
    if(id==0)call system('mkdir out1.output')
    if(id==0)call system('mkdir out1.output/profile')
    if(id==0)call system('mkdir out3.field_data')
    if(id==0)print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    if(id==0)print*,''


    if(id==0)print*,'**************************** Processing geometry *****************************'
    if(id==0)  print*,'.........'
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !initializing walls
    call MemAllocate_geometry(1)
    call set_walls
    
    !~~~~~~~~~~~~~~~~~~ dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! effective dimensions for an open channel (considering half-way bounceback)
    ! those parameters are important in inlet velocity distribution calculation
    la_z =  nzGlobal - 1   
    la_y =  nyGlobal - 1 - 0.5 - 0.5    ! half-way bounceback 0.5 + 0.5
    la_x =  nxGlobal - 1 - 0.5 - 0.5
    A_xy = la_x*la_y  ! cross-sectional area for an open duct
    !~~~~~~~~~~~~~~~~~~~~~~~~~

    if(id==0)then
        open(unit=11,file='out1.output/info.txt',status='replace')
        write(11,*)'Grid information:'
        write(11,"('nxGlobal = ', I4, ', nyGlobal = ', I4, ', nzGlobal = ', I4)") nxGlobal, nyGlobal, nzGlobal
        write(11,"('Inlet open cross sectional area = ', F14.2)") A_xy

        write(11,*)' Pore information:'
        write(11,"('Total number of pore nodes = ', I14)") pore_sum  

        write(11,"('Total number of effective pore nodes (middle section for analysis) = ', I14)") pore_sum_effective
    
        porosity_full = dble(pore_sum)/(la_x*la_y*la_z)
        porosity_effective = dble(pore_sum_effective)/(la_x*la_y)/dble(nzglobal - n_exclude_outlet - n_exclude_inlet)
        write(11,"('Full domain porosity = ', F6.4)") porosity_full
        write(11,"('Effective domain porosity (eclude inlet/outlet) = ', F6.4)") porosity_effective
        close(11)

        write(*,"(' Total number of pore nodes (excluding buffer layers) = ', I14)") pore_sum 
        write(*,"(' Inlet cross sectional area = ', F14.2)") A_xy
        print*, 'Total number of effective pore nodes (middle section for analysis) = '
        print*, pore_sum_effective
        write(*,"(' Full domain porosity = ', F6.4)") porosity_full
        write(*,"(' Effective domain porosity (eclude inlet/outlet)  = ', F6.4)") porosity_effective
    endif
    if(id==0)print*,'************************** End Processing geometry ***************************'
    if(id==0)print*,''


    if(id==0)print*,'************************ Processing fluid and flow info **********************'
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !allocate fluid memory
    call MemAllocate(1)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !parameters
    rt=3d0*la_nu+0.5d0      !relaxation time in collision, related to viscosity
    rti=1d0/rt
    la_nui=1d0/la_nu
    if(id==0)THEN
        write(*,"(' Fluid relaxation time = ', F8.6)") rt
    endif

    !MRT PARAMETERS, assuming viscosity is an constant
    omega = 1d0/(3d0*la_nu+0.5d0)
    s_nu =  omega

    if(mrt_para_preset==1)then
      !************ optimized for bounceback ************
      s_e =  omega
      s_e2 = omega
      s_pi = omega
      s_q =  8.0d0*(2.0d0-omega)/(8.0d0-omega)
      s_t = s_q
    elseif(mrt_para_preset==1)then
      !************ original preset, more stable ************
      s_e =  1.19d0
      s_e2 = 1.4d0
      s_pi = 1.4d0
      s_q= 1.2d0
      s_t = 1.98d0
    else
      !************ single relaxation time ************
      s_e =  omega
      s_e2 = omega
      s_pi = omega
      s_q =  omega
      s_t = omega
    endif

    ! body force/pressure gradient
    ! if pressure or velocity BC is used to drive the flow, force_z should 0
    force_Z=force_z0
    
    !outlet pressure (density) only used when outlet_BC==2
    rho_out = 1d0

    !open pressure or velocity inlet BC
    if(kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then    !non-periodic BC along flow direction (z)
        if(inlet_BC==1)then  !velocity inlet bc
            call initialization_open_velocity_inlet_bc
        elseif(inlet_BC==2)then   !pressure inlet bc           
            rho_in = rho_out + rho_drop
            if(id==0)write(*,"(' Inlet density (pressure inlet BC) = ', F6.4)") rho_in
        endif
    endif
    if(id==0)print*,'******************** End processing fluid and flow info **********************'
    if(id==0)print*,''

    if(d_vol_detail>0d0.and.inlet_BC==1.and.kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then 
        ntime_visual = dble(d_vol_detail * pore_sum) / flowrate
        ! due to AA pattern streaming, PDFs after odd step is stored in oppotite way.
        ! only use even step for outputs 
        if (MOD(ntime_visual,2) .eq. 1) then  
            ntime_visual = ntime_visual + 1   
        endif 
        if(id==0)print*,'Full VTK files timer is modified based on injected volume:', ntime_visual
    endif

    if(d_vol_monitor>0d0.and.inlet_BC==1.and.kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then
        ! due to AA pattern streaming, PDFs after odd step is stored in oppotite way.
        ! only use even step for outputs
        ntime_monitor = dble(d_vol_monitor * pore_sum) / flowrate
        if (MOD(ntime_monitor,2) .eq. 1) then  
            ntime_monitor = ntime_monitor + 1   
        endif 
        ntime_monitor_profile =  ntime_monitor_profile_ratio * ntime_monitor 
        if(id==0)print*,'Monitor timer is modified based on injected volume:', ntime_monitor
        if(id==0)print*,'Monitor_profile timer is modified based on injected volume:', ntime_monitor_profile
    endif

    monitor_previous_value=0d0
    monitor_current_value=0d0

    return
end subroutine initialization_basic

!~~~~~~~~~~~~~~~~~~~~~open velocity inlet BC initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine initialization_open_velocity_inlet_BC
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    real(kind=8) :: temp,temp3,x,y
    integer :: i,j, n_vin

    uin_avg_0 = Re * la_nu / char_length
    flowrate = uin_avg_0 * A_xy
    uin_avg = uin_avg_0
    if(id==0)write(*,"(' Inlet average velocity = ', F13.6)") uin_avg_0
    if(id==0)write(*,"(' Inlet flowrate = ', E13.6)") flowrate

    ! default uniform inlet velocity profile
    !$OMP PARALLEL DO private(i,x,y)
    do j=1,ny
        do i=1,nx
          x = idx*nx + i
          y = idy*ny + j
          w_in(i,j) = 0d0
          if(x>1.and.x<nxGlobal.and.y>1.and.y<nyGlobal)then  
            w_in(i,j) = uin_avg    
          endif
        enddo
    enddo

    n_vin = 1000
    ! analytical velocity profile
    call inlet_vel_profile_rectangular(uin_avg_0, n_vin)

    if(target_inject_pore_volume>0)then  ! stop simulation based on injected volume is enabled
        ntime_max = dble(target_inject_pore_volume * pore_sum) / flowrate
        if (MOD(ntime_max,2) .eq. 1) then  
            ntime_max = ntime_max + 1   
        endif 
        if(id==0)print*,'Maximum time step is modified based on target inject volume! ntime_max = '
        if(id==0)print*, ntime_max
    endif

return
end subroutine initialization_open_velocity_inlet_BC

!=====================================================================================================================================
!---------------------- initialization for new simulation - field variables ----------------------
!=====================================================================================================================================
subroutine initialization_new
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer:: i,j,k ,m,n
    real(kind=8) :: temp,x,y,z,rho1,ux1,uy1,uz1

    ntime0=1
    !initializing fluid
    !$OMP parallel DO private(i,j,x,y,z) collapse(2)
    do k=1-overlap,nz+overlap
        do j=1-overlap,ny+overlap
            do i=1-overlap,nx+overlap
                x = idx*nx + i
                y = idy*ny + j
                z = idz*nz + k
                u(i,j,k) = 0d0
                v(i,j,k) = 0d0
                w(i,j,k) = 0d0
                rho(i,j,k) = 1d0
            enddo
        enddo
    enddo

    call initialization_new_pdf

    return
end subroutine initialization_new

! initial particle distribution functions
subroutine initialization_new_pdf
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    integer:: i,j,k 
    real(kind=8) :: usqrt,udotc,rho1

    !$OMP parallel DO private(i,j,rho1,usqrt) collapse(2)
    do k=1-overlap,nz+overlap
        do j=1-overlap,ny+overlap
            do i=1-overlap,nx+overlap
                usqrt = u(i,j,k) * u(i,j,k) + v(i,j,k) * v(i,j,k) +  w(i,j,k) * w(i,j,k)
                rho1 = rho(i,j,k)
    
                f0(i,j,k) =   rho1 * w_equ_0 + rho1* w_equ_0*(- 1.5d0 * usqrt)
                f1(i,j,k) =   rho1 * w_equ_1 + rho1* w_equ_1*(3.0d0 *    u(i,j,k) + 4.5d0 * u(i,j,k) * u(i,j,k) - 1.5d0 * usqrt)
                f2(i,j,k) =   rho1 * w_equ_1 + rho1* w_equ_1*(3.0d0 * (-u(i,j,k)) + 4.5d0 * u(i,j,k) * u(i,j,k) - 1.5d0 * usqrt)
                f3(i,j,k) =   rho1 * w_equ_1 + rho1* w_equ_1*(3.0d0 *    v(i,j,k) + 4.5d0 * v(i,j,k) * v(i,j,k) - 1.5d0 * usqrt)
                f4(i,j,k) =   rho1 * w_equ_1 + rho1* w_equ_1*(3.0d0 * (-v(i,j,k)) + 4.5d0 * v(i,j,k) * v(i,j,k) - 1.5d0 * usqrt)
                f5(i,j,k) =   rho1 * w_equ_1 + rho1* w_equ_1*(3.0d0 *    w(i,j,k) + 4.5d0 * w(i,j,k) * w(i,j,k) - 1.5d0 * usqrt)
                f6(i,j,k) =   rho1 * w_equ_1 + rho1* w_equ_1*(3.0d0 * (-w(i,j,k)) + 4.5d0 * w(i,j,k) * w(i,j,k) - 1.5d0 * usqrt)
                f7(i,j,k)  =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * ( u(i,j,k) + v(i,j,k))  + 4.5d0 * ( u(i,j,k) + v(i,j,k)) * ( u(i,j,k) + v(i,j,k))    - 1.5d0 * usqrt)
                f8(i,j,k)  =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * (-u(i,j,k) + v(i,j,k))  + 4.5d0 * (-u(i,j,k) + v(i,j,k)) * (-u(i,j,k) + v(i,j,k))    - 1.5d0 * usqrt)
                f9(i,j,k)  =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * ( u(i,j,k) - v(i,j,k))  + 4.5d0 * ( u(i,j,k) - v(i,j,k)) * ( u(i,j,k) - v(i,j,k))    - 1.5d0 * usqrt)
                f10(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * (-u(i,j,k) - v(i,j,k))  + 4.5d0 * (-u(i,j,k) - v(i,j,k)) * (-u(i,j,k) - v(i,j,k))    - 1.5d0 * usqrt)
                f11(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * ( u(i,j,k) + w(i,j,k))  + 4.5d0 * ( u(i,j,k) + w(i,j,k)) * ( u(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                f12(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * (-u(i,j,k) + w(i,j,k))  + 4.5d0 * (-u(i,j,k) + w(i,j,k)) * (-u(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                f13(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * ( u(i,j,k) - w(i,j,k))  + 4.5d0 * ( u(i,j,k) - w(i,j,k)) * ( u(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
                f14(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * (-u(i,j,k) - w(i,j,k))  + 4.5d0 * (-u(i,j,k) - w(i,j,k)) * (-u(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
                f15(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * ( v(i,j,k) + w(i,j,k))  + 4.5d0 * ( v(i,j,k) + w(i,j,k)) * ( v(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                f16(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * (-v(i,j,k) + w(i,j,k))  + 4.5d0 * (-v(i,j,k) + w(i,j,k)) * (-v(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                f17(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * ( v(i,j,k) - w(i,j,k))  + 4.5d0 * ( v(i,j,k) - w(i,j,k)) * ( v(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
                f18(i,j,k) =  rho1 * w_equ_2 + rho1* w_equ_2*(3.0d0 * (-v(i,j,k) - w(i,j,k))  + 4.5d0 * (-v(i,j,k) - w(i,j,k)) * (-v(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
            enddo
        enddo
    enddo

    !~~~~~~~~~~~~~~~~~~~convective outflow BC ~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(outlet_BC==1)then 
        if(idz==npz-1)then
            !$OMP parallel DO private(i)
            do j=1-overlap,ny+overlap
                do i=1-overlap,nx+overlap
                    f_convec_bc(i,j,0) = f0(i,j,nz)
                    f_convec_bc(i,j,1) = f1(i,j,nz)
                    f_convec_bc(i,j,2) = f2(i,j,nz)
                    f_convec_bc(i,j,3) = f3(i,j,nz)
                    f_convec_bc(i,j,4) = f4(i,j,nz)
                    f_convec_bc(i,j,5) = f5(i,j,nz)
                    f_convec_bc(i,j,6) = f6(i,j,nz)
                    f_convec_bc(i,j,7) = f7(i,j,nz)
                    f_convec_bc(i,j,8) = f8(i,j,nz)
                    f_convec_bc(i,j,9) = f9(i,j,nz)
                    f_convec_bc(i,j,10) = f10(i,j,nz)
                    f_convec_bc(i,j,11) = f11(i,j,nz)
                    f_convec_bc(i,j,12) = f12(i,j,nz)
                    f_convec_bc(i,j,13) = f13(i,j,nz)
                    f_convec_bc(i,j,14) = f14(i,j,nz)
                    f_convec_bc(i,j,15) = f15(i,j,nz)
                    f_convec_bc(i,j,16) = f16(i,j,nz)
                    f_convec_bc(i,j,17) = f17(i,j,nz)
                    f_convec_bc(i,j,18) = f18(i,j,nz)
                enddo
            enddo
        endif
    endif

    return
end subroutine initialization_new_pdf


!=====================================================================================================================================
!---------------------- initialization for old simulation - field variables ----------------------
!=====================================================================================================================================
subroutine initialization_old
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    LOGICAL ALIVE
    character (len=30) :: flnm ! file name

    write(flnm,"('id',i4.4)")id

    INQUIRE(FILE='out2.checkpoint/'//trim(flnm),EXIST=ALIVE)

    if(alive)then

        if(id==0)print*,'loading checkpoint data!'

        open(unit=9+id, file='out2.checkpoint/'//trim(flnm), FORM='unformatted', status='old',access='stream')

        rewind(9+id)
        read(9+id)ntime0,force_z,rho_in
        !fluid PDF
        read(9+id)f0
        read(9+id)f1
        read(9+id)f2
        read(9+id)f3
        read(9+id)f4
        read(9+id)f5
        read(9+id)f6
        read(9+id)f7
        read(9+id)f8
        read(9+id)f9
        read(9+id)f10
        read(9+id)f11
        read(9+id)f12
        read(9+id)f13
        read(9+id)f14
        read(9+id)f15
        read(9+id)f16
        read(9+id)f17
        read(9+id)f18
        !~~~~~~~~~~~~~~~~~~~convective outflow BC ~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(outlet_BC==1)then 
            read(9+id)f_convec_bc
        endif
        if(id==0)print*,'Checkpoint data loaded!'
    else
        if(id==0)write(*,*)'Checkpoint data not found! Exiting program!'
        call MPI_Barrier(MPI_COMM_vgrid,ierr)
        call mpi_abort(MPI_COMM_vgrid,ierr)
    endif
    close(9+id)
    
    ntime = ntime0

    return
end subroutine initialization_old


!=====================================================================================================================================
!---------------------- memory allocate/deallocate ----------------------
!=====================================================================================================================================
!*************** geometry related memory allocate/deallocate ********************************
subroutine MemAllocate_geometry(flag)
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: FLAG

    IF(FLAG == 1) THEN
        allocate(walls(1-overlap_walls:nx+overlap_walls,1-overlap_walls:ny+overlap_walls,1-overlap_walls:nz+overlap_walls))
        allocate(walls_global(1:nxGlobal,1:nyGlobal,1:nzGlobal))
    else
        deallocate(walls)
        deallocate (walls_global)
    endif

    return
end subroutine MemAllocate_geometry

!************* fluid flow related memory allocate/deallocate ******************************
subroutine MemAllocate(flag)
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER, INTENT (IN) :: FLAG
    integer:: lx,ly,lz,n,lap_t,overlap_temp

    IF(FLAG == 1) THEN
        allocate( &
            u(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            v(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            w(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            rho(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&                     
            w_in(1-overlap:nx+overlap,1-overlap:ny+overlap),&
            f0(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f1(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f2(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f3(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f4(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f5(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f6(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f7(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f8(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f9(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f10(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f11(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f12(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f13(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f14(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f15(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f16(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f17(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            f18(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap)&
            )      

        !convective BC
        if(outlet_BC==1)then 
            allocate(f_convec_bc(1-overlap:nx+overlap,1-overlap:ny+overlap,0:18))
        endif

        !-MPI BUFFER
        allocate( &
            send_pdf_xP(1 :ny  ,1 :nz  ,5),&
            send_pdf_xM(1 :ny  ,1 :nz  ,5),&
            recv_pdf_xP(1 :ny  ,1 :nz  ,5),&
            recv_pdf_xM(1 :ny  ,1 :nz  ,5),&
            send_pdf_yP(1 :nx  ,1 :nz  ,5),&
            send_pdf_yM(1 :nx  ,1 :nz  ,5),&
            recv_pdf_yP(1 :nx  ,1 :nz  ,5),&
            recv_pdf_yM(1 :nx  ,1 :nz  ,5),&
            send_pdf_zP(1 :nx  ,1 :ny  ,5),&
            send_pdf_zM(1 :nx  ,1 :ny  ,5),&
            recv_pdf_zP(1 :nx  ,1 :ny  ,5),&
            recv_pdf_zM(1 :nx  ,1 :ny  ,5))

        allocate( &
            send_pdf_yPzP(1 :nx  ),&
            send_pdf_yMzP(1 :nx  ),&
            send_pdf_yPzM(1 :nx  ),&
            send_pdf_yMzM(1 :nx  ),&
            recv_pdf_yPzP(1 :nx  ),&
            recv_pdf_yMzP(1 :nx  ),&
            recv_pdf_yPzM(1 :nx  ),&
            recv_pdf_yMzM(1 :nx  ),&

            send_pdf_xPzP(1 :ny  ),&
            send_pdf_xMzP(1 :ny  ),&
            send_pdf_xPzM(1 :ny  ),&
            send_pdf_xMzM(1 :ny  ),&
            recv_pdf_xPzP(1 :ny  ),&
            recv_pdf_xMzP(1 :ny  ),&
            recv_pdf_xPzM(1 :ny  ),&
            recv_pdf_xMzM(1 :ny  ),&

            send_pdf_xPyP(1 :nz  ),&
            send_pdf_xMyP(1 :nz  ),&
            send_pdf_xPyM(1 :nz  ),&
            send_pdf_xMyM(1 :nz  ),&
            recv_pdf_xPyP(1 :nz  ),&
            recv_pdf_xMyP(1 :nz  ),&
            recv_pdf_xPyM(1 :nz  ),&
            recv_pdf_xMyM(1 :nz  ))

        !mpi buffer size
        isize_pdf_x = 5*1*(ny)*(nz)   ! 10 varaibles, 1 layer
        isize_pdf_y = 5*1*(nx)*(nz)
        isize_pdf_z = 5*1*(ny)*(nx)
        isize_pdf_ex = nx !x edge
        isize_pdf_ey = ny
        isize_pdf_ez = nz

        !MPI status
        allocate(MPI_STAT(MPI_STATUS_SIZE,4), MPI_ESTAT(MPI_STATUS_SIZE,8))

        !temporary arrays
        allocate(fl(nz),pre(nz))
        tk_isize=2*nz+1
        allocate(tk(tk_isize))  ! isize = 5*nz+1          !fl, pre  +  umax 

    else
        deallocate(u,v,w,w_in,rho)
        deallocate(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18)
        
        
        if(outlet_bc==1)then
            deallocate(f_convec_bc)
        endif
        
        deallocate(send_pdf_xP,send_pdf_xM,send_pdf_yP,send_pdf_yM,send_pdf_zP,send_pdf_zM,&
        recv_pdf_xM,recv_pdf_xP,recv_pdf_yM,recv_pdf_yP,recv_pdf_zM,recv_pdf_zP)

        deallocate(send_pdf_yPzP,send_pdf_yMzP,send_pdf_yPzM,send_pdf_yMzM,send_pdf_xPzP,&
        send_pdf_xMzP,send_pdf_xPzM,send_pdf_xMzM,send_pdf_xPyP,send_pdf_xMyP,send_pdf_xPyM,send_pdf_xMyM)

        deallocate(recv_pdf_yPzP,recv_pdf_yMzP,recv_pdf_yPzM,recv_pdf_yMzM,recv_pdf_xPzP,&
        recv_pdf_xMzP,recv_pdf_xPzM,recv_pdf_xMzM,recv_pdf_xPyP,recv_pdf_xMyP,recv_pdf_xPyM,recv_pdf_xMyM)
        
        deallocate(fl,pre,tk)

    endif

    return
end subroutine MemAllocate
