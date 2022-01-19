#include "./preprocessor.h"
!=====================================================================================================================================
!---------------------- initialization basic ----------------------
!=====================================================================================================================================
subroutine initialization_basic_multi 
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer:: i,j,k, lx,ly,lz,n,z,n_vin, n_small,n_large    !n_vin number of terms in inlet velocity profile
    real(kind=8) :: temp,TEMP3,x,y, dp_small, dp_large, A_xy_full
    LOGICAL :: ALIVE

    call read_parameter_multi
    
    call set_MPI

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !create folders
    if(id==0)print*,''
    if(id==0)print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    if(id==0)print*,'Creating directories if not exist'
    if(id==0)call system('mkdir out2.checkpoint')
    if(id==0)call system('mkdir out2.checkpoint/2rd_backup')
    if(id==0)call system('mkdir out1.output')
    if(id==0)call system('mkdir out1.output/profile')
    if(id==0)call system('mkdir out3.field_data')
    if(id==0)call system('mkdir out3.field_data/phase_distribution')
    if(id==0)call system('mkdir out3.field_data/full_flow_field') 
    if(id==0)print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    if(id==0)print*,''


    if(id==0)print*,'**************************** Processing geometry *****************************'
    if(id==0)  print*,'.........'
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !initializing walls
    call MemAllocate_geometry(1)
    call set_walls
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !geometry related preprocessing
    if(geometry_preprocess_cmd==0)then
        if(id==0)print*,'------ Start processing boundary nodes info ------'
        call geometry_preprocessing_new
        if(id==0)print*,'------ End processing boundary nodes info --------'
    elseif(geometry_preprocess_cmd==1)then
        INQUIRE(FILE=trim(geo_boundary_file_path),EXIST=ALIVE)
        if(alive)then
            if(id==0)print*,'------ Start loading boundary nodes info ---------'
            call geometry_preprocessing_load
            if(id==0)print*,'------ End loading boundary nodes info -----------'
        else
            if(id==0)print*,'Error! No precomputed boundary info file found! Exiting program!'
            call MPI_Barrier(MPI_COMM_vgrid,ierr)
            call mpi_abort(MPI_COMM_vgrid,ierr)
        endif
    else
        if(id==0)then
            print*, 'Wrong value of geometry_preprocess_cmd! Stop program!'
        endif
        call MPI_Barrier(MPI_COMM_vgrid,ierr)
        call mpi_abort(MPI_COMM_vgrid,ierr)
    endif
    !call check_geometry_linked_data

    !~~~~~~~~~~~~~~~~~~ dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! la_x,y,z only used in fluid displacement simulation
    ! channel walls are considered below for effective sample volume and cross sectional area
    la_z =  nzGlobal - 1       
    la_y =  nyGlobal - 1 - 0.5 - 0.5    ! half-way bounceback 0.5 + 0.5
    la_x =  nxGlobal - 1 - 0.5 - 0.5
    A_xy = la_x*la_y  !cross-sectional area for an open duct
    volume_sample = A_xy*la_z     !domain volume
    !~~~~~~~~~~~~~~~~~~~~~~~~~

    if(id==0)then
        open(unit=11,file='out1.output/info.txt',status='replace')
        write(11,*)'Grid information:'
        write(11,"('nxGlobal = ', I4, ', nyGlobal = ', I4, ', nzGlobal = ', I4)") nxGlobal, nyGlobal, nzGlobal
        write(11,"('Inlet open cross sectional area = ', F14.2)") A_xy

        write(11,*)' Pore information:'
        write(11,"('Total number of pore nodes = ', I14)") pore_sum  
        write(11,"('Total number of solid boundary nodes = ', I14)") num_solid_boundary_global
        write(11,"('Total number of fluid boundary nodes = ', I14)") num_fluid_boundary_global

        write(11,"('Total number of effective pore nodes (excluding inlet/outlet) = ', I14)") pore_sum_effective
    
        porosity_full = dble(pore_sum)/((nxGlobal-2)*(nyGlobal-2)*(nzGlobal))
        porosity_effective = dble(pore_sum_effective)/((nxGlobal-2)*(nyGlobal-2))/dble(nzglobal - n_exclude_outlet - n_exclude_inlet)
        write(11,"('Full domain porosity = ', F6.4)") porosity_full
        write(11,"('Effective domain porosity  (excluding inlet/outlet) = ', F6.4)") porosity_effective
        close(11)

        write(*,"(' Total number of pore nodes = ', I14)") pore_sum 
        write(*,"(' Inlet open cross sectional area = ', F14.2)") A_xy
        print*, 'Total number of effective pore nodes (excluding inlet/outlet) = '
        print*, pore_sum_effective
        write(*,"(' Full domain porosity = ', F6.4)") porosity_full
        write(*,"(' Effective domain porosity  (excluding inlet/outlet) = ', F6.4)") porosity_effective
    endif
    if(id==0)print*,'************************** End Processing geometry ***************************'
    if(id==0)print*,''




    if(id==0)print*,'************************ Processing fluid and flow info **********************'
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !allocate fluid memory
    call MemAllocate_multi(1)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !parameters
    rt1=3d0*la_nu1+0.5d0      !relaxation time in collision, related to viscosity
    rti1=1d0/rt1
    rt2=3d0*la_nu2+0.5d0      
    rti2=2d0/rt2
    la_nui1=1d0/la_nu1
    la_nui2=1d0/la_nu2
    if(id==0)THEN
        write(*,"(' Fluid 1 relaxation time = ', F8.6)") rt1
        write(*,"(' Fluid 2 relaxation time = ', F8.6)") rt2
    endif

    !injection parameters
    phi_inlet = 2d0*sa_inject-1d0     !inlet BC order parameter, consistent with injection fluid

    ! body force, pressure gradient
    ! if pressure or velocity BC is used to drive the flow, force_z should be set to 0
    force_Z=force_z0
    
    !outlet pressure (density) only used when outlet_BC==2
    rho_out = 1d0
    if(id==0)THEN
        write(*,"(' Inlet boundary nodes order parameter = ', F6.4)") phi_inlet
    endif

    !open pressure or velocity inlet BC
    if(kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then    !non-periodic BC along flow direction (z)
        if(inlet_BC==1)then  !velocity inlet bc
            call initialization_open_velocity_inlet_bc
        elseif(inlet_BC==2)then   !pressure inlet bc
            p_gradient = -force_z0/3d0         !pressure gradient
            rho_in = rho_out- p_gradient*nzglobal
            if(id==0)write(*,"(' Inlet density (pressure inlet BC) = ', F6.4)") rho_in
        endif
    endif
    if(id==0)print*,'******************** End processing fluid and flow info **********************'
    if(id==0)print*,''

    if(d_vol_animation>0d0.and.inlet_BC==1.and.kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then 
        ntime_animation = dble(d_vol_animation * pore_sum) / flowrate
        ! due to AA pattern streaming, PDFs after odd step is stored in oppotite way.
        ! only use even step for outputs 
        if (MOD(ntime_animation,2) .eq. 1) then  
            ntime_animation = ntime_animation + 1   
        endif 
        if(id==0)print*,'Animation VTK files timer is modified based on injected volume:', ntime_animation
    endif

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
    kinetic_energy(:)=0d0

    return
end subroutine initialization_basic_multi

!~~~~~~~~~~~~~~~~~~~~~open velocity inlet BC initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine initialization_open_velocity_inlet_BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    real(kind=8) :: x, y
    integer :: i,j,k,n_vin

    uin_avg_0 = ca_0*gamma/La_nu1            !common definition
    uin_avg = uin_avg_0
    flowrate = uin_avg_0 * A_xy
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
subroutine initialization_new_multi
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer:: i,j,k ,m,n
    real(kind=8) :: temp,x,y,z,tempa,tempc,usqrt,udotc,random,rho1,rho2,ux1,uy1,uz1,z1,z2

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
                cn_x(i,j,k) = 0d0
                cn_y(i,j,k) = 0d0
                cn_z(i,j,k) = 0d0
                c_norm(i,j,k) = 0d0
                curv(i,j,k) = 0d0
            enddo
        enddo
    enddo

    call random_seed()
    ! Initial fluid phase distribution
    !$OMP parallel DO private(i,j,x,y,z,random) collapse(2)
    do k=1-overlap_phi,nz+overlap_phi
        do j=1-overlap_phi,ny+overlap_phi
            do i=1-overlap_phi,nx+overlap_phi              
                x = idx*nx + i
                y = idy*ny + j
                z = idz*nz + k
                if(initial_fluid_distribution_option==1)then ! drainage
                    phi(i,j,k) = -1d0
                    if(z<=interface_z0)then
                        phi(i,j,k) = 1d0
                    endif
                elseif(initial_fluid_distribution_option==2)then    ! imbibition
                    phi(i,j,k) = 1d0
                    if(z<=interface_z0)then
                        phi(i,j,k) = -1d0
                    endif
                elseif(initial_fluid_distribution_option==3)then   ! contact angle measurement
                    phi(i,j,k) = -1d0
                    if((x-(nxglobal+1)*0.0d0)**2+(z-(nzglobal+1)*0.5d0)**2+(y-(nyglobal+1)*0.5d0)**2<=interface_Z0**2)then
                        phi(i,j,k) = 1d0
                    endif
                elseif(initial_fluid_distribution_option==4)then   ! contact angle measurement
                    phi(i,j,k) = 1d0
                    if((x-(nxglobal+1)*0.0d0)**2+(z-(nzglobal+1)*0.5d0)**2+(y-(nyglobal+1)*0.5d0)**2<=interface_Z0**2)then
                        phi(i,j,k) = -1d0
                    endif
                elseif(initial_fluid_distribution_option==5)then   ! surface tension measurement
                    phi(i,j,k) = -1d0
                    if((x-(nxglobal+1)*0.5d0)**2+(z-(nzglobal+1)*0.5d0)**2+(y-(nyglobal+1)*0.5d0)**2<=interface_Z0**2)then
                        phi(i,j,k) = 1d0
                    endif
                elseif(initial_fluid_distribution_option==6)then   !randomly distributed: steady state relative perm measurement
                    phi(i,j,k) = 1d0
                    call RANDOM_NUMBER(random)
                    if(random>Sa_target)then
                        phi(i,j,k) = -1d0
                    endif
                else   
                    if(id==0)print*,'Input parameter initial_fluid_distribution_option error! Stop program!!!'
                    call MPI_Barrier(MPI_COMM_vgrid,ierr)
                    call mpi_abort(MPI_COMM_vgrid,ierr)                 
                endif     
            enddo
        enddo
    enddo
    if(kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then    !non-periodic BC along flow direction (z)
        !$OMP parallel DO private(i,j,x,y,z) collapse(2)
        do k=1-overlap_phi,nz+overlap_phi
            do j=1-overlap_phi,ny+overlap_phi
                do i=1-overlap_phi,nx+overlap_phi
                    x = idx*nx + i
                    y = idy*ny + j
                    z = idz*nz + k
                    if(z<=0)then
                        phi(i,j,k) = phi_inlet   ! phi_inlet = 2d0*sa_inject-1d0, inlet BC order parameter, consistent with injecting fluid saturation
                    endif                     
                enddo
            enddo
        enddo
    endif

    call initialization_new_multi_pdf

    if(steady_state_option==2)then
        !using the phase field difference to indicate steady state
        !$OMP parallel DO private(i,j) collapse(2)
        do k=1-overlap_phi,nz+overlap_phi
            do j=1-overlap_phi,ny+overlap_phi
                do i=1-overlap_phi,nx+overlap_phi
                    phi_old(i,j,k) = phi(i,j,k)
                enddo
            enddo
        enddo
    endif

    return
end subroutine initialization_new_multi

! initial particle distribution functions
subroutine initialization_new_multi_pdf
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer:: i,j,k ,m
    real(kind=8) :: temp,x,y,z,tempa,tempc,usqrt,udotc,random,rho1,rho2,ux1,uy1,uz1,z1,z2

    !$OMP parallel DO private(i,j,rho1,rho2,usqrt) collapse(2)
    do k=1-overlap,nz+overlap
        do j=1-overlap,ny+overlap
            do i=1-overlap,nx+overlap
                usqrt = u(i,j,k) * u(i,j,k) + v(i,j,k) * v(i,j,k) +  w(i,j,k) * w(i,j,k)
                rho1 = rho(i,j,k)*(1.0d0+phi(i,j,k))*0.5d0
                rho2 = rho(i,j,k)*(1.0d0-phi(i,j,k))*0.5d0

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

                g0(i,j,k) =   rho2 * w_equ_0 + rho2* w_equ_0*(- 1.5d0 * usqrt)
                g1(i,j,k) =   rho2 * w_equ_1 + rho2* w_equ_1*(3.0d0 *    u(i,j,k) + 4.5d0 * u(i,j,k) * u(i,j,k) - 1.5d0 * usqrt)
                g2(i,j,k) =   rho2 * w_equ_1 + rho2* w_equ_1*(3.0d0 * (-u(i,j,k)) + 4.5d0 * u(i,j,k) * u(i,j,k) - 1.5d0 * usqrt)
                g3(i,j,k) =   rho2 * w_equ_1 + rho2* w_equ_1*(3.0d0 *    v(i,j,k) + 4.5d0 * v(i,j,k) * v(i,j,k) - 1.5d0 * usqrt)
                g4(i,j,k) =   rho2 * w_equ_1 + rho2* w_equ_1*(3.0d0 * (-v(i,j,k)) + 4.5d0 * v(i,j,k) * v(i,j,k) - 1.5d0 * usqrt)
                g5(i,j,k) =   rho2 * w_equ_1 + rho2* w_equ_1*(3.0d0 *    w(i,j,k) + 4.5d0 * w(i,j,k) * w(i,j,k) - 1.5d0 * usqrt)
                g6(i,j,k) =   rho2 * w_equ_1 + rho2* w_equ_1*(3.0d0 * (-w(i,j,k)) + 4.5d0 * w(i,j,k) * w(i,j,k) - 1.5d0 * usqrt)
                g7(i,j,k)  =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * ( u(i,j,k) + v(i,j,k))  + 4.5d0 * ( u(i,j,k) + v(i,j,k)) * ( u(i,j,k) + v(i,j,k))    - 1.5d0 * usqrt)
                g8(i,j,k)  =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * (-u(i,j,k) + v(i,j,k))  + 4.5d0 * (-u(i,j,k) + v(i,j,k)) * (-u(i,j,k) + v(i,j,k))    - 1.5d0 * usqrt)
                g9(i,j,k)  =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * ( u(i,j,k) - v(i,j,k))  + 4.5d0 * ( u(i,j,k) - v(i,j,k)) * ( u(i,j,k) - v(i,j,k))    - 1.5d0 * usqrt)
                g10(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * (-u(i,j,k) - v(i,j,k))  + 4.5d0 * (-u(i,j,k) - v(i,j,k)) * (-u(i,j,k) - v(i,j,k))    - 1.5d0 * usqrt)
                g11(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * ( u(i,j,k) + w(i,j,k))  + 4.5d0 * ( u(i,j,k) + w(i,j,k)) * ( u(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                g12(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * (-u(i,j,k) + w(i,j,k))  + 4.5d0 * (-u(i,j,k) + w(i,j,k)) * (-u(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                g13(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * ( u(i,j,k) - w(i,j,k))  + 4.5d0 * ( u(i,j,k) - w(i,j,k)) * ( u(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
                g14(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * (-u(i,j,k) - w(i,j,k))  + 4.5d0 * (-u(i,j,k) - w(i,j,k)) * (-u(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
                g15(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * ( v(i,j,k) + w(i,j,k))  + 4.5d0 * ( v(i,j,k) + w(i,j,k)) * ( v(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                g16(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * (-v(i,j,k) + w(i,j,k))  + 4.5d0 * (-v(i,j,k) + w(i,j,k)) * (-v(i,j,k) + w(i,j,k))    - 1.5d0 * usqrt)
                g17(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * ( v(i,j,k) - w(i,j,k))  + 4.5d0 * ( v(i,j,k) - w(i,j,k)) * ( v(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
                g18(i,j,k) =  rho2 * w_equ_2 + rho2* w_equ_2*(3.0d0 * (-v(i,j,k) - w(i,j,k))  + 4.5d0 * (-v(i,j,k) - w(i,j,k)) * (-v(i,j,k) - w(i,j,k))    - 1.5d0 * usqrt)
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

                    g_convec_bc(i,j,0) = g0(i,j,nz)
                    g_convec_bc(i,j,1) = g1(i,j,nz)
                    g_convec_bc(i,j,2) = g2(i,j,nz)
                    g_convec_bc(i,j,3) = g3(i,j,nz)
                    g_convec_bc(i,j,4) = g4(i,j,nz)
                    g_convec_bc(i,j,5) = g5(i,j,nz)
                    g_convec_bc(i,j,6) = g6(i,j,nz)
                    g_convec_bc(i,j,7) = g7(i,j,nz)
                    g_convec_bc(i,j,8) = g8(i,j,nz)
                    g_convec_bc(i,j,9) = g9(i,j,nz)
                    g_convec_bc(i,j,10) = g10(i,j,nz)
                    g_convec_bc(i,j,11) = g11(i,j,nz)
                    g_convec_bc(i,j,12) = g12(i,j,nz)
                    g_convec_bc(i,j,13) = g13(i,j,nz)
                    g_convec_bc(i,j,14) = g14(i,j,nz)
                    g_convec_bc(i,j,15) = g15(i,j,nz)
                    g_convec_bc(i,j,16) = g16(i,j,nz)
                    g_convec_bc(i,j,17) = g17(i,j,nz)
                    g_convec_bc(i,j,18) = g18(i,j,nz)

                    phi_convec_bc(i,j) = phi(i,j,nz)
                enddo
            enddo
        endif
    endif

    return
end subroutine initialization_new_multi_pdf



!=====================================================================================================================================
!---------------------- initialization for old simulation - field variables ----------------------
!=====================================================================================================================================
subroutine initialization_old_multi
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
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
        read(9+id)g0
        read(9+id)g1
        read(9+id)g2
        read(9+id)g3
        read(9+id)g4
        read(9+id)g5
        read(9+id)g6
        read(9+id)g7
        read(9+id)g8
        read(9+id)g9
        read(9+id)g10
        read(9+id)g11
        read(9+id)g12
        read(9+id)g13
        read(9+id)g14
        read(9+id)g15
        read(9+id)g16
        read(9+id)g17
        read(9+id)g18
        read(9+id)phi
        !~~~~~~~~~~~~~~~~~~~convective outflow BC ~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(outlet_BC==1)then 
            read(9+id)f_convec_bc
            read(9+id)g_convec_bc
            read(9+id)phi_convec_bc
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
end subroutine initialization_old_multi



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
subroutine MemAllocate_multi(flag)
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
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
            curv(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
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
            f18(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g0(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g1(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g2(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g3(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g4(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g5(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g6(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g7(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g8(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g9(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g10(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g11(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g12(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g13(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g14(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g15(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g16(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g17(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap),&
            g18(1-overlap:nx+overlap,1-overlap:ny+overlap,1-overlap:nz+overlap)&
            )

        overlap_temp = 2  ! used in CSF calculation
        allocate( &            
            cn_x(1-overlap_temp:nx+overlap_temp,1-overlap_temp:ny+overlap_temp,1-overlap_temp:nz+overlap_temp),&
            cn_y(1-overlap_temp:nx+overlap_temp,1-overlap_temp:ny+overlap_temp,1-overlap_temp:nz+overlap_temp),&
            cn_z(1-overlap_temp:nx+overlap_temp,1-overlap_temp:ny+overlap_temp,1-overlap_temp:nz+overlap_temp),&     
            c_norm(1-overlap_temp:nx+overlap_temp,1-overlap_temp:ny+overlap_temp,1-overlap_temp:nz+overlap_temp)& 
            )

        !convective BC
        if(outlet_BC==1)then 
            allocate( &
                phi_convec_bc(1-overlap:nx+overlap,1-overlap:ny+overlap),&
                g_convec_bc(1-overlap:nx+overlap,1-overlap:ny+overlap,0:18),&
                f_convec_bc(1-overlap:nx+overlap,1-overlap:ny+overlap,0:18))
        endif

        !three layers used in cn_x,y,z calculation (without MPI transfer of cn_x,y,z)    
        allocate(phi(1-overlap_phi:nx+overlap_phi,1-overlap_phi:ny+overlap_phi,1-overlap_phi:nz+overlap_phi))

        !phase field steady state
        if(steady_state_option==2)then
            allocate(phi_old(1-overlap_phi:nx+overlap_phi,1-overlap_phi:ny+overlap_phi,1-overlap_phi:nz+overlap_phi))   !using the phase field difference to indicate steady state
        endif

        !-MPI BUFFER
        allocate( &
            send_pdf_xP(1 :ny  ,1 :nz  ,10),&
            send_pdf_xM(1 :ny  ,1 :nz  ,10),&
            recv_pdf_xP(1 :ny  ,1 :nz  ,10),&
            recv_pdf_xM(1 :ny  ,1 :nz  ,10),&
            send_pdf_yP(1 :nx  ,1 :nz  ,10),&
            send_pdf_yM(1 :nx  ,1 :nz  ,10),&
            recv_pdf_yP(1 :nx  ,1 :nz  ,10),&
            recv_pdf_yM(1 :nx  ,1 :nz  ,10),&
            send_pdf_zP(1 :nx  ,1 :ny  ,10),&
            send_pdf_zM(1 :nx  ,1 :ny  ,10),&
            recv_pdf_zP(1 :nx  ,1 :ny  ,10),&
            recv_pdf_zM(1 :nx  ,1 :ny  ,10))

        allocate( &
            send_pdf_yPzP(1 :2*nx  ),&
            send_pdf_yMzP(1 :2*nx  ),&
            send_pdf_yPzM(1 :2*nx  ),&
            send_pdf_yMzM(1 :2*nx  ),&
            recv_pdf_yPzP(1 :2*nx  ),&
            recv_pdf_yMzP(1 :2*nx  ),&
            recv_pdf_yPzM(1 :2*nx  ),&
            recv_pdf_yMzM(1 :2*nx  ),&

            send_pdf_xPzP(1 :2*ny  ),&
            send_pdf_xMzP(1 :2*ny  ),&
            send_pdf_xPzM(1 :2*ny  ),&
            send_pdf_xMzM(1 :2*ny  ),&
            recv_pdf_xPzP(1 :2*ny  ),&
            recv_pdf_xMzP(1 :2*ny  ),&
            recv_pdf_xPzM(1 :2*ny  ),&
            recv_pdf_xMzM(1 :2*ny  ),&

            send_pdf_xPyP(1 :2*nz  ),&
            send_pdf_xMyP(1 :2*nz  ),&
            send_pdf_xPyM(1 :2*nz  ),&
            send_pdf_xMyM(1 :2*nz  ),&
            recv_pdf_xPyP(1 :2*nz  ),&
            recv_pdf_xMyP(1 :2*nz  ),&
            recv_pdf_xPyM(1 :2*nz  ),&
            recv_pdf_xMyM(1 :2*nz  ))

        !mpi buffer size
        isize_pdf_x = 10*1*(ny)*(nz)   ! 10 varaibles, 1 layer
        isize_pdf_y = 10*1*(nx)*(nz)
        isize_pdf_z = 10*1*(ny)*(nx)
        isize_pdf_ex = 2*nx !x edge
        isize_pdf_ey = 2*ny
        isize_pdf_ez = 2*nz
        isize_phi_x = (nz+overlap_phi+overlap_phi)*(ny+overlap_phi+overlap_phi)*overlap_phi  ! 
        isize_phi_y = (nx+overlap_phi+overlap_phi)*(nz+overlap_phi+overlap_phi)*overlap_phi
        isize_phi_z = (nx+overlap_phi+overlap_phi)*(ny+overlap_phi+overlap_phi)*overlap_phi

        !MPI status
        allocate(MPI_STAT(MPI_STATUS_SIZE,4), MPI_ESTAT(MPI_STATUS_SIZE,8))

        !for phi MPI communication, three layers required for cn_x,y,z calculation

        allocate( &     !face
            send_phi_xM(overlap_phi,1 :ny  ,1 :nz  ),&
            send_phi_xP(overlap_phi,1 :ny  ,1 :nz  ),&
            recv_phi_xM(overlap_phi,1 :ny  ,1 :nz  ),&
            recv_phi_xP(overlap_phi,1 :ny  ,1 :nz  ),&
            send_phi_yM(1 :nx  ,overlap_phi,1 :nz  ),&
            send_phi_yP(1 :nx  ,overlap_phi,1 :nz  ),&
            recv_phi_yM(1 :nx  ,overlap_phi,1 :nz  ),&
            recv_phi_yP(1 :nx  ,overlap_phi,1 :nz  ),&
            send_phi_zM(1 :nx  ,1 :ny  ,overlap_phi),&
            send_phi_zP(1 :nx  ,1 :ny  ,overlap_phi),&
            recv_phi_zM(1 :nx  ,1 :ny  ,overlap_phi),&
            recv_phi_zP(1 :nx  ,1 :ny  ,overlap_phi))

        !buffer length: 1 varaibles, overlap_phi layer
        isize_phi_x = overlap_phi*(ny)*(nz)
        isize_phi_y = overlap_phi*(nx)*(nz)
        isize_phi_z = overlap_phi*(ny)*(nx)

        allocate( &    !edge
            send_phi_yPzP(1:nx, overlap_phi, overlap_phi),&
            send_phi_yPzM(1:nx, overlap_phi, overlap_phi),&
            send_phi_yMzP(1:nx, overlap_phi, overlap_phi),&
            send_phi_yMzM(1:nx, overlap_phi, overlap_phi),&
            send_phi_xPzP(overlap_phi, 1:ny, overlap_phi),&
            send_phi_xPzM(overlap_phi, 1:ny, overlap_phi),&
            send_phi_xMzP(overlap_phi, 1:ny, overlap_phi),&
            send_phi_xMzM(overlap_phi, 1:ny, overlap_phi),&
            send_phi_xPyP(overlap_phi, overlap_phi, 1:nz),&
            send_phi_xPyM(overlap_phi, overlap_phi, 1:nz),&
            send_phi_xMyP(overlap_phi, overlap_phi, 1:nz),&
            send_phi_xMyM(overlap_phi, overlap_phi, 1:nz),&

            recv_phi_yPzP(1:nx, overlap_phi, overlap_phi),&
            recv_phi_yPzM(1:nx, overlap_phi, overlap_phi),&
            recv_phi_yMzP(1:nx, overlap_phi, overlap_phi),&
            recv_phi_yMzM(1:nx, overlap_phi, overlap_phi),&
            recv_phi_xPzP(overlap_phi, 1:ny, overlap_phi),&
            recv_phi_xPzM(overlap_phi, 1:ny, overlap_phi),&
            recv_phi_xMzP(overlap_phi, 1:ny, overlap_phi),&
            recv_phi_xMzM(overlap_phi, 1:ny, overlap_phi),&
            recv_phi_xPyP(overlap_phi, overlap_phi, 1:nz),&
            recv_phi_xPyM(overlap_phi, overlap_phi, 1:nz),&
            recv_phi_xMyP(overlap_phi, overlap_phi, 1:nz),&
            recv_phi_xMyM(overlap_phi, overlap_phi, 1:nz))

        !edge buffer length: 1 varaibles
        isize_phi_ex = overlap_phi*overlap_phi*(nx)
        isize_phi_ey = overlap_phi*overlap_phi*(ny)
        isize_phi_ez = overlap_phi*overlap_phi*(nz)

        allocate( &      !corner
            send_phi_xPyPzP(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xPyPzM(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xPyMzP(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xPyMzM(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xMyPzP(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xMyPzM(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xMyMzP(overlap_phi, overlap_phi, overlap_phi),&
            send_phi_xMyMzM(overlap_phi, overlap_phi, overlap_phi),&

            recv_phi_xPyPzP(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xPyPzM(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xPyMzP(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xPyMzM(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xMyPzP(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xMyPzM(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xMyMzP(overlap_phi, overlap_phi, overlap_phi),&
            recv_phi_xMyMzM(overlap_phi, overlap_phi, overlap_phi))

        !corner buffer size
        isize_phi_corner = overlap_phi*overlap_phi*overlap_phi

        !temporary arrays
        allocate(fl1(nz),fl2(nz),pre(nz),mass1(nz),mass2(nz),vol1(nz),vol2(nz))
        tk_isize=7*nz+3
        allocate(tk(tk_isize))  ! isize = 5*nz+1          !fl1, fl2, pre,mass1,mass2,vol1,vol2  +  umax + usq1 + usq2 (kinetic energy)

    else
        deallocate(u,v,w,w_in,rho,cn_x,cn_y,cn_z,phi,curv,c_norm)
        deallocate(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
        
        if(steady_state_option==2)then
            deallocate(phi_old)   !steady_state_option based on phase field
        endif
        
        if(outlet_bc==1)then
            deallocate(f_convec_bc,g_convec_bc,phi_convec_bc)
        endif
        
        deallocate( send_pdf_xP,send_pdf_xM,send_pdf_yP,send_pdf_yM,send_pdf_zP,send_pdf_zM,recv_pdf_xM,recv_pdf_xP,recv_pdf_yM,recv_pdf_yP,recv_pdf_zM,recv_pdf_zP)
        deallocate( send_phi_xM,send_phi_xP,send_phi_yM,send_phi_yP,send_phi_zM,send_phi_zP,recv_phi_xM,recv_phi_xP,recv_phi_yM,recv_phi_yP,recv_phi_zM,recv_phi_zP)
        deallocate(send_pdf_yPzP,send_pdf_yMzP,send_pdf_yPzM,send_pdf_yMzM,send_pdf_xPzP,send_pdf_xMzP,send_pdf_xPzM,send_pdf_xMzM,send_pdf_xPyP,send_pdf_xMyP,send_pdf_xPyM,send_pdf_xMyM)
        deallocate(recv_pdf_yPzP,recv_pdf_yMzP,recv_pdf_yPzM,recv_pdf_yMzM,recv_pdf_xPzP,recv_pdf_xMzP,recv_pdf_xPzM,recv_pdf_xMzM,recv_pdf_xPyP,recv_pdf_xMyP,recv_pdf_xPyM,recv_pdf_xMyM)
        deallocate(send_phi_yPzP,send_phi_yMzP,send_phi_yPzM,send_phi_yMzM,send_phi_xPzP,send_phi_xMzP,send_phi_xPzM,send_phi_xMzM,send_phi_xPyP,send_phi_xMyP,send_phi_xPyM,send_phi_xMyM)
        deallocate(recv_phi_yPzP,recv_phi_yMzP,recv_phi_yPzM,recv_phi_yMzM,recv_phi_xPzP,recv_phi_xMzP,recv_phi_xPzM,recv_phi_xMzM,recv_phi_xPyP,recv_phi_xMyP,recv_phi_xPyM,recv_phi_xMyM)
        deallocate(send_phi_xPyPzP,send_phi_xPyMzP,send_phi_xPyPzM,send_phi_xPyMzM,send_phi_xMyPzP,send_phi_xMyMzP,send_phi_xMyPzM,send_phi_xMyMzM)
        deallocate(recv_phi_xPyPzP,recv_phi_xPyMzP,recv_phi_xPyPzM,recv_phi_xPyMzM,recv_phi_xMyPzP,recv_phi_xMyMzP,recv_phi_xMyPzM,recv_phi_xMyMzM)

        deallocate(fl1, fl2, mass1,mass2,vol1,vol2, pre,tk)

    endif

    return
end subroutine MemAllocate_multi
