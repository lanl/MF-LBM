!==============================================================================
!                              PROGRAM: MF-LBM                                 
!~~~~~~~~~~~~~~~~~~~~~~~  April 1st, 2021.  C20125  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! This program was prepared at Los Alamos National Laboratory (LANL), Earth and 
! Environmental Sciences Division, Computational Earth Science Group (EES-16),  
! Subsurface Flow and Transport Team.
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! https://github.com/lanl/MF-LBM

! Authors: Yu Chen (yu_chen_007@outlook.com), 
! Qinjun Kang, Albert Valocchi (UIUC), Hari Viswanathan
!==============================================================================

#include "./preprocessor.h"
!===============================================================================================================================
!---------------------- Main Starts ----------------------
!===============================================================================================================================
program main_multiphase
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    USE,INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    include 'mpif.h'
    !indicator used to save extra backup checkpoint(pdf) data
    integer :: save_checkpoint_data_indicator, save_2rd_checkpoint_data_indicator   
    integer :: counter_checkpoint_save, counter_2rd_checkpoint_save
    real (kind = 8) :: t_all_sum, t_all, t_run_time, t_start1, t_end1, t_start2, t_end2
#ifdef _openacc
    integer, external :: setDevice
#endif

    !-set mpi environment-
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,id0,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,np,ierr)

    if(id0==0)then
        print*,"=============================================================================="
        print*,"                              PROGRAM: MF-LBM                                 "
        print*,"~~~~~~~~~~~~~~~~~~~~~~~  April 1st, 2021.  C20125  ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*,''
        print*,"This program was prepared at Los Alamos National Laboratory (LANL), Earth and "
        print*,"Environmental Sciences Division, Computational Earth Science Group (EES-16),  "
        print*,"Subsurface Flow and Transport Team."
        print*,''
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*,"https://github.com/lanl/MF-LBM"
        print*,''
        print*,"Authors: Yu Chen (yu_chen_007@outlook.com), "
        print*,"Qinjun Kang, Albert Valocchi (UIUC), Hari Viswanathan"
        print*,"=============================================================================="
        print*,''
    endif

    !################################################################################################################################
    !                                               Preparation
    !################################################################################################################################
    if(id0==0.and.(.not.ieee_support_nan()))then
        print*, 'ieee_support_nan is not True! Stop program!' 
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call mpi_abort(MPI_COMM_WORLD,ierr)
    endif

#ifdef _openacc
    !identify GPU device
    devnum = setDevice(np,id0)
    !print*,'MPI process id:', id0, 'GPU device number:', devnum
    call MPI_Barrier(MPI_COMM_WORLD,ierr )
#endif

    t_start2 = MPI_WTIME()

    simulation_end_indicator = 0
    save_checkpoint_data_indicator = 0          !default 0, saving data 1, after saving data 0
    save_2rd_checkpoint_data_indicator = 0       !default 0, saving data 1, after saving data 0

    ! initial value 1; after each checkpoint data saving, counter_checkpoint_save = counter_checkpoint_save + 1
    counter_checkpoint_save = 1.0  
    ! initial value 1; after each checkpoint data saving, counter_2rd_checkpoint_save = counter_2rd_checkpoint_save + 1
    counter_2rd_checkpoint_save = 1.0 
    
    relaxation = 1d0   !no relaxation if 1; relaxation feature temporary disabled, should keep 1

    if(id0==0)print*, ''
    if(id0==0)print*,'***************************** Initialization **********************************'
    call initialization_basic_multi
    if(trim(job_status)=='new_simulation')then
        call initialization_new_multi
    elseif(trim(job_status)=='continue_simulation')then
        call initialization_old_multi
    endif
    if(change_inlet_fluid_phase_cmd/=0)then
        call change_inlet_fluid_phase
        if(id==0)write(*,"(1X,'Inlet fluid phase was changed (1 to fluid1; 2 to fluid2): ', I1)")change_inlet_fluid_phase_cmd
    endif
    if(id0==0)print*, ''
    if(id0==0)print*,'************************** Initialization ends ********************************'
    if(id0==0)print*, ''

    !$acc data copy(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,u,v,w,rho,phi) &
    !$acc copyin(walls,w_in) copy(cn_x,cn_y,cn_z,c_norm,curv) create(fl1,fl2,sa1,pre,mass1,mass2,vol1,vol2,tk) &
    !$acc copyin(solid_boundary_nodes,fluid_boundary_nodes) &
    !$acc copy(f_convec_bc,g_convec_bc,phi_convec_bc) &
    !$acc copy(phi_old) &
    !$acc create(send_pdf_xP,send_pdf_xM,send_pdf_yP,send_pdf_yM,send_pdf_zP,send_pdf_zM, recv_pdf_xM,recv_pdf_xP,recv_pdf_yM,recv_pdf_yP,recv_pdf_zM,recv_pdf_zP, &
    !$acc send_pdf_yPzP,send_pdf_yMzP,send_pdf_yPzM,send_pdf_yMzM, send_pdf_xPzP,send_pdf_xMzP,send_pdf_xPzM,send_pdf_xMzM,send_pdf_xPyP,send_pdf_xMyP,send_pdf_xPyM,send_pdf_xMyM,&
    !$acc recv_pdf_yPzP,recv_pdf_yMzP,recv_pdf_yPzM,recv_pdf_yMzM, recv_pdf_xPzP,recv_pdf_xMzP,recv_pdf_xPzM,recv_pdf_xMzM,recv_pdf_xPyP,recv_pdf_xMyP,recv_pdf_xPyM,recv_pdf_xMyM,&
    !$acc send_phi_xM,send_phi_xP,send_phi_yM,send_phi_yP,send_phi_zM,send_phi_zP,recv_phi_xM,recv_phi_xP,recv_phi_yM,recv_phi_yP,recv_phi_zM,recv_phi_zP,&
    !$acc send_phi_yPzP,send_phi_yMzP,send_phi_yPzM,send_phi_yMzM, send_phi_xPzP,send_phi_xMzP,send_phi_xPzM,send_phi_xMzM,send_phi_xPyP,send_phi_xMyP,send_phi_xPyM,send_phi_xMyM,&
    !$acc recv_phi_yPzP,recv_phi_yMzP,recv_phi_yPzM,recv_phi_yMzM, recv_phi_xPzP,recv_phi_xMzP,recv_phi_xPzM,recv_phi_xMzM,recv_phi_xPyP,recv_phi_xMyP,recv_phi_xPyM,recv_phi_xMyM,&
    !$acc send_phi_xPyPzP,send_phi_xPyMzP,send_phi_xPyPzM,send_phi_xPyMzM,send_phi_xMyPzP,send_phi_xMyMzP,send_phi_xMyPzM,send_phi_xMyMzM)

    t_all_sum=0.0d0
    ntime=ntime0

    call color_gradient

    call cal_saturation

    if(trim(job_status)=='new_simulation')then
        saturation_old = saturation
        if(benchmark_cmd==0)then
            if(extreme_large_sim_cmd==0)then  ! initial distribution
                ! vtk_type 1: full flow field info for detailed analysis 
                ! vtk_type 2: phase field info, with single precision to save space
                ! vtk_type 3: force vectors from the CSF model
                call VTK_legacy_writer_3D(ntime, 2)
                call VTK_walls_bin
            else
                ! parallel I/O, distributed files, require post processing
                ! should be replaced by parallel version of VTK in the future
                call save_phi(ntime)   
            endif
        endif
    elseif(trim(job_status)=='continue_simulation')then
        saturation_old = -1d0    ! a negative value, in order to avoid convergence check at the first step after reading old data
    endif

    if(id==0)then
        write(*,"(1X,'Initial saturation: ', F6.4)") saturation_full_domain
    endif

    if(id==0)print*,''
    !################################################################################################################################
    !                                              Performance benchmarking
    !################################################################################################################################
    if(benchmark_cmd==1)then
        if(id==0)print*,'********************** Performance benchmarking *******************************'
        call benchmark
        if(id==0)print*,'******************** Performance benchmarking ends ****************************'

    !################################################################################################################################
    !                                                   Main loop Starts 
    !################################################################################################################################
    else
        if(id==0)print*,'************************** Entering main loop *********************************'

        do ntime=ntime0, ntime_max

            t_start1 = MPI_WTIME()

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            call main_iteration_kernel
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
            ! ----------- computation time ----------- 
            t_end1 = MPI_WTIME()
            t_all = t_end1 - t_start1 
            if(id==0)then
                IF(MOD(ntime,ntime_clock_sum)==0)THEN
                    open(unit=78,file='out1.output/time.dat',status='unknown',position='append')
                    write(78,"(I10,1x,F10.2, ' sec')")ntime,t_all_sum
                    close(78)
                    t_all_sum = 0d0
                ENDIF
                t_all_sum = t_all_sum + t_all
            endif
            ! ----------- monitors -----------
            if(MOD(ntime,ntime_monitor)==0)then
                if(steady_state_option==0.or.steady_state_option==3)then   
                    !displacement simulation: finished based on saturation reaches steady state or injected volume
                    call monitor
                    if(breakthrough_check==1)then
                        call monitor_breakthrough
                    endif
                elseif(steady_state_option==1)then !steady state simulation based on capillary pressure
                    call monitor_multiphase_steady_capillarypressure
                elseif(steady_state_option==2)then !steady state simulation based on phase field
                    call monitor_multiphase_steady_phasefield
                endif
            endif
            ! ----------- simulation progress - time steps -----------
            if(mod(ntime,ntime_display_steps)==0)then
                    if(id==0)print*,'ntime = ',ntime
            endif
            ! ----------- phase field (vtk) files for visualization and other analysis -----------
            if(MOD(ntime,ntime_animation)==0)then
                if(extreme_large_sim_cmd==0)then
                    call VTK_legacy_writer_3D(ntime, 2)
                    !call VTK_small_bin_2d_micromodel(ntime)
                else
                    ! parallel I/O, distributed files, require post processing
                    ! should be replaced by parallel version of VTK in the future
                    call save_phi(ntime)   
                endif
            endif
            ! ----------- full flow field (VTK) files for further analysis -----------
            if(MOD(ntime,ntime_visual)==0)then
                if(extreme_large_sim_cmd==0)then
                    call VTK_legacy_writer_3D(ntime, 1)
                else
                    call save_macro(ntime)    ! parallel I/O, distributed files, require post processing
                endif
            endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
                
            !~~~~~~~~~~~~~~~~~~ SAVE checkpoint PDF DATA for restarting simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF(MOD(ntime,ntime_clock_sum)==0)THEN   ! frequency to check checkpoint data saving
                if(id==0)then
                    t_end2 = MPI_WTIME()                   
                    t_run_time = (t_end2 - t_start2 )*2.77777777777d-4   !  second to hour  1/3600           
                    write(*, "(A20,1X,F6.3,A8)")'simulation has run ', t_run_time,' hours'

                    if(t_run_time>=simulation_duration_timer)then                   
                        print*,'Time to save checkpoint data and exit the program!'
                        simulation_end_indicator = 2 
                    endif

                    if(t_run_time>= counter_checkpoint_save * checkpoint_save_timer)then                   
                        print*,'Time to save checkpoint data!'
                        save_checkpoint_data_indicator = 1 
                        counter_checkpoint_save = counter_checkpoint_save + 1
                    endif

                    if(t_run_time>=counter_2rd_checkpoint_save * checkpoint_2rd_save_timer)then                   
                        print*,'Time to save secondary checkpoint data!'
                        save_2rd_checkpoint_data_indicator = 1 
                        counter_2rd_checkpoint_save = counter_2rd_checkpoint_save + 1
                    endif
                endif 

                call MPI_Bcast(simulation_end_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr) 
                call MPI_Bcast(save_checkpoint_data_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)    
                call MPI_Bcast(save_2rd_checkpoint_data_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)    

                if(save_checkpoint_data_indicator==1)then
                    call save_checkpoint(0)    ! save pdf to the default location when option is 0
                    save_checkpoint_data_indicator = 0   !reset status
                endif 
                if(save_2rd_checkpoint_data_indicator==1)then
                    call save_checkpoint(1)    ! save pdf to the secondary backup location when option is 1
                    save_2rd_checkpoint_data_indicator = 0   !reset status
                endif
            ENDIF
            !~~~~~~~~~~~~~~~~~~ SAVE checkpoint PDF DATA for restarting simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if(simulation_end_indicator>0)exit

        ENDDO
        if(id==0)print*,'************************ Exiting main iteration ******************************'
    !###############################################################################################################################################
    !                                                   Main loop ends 
    !###############################################################################################################################################
        
        if(simulation_end_indicator==0)then
            ntime = ntime - 1  ! dial back ntime
            call save_checkpoint(0)
            if(extreme_large_sim_cmd==0)then
                call VTK_legacy_writer_3D(ntime, 2)
            else
                call save_phi(ntime)   ! parallel I/O, distributed files, require post processing
            endif         
            if(id==0)then
                if(simulation_end_indicator==0)print*, 'Simulation ended after ', ntime, 'iterations which reached the maximum time step!'
                OPEN(UNIT=9,FILE='./job_status.txt',form='formatted',status='replace')
                write(9,'(a)')"simulation_reached_max_step"
                close(9)
            endif

        elseif(simulation_end_indicator==1)then
            call save_checkpoint(0)
            if(extreme_large_sim_cmd==0)then
                call VTK_legacy_writer_3D(ntime, 2)
            else
                call save_phi(ntime)   ! parallel I/O, distributed files, require post processing
            endif           
            if(id==0)then
                print*, 'Simulation ended successfully after ', ntime, 'iterations!'               
                OPEN(UNIT=9,FILE='./job_status.txt',form='formatted',status='replace')
                write(9,'(a)')'simulation_done'
                close(9)
            endif

        elseif(simulation_end_indicator==2)then
            call save_checkpoint(0)
            if(id==0)then
                print*, 'Simulation data saved but not finished after ', ntime, 'iterations!'
                OPEN(UNIT=9,FILE='./job_status.txt',form='formatted',status='replace')
                write(9,'(a)')'continue_simulation'
                close(9)
            endif

        elseif(simulation_end_indicator==3)then
            if(extreme_large_sim_cmd==0)then
                call VTK_legacy_writer_3D(ntime, 1)
            else
                call save_phi(ntime)   ! parallel I/O, distributed files, require post processing
            endif 
            if(id==0)then       
                print*, 'Simulation failed after ', ntime, 'iterations!'  
                OPEN(UNIT=9,FILE='./job_status.txt',form='formatted',status='replace')
                write(9,'(a)')'simulation_failed'
                close(9)
            endif 
        endif
        call MPI_Barrier( MPI_COMM_vgrid,ierr)

    endif  ! performance benchmarking conditional branch

    !$acc end data
    
    call MemAllocate_geometry(2)
    call MemAllocate_multi(2)
    call mpi_finalize(IERR)
End program main_multiphase
!================================================================================================================================================
!---------------------- Main Ends ----------------------
!================================================================================================================================================




!================================================================================================================================================
!---------------------- iteration kernel starts ----------------------
!================================================================================================================================================
subroutine main_iteration_kernel
    use mpi_variable
    use Misc_module
    use Fluid_multiphase
    IMPLICIT NONE
    include 'mpif.h'

    ! MPI enabled
    if(mpi_x.or.mpi_y.or.mpi_z)then
        call mpi_irecv_initialization 
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ kernel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(mod(ntime,2)==0)then
        !************************** even step *****************************************
        !~~~~~~~~~~~~~~~ overlapped communication and computation ~~~~~~~~~~~~~~~~~

        if(mpi_z)then
            call kernel_even_color(1,nx, 1,ny, 1,            iz_async, LBM_async_z)
            call kernel_even_color(1,nx, 1,ny, nz-iz_async+1,nz      , LBM_async_z)
        endif

        if(mpi_y)then 
            call kernel_even_color(1,nx, 1,             iy_async, iz_async+1,nz-iz_async, LBM_async_y)
            call kernel_even_color(1,nx, ny-iy_async+1, ny,       iz_async+1,nz-iz_async, LBM_async_y)
        endif  
        ! if(mpi_x)then 
        !     call kernel_even_color(1,            ix_async, iy_async+1,ny-iy_async, iz_async+1,nz-iz_async, LBM_async_x)
        !     call kernel_even_color(nx-ix_async+1,nx,       iy_async+1,ny-iy_async, iz_async+1,nz-iz_async, LBM_async_x)
        ! endif    


        ! MPI enabled
        if(mpi_x.or.mpi_y.or.mpi_z)then  
            call mpi_pdf_pull_async_pack
            call mpi_phi_async_pack
            !$acc wait(LBM_async_z,LBM_async_y,LBM_async_x)
#ifdef _openacc
            call kernel_even_color(ix_async+1,nx-ix_async,iy_async+1,ny-iy_async,iz_async+1,nz-iz_async, LBM_async)
            call mpi_send_req
#else
            call mpi_send_req
            call kernel_even_color(ix_async+1,nx-ix_async,iy_async+1,ny-iy_async,iz_async+1,nz-iz_async, LBM_async)
#endif
            call mpi_pdf_pull_sync_update
            call mpi_phi_sync_update
            !$acc wait

        ! MPI disabled
        else

            call kernel_even_color(1,nx,1,ny,1,nz, LBM_async)
            !$acc wait
        endif

        !~~~~~~~~~~~~~~~ overlapped communication and computation ~~~~~~~~~~~~~~~~~

        !%%%%%%%%%%%%%%%% boundary conditions %%%%%%%%%%%%%%%%%%%
        if(kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then    !non-periodic BC along flow direction (z)
            if(inlet_BC==1)then
                call inlet_bounce_back_velocity_BC_before_odd   !velocity inlet bc
            elseif(inlet_BC==2)then
                call inlet_Zou_He_pressure_BC_before_odd   !pressure inlet bc
            endif
            if(outlet_BC==1)then
                call outlet_convective_BC_before_odd   !convective outlet bc
            elseif(outlet_BC==2)then
                call outlet_Zou_He_pressure_BC_before_odd   !pressure outlet bc
            endif
        endif
        if(porous_plate_cmd/=0)then      !place a porous plate       
            call porous_plate_BC_before_odd
        endif
        !%%%%%%%%%%%%%%%% boundary conditions %%%%%%%%%%%%%%%%%%%

        !************************** even step *****************************************

    else

        !************************** odd step *****************************************

        !~~~~~~~~~~~~~~~ overlapped communication and computation ~~~~~~~~~~~~~~~~~
        if(mpi_z)then
            call kernel_odd_color(1,nx, 1,ny, 1,            iz_async, LBM_async_z)
            call kernel_odd_color(1,nx, 1,ny, nz-iz_async+1,nz      , LBM_async_z)
        endif
        if(mpi_y)then
            call kernel_odd_color(1,nx, 1,             iy_async, iz_async+1,nz-iz_async, LBM_async_y)
            call kernel_odd_color(1,nx, ny-iy_async+1, ny,       iz_async+1,nz-iz_async, LBM_async_y)
        endif  
        !if(mpi_x)then 
        !   call kernel_odd_color(1,            ix_async, iy_async+1,ny-iy_async, iz_async+1,nz-iz_async, LBM_async_x)
        !   call kernel_odd_color(nx-ix_async+1,nx,       iy_async+1,ny-iy_async, iz_async+1,nz-iz_async, LBM_async_x)
        !endif 

    ! MPI enabled
        if(mpi_x.or.mpi_y.or.mpi_z)then  
            call mpi_pdf_push_async_pack
            call mpi_phi_async_pack
            !$acc wait(LBM_async_z,LBM_async_y,LBM_async_x)
#ifdef _openacc
            call kernel_odd_color(ix_async+1,nx-ix_async,iy_async+1,ny-iy_async,iz_async+1,nz-iz_async, LBM_async)
            call mpi_send_req
#else
            call mpi_send_req
            call kernel_odd_color(ix_async+1,nx-ix_async,iy_async+1,ny-iy_async,iz_async+1,nz-iz_async, LBM_async)
#endif

            call mpi_pdf_push_sync_update
            call mpi_phi_sync_update       
            !$acc wait

        ! MPI disabled
        else

            call kernel_odd_color(1,nx,1,ny,1,nz, LBM_async)
            !$acc wait
        endif
        !~~~~~~~~~~~~~~~ overlapped communication and computation ~~~~~~~~~~~~~~~~~

        !%%%%%%%%%%%%%%%% boundary conditions %%%%%%%%%%%%%%%%%%%
        if(kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then    !non-periodic BC along flow direction (z)
            if(inlet_BC==1)then
                call inlet_bounce_back_velocity_BC_after_odd   !velocity inlet bc
            elseif(inlet_BC==2)then
                call inlet_Zou_He_pressure_BC_after_odd   !pressure inlet bc
            endif
            if(outlet_BC==1)then
                call outlet_convective_BC_after_odd   !convective outlet bc
            elseif(outlet_BC==2)then
                call outlet_Zou_He_pressure_BC_after_odd   !pressure outlet bc
            endif
        endif
        if(porous_plate_cmd/=0)then      !place a porous plate           
            call porous_plate_BC_after_odd
        endif
        !%%%%%%%%%%%%%%%% boundary conditions %%%%%%%%%%%%%%%%%%%
        
        !************************** odd step *****************************************
    endif
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ kernel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call color_gradient       !calculate color gradient

    return
end subroutine main_iteration_kernel
!================================================================================================================================================
!---------------------- iteration kernel ends ----------------------
!================================================================================================================================================





!================================================================================================================================================
!---------------------- performance benchmarking starts ----------------------
!================================================================================================================================================
subroutine benchmark
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
#ifdef _openacc
    use cudafor
#endif
IMPLICIT NONE
  integer :: clock2, clock1, clockmax, clockrate, ticks, clock01, clock02, round, maxRound
  real(kind=8) :: t_all

  maxRound = 3

  call system_clock(count_max=clockmax, count_rate=clockrate)

  if(id==0)print*,'------------ Warm-up section starts --------------------'
  do ntime=1, 20
      call main_iteration_kernel
  ENDDO
  if(id==0)print*,'------------ Warm-up section ends ----------------------'

  if(id==0)print*,'----------- Main benchmarking starts -------------------'
  if (id == 0) open (unit=78, file='out1.output/benchmark_time.dat', status='unknown', position='append')

  do round = 1, maxRound
    call system_clock(clock01)
#ifdef gpu_profiling
    call cudaProfilerStart()
#endif
    do ntime = 1, ntime_max_benchmark
        call main_iteration_kernel
    END DO
#ifdef gpu_profiling
    call cudaProfilerStop()
#endif
    call MPI_Barrier(MPI_COMM_vgrid, ierr)
    call system_clock(clock02)
    ntime = ntime - 1
    ticks = clock02 - clock01
    t_all = float(ticks)/float(clockrate)
    if (id == 0) then
      write (78, "('Code performance: ', F12.4, ' MLUPS')") pore_sum*dble(ntime)/(t_all*1000000d0)
      write (*, "(1X,'Code performance: ', F12.4, ' MLUPS')") pore_sum*dble(ntime)/(t_all*1000000d0)
    end if
  enddo
  if(id==0)print*,'----------- Main benchmarking ends ---------------------'

  call monitor

  if(id==0)then
      write(78,*)'time step  ','wallclock time  ','saturation ', 'capillary number'
      write(78,"(I8,3(1x,E14.6))")ntime,t_all,saturation,ca     ! sa and ca used to check if simulation doesnot work properly
      write(*,"(1X,'Benchmarking ended successfully after ', I6, ' iterations')")ntime*maxRound
  endif
    !call VTK_legacy_writer_3D(ntime, 2)
    
    return
end subroutine benchmark
!====================================================================================================================================================
!---------------------- performance benchmarking ends ----------------------
!====================================================================================================================================================

