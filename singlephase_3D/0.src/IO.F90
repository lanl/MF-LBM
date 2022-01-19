!============================================================================================================================
!---------------------- Read input parameters ----------------------
!============================================================================================================================
subroutine read_parameter
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    LOGICAL :: ALIVE
    integer :: N(100)
    real(kind=8)  :: A(100)
    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos, i
    integer :: ios = 0
    integer :: error_signal

    if (id0 == 0) print *, 'Checking simulation status ...'
    INQUIRE (FILE='./job_status.txt', EXIST=ALIVE)
    if (alive) then
        OPEN (UNIT=10, FILE='./job_status.txt', form='formatted', action='READ')
        read (10, *) job_status
        close (10)
        if (trim(job_status) == 'new_simulation') then
            if (id0 == 0) print *, 'Simulation status is: New Simulation!'

        elseif (trim(job_status) == 'continue_simulation') then
            if (id0 == 0) print *, 'Simulation status is: Continue existing simulation!'

        elseif (trim(job_status) == 'simulation_done') then
            if (id0 == 0) print *, 'Previous simulation is already finished! Exiting program!'
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call mpi_abort(MPI_COMM_WORLD, ierr)

        elseif (trim(job_status) == 'simulation_reached_max_step') then
            if (id0 == 0) print *, 'Previous simulation reached maximum step! Exiting program!'
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call mpi_abort(MPI_COMM_WORLD, ierr)

        elseif (trim(job_status) == 'simulation_failed') then
            if (id0 == 0) print *, 'Previous simulation failed! Exiting program!'
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call mpi_abort(MPI_COMM_WORLD, ierr)

        else
            if (id0 == 0) print *, 'Wrong simlation status! Exiting program!'
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call mpi_abort(MPI_COMM_WORLD, ierr)
        end if
    else
        if (id0 == 0) print *, 'Missing job status file! Exiting program!'
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call mpi_abort(MPI_COMM_WORLD, ierr)
    end if

    INQUIRE (FILE='./simulation_control.txt', EXIST=ALIVE)
    if (.not. alive) then
        if (id0 == 0) print *, 'Missing simulation control file! Exiting program!'
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call mpi_abort(MPI_COMM_WORLD, ierr)
    end if

    if (id0 .eq. 0) then   !here id0 is used before setMPI
        open (90, file='./simulation_control.txt', form='formatted', action='READ')
        print *, ''
        print *, '****************** Reading in parameters from control file *******************'
        print *, '.........'
        ! ios is negative if an end of record condition is encountered or if an endfile condition was detected.
        ! It is positive if an error was detected.  ios is zero otherwise.
        do while (ios == 0)
            read (90, '(A)', iostat=ios) buffer
            if (ios == 0) then
                ! Find the first instance of whitespace.  Split label and data.
                pos = scan(buffer, '    ')
                label = buffer(1:pos)
                buffer = buffer(pos + 1:)

                select case (label)
                case ('steady_state_option')
                    read (buffer, *, iostat=ios) steady_state_option
                    write (*, "(1X,'Steady_state_option: ', I2)") steady_state_option
                    print *, '---------------------------'
                case ('convergence_criteria')
                    read (buffer, *, iostat=ios) convergence_criteria
                    write (*, "(1X,'Convergence_criteria: ', ES13.6)") convergence_criteria
                    print *, '---------------------------'

                case ('MRT_collision_parameter_preset')
                    read (buffer, *, iostat=ios) mrt_para_preset
                    write (*, "(1X,'MRT_collision_parameter_preset: ', I2)") mrt_para_preset
                    print *, '---------------------------'

                case ('benchmark_cmd')
                    read (buffer, *, iostat=ios) benchmark_cmd
                    write (*, "(1X,'Benchmark_cmd: ', I2)") benchmark_cmd
                    print *, '---------------------------'

                case ('output_fieldData_precision_cmd')
                    read (buffer, *, iostat=ios) output_fieldData_precision_cmd
                    write (*, "(1X,'output_fieldData_precision_cmd: ', I2)") output_fieldData_precision_cmd
                    print *, '---------------------------'

                case ('extreme_large_sim_cmd')
                    read (buffer, *, iostat=ios) extreme_large_sim_cmd
                    write (*, "(1X,'Extreme_large_sim_cmd: ', I2)") extreme_large_sim_cmd
                    print *, '---------------------------'

                case ('modify_geometry_cmd')
                    read (buffer, *, iostat=ios) modify_geometry_cmd
                    write (*, "(1X,'Modify_geometry_cmd: ', I2)") modify_geometry_cmd
                    print *, '---------------------------'

                case ('external_geometry_read_cmd')
                    read (buffer, *, iostat=ios) external_geometry_read_cmd
                    write (*, "(1X,'External_geometry_read_cmd: ', I2)") external_geometry_read_cmd
                    print *, '---------------------------'

                case ('lattice_dimensions') ! total nodes for the whole domain
                    read (buffer, *, iostat=ios) nxGlobal, nyGlobal, nzGlobal
                    write (*, "(1X,'nxGlobal = ', I5)") nxGlobal
                    write (*, "(1X,'nYGlobal = ', I5)") nYGlobal
                    write (*, "(1X,'nZGlobal = ', I5)") nZGlobal
                case ('excluded_layers') !excluded_layers
                    read (buffer, *, iostat=ios) n_exclude_inlet, n_exclude_outlet
                    write (*, "(1X,'Inlet_excluded_layers = ', I3)") n_exclude_inlet
                    write (*, "(1X,'Outlet_excluded_layers = ', I3)") n_exclude_outlet
                case ('char_length') ! characteristic length of the flow domain (used in the definition of Re number)
                    read (buffer, *, iostat=ios) char_length
                    write (*, "(1X,'char_length = ', F10.3)") char_length

                case ('domain_wall_status_x') !domain_wall_status
                    read (buffer, *, iostat=ios) domain_wall_status_x_min, domain_wall_status_x_max
                    write (*, "(1X,'Apply nonslip boundary condition at x = 1       : ', I2)") domain_wall_status_x_min
                    write (*, "(1X,'Apply nonslip boundary condition at x = nxGlobal: ', I2)") domain_wall_status_x_max
                case ('domain_wall_status_y') !domain_wall_status
                    read (buffer, *, iostat=ios) domain_wall_status_y_min, domain_wall_status_y_max
                    write (*, "(1X,'Apply nonslip boundary condition at y = 1       : ', I2)") domain_wall_status_y_min
                    write (*, "(1X,'Apply nonslip boundary condition at y = nyGlobal: ', I2)") domain_wall_status_y_max
                case ('domain_wall_status_z') !domain_wall_status
                    read (buffer, *, iostat=ios) domain_wall_status_z_min, domain_wall_status_z_max
                    write (*, "(1X,'Apply nonslip boundary condition at z = 1       : ', I2)") domain_wall_status_z_min
                    write (*, "(1X,'Apply nonslip boundary condition at z = nzGlobal: ', I2)") domain_wall_status_z_max
                case ('periodic_indicator') ! periodic BC indicator, 1 periodic, 0 non-periodic
                    read (buffer, *, iostat=ios) iper, jper, kper
                    write (*, "(1X,'X_periodic_option = ', I2)") iper
                    write (*, "(1X,'Y_periodic_option = ', I2)") jper
                    write (*, "(1X,'Z_periodic_option = ', I2)") kper
                case ('MPI_async_layers_num') !layers of nodes used in overlaping computation and communication
                    read (buffer, *, iostat=ios) ix_async, iy_async, iz_async
                    write (*, "(1X,'MPI_async_layers_num_X = ', I2)") ix_async
                    write (*, "(1X,'MPI_async_layers_num_Y = ', I2)") iy_async
                    write (*, "(1X,'MPI_async_layers_num_Z = ', I2)") iz_async
                case ('MPI_process_num') ! mpi process numbers along each axis
                    read (buffer, *, iostat=ios) npx, npy, npz
                    write (*, "(1X,'MPI_process_num_X = ', I4)") npx
                    write (*, "(1X,'MPI_process_num_Y = ', I4)") npy
                    write (*, "(1X,'MPI_process_num_Z = ', I4)") npz
                    print *, '---------------------------'

                case ('fluid_viscosity') ! fluid1_viscosity
                    read (buffer, *, iostat=ios) la_nu
                    write (*, "(1X,'Fluid_viscosity = ', F8.6)") la_nu

                case ('inlet_BC') !inlet_BC
                    read (buffer, *, iostat=ios) inlet_BC
                    write (*, "(1X,'Inlet_BC = ', I2)") inlet_BC
                case ('outlet_BC') !outlet_BC
                    read (buffer, *, iostat=ios) outlet_BC
                    write (*, "(1X,'Outlet_BC = ', I2)") outlet_BC

                case ('body_force_0') ! initial value of body force Z or pressure gradient
                    read (buffer, *, iostat=ios) force_Z0
                    write (*, "(1X,'Body_force_Z0 = ', ES13.6)") force_Z0
                case ('body_force_increment') ! body force increment
                    read (buffer, *, iostat=ios) D_force_Z
                    write (*, "(1X,'Body_force_Z0 = ', ES13.6)") D_force_Z

                case ('Reynolds_number')
                    read (buffer, *, iostat=ios) Re
                    write (*, "(1X,'Reynolds_number = ', F10.3)") Re

                case ('rho_drop')
                    read (buffer, *, iostat=ios) rho_drop
                    write (*, "(1X,'inlet/outlet density drop = ', F6.4)") rho_drop

                case ('target_inject_pore_volume')
                    read (buffer, *, iostat=ios) target_inject_pore_volume
                    write (*, "(1X,'Target_inject_pore_volume = ', F5.2)") target_inject_pore_volume

                case ('max_time_step') ! timer: max iterations
                    read (buffer, *, iostat=ios) ntime_max
                    write (*, "(1X,'Max_iterations = ', I10)") ntime_max
                case ('max_time_step_benchmark') ! timer: max iterations for the benchmark case
                    read (buffer, *, iostat=ios) ntime_max_benchmark
                    write (*, "(1X,'Max_iterations_benchmark = ', I8)") ntime_max_benchmark
                case ('ntime_visual') ! timer: when to output detailed visualization data
                    read (buffer, *, iostat=ios) ntime_visual
                    write (*, "(1X,'Full_VTK_timer = ', I8)") ntime_visual

                case ('monitor_timer') ! timer: when to output monitor data, profile along flow direction
                    read (buffer, *, iostat=ios) ntime_monitor
                    write (*, "(1X,'Monitor_timer = ', I8)") ntime_monitor
                case ('monitor_profile_timer_ratio')
                    read (buffer, *, iostat=ios) ntime_monitor_profile_ratio
                    write (*, "(1X,'Monitor_profile_timer_ratio = ', I4)") ntime_monitor_profile_ratio
                case ('computation_time_timer') ! timer: gaps used to record computation time
                    read (buffer, *, iostat=ios) ntime_clock_sum
                    write (*, "(1X,'Computation_time_timer = ', I8)") ntime_clock_sum
                case ('display_steps_timer') ! timer: when to display time steps
                    read (buffer, *, iostat=ios) ntime_display_steps
                    write (*, "(1X,'Display_time_steps_timer = ', I8)") ntime_display_steps
                    print *, '---------------------------'

                case ('checkpoint_save_timer') ! save PDF data for restart simulation based on wall clock timer  (unit hours)
                    read (buffer, *, iostat=ios) checkpoint_save_timer
                    write (*, "(1X,'Checkpoint_save_timer (wall clock time, hours) = ', F6.2)") checkpoint_save_timer
                case ('checkpoint_2rd_save_timer') ! save secondary PDF data for restart simulation based on wall clock timer  (unit hours)
                    read (buffer, *, iostat=ios) checkpoint_2rd_save_timer
                    write (*, "(1X,'Checkpoint_2rd_save_timer (wall clock time, hours) = ', F6.2)") checkpoint_2rd_save_timer
                case ('simulation_duration_timer') ! simulation duration in hours, exit and save simulation afterwards
                    read (buffer, *, iostat=ios) simulation_duration_timer
                    write (*, "(1X,'Simulation_duration (wall clock time, hours) = ', F6.2)") simulation_duration_timer
                    print *, '---------------------------'

                case ('d_vol_detail') ! when to output detailed visulization data - based on injected volume
                    read (buffer, *, iostat=ios) d_vol_detail
                    write (*, "(1X,'Full_VTK_by_injected_vol = ', F6.2)") d_vol_detail
                case ('d_vol_monitor') ! when to output bulk property data - based on injected volume
                    read (buffer, *, iostat=ios) d_vol_monitor
                    write (*, "(1X,'Monitor_by_injected_vol = ', F6.2)") d_vol_monitor
                    print *, '---------------------------'

                end select
            end if
        end do
        close (90)

        ntime_monitor_profile = ntime_monitor_profile_ratio*ntime_monitor
        if (id0 == 0) write (*, "(1X,'monitor_profile_timer = ', I8)") ntime_monitor_profile

        print *, '************** End reading in parameters from control file *******************'
        print *, ''

        N(1) = output_fieldData_precision_cmd
        N(2) = mrt_para_preset
        N(3) = steady_state_option
        N(4) = benchmark_cmd

        N(5) = nxGlobal
        N(6) = nyGlobal
        N(7) = nzGlobal
        N(10) = n_exclude_inlet
        N(11) = n_exclude_outlet

        N(12) = iper
        N(13) = jper
        N(14) = kper
        N(15) = inlet_BC
        N(16) = outlet_BC
        N(17) = domain_wall_status_x_min
        N(18) = domain_wall_status_x_max
        N(19) = domain_wall_status_y_min
        N(20) = domain_wall_status_y_max
        N(21) = domain_wall_status_z_min
        N(22) = domain_wall_status_z_max

        N(23) = modify_geometry_cmd
        N(24) = external_geometry_read_cmd
        N(25) = geometry_preprocess_cmd

        N(26) = ntime0
        N(27) = ntime_animation
        N(28) = ntime_macro
        N(29) = ntime_visual
        N(30) = ntime_pdf
        N(31) = ntime_max
        N(32) = ntime_max_benchmark
        N(33) = ntime_clock_sum
        N(34) = ntime_monitor
        N(35) = ntime_monitor_profile
        N(36) = ntime_relaxation
        N(37) = ntime_display_steps

        N(38) = npx
        N(39) = npy
        N(40) = npz
        N(41) = ix_async
        N(42) = iy_async
        N(43) = iz_async

        N(44) = double_bak_checkpoint_pdf_cmd

        N(47) = extreme_large_sim_cmd

        A(1) = char_length
        A(2) = force_z0
        A(3) = D_force_Z
        A(4) = Re
        A(5) = rho_drop

        A(10) = la_nu
        A(12) = d_vol_monitor

        A(14) = d_vol_detail

        A(15) = checkpoint_save_timer
        A(16) = simulation_duration_timer

        A(17) = convergence_criteria

        A(18) = target_inject_pore_volume

        A(19) = checkpoint_2rd_save_timer

    end if

    call mpi_bcast(N, 100, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(A, 100, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    output_fieldData_precision_cmd = N(1)
    mrt_para_preset = N(2)
    steady_state_option = N(3)
    benchmark_cmd = N(4)

    nxGlobal = N(5)
    nyGlobal = N(6)
    nzGlobal = N(7)
    n_exclude_inlet = N(10)
    n_exclude_outlet = N(11)

    iper = N(12)
    jper = N(13)
    kper = N(14)
    inlet_BC = N(15)
    outlet_BC = N(16)
    domain_wall_status_x_min = N(17)
    domain_wall_status_x_max = N(18)
    domain_wall_status_y_min = N(19)
    domain_wall_status_y_max = N(20)
    domain_wall_status_z_min = N(21)
    domain_wall_status_z_max = N(22)

    modify_geometry_cmd = N(23)
    external_geometry_read_cmd = N(24)
    geometry_preprocess_cmd = N(25)

    ntime0 = N(26)

    ntime_macro = N(28)
    ntime_visual = N(29)
    ntime_pdf = N(30)
    ntime_max = N(31)
    ntime_max_benchmark = N(32)
    ntime_clock_sum = N(33)
    ntime_monitor = N(34)
    ntime_monitor_profile = N(35)

    ntime_relaxation = N(36)
    ntime_display_steps = N(37)

    npx = N(38)
    npy = N(39)
    npz = N(40)
    ix_async = N(41)
    iy_async = N(42)
    iz_async = N(43)

    double_bak_checkpoint_pdf_cmd = N(44)
    extreme_large_sim_cmd = N(47)

    char_length = A(1)
    force_z0 = A(2)
    D_force_Z = A(3)
    Re = A(4)
    rho_drop = A(5)
    la_nu = A(10)

    d_vol_monitor = A(12)

    d_vol_detail = A(14)

    checkpoint_save_timer = A(15)
    simulation_duration_timer = A(16)

    convergence_criteria = A(17)

    target_inject_pore_volume = A(18)

    checkpoint_2rd_save_timer = A(19)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ check correctness of input parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (id0 == 0) then
        print *, '************ Start checking correctness of input parameters ******************'
    end if
    error_signal = 0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! MPI parameters
    if (iper == 1 .or. domain_wall_status_x_max == 0 .or. domain_wall_status_x_min == 0) then
if (id0 == 0) print *, 'Error: X direction periodic BC enabled or non-slip BC not applied at x = xmin or x = xmax! Exiting program!'
        error_signal = 1
    end if
    if (jper == 0 .and. (domain_wall_status_y_max == 0 .or. domain_wall_status_y_min == 0)) then
        if(id0==0)print*,'Error: non-slip BC not applied at y = ymin or y = ymax while y direction periodic BC not enabled! Exiting program!'
        error_signal = 1
    end if
    if (jper == 1 .and. (domain_wall_status_y_max == 1 .or. domain_wall_status_y_min == 1)) then
 if (id0 == 0) print *, 'Error: non-slip BC applied at y = ymin or y = ymax while y direction periodic BC enabled! Exiting program!'
        error_signal = 1
    end if
    if (kper == 1 .and. (domain_wall_status_z_max == 1 .or. domain_wall_status_z_min == 1)) then
 if (id0 == 0) print *, 'Error: non-slip BC applied at z = zmin or z = zmax while z direction periodic BC enabled! Exiting program!'
        error_signal = 1
    end if

    if (npx*npy*npz .ne. np) then
        if (id0 == 0) print *, 'MPI error: npx*npy*npz is not equal to np! Exiting program!'
        error_signal = 1
    end if

    ! this code has temporary disabled all X direction MPI communication
    if (npx /= 1) then
        if (id0 == 0) print *, 'MPI error: MPI_process_num_X is not equal to 1! Exiting program!'
        error_signal = 1
    end if

    ! this code has temporary disabled all X direction MPI communication
    mpi_x = .false.
    if (ix_async /= 0) then
        ix_async = 0
    end if
    ! npy > 1 or y direction periodic BC
    if (npy > 1 .or. jper == 1) then
        mpi_y = .true.
        if (iy_async == 0) then
          if (id0 == 0) print *, 'MPI error: iy_async is zero when MPI communication along y direction is enabled! Exiting program!'
            error_signal = 1
        else
            if (id0 == 0) print *, 'MPI communication along y direction is enabled.'
        end if
    else
        mpi_y = .false.
        if (iy_async /= 0) then
            iy_async = 0
        end if
    end if
    ! npz > 1 or z direction periodic BC
    if (npz > 1 .or. kper == 1) then
        mpi_z = .true.
        if (iz_async == 0) then
          if (id0 == 0) print *, 'MPI error: iz_async is zero when MPI communication along z direction is enabled! Exiting program!'
            error_signal = 1
        else
            if (id0 == 0) print *, 'MPI communication along z direction is enabled.'
        end if
    else
        mpi_z = .false.
        if (iz_async /= 0) then
            iz_async = 0
        end if
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Boundary conditions
    if (outlet_BC == 1 .and. inlet_BC == 2) then
  if(id0==0)print*,'Inlet/outlet boundary condition error: Inlet pressure + outlet convective BC is not supported! Exiting program!'
        error_signal = 1
    end if

    if (error_signal == 1) then
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call mpi_abort(MPI_COMM_WORLD, ierr)
    else
        if (id0 == 0) then
            print *, 'Everything looks good!'
            print *, '************** End checking correctness of input parameters ******************'
            print *, ''
        end if
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ check correctness of input parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    return
end subroutine read_parameter

!=======================================================================================================================================================================================================================
!---------------------- save data ----------------------
!=======================================================================================================================================================================================================================
!******************************* save checkpoint data *************************************
subroutine save_checkpoint(save_option)   ! option 0 - default location; option 1 - secondary location
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    integer :: save_option
    character(len=30) :: flnm   !file name
    character(len=255) :: cwd

    !$acc update host(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18)
    if (outlet_bc == 1) then    !convective BC
        !$acc update host(f_convec_bc)
    end if

    write (flnm, "('id',i4.4)") id
    if (save_option == 0) then
        open (unit=9 + id, file='out2.checkpoint/'//trim(flnm), FORM='unformatted', status='replace', access='stream')
    else
        open (unit=9 + id, file='out2.checkpoint/2rd_backup/'//trim(flnm), FORM='unformatted', status='replace', access='stream')
    end if

    rewind (9 + id)
    !time
    write (9 + id) ntime + 1, force_z, rho_in
    !fluid PDF
    write (9 + id) f0
    write (9 + id) f1
    write (9 + id) f2
    write (9 + id) f3
    write (9 + id) f4
    write (9 + id) f5
    write (9 + id) f6
    write (9 + id) f7
    write (9 + id) f8
    write (9 + id) f9
    write (9 + id) f10
    write (9 + id) f11
    write (9 + id) f12
    write (9 + id) f13
    write (9 + id) f14
    write (9 + id) f15
    write (9 + id) f16
    write (9 + id) f17
    write (9 + id) f18
    if (outlet_bc == 1) then    !convective BC
        write (9 + id) f_convec_bc
    end if
    close (9 + id)
    if (id == 0 .and. save_option == 0) then
        print *, 'Saving checkpoint data completed!'
        OPEN (UNIT=9, FILE='./job_status.txt', form='formatted', status='replace')
        write (9, '(a)') 'continue_simulation'
        close (9)
    end if
    if (id == 0 .and. save_option == 1) print *, 'Saving secondary checkpoint data completed!'

    return
end subroutine save_checkpoint

!******************************* save data - macro variables *************************************
subroutine save_macro(nt)
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    integer :: nt
    integer :: i, j, k
    character(len=30) :: flnm  !file name

    call compute_macro_vars

    !$acc update host(u,v,w,rho)

    write (flnm, "('full_nt',i9.9,'_id',i5.5)") nt, id

    if (id == 0) print *, 'Start to save full flow field data.'

    open (unit=9 + id, file='out3.field_data/'//trim(flnm), FORM='unformatted', status='replace', access='stream')
    write (9 + id) idx, idy, idz, nx, ny, nz
    if (output_fieldData_precision_cmd == 0) then
        write (9 + id) (((real(u(i, j, k)), i=1, nx), j=1, ny), k=1, nz)
        write (9 + id) (((real(v(i, j, k)), i=1, nx), j=1, ny), k=1, nz)
        write (9 + id) (((real(w(i, j, k)), i=1, nx), j=1, ny), k=1, nz)
        write (9 + id) (((real(rho(i, j, k)), i=1, nx), j=1, ny), k=1, nz)
    else
        write (9 + id) (((u(i, j, k), i=1, nx), j=1, ny), k=1, nz)
        write (9 + id) (((v(i, j, k), i=1, nx), j=1, ny), k=1, nz)
        write (9 + id) (((w(i, j, k), i=1, nx), j=1, ny), k=1, nz)
        write (9 + id) (((rho(i, j, k), i=1, nx), j=1, ny), k=1, nz)
    end if
    close (9 + id)

    if (id == 0) print *, 'Full flow field data saved!'

    return
end subroutine save_macro

!======================================================================================================================================
!---------------------- save data - vtk ----------------------
!======================================================================================================================================
! *********************** 3D legacy VTK writer ***********************
subroutine VTK_legacy_writer_3D(nt)
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    character*30 :: flnm, fmt
    integer :: vtk_type
    integer :: i, j, k, i1, i2, j1, j2, k1, k2, rg, l1, l2, m1, m2, n1, n2, rf, m, nt, wall_indicator
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10, str4*14
    integer   :: ivtk = 9, int
    real(kind=8), allocatable, dimension(:, :, :)::dd, utt, vtt, wtt

    rg = 0; rf = 0; 
    i1 = 1; i2 = nx; 
    j1 = 1; j2 = ny; 
    k1 = 1; k2 = nz; 
    l1 = 1; l2 = nxGlobal; 
    m1 = 1; m2 = nyGlobal; 
    n1 = 1; n2 = nzGlobal; 
    if (id .eq. 0) then
        allocate (dd(l1:l2, m1:m2, n1:n2), utt(l1:l2, m1:m2, n1:n2), vtt(l1:l2, m1:m2, n1:n2), wtt(l1:l2, m1:m2, n1:n2))
    end if
    call compute_macro_vars
    !$acc update host(u,v,w,rho)
    call AllGather(u(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, utt, l1, l2, m1, m2, n1, n2, rf)
    call AllGather(v(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, vtt, l1, l2, m1, m2, n1, n2, rf)
    call AllGather(w(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, wtt, l1, l2, m1, m2, n1, n2, rf)
    call AllGather(rho(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, dd, l1, l2, m1, m2, n1, n2, rf)

    if (id .eq. 0) then
        if (output_fieldData_precision_cmd == 0) then
            fmt = 'float'
        else
            fmt = 'double'
        end if

        write (flnm, '(i10.10,".vtk")') nt

        OPEN(UNIT = ivtk, FILE ='out3.field_data/full_flow_field_'//flnm, FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

        lf = char(10) ! line feed character
        buffer = '# vtk DataFile Version 3.0'//lf
        write (ivtk) trim(buffer)
        buffer = 'vtk output'//lf
        write (ivtk) trim(buffer)
        buffer = 'BINARY'//lf
        write (ivtk) trim(buffer)
        buffer = 'DATASET STRUCTURED_POINTS '//lf
        write (ivtk) trim(buffer)
        write (str1(1:10), '(i10)') nxglobal
        write (str2(1:10), '(i10)') nyglobal
        write (str3(1:10), '(i10)') nzglobal
        buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf
        write (ivtk) trim(buffer)
        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf
        write (ivtk) trim(buffer)
        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'SPACING '//str1//' '//str2//' '//str3//lf
        write (ivtk) trim(buffer)
        num = int(nxGlobal, kind=8)*int(nyGlobal, kind=8)*int(nzGlobal, kind=8)
        write (str4(1:14), '(i14)') num
        buffer = 'POINT_DATA '//str4//lf
        write (ivtk) trim(buffer)

        if (output_fieldData_precision_cmd == 0) then
            fmt = 'float'
            buffer = 'SCALARS density '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((real(dd(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
            buffer = 'SCALARS velocity_X '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((real(utt(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
            buffer = 'SCALARS velocity_Y '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((real(vtt(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
            buffer = 'SCALARS velocity_Z '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((real(wtt(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        else
            fmt = 'double'
            buffer = 'SCALARS density '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((dd(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
            buffer = 'SCALARS velocity_X '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((utt(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
            buffer = 'SCALARS velocity_Y '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((vtt(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
            buffer = 'SCALARS velocity_Z '//fmt//lf
            write (ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write (ivtk) trim(buffer)
            write (ivtk) (((wtt(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        end if

        deallocate (dd, utt, vtt, wtt)

        close (ivtk)
    end if

    return
end subroutine VTK_legacy_writer_3D

subroutine VTK_legacy_writer_3D_half(nt)
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    character*30 :: flnm, fmt
    integer :: vtk_type
    integer :: i, j, k, i1, i2, j1, j2, k1, k2, rg, l1, l2, m1, m2, n1, n2, rf, m, nt, wall_indicator
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10, str4*14
    integer   :: ivtk = 9, int
    real(kind=8), allocatable, dimension(:, :, :)::dd, utt, vtt, wtt

    rg = 0; rf = 0; 
    i1 = 1; i2 = nx; 
    j1 = 1; j2 = ny; 
    k1 = 1; k2 = nz; 
    l1 = 1; l2 = nxGlobal; 
    m1 = 1; m2 = nyGlobal; 
    n1 = 1; n2 = nzGlobal; 
    fmt = 'float'
    if (id .eq. 0) then
        allocate (dd(l1:l2, m1:m2, n1:n2), utt(l1:l2, m1:m2, n1:n2), vtt(l1:l2, m1:m2, n1:n2), wtt(l1:l2, m1:m2, n1:n2))
    end if
    call compute_macro_vars
    !$acc update host(u,v,w,rho)
    call AllGather(u(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, utt, l1, l2, m1, m2, n1, n2, rf)
    call AllGather(v(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, vtt, l1, l2, m1, m2, n1, n2, rf)
    call AllGather(w(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, wtt, l1, l2, m1, m2, n1, n2, rf)
    call AllGather(rho(1:nx, 1:ny, 1:nz), i1, i2, j1, j2, k1, k2, rg, dd, l1, l2, m1, m2, n1, n2, rf)

    if (id .eq. 0) then
        write (flnm, '(i10.10,".vtk")') nt

      OPEN(UNIT = ivtk, FILE ='out3.field_data/full_flow_field_half_'//flnm, FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

        lf = char(10) ! line feed character
        buffer = '# vtk DataFile Version 3.0'//lf
        write (ivtk) trim(buffer)
        buffer = 'vtk output'//lf
        write (ivtk) trim(buffer)
        buffer = 'BINARY'//lf
        write (ivtk) trim(buffer)
        buffer = 'DATASET STRUCTURED_POINTS '//lf
        write (ivtk) trim(buffer)
        write (str1(1:10), '(i10)') nxglobal/2
        write (str2(1:10), '(i10)') nyglobal/2
        write (str3(1:10), '(i10)') nzglobal/2
        buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf
        write (ivtk) trim(buffer)
        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf
        write (ivtk) trim(buffer)
        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'SPACING '//str1//' '//str2//' '//str3//lf
        write (ivtk) trim(buffer)
        num = int(nxGlobal, kind=8)*int(nyGlobal, kind=8)*int(nzGlobal, kind=8)/8
        write (str4(1:14), '(i14)') num
        buffer = 'POINT_DATA '//str4//lf
        write (ivtk) trim(buffer)

        buffer = 'SCALARS density '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(dd(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)
        buffer = 'SCALARS velocity_X '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(utt(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)
        buffer = 'SCALARS velocity_Y '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(vtt(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)
        buffer = 'SCALARS velocity_Z '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(wtt(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)
        deallocate (dd, utt, vtt, wtt)

        close (ivtk)
    end if

    return
end subroutine VTK_legacy_writer_3D_half

!****************** save geometry VTK ***********************
subroutine VTK_walls_bin   !solid geometry
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i, j, k, m, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10, str4*14
    integer   :: ivtk = 9, int

    if (id .eq. 0) then
 open (unit=ivtk, file='out3.field_data/walls_bin.vtk', FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')
        lf = char(10) ! line feed character
        buffer = '# vtk DataFile Version 3.0'//lf; write (ivtk) trim(buffer)
        buffer = 'vtk output'//lf; write (ivtk) trim(buffer)
        buffer = 'BINARY'//lf; write (ivtk) trim(buffer)
        buffer = 'DATASET STRUCTURED_POINTS '//lf; write (ivtk) trim(buffer)

        write (str1(1:10), '(i10)') nxglobal
        write (str2(1:10), '(i10)') nyglobal
        write (str3(1:10), '(i10)') nzglobal
        buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'SPACING '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)
        num = int(nxGlobal, kind=8)*int(nyGlobal, kind=8)*int(nzGlobal, kind=8)
        write (str4(1:14), '(i14)') num
        buffer = 'POINT_DATA '//str4//lf; write (ivtk) trim(buffer)

        !scalar - walls
        buffer = 'SCALARS walls int'//lf; write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf; write (ivtk) trim(buffer)
        write (ivtk) (((int(walls_global(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)

        close (ivtk)

    end if
    return
end subroutine VTK_walls_bin

subroutine VTK_walls_bin_half   !solid geometry
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i, j, k, m, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10, str4*14
    integer   :: ivtk = 9, int

    if (id .eq. 0) then
 open (unit=ivtk, file='out3.field_data/walls_bin.vtk', FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')
        lf = char(10) ! line feed character
        buffer = '# vtk DataFile Version 3.0'//lf; write (ivtk) trim(buffer)
        buffer = 'vtk output'//lf; write (ivtk) trim(buffer)
        buffer = 'BINARY'//lf; write (ivtk) trim(buffer)
        buffer = 'DATASET STRUCTURED_POINTS '//lf; write (ivtk) trim(buffer)

        write (str1(1:10), '(i10)') nxglobal/2
        write (str2(1:10), '(i10)') nyglobal/2
        write (str3(1:10), '(i10)') nzglobal/2
        buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

        write (str1(1:10), '(i10)') 1
        write (str2(1:10), '(i10)') 1
        write (str3(1:10), '(i10)') 1
        buffer = 'SPACING '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

        num = int(nxGlobal, kind=8)*int(nyGlobal, kind=8)*int(nzGlobal, kind=8)/8
        write (str4(1:14), '(i14)') num
        buffer = 'POINT_DATA '//str4//lf; write (ivtk) trim(buffer)

        !scalar - walls
        buffer = 'SCALARS walls int'//lf; write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf; write (ivtk) trim(buffer)
        write (ivtk) (((int(walls_global(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)

        close (ivtk)

    end if
    return
end subroutine VTK_walls_bin_half
