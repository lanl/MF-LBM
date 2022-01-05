!*************************************************************************** Misc_module ****************************************************************************
module Misc_module
    IMPLICIT none
    save
    real(kind=8), parameter :: Pi = 3.14159265358979323846d0
    real(kind=8), parameter :: eps = 1d-14
    real(kind=8) :: convergence_criteria
    integer(kind=1), allocatable, dimension(:, :, :) :: walls, walls_global

    ! job status: new_simulation; continue_simulation; simulation_done; simulation_failed
    character(len=100) :: job_status
    character(len=200) :: geo_file_path, geo_boundary_file_path

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ input commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1 - exit as simulation finished ; 0 - continue simulation ; 2 - save checkpoint data and exit; 3 - exit as simulation failed
    integer :: simulation_end_indicator

    ! modify geometry command (not modify 1; modify 2)
    integer :: modify_geometry_cmd

    ! steady state option: 0 - unsteady state simulation; 1 - steady state simulation based on flow rate
    integer :: steady_state_option

    ! read external geometry cmd:  0 - no; 1 - yes, default name: ../../walldata/walls.dat
    integer :: external_geometry_read_cmd

    ! geometry preprocess cmd:  0 - process the geometry during the simulation; 1 - load external preprocessed geometry data
    integer :: geometry_preprocess_cmd

    ! performance benchmarking = 1; regular simulation = 0
    integer :: benchmark_cmd

    ! double back-up checkpoint data: 0 - no; 1 - yes (will increase storage space)
    integer :: double_bak_checkpoint_pdf_cmd

    ! necesssary modifications for extreme large simulations including number limit and I/O related issues: 0 - no; 1 - yes
    integer :: extreme_large_sim_cmd

    ! output field data precision (simulation is always double precision): 0 - single precision; 1 - double precision
    integer :: output_fieldData_precision_cmd
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ input commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ domain and geometry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    real(kind=8) ::  la_x, la_y, la_z    !effective simulation domain dimensions, in lattice units
    integer :: nxGlobal, nyGlobal, nzGlobal        !full grid: 1 to n?Global

    integer :: iper, jper, kper   !periodic BC indicators

    !domain_wall_status of the domain boundaries: 1 - solid wall; 0 - fluid
    integer :: domain_wall_status_x_min, domain_wall_status_x_max
    integer :: domain_wall_status_y_min, domain_wall_status_y_max
    integer :: domain_wall_status_z_min, domain_wall_status_z_max

    integer :: nx_sample, ny_sample, nz_sample     !rock sample size

    ! inlet/outlet zone excluded in the bulk properties calculation
    ! 1<=k<=n_exclude_inlet and nzglobal-n_exclude_outlet+1<=k<=nzglobal
    ! does not affect simulation
    integer :: n_exclude_inlet, n_exclude_outlet

    integer :: nx, ny, nz   !local MPI domain size
    real(kind=8) ::  A_xy, A_xy_effective, volume_sample   !Axy: cross section area
    real(kind=8) :: char_length   ! characteristic length

    ! inlet_BC selection (overridden by x direction periodic BC): 1-velocity; 2-pressure
    ! outlet_BC selection (overridden by x direction periodic BC): 1-convective; 2-pressure
    integer :: inlet_BC, outlet_BC

    integer :: n_fluid_node_local  !total fluid nodes of the partitioned domain

    integer, allocatable, dimension(:) :: pore_profile_z   !used in IO for easier data processing. Only available on rank id0
    integer(kind=8) :: pore_sum, pore_sum_effective  !effeictive pore sum excludes inlet and outlet portion
    real(kind=8) :: porosity_full, porosity_effective

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ lattice ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !streaming direction - D3Q19
    integer, dimension(:), parameter ::  ex(0:18) = (/0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0/)
    integer, dimension(:), parameter ::  ey(0:18) = (/0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1/)
    integer, dimension(:), parameter ::  ez(0:18) = (/0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1/)
    integer, dimension(:), parameter :: opc(0:18) = (/0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15/)        !oppsite directions
    !MRT
    real(kind=8), parameter :: mrt_coef1 = 1d0/19d0
    real(kind=8), parameter :: mrt_coef2 = 1d0/2394d0
    real(kind=8), parameter :: mrt_coef3 = 1d0/252d0
    real(kind=8), parameter :: mrt_coef4 = 1d0/72d0
    real(kind=8), parameter :: mrt_e2_coef1 = 0d0
    real(kind=8), parameter :: mrt_e2_coef2 = -475d0/63d0
    real(kind=8), parameter :: mrt_omega_xx = 0d0
    ! real(kind=8), parameter :: mrt_e2_coef1=3d0
    ! real(kind=8), parameter :: mrt_e2_coef2=-11d0/2d0
    ! real(kind=8), parameter :: mrt_omega_xx=-0.5d0
    !D3Q19 MODEL
   real(kind=8), dimension(:), parameter :: w_equ(0:18)=(/1d0/3d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,&
     1.0d0/18.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
                                                           1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0/)
    real(kind=8), parameter :: w_equ_0 = 1.0d0/3.0d0
    real(kind=8), parameter :: w_equ_1 = 1.0d0/18.0d0
    real(kind=8), parameter :: w_equ_2 = 1.0d0/36.0d0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ unsorted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !timers
   integer :: ntime0, ntime, ntime_max, ntime_max_benchmark, ntime_macro, ntime_visual, ntime_pdf, ntime_relaxation, ntime_animation
    integer :: ntime_clock_sum, ntime_display_steps, ntime_monitor, ntime_monitor_profile, ntime_monitor_profile_ratio
    real(kind=8) :: d_vol_animation, d_vol_detail, d_vol_monitor, d_vol_monitor_prof
    real(kind=8) :: checkpoint_save_timer, checkpoint_2rd_save_timer, simulation_duration_timer

    !monitors
    real(kind=8) :: monitor_previous_value, monitor_current_value

#ifdef _openacc
    integer :: devnum
#endif

    ! openacc async queue
    integer, parameter :: LBM_kernels = 0
    integer, parameter :: LBM_sync = 0, LBM_async = 1, LBM_async_z = 2, LBM_async_y = 3, LBM_async_x = 4
    integer, PARAMETER :: z_pdf_update_async = 5, y_pdf_update_async = 6, x_pdf_update_async = 7
    integer, PARAMETER :: edge_pdf_async = 11, edge_pdf_update_async = 12

end module Misc_module
!*************************************************************************** Misc_module ****************************************************************************

!********************************************************************** Fluid - singlephase *************************************************************************
MODULE Fluid_singlephase
    IMPLICIT NONE
    SAVE
    real(kind=8), ALLOCATABLE, DIMENSION(:, :) :: W_in   !inlet velocity distribution
    real(kind=8), ALLOCATABLE, DIMENSION(:, :, :) :: u, v, w, rho
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18      !fluid PDF
    real(kind=8) ::  rt, rti, la_nu, la_nui !collision relaxation time and viscosity
    real(kind=8) ::  force_x, force_y, force_z, D_force_z, force_z0   !body force, default flow direction: z
    real(kind=8) ::  uin_max, uin_avg, uin_avg_0, flowrate
    real(kind=8) ::  rho_out, rho_in, rho_in_avg, rho_avg_inlet, rho_avg_outlet, rho_drop
    real(kind=8) ::  relaxation, uin_avg_convec, p_max, rho_in_max, umax_global, target_inject_pore_volume

    !MRT relaxation parameters
    real(kind=8) ::   s_e, s_e2, s_q, s_nu, s_pi, s_t
    integer :: mrt_para_preset

    !fluid distribution of previous step for convective outlet BC
    real(kind=8), ALLOCATABLE, DIMENSION(:, :, :) :: f_convec_bc

    ! flow condition
    real(kind=8) :: Re

    !temporary arrays
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: fl, pre   !flowrate, pressure

END MODULE Fluid_singlephase
!********************************************************************** Fluid - singlephase *************************************************************************

!********************************************************************** MPI and OpenACC *****************************************************************************
module mpi_variable
    IMPLICIT NONE
    SAVE
    INTEGER, PARAMETER :: mpi_dim = 3
    INTEGER, DIMENSION(1:mpi_dim) :: mpi_coords
    INTEGER, PARAMETER :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4, TAG5 = 5, TAG6 = 6, TAG7 = 7, TAG8 = 8, TAG9 = 9, &
                          TAG10 = 10, TAG11 = 11, TAG12 = 12, TAG13 = 13, TAG14 = 14, TAG15 = 15, TAG16 = 16, TAG17 = 17, TAG18 = 18

    INTEGER, DIMENSION(:) :: MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
    INTEGER, DIMENSION(:) :: MPI_REQ_EX(8), MPI_REQ_EY(8), MPI_REQ_EZ(8)  !edge

    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: MPI_STAT, MPI_ESTAT
    logical :: mpi_x, mpi_y, mpi_z   ! MPI communication flags along different directions: 0 - no MPI communication; 1 - with MPI communication
    integer :: id0, id, idx, idy, idz, np, npx, npy, npz, MPI_COMM_VGRID   !id0: processor id in comm_world, id: processor id in vgrid
    integer, parameter :: overlap = 1, overlap_walls = 2, overlap_phi = 4  !overlaps for MPI communications=

    integer :: idxM, idxP, idyM, idyP, idzM, idzP, ierr   !node connections, x plus, x minus, y+-, z+-
    integer :: idyPzP, idyMzP, idyPzM, idyMzM  !diagonal node connections, x edge
    integer :: idxPzP, idxMzP, idxPzM, idxMzM  !diagonal node connections, y edge
    integer :: idxPyP, idxMyP, idxPyM, idxMyM  !diagonal node connections, z edge
    integer :: idxPyPzP, idxPyMzP, idxPyPzM, idxPyMzM, idxMyPzP, idxMyMzP, idxMyPzM, idxMyMzM  !diagonal node connections, corners

    integer :: isize_pdf_x  !buffer size for PDF
    integer :: isize_pdf_y
    integer :: isize_pdf_z
    integer :: isize_pdf_ex !buffer size for pdf, edge
    integer :: isize_pdf_ey
    integer :: isize_pdf_ez

    integer :: tk_isize   !!buffer size used in monitor subroutine
    integer :: iz_async, iy_async, ix_async ! MPI boundary layers number for overlaping communication and computation
    !MPI BUFFER, to overcome the allocating memory overhead of Intel Phi
    real(kind=8), ALLOCATABLE, DIMENSION(:, :, :) :: send_pdf_xP, send_pdf_xM, send_pdf_yP, send_pdf_yM, send_pdf_zP, send_pdf_zM
    real(kind=8), ALLOCATABLE, DIMENSION(:, :, :) :: recv_pdf_xM, recv_pdf_xP, recv_pdf_yM, recv_pdf_yP, recv_pdf_zM, recv_pdf_zP
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: send_pdf_yPzP, send_pdf_yMzP, send_pdf_yPzM, send_pdf_yMzM    !edge send
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: send_pdf_xPzP, send_pdf_xMzP, send_pdf_xPzM, send_pdf_xMzM
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: send_pdf_xPyP, send_pdf_xMyP, send_pdf_xPyM, send_pdf_xMyM
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: recv_pdf_yPzP, recv_pdf_yMzP, recv_pdf_yPzM, recv_pdf_yMzM    !edge recv
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: recv_pdf_xPzP, recv_pdf_xMzP, recv_pdf_xPzM, recv_pdf_xMzM
    real(kind=8), ALLOCATABLE, DIMENSION(:) :: recv_pdf_xPyP, recv_pdf_xMyP, recv_pdf_xPyM, recv_pdf_xMyM

    real(kind=8), ALLOCATABLE, DIMENSION(:) :: tk   !buffer used in monitor subroutine

end module mpi_variable
!********************************************************************** MPI and OpenACC *****************************************************************************
