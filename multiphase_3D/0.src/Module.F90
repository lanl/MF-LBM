#include "./preprocessor.h"

!*************************************************************************** Misc_module ****************************************************************************
module Misc_module
    IMPLICIT none
    save
    real(kind=8),parameter :: Pi = 3.14159265358979323846d0
    real(kind=8),parameter :: eps=1.110223025d-16
    real(kind=8) :: convergence_criteria
    integer(kind=1) ,  allocatable, dimension(:,:,:) :: walls, walls_global

    ! job status: new_simulation; continue_simulation; simulation_done; simulation_failed
    character(len=100) :: job_status
    character(len=200) :: geo_file_path, geo_boundary_file_path

    ! error signal to exit the program
    integer :: error_signal   

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ input commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    !1 - simulation ends ; 0 - initial value (exceed max time step) ; 2 - simulation not finished, but save data and exit program; 3 - simulation failed
    integer :: simulation_end_indicator  

    ! Modify geometry command (not modify 1; modify 2)  
    integer :: modify_geometry_cmd  

    ! Initial fluid distribution option: 
    ! 1 - fluid1 (nonwetting phase) at inlet: drainage 
    ! 2 - fluid2 (wetting phase) at inlet: imbibition 
    ! 3 - a fluid1 (nonwetting phase) drop attached to y=1: contact angle measurement
    ! 4 - a fluid2 (wetting phase) drop attached to y=1: contact angle measurement
    ! 5 - a fluid1 drop located in the center of the domain: Laplace law test
    ! 6 - fluid1 and fluid2 raondomely distributed in the pore space with desired saturation: steady state relative perm measurement
    ! 0 or else - error!
    integer :: initial_fluid_distribution_option  

    ! steady state option: 0 - unsteady state simulation; 1 - steady state simulation based on capillary pressure; 
    ! 2 - steady state simulation based on phase field; 3 - steady state simulation based on fractional flowrate
    integer :: steady_state_option

    ! read external geometry cmd:  0 - no; 1 - yes, default name: ../../walldata/walls.dat 
    integer :: external_geometry_read_cmd

    ! geometry preprocess cmd:  0 - process the geometry during the simulation; 1 - load external preprocessed geometry data 
    integer :: geometry_preprocess_cmd 

    ! benchmarking simulation = 1; regular simulation = 0
    integer :: benchmark_cmd    

    ! double back up checkpoint data: 0 - no; 1 - yes (increase storage space)
    integer :: double_bak_checkpoint_pdf_cmd    

    ! place porous plate in the domain (usually near inlet or outlet): 0 - nothing; 1 - block fluid 1; 2 - block fluid 2
    integer :: porous_plate_cmd 

    ! Change inlet fluid phase before initial_interface_position: 0 - nothing; 1 - fluid1; 2- fluid 2
    integer :: change_inlet_fluid_phase_cmd

    ! necesssary modifications for extreme large simulations including number limit and I/O related issues: 0 - no; 1 - yes
    integer :: extreme_large_sim_cmd
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ input commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ domain and geometry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    real(kind=8) ::  la_x,la_y,la_z    !effective simulation domain dimension, in lattice units
    integer :: nxGlobal,nyGlobal,nzGlobal        !full grid: 1 to n?Global

    integer :: iper,jper,kper   !periodic BC indicator

    !domain_wall_status of the domain boundaries: 1 - solid wall; 0 - fluid 
    integer :: domain_wall_status_x_min,domain_wall_status_x_max
    integer :: domain_wall_status_y_min,domain_wall_status_y_max
    integer :: domain_wall_status_z_min,domain_wall_status_z_max

    integer :: nx_sample,ny_sample,nz_sample     !rock sample size
    ! inlet/outlet zone excluded in the bulk properties calculation 
    ! 1<=k<=n_exclude_inlet and nzglobal-n_exclude_outlet+1<=k<=nzglobal
    ! does not affect simulation
    integer :: n_exclude_inlet, n_exclude_outlet      

    integer :: nx,ny,nz   !local MPI domain size
    real(kind=8) ::  A_xy, A_xy_effective, volume_sample   !Axy: cross section area

    ! inlet_BC selection (overridden by x direction periodic BC): 1-velocity; 2-pressure 
    ! outlet_BC selection (overridden by x direction periodic BC): 1-convective; 2-pressure 
    integer :: inlet_BC,outlet_BC 
    real(kind=8) :: target_inject_pore_volume   

    integer :: n_fluid_node_local  !total fluid nodes of local domain

    type indirect_solid_boundary_nodes
        integer :: ix,iy,iz,i_fluid_num   !number of fluid nodes connected to the boundary node
        integer, dimension (1:18) :: neighbor_list
        real(kind=8) :: la_weight              !weight factor used in extrapolation
    end type indirect_solid_boundary_nodes
    type(indirect_solid_boundary_nodes), allocatable, dimension(:) :: solid_boundary_nodes
    integer :: num_solid_boundary_global,num_solid_boundary     ! number of solid boundary nodes
    
    type indirect_fluid_boundary_nodes
        integer :: ix,iy,iz
        real(kind=8) :: nwx,nwy,nwz,cos_theta  !normal direction of the solid surface, local contact angle
    end type indirect_fluid_boundary_nodes
    type(indirect_fluid_boundary_nodes), allocatable, dimension(:) :: fluid_boundary_nodes
    integer :: num_fluid_boundary_global,num_fluid_boundary     ! number of fluid boundary nodes

    integer,  allocatable, dimension(:) :: pore_profile_z   !used in IO for easier data process. Only available in MPI process id0
    integer :: pore_sum, pore_sum_effective  !effeictive pore sum excludes inlet and outlet portion
    real(kind=8) :: porosity_full, porosity_effective

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ lattice ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    !streaming direction - D3Q19
    integer, dimension(:), parameter ::  ex(0:18)=(/0, 1, -1,  0,  0,  0,  0, 1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 /)
    integer, dimension(:), parameter ::  ey(0:18)=(/0, 0,  0,  1, -1,  0,  0, 1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1 /)
    integer, dimension(:), parameter ::  ez(0:18)=(/0, 0,  0,  0,  0,  1, -1, 0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    integer, dimension(:), parameter :: opc(0:18)=(/0, 2,  1,  4,  3,  6,  5, 10, 9,  8,  7, 14, 13, 12, 11, 18, 17, 16, 15 /)        !oppsite directions
    !MRT
    real(kind=8), parameter :: mrt_coef1=1d0/19d0
    real(kind=8), parameter :: mrt_coef2=1d0/2394d0
    real(kind=8), parameter :: mrt_coef3=1d0/252d0
    real(kind=8), parameter :: mrt_coef4=1d0/72d0
    !D3Q19 MODEL
    real(kind=8), dimension(:), parameter :: w_equ(0:18)=(/1d0/3d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,&
        & 1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
    real(kind=8), parameter :: w_equ_0=1.0d0/3.0d0
    real(kind=8), parameter :: w_equ_1=1.0d0/18.0d0
    real(kind=8), parameter :: w_equ_2=1.0d0/36.0d0
    real(kind=8), parameter :: la_vel_norm_i=1.0d0/sqrt(2d0)   !account for the diagonal length


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ unsorted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    !timers
    integer :: ntime0, ntime, ntime_max, ntime_max_benchmark,ntime_macro,ntime_visual,ntime_pdf,ntime_relaxation,ntime_animation
    integer :: ntime_clock_sum, ntime_display_steps,ntime_monitor,ntime_monitor_profile,ntime_monitor_profile_ratio
    integer :: num_slice   !how many monitoring slices (evenly divided through y axis)
    real(kind=8) :: d_sa_vol, sa_vol, d_vol_animation, d_vol_detail, d_vol_monitor, d_vol_monitor_prof
    real(kind=8) :: checkpoint_save_timer,checkpoint_2rd_save_timer, simulation_duration_timer  
    integer :: wallclock_timer_status_sum   !wall clock time reach desired time for each MPI process
    double precision :: wallclock_pdfsave_timer  ! save PDF data for restart simulation based on wall clock timer

    !monitors
    real(kind=8) :: monitor_previous_value , monitor_current_value
    integer, parameter :: time_pt_max = 10      ! maximum 10 points in the history    
    type monitor_time_pt
        real(kind=8) :: v1,v2,v3,it  !  ntime / ntime_monitor
    end type monitor_time_pt
    type(monitor_time_pt),dimension(time_pt_max) :: monitor_pt          
    integer :: n_current_monitor  !current index in the monitor_time_pt array

#ifdef _openacc
    integer :: devnum
#endif

    ! openacc async queue
    integer, parameter :: LBM_kernels=0
    integer, parameter :: LBM_sync=0, LBM_async=1, LBM_async_z=2, LBM_async_y=3, LBM_async_x=4
    integer, PARAMETER :: z_pdf_update_async=5, y_pdf_update_async=6, x_pdf_update_async=7
    integer, PARAMETER :: z_phi_update_async=8, y_phi_update_async=9, x_phi_update_async=10
    integer, PARAMETER :: edge_pdf_async=11, edge_pdf_update_async=12
    integer, PARAMETER :: edge_phi_async=13, edge_phi_update_async=14

end module Misc_module
!*************************************************************************** Misc_module ****************************************************************************



!********************************************************************** Fluid - singlephase *************************************************************************
MODULE Fluid_singlephase
    IMPLICIT NONE
    SAVE
    real(kind=8),ALLOCATABLE,DIMENSION(:,:) :: W_in   !inlet velocity profile
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: u,v,w,rho
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18      !fluid PDF1
    real(kind=8) ::  rt1,rti1,la_nu1,la_nui1 !viscosity for singlephase fluid
    real(kind=8) ::  force_x,force_y,force_z, D_force_z, force_z0   !body force, default flow direction: z
    real(kind=8) ::  uin_max,uin_avg, uin_avg_0,flowrate, p_gradient,rho_out ,rho_in,rho_in_avg,rho_avg_inlet,rho_avg_outlet,rho_drop
    real(kind=8) ::  relaxation,uin_avg_convec,p_max, rho_in_max, umax_global  
    !temporary arrays
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: fl,pre   !flowrate, pressure, saturation vs z
END MODULE Fluid_singlephase
!********************************************************************** Fluid - singlephase *************************************************************************



!********************************************************************** Fluid - multiphase **************************************************************************
MODULE Fluid_multiphase
    IMPLICIT NONE
    SAVE

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ command input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    ! breakthrough check (only works when use open BCs and Ca>0): whether exiting the simulations when the invading phase reaches the outlet
    ! 0 - not check; 1 - check
    integer :: breakthrough_check
    integer :: phase_dist_cmd   !initialization of the phase field: 0 - interface location at z axis; 1 - random distribution according to target saturation; 2 - other distribution
    integer :: fluid_injection_indicator   ! 1-inject fluid 1; 2-inject fluid 2; to be used with buffer zone injection for rel-perm measurement
 

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ fluid ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18      !fluid PDF2
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: phi,cn_x,cn_y,cn_z,c_norm,curv          !order parameter, phase gradient, curvature, norm of color gradient
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: phi_old   !using the phase field difference to indicate steady state

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ flow condition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    !fluid distribution of previous step for convective outlet BC
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: f_convec_bc
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: g_convec_bc
    real(kind=8),ALLOCATABLE,DIMENSION(:,:) :: phi_convec_bc

    real(kind=8), dimension (3) :: w_darcy   !darcy velocity

    real(kind=8), dimension (2) :: kinetic_energy  !kinetic enegy of each phase

    integer :: Z_porous_plate  !porous plate or membrane used to prevent one phase from passing through

    real(kind=8) :: rt2,rti2,la_nu2,la_nui2, phys_nu2
    real(kind=8) :: phi_inlet   !injecting fluid selection
    integer :: outlet_phase1_sum   !total fluid 1 nodes at the outlet (indicating breakthrough)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ multiphase model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    real(kind=8), parameter :: RK_weight2=1.0d0/sqrt(2d0)/36d0   !R-K recolor scheme weight for diagonal directions
    !phase gradient calculation weights
    real(kind=8), parameter, dimension(:) :: ISO4(1:2) = (/1.0d0/6.0d0, 1.0d0/12.0d0/)
    real(kind=8), parameter, dimension(:) :: ISO8(1:7) = (/4.0d0/45.0d0, 1.0d0/21.0d0, 2.0d0/105.0d0, 5.0d0/504.0d0, 1.0d0/315.0d0, 1.0d0/630.0d0, 1.0d0/5040.0d0/)
    !surface tension, contact angle
    real(kind=8) :: gamma,phys_gamma,theta,beta,ca_0, u_filter, ca
    !initial interface position, saturation
    real(kind=8) ::  interface_z0, Sa0,saturation,saturation_full_domain, sa_inject, sa_target, mass1_sum, mass2_sum,saturation_old,vol1_sum, vol2_sum
    integer :: lap_mr   !multirange color gradient model overlaping area
    !temporary arrays
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: fl1,fl2,sa1,mass1,mass2,vol1,vol2   !fluid 1 and fluid 2 fractional flowrate, pressure, saturation vs z
END MODULE Fluid_multiphase
!********************************************************************** Fluid - multiphase **************************************************************************



!********************************************************************** MPI and OpenACC *****************************************************************************
module mpi_variable
    IMPLICIT NONE
    SAVE
    INTEGER, PARAMETER :: mpi_dim = 3
    INTEGER, DIMENSION(1:mpi_dim) :: mpi_coords
    INTEGER, PARAMETER :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4, TAG5 = 5, TAG6 = 6, TAG7 = 7, TAG8 = 8, TAG9 = 9, &
        TAG10 = 10, TAG11 = 11, TAG12 = 12, TAG13 = 13, TAG14 = 14, TAG15 = 15, TAG16 = 16, TAG17 = 17, TAG18 = 18
    INTEGER, PARAMETER :: TAG_phi_zM =19,TAG_phi_zP =20, TAG_phi_yM =21,TAG_phi_yP =22, TAG_phi_xM =23,TAG_phi_xP =24, &
        TAG_phi_xM_yM =25, TAG_phi_xP_yM =26, TAG_phi_xM_yP =27, TAG_phi_xP_yP =28, &
        TAG_phi_xM_zM =29, TAG_phi_xP_zM =30, TAG_phi_xM_zP =31, TAG_phi_xP_zP =32, &
        TAG_phi_yM_zM =33, TAG_phi_yP_zM =34, TAG_phi_yM_zP =35, TAG_phi_yP_zP =36, &
        TAG_phi_xM_yM_zM =37 ,  TAG_phi_xM_yP_zM =38 ,  TAG_phi_xM_yM_zP =39 ,  TAG_phi_xM_yP_zP =40 ,  TAG_phi_xP_yM_zM =41 ,  TAG_phi_xP_yP_zM =42 ,  TAG_phi_xP_yM_zP =43 ,  TAG_phi_xP_yP_zP =44

    INTEGER, DIMENSION(:) :: MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
    INTEGER, DIMENSION(:) :: MPI_REQ_phi_X(4), MPI_REQ_phi_Y(4), MPI_REQ_phi_Z(4)  !macro varaibles such as phi
    INTEGER, DIMENSION(:) :: MPI_REQ_EX(8), MPI_REQ_EY(8), MPI_REQ_EZ(8)  !edge
    INTEGER, DIMENSION(:) :: MPI_REQ_phi_EX(8), MPI_REQ_phi_EY(8), MPI_REQ_phi_EZ(8)  !macro varaibles such as phi, edge

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: MPI_STAT, MPI_ESTAT
    integer :: id0,id,idx,idy,idz,np,npx,npy,npz,MPI_COMM_VGRID   !id0: processor id in comm_world, id: processor id in vgrid    
    integer, parameter :: overlap = 1, overlap_walls = 2, overlap_phi=4  !overlaps for MPI communications=
    
    integer :: idxM,idxP,idyM,idyP,idzM,idzP,ierr   !node connections, x plus, x minus, y+-, z+-
    integer :: idyPzP,idyMzP,idyPzM,idyMzM  !diagonal node connections, x edge
    integer :: idxPzP,idxMzP,idxPzM,idxMzM  !diagonal node connections, y edge
    integer :: idxPyP,idxMyP,idxPyM,idxMyM  !diagonal node connections, z edge
    integer :: idxPyPzP,idxPyMzP,idxPyPzM,idxPyMzM,idxMyPzP,idxMyMzP,idxMyPzM,idxMyMzM  !diagonal node connections, corners
  
    integer :: isize_pdf_x  !buffer size for PDF
    integer :: isize_pdf_y
    integer :: isize_pdf_z
    integer :: isize_pdf_ex !buffer size for pdf, edge
    integer :: isize_pdf_ey
    integer :: isize_pdf_ez
    integer :: isize_phi_x  !buffer size for phi
    integer :: isize_phi_y
    integer :: isize_phi_z
    integer :: isize_phi_ex  !buffer size for phi, edge
    integer :: isize_phi_ey
    integer :: isize_phi_ez
    integer :: isize_phi_corner  !buffer size for phi, corner

    integer :: tk_isize   !!buffer size used in monitor subroutine
    integer :: iz_async, iy_async, ix_async ! MPI boundary layers number for overlaping communication and computation
    !MPI BUFFER, to overcome the allocating memory overhead of Intel Phi
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: send_pdf_xP,send_pdf_xM,send_pdf_yP,send_pdf_yM,send_pdf_zP,send_pdf_zM
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: recv_pdf_xM,recv_pdf_xP,recv_pdf_yM,recv_pdf_yP,recv_pdf_zM,recv_pdf_zP
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: send_pdf_yPzP,send_pdf_yMzP,send_pdf_yPzM,send_pdf_yMzM    !edge send
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: send_pdf_xPzP,send_pdf_xMzP,send_pdf_xPzM,send_pdf_xMzM
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: send_pdf_xPyP,send_pdf_xMyP,send_pdf_xPyM,send_pdf_xMyM
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: recv_pdf_yPzP,recv_pdf_yMzP,recv_pdf_yPzM,recv_pdf_yMzM    !edge recv
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: recv_pdf_xPzP,recv_pdf_xMzP,recv_pdf_xPzM,recv_pdf_xMzM
    real(kind=8),ALLOCATABLE,DIMENSION(:) :: recv_pdf_xPyP,recv_pdf_xMyP,recv_pdf_xPyM,recv_pdf_xMyM

    real(kind=8),ALLOCATABLE,DIMENSION(:) :: tk   !buffer used in monitor subroutine

    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: send_phi_xM,send_phi_xP,send_phi_yM,send_phi_yP,send_phi_zM,send_phi_zP
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: recv_phi_xM,recv_phi_xP,recv_phi_yM,recv_phi_yP,recv_phi_zM,recv_phi_zP
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: send_phi_yPzP,send_phi_yMzP,send_phi_yPzM,send_phi_yMzM    !edge send
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: send_phi_xPzP,send_phi_xMzP,send_phi_xPzM,send_phi_xMzM
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: send_phi_xPyP,send_phi_xMyP,send_phi_xPyM,send_phi_xMyM
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: recv_phi_yPzP,recv_phi_yMzP,recv_phi_yPzM,recv_phi_yMzM    !edge recv
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: recv_phi_xPzP,recv_phi_xMzP,recv_phi_xPzM,recv_phi_xMzM
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: recv_phi_xPyP,recv_phi_xMyP,recv_phi_xPyM,recv_phi_xMyM
    !corner
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: send_phi_xPyPzP,send_phi_xPyMzP,send_phi_xPyPzM,send_phi_xPyMzM,send_phi_xMyPzP,send_phi_xMyMzP,send_phi_xMyPzM,send_phi_xMyMzM
    real(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: recv_phi_xPyPzP,recv_phi_xPyMzP,recv_phi_xPyPzM,recv_phi_xPyMzM,recv_phi_xMyPzP,recv_phi_xMyMzP,recv_phi_xMyPzM,recv_phi_xMyMzM

end module mpi_variable
!********************************************************************** MPI and OpenACC *****************************************************************************
