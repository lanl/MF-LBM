#*********************************************************************************#
#                      Simulation Configuration Script
#*********************************************************************************#
#
#--------------------------------------- PATH --------------------------------------------- 
# absolute or relative paths of corresponding files or directories 
template_directory="../../../multiphase_3D/run_template/"
exec_location="../../../multiphase_3D/1.exec/MF_LBM.cpu"

# if external geometry is used (external_geometry_read_cmd = 1)
geometry_file="placeholder"
# if pre-computed boundary info is used (geometry_preprocess_cmd = 2)
geometry_boundary_info_file="placeholder"

#------------------------ choose Linux or OSX system --------------------------------------
# sed command is slightly different between Linux and OSX. 
# Uncomment OSX below if you are using a MAC
# OS="OSX"

# --- Linux
sed_option="-i"
# --- OSX           sed -i.bak should be used. One can delete .bak files afterwards 
if [[ "$OS" == "OSX" ]] 
then
    sed_option="-i.bak"
fi



# -------------------------------------- modify path --------------------------------------
#  write above paths to path_info.txt
cp $template_directory/template-path_info.txt path_info.txt
sed $sed_option "s|exec_path_holder|$exec_location|" path_info.txt
sed $sed_option "s|geometry_file_path_holder|$geometry_file|" path_info.txt
sed $sed_option "s|geometry_boundary_file_path_holder|$geometry_boundary_info_file|" path_info.txt
cp $template_directory/clean.sh ./


#----------------- parameters to be modified for the present simulation ------------------- 
# Only the most common parameters that require modification are listed here with default value.
# Comment out the ones you want to change, and modify the values.
# One can further change the script to modify the control file or edit the control file directly.

cp $template_directory/template-simulation_control.txt simulation_control.txt

# initial_fluid_distribution_option=1
# sed $sed_option "s|initial_fluid_distribution_option .*|initial_fluid_distribution_option $initial_fluid_distribution_option|g" ./simulation_control.txt
# benchmark_cmd=0
# sed $sed_option "s|benchmark_cmd .*|benchmark_cmd $benchmark_cmd|g" ./simulation_control.txt
# breakthrough_check=0
# sed $sed_option "s|breakthrough_check .*|breakthrough_check $breakthrough_check|g" ./simulation_control.txt
# steady_state_option=0
# sed $sed_option "s|steady_state_option .*|steady_state_option $steady_state_option|g" ./simulation_control.txt
# convergence_criteria=1d-6
# sed $sed_option "s|convergence_criteria .*|convergence_criteria $convergence_criteria|g" ./simulation_control.txt
# external_geometry_read_cmd=0
# sed $sed_option "s|external_geometry_read_cmd .*|external_geometry_read_cmd $external_geometry_read_cmd|g" ./simulation_control.txt
# geometry_preprocess_cmd=0
# sed $sed_option "s|geometry_preprocess_cmd .*|geometry_preprocess_cmd $geometry_preprocess_cmd|g" ./simulation_control.txt
# extreme_large_sim_cmd=0
# sed $sed_option "s|extreme_large_sim_cmd .*|extreme_large_sim_cmd $extreme_large_sim_cmd|g" ./simulation_control.txt

# domain info
# nx=40
# ny=40
# nz=60
# sed $sed_option "s|lattice_dimensions .*|lattice_dimensions $nx,$ny,$nz|g" ./simulation_control.txt
# exclude_layers_inlet=10
# exclude_layers_outlet=10
# sed $sed_option "s|excluded_layers .*|excluded_layers $exclude_layers_inlet,$exclude_layers_outlet|g" ./simulation_control.txt
# periodic_x=0
# periodic_y=0
# periodic_z=0
# sed $sed_option "s|periodic_indicator .*|periodic_indicator $periodic_x,$periodic_y,$periodic_z|g" ./simulation_control.txt
# mpi_npx=1  # x-axis domain decomposition currently disabled, keep using mpi_npx=1 
# mpi_npy=1
# mpi_npz=2
# sed $sed_option "s|MPI_process_num .*|MPI_process_num $mpi_npx,$mpi_npy,$mpi_npz|g" ./simulation_control.txt

# fluid properties
# fluid1_viscosity=0.004
# fluid2_viscosity=0.06
# sed $sed_option "s|fluid1_viscosity .*|fluid1_viscosity $fluid1_viscosity|g" ./simulation_control.txt
# sed $sed_option "s|fluid2_viscosity .*|fluid2_viscosity $fluid2_viscosity|g" ./simulation_control.txt
# surface_tension=0.03
# sed $sed_option "s|surface_tension .*|surface_tension $surface_tension|g" ./simulation_control.txt
# theta=30
# sed $sed_option "s|theta .*|theta $theta|g" ./simulation_control.txt

# flow condition
# target_inject_pore_volume=1.0
# sed $sed_option "s|target_inject_pore_volume .*|target_inject_pore_volume $target_inject_pore_volume|g" ./simulation_control.txt
# initial_interface_position=10.0
# sed $sed_option "s|initial_interface_position .*|initial_interface_position $initial_interface_position|g" ./simulation_control.txt
# capillary_number=100d-6
# sed $sed_option "s|capillary_number .*|capillary_number $capillary_number|g" ./simulation_control.txt
# saturation_injection=1.0
# sed $sed_option "s|saturation_injection .*|saturation_injection $saturation_injection|g" ./simulation_control.txt
# body_force_0=0d-6
# sed $sed_option "s|body_force_0 .*|body_force_0 $body_force_0|g" ./simulation_control.txt

# timers
# ntime_visual=10000000
# sed $sed_option "s|ntime_visual .*|ntime_visual $ntime_visual|g" ./simulation_control.txt
# ntime_animation=2000
# sed $sed_option "s|ntime_animation .*|ntime_animation $ntime_animation|g" ./simulation_control.txt
# monitor_timer=1000
# sed $sed_option "s|monitor_timer .*|monitor_timer $monitor_timer|g" ./simulation_control.txt
# display_steps_timer=1000
# sed $sed_option "s|display_steps_timer .*|display_steps_timer $display_steps_timer|g" ./simulation_control.txt
# checkpoint_save_timer=2.0
# sed $sed_option "s|checkpoint_save_timer .*|checkpoint_save_timer $checkpoint_save_timer|g" ./simulation_control.txt
# simulation_duration_timer=15.7
# sed $sed_option "s|simulation_duration_timer .*|simulation_duration_timer $simulation_duration_timer|g" ./simulation_control.txt

# d_vol_animation=0.05
# sed $sed_option "s|d_vol_animation .*|d_vol_animation $d_vol_animation|g" ./simulation_control.txt
# d_vol_monitor=0.01
# sed $sed_option "s|d_vol_monitor .*|d_vol_monitor $d_vol_monitor|g" ./simulation_control.txt



# -------------------------- modify interactive run script ----------------------------------

MPI_process_num=2   # MPI_process_num must be equal to mpi_npx * mpi_npy * mpi_npz

# Below is the run command for typical CPU computing node with two sockets
run_command="mpirun -n $MPI_process_num --map-by ppr:2:node --bind-to numa $exec_location"
# check out other run commands below for different platforms

################# sample run command #################

# ~~~~~ typical CPU computing nodes with 2 CPUs per node
# mpirun -n $MPI_process_num --map-by ppr:2:node --bind-to numa $exec_location

# ~~~~~ local consumer machine
# mpirun -n $MPI_process_num $exec_location

# ~~~~~ LANL Darwin IBM POWER9 GPU nodes with 4 GPUs per node
# mpirun -n $MPI_process_num --map-by ppr:4:node --bind-to socket --mca btl ^openib $exec_location
# ~~~~~ nvprof
# mpirun -n $MPI_process_num --map-by ppr:4:node --bind-to socket --mca btl ^openib nvprof --profile-from-start off  --print-gpu-trace -o output.%h.%p.%q{OMPI_COMM_WORLD_RANK}.nvvp $exec_location

# ~~~~~ LANL Kodiak GPU nodes with 4 GPUs per node
# mpirun -n $MPI_process_num --map-by ppr:4:node --bind-to socket $exec_location

# ~~~~~ use srun
# srun -n $MPI_process_num $exec_location

# Intel Xeon Phi 'KNL' node
# export KMP_AFFINITY=noverbose,warnings,respect,granularity=fine,duplicates,scatter,0,0
# mpirun -n $MPI_process_num numactl --membind=1 $exec_location

################# sample run command #################

# Interactive run script. Slurm job scripts can be generated easily based on irun.sh
cp $template_directory/template-irun.sh irun.sh
sed $sed_option "s|run_command|$run_command|" irun.sh

# --------------------------    clean up    ----------------------------------
if [[ "$OS" == "OSX" ]] 
then
    rm *.bak
fi