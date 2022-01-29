#!/bin/bash
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
# write above paths to path_info.txt
cp $template_directory/template-path_info.txt path_info.txt
sed $sed_option "s|exec_path_holder|$exec_location|" path_info.txt
sed $sed_option "s|geometry_file_path_holder|$geometry_file|" path_info.txt
sed $sed_option "s|geometry_boundary_file_path_holder|$geometry_boundary_info_file|" path_info.txt
cp $template_directory/clean.sh ./


#---------------------- copy and modify the control file template -------------------------- 
cp $template_directory/template-simulation_control.txt simulation_control.txt

initial_fluid_distribution_option=3
sed $sed_option "s|initial_fluid_distribution_option .*|initial_fluid_distribution_option $initial_fluid_distribution_option|g" ./simulation_control.txt
steady_state_option=1
sed $sed_option "s|steady_state_option .*|steady_state_option $steady_state_option|g" ./simulation_control.txt
nx=50
ny=80
nz=80
sed $sed_option "s|lattice_dimensions .*|lattice_dimensions $nx,$ny,$nz|g" ./simulation_control.txt
periodic_x=0
periodic_y=1
periodic_z=1
sed $sed_option "s|periodic_indicator .*|periodic_indicator $periodic_x,$periodic_y,$periodic_z|g" ./simulation_control.txt
theta=45
sed $sed_option "s|theta .*|theta $theta|g" ./simulation_control.txt
initial_interface_position=18.0
sed $sed_option "s|initial_interface_position .*|initial_interface_position $initial_interface_position|g" ./simulation_control.txt
display_steps_timer=100
sed $sed_option "s|display_steps_timer .*|display_steps_timer $display_steps_timer|g" ./simulation_control.txt
domain_wall_status_y=0
sed $sed_option "s|domain_wall_status_y .*|domain_wall_status_y $domain_wall_status_y, $domain_wall_status_y|g" ./simulation_control.txt


# -------------------------- modify interactive run script ----------------------------------
MPI_process_num=2   # MPI_process_num must be equal to mpi_npx * mpi_npy * mpi_npz

# Below is the run command for typical CPU computing node with two NUMAs
run_command="mpirun -n $MPI_process_num --map-by ppr:2:node --bind-to numa $exec_location"
# check out other run commands below for different platforms

################# sample run command #################

# ~~~~~ typical CPU computing nodes with two NUMAs per node
# mpirun -n $MPI_process_num --map-by ppr:2:node --bind-to numa $exec_location

# ~~~~~ consumer-grade desktop or laptop
# mpirun -n $MPI_process_num $exec_location

# ~~~~~ GPU nodes with 4 GPUs per node
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
