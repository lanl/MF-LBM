#!/bin/bash
#*********************************************************************************#
#                      Simulation Configuration Script
#*********************************************************************************#
#
#--------------------------------------- PATH --------------------------------------------- 
# absolute or relative paths of corresponding files or directories 
template_directory="../../../multiphase_3D/run_template/"
exec_location="../../../multiphase_3D/1.exec/MF_LBM.gpu"

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

breakthrough_check=1
sed $sed_option "s|breakthrough_check .*|breakthrough_check $breakthrough_check|g" ./simulation_control.txt
modify_geometry_cmd=1
sed $sed_option "s|modify_geometry_cmd .*|modify_geometry_cmd $modify_geometry_cmd|g" ./simulation_control.txt

mpi_npx=1  # x-axis domain decomposition currently disabled, keep using mpi_npx=1 
mpi_npy=1
mpi_npz=1
sed $sed_option "s|MPI_process_num .*|MPI_process_num $mpi_npx,$mpi_npy,$mpi_npz|g" ./simulation_control.txt

# -------------------------- modify interactive run script ----------------------------------
MPI_process_num=1   # MPI_process_num must be equal to mpi_npx * mpi_npy * mpi_npz

run_command="mpirun -n $MPI_process_num $exec_location"

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
