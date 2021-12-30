#!/bin/bash
#*********************************************************************************#
#                      Simulation Configuration Script
#*********************************************************************************#
#
#--------------------------------------- PATH --------------------------------------------- 
# absolute or relative paths of corresponding files or directories 
template_directory="../run_template/"
exec_location="../1.exec/MF_LBM_singlephase.cpu"

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


#----------------- simulation control ------------------- 
# copy the sample simulation control file and further edit the parameters
cp $template_directory/template-simulation_control.txt simulation_control.txt


# -------------------------- modify interactive run script ----------------------------------

MPI_process_num=2   # MPI_process_num must be equal to mpi_npx * mpi_npy * mpi_npz

# Below is the run command for a typical CPU computing node with two sockets
run_command="mpirun -n $MPI_process_num --map-by ppr:2:node --bind-to numa $exec_location"
# check out other run commands below for different platforms

################# sample run command #################

# ~~~~~ typical CPU computing nodes with 2 CPUs per node
# mpirun -n $MPI_process_num --map-by ppr:2:node --bind-to numa $exec_location

# ~~~~~ single socket consumer-grade machine
# mpirun -n $MPI_process_num $exec_location

# ~~~~~ LANL Darwin IBM POWER9 GPU nodes with 4 GPUs per node
# mpirun -n $MPI_process_num --map-by ppr:4:node --bind-to socket --mca btl ^openib $exec_location
# ~~~~~ nvprof
# mpirun -n $MPI_process_num --map-by ppr:4:node --bind-to socket --mca btl ^openib nvprof --profile-from-start off  --print-gpu-trace -o output.%h.%p.%q{OMPI_COMM_WORLD_RANK}.nvvp $exec_location

# ~~~~~ LANL Kodiak GPU nodes with 4 GPUs per node
# mpirun -n $MPI_process_num --map-by ppr:4:node --bind-to socket $exec_location

# ~~~~~ slurm
# srun -n $MPI_process_num $exec_location

# Intel Xeon Phi 'KNL' node
# export KMP_AFFINITY=noverbose,warnings,respect,granularity=fine,duplicates,scatter,0,0
# mpirun -n $MPI_process_num numactl --membind=1 $exec_location

################# sample run command #################

# Interactive run script (Slurm or other job scripts can be generated easily based on irun.sh)
cp $template_directory/template-irun.sh irun.sh
sed $sed_option "s|run_command|$run_command|" irun.sh

# --------------------------    clean up    ----------------------------------
if [[ "$OS" == "OSX" ]] 
then
    rm *.bak
fi
