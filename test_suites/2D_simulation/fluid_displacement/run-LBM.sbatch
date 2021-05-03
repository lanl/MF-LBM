#!/bin/bash -l
#SBATCH --job-name="MF-LBM"
#SBATCH --output="%j.out"
#SBATCH --qos=normal
#!SBATCH --qos=long
#SBATCH -N 1
#SBTACH -n 2
#SBATCH -t 10:00:00
#!SBATCH -p power9
#!SBATCH --constraint=gpu_count:4 
#!SBATCH -p power8+
#!SBATCH -p scaling
#!SBATCH -p knl-quad_flat
#!SBATCH -p amd-epyc
#!SBATCH -p skylake-platinum
#!SBATCH -p arm
#SBATCH -p skylake-gold

exec_location="/scratch/users/ychen2018/2d-micromodel/mf-lbm-multiphase/2d-code-multiphase/0.exec/MF_LBM.cpu"

module purge

#power9
#module load pgi/20.1; module load openmpi/p9/4.0.4-pgi_20.1; module load cuda/10.1
#power8+
#module load pgi/18.10; module load openmpi/3.1.3-pgi_18.10;module load cuda/9.2
#Intel-MIC
#module load intel/19.0.2;impi/5.1.1.109-mic
#Intel-CPU
module load intel/19.0.5;module load openmpi/4.0.2-intel_19.0.5
#ARM
#module load openmpi/3.1.3-arm_hpc-compiler-19.0.0

#GPU
#mpirun -np 16 --map-by ppr:4:node --bind-to socket --mca btl ^openib nvprof --profile-from-start off  --print-gpu-trace -o output.%h.%p.%q{OMPI_COMM_WORLD_RANK}.nvvp ../0.exec/MF_LBM.gpu
#mpirun -n 4 --map-by ppr:4:node --bind-to socket --mca btl ^openib ../0.exec/MF_LBM.gpu
#mpirun -n 6 --map-by ppr:2:node --bind-to socket ../0.exec/MF_LBM.gpu

#KNL
#export OMP_NUM_THREADS=68
#export KMP_AFFINITY=noverbose,warnings,respect,granularity=fine,duplicates,scatter,0,0
#mpirun -n 1 numactl --membind=1 ../0.exec/MF_LBM.mic

#CPU
export OMP_NUM_THREADS=22
mpirun -n 2 --map-by ppr:2:node --bind-to numa $exec_location
