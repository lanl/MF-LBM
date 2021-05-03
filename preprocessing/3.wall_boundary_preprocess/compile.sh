# due to different implementation of MPI_Type_struct (padding) of different compilers
# the same compiler should be used to compile this preprocessing code and the main simulation code

# gcc
gfortran -ffree-line-length-300 -fopenmp -O2 3d_wall_boundary_process.f90 -o a.out
# PGI
# pgf90 -mp -O2 3d_wall_boundary_process.f90 -o a.out 
# Intel
# ifort -qopenmp -O2 3d_wall_boundary_process.f90 -o a.out



