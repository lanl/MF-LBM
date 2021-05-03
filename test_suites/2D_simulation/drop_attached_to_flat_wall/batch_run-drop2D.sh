function simulation_kernel () {
    code_location=$1
    viscosity_r=$2
    viscosity_b=$3
    surface_tension=$4
    theta=$5   #theta mesured through fluid b(blue)
    nx=$6
    ny=$7
    drop_radius=$8

    folder_name="vis_r"$viscosity_r"vis_b"$viscosity_b"sigma"$surface_tension"theta"$theta
    mkdir $folder_name
    cd $folder_name
    mkdir run_script
    cp ../../run-LBM.sbatch ./run_script  # job script should already be properly configured
    cp ../../multiphase_control.txt ./run_script
    cd run_script
    echo "new_simulation" > job_status.txt


    sed -i "s/lattice_dimensions .*/lattice_dimensions $nx,$ny/g" ./multiphase_control.txt
    sed -i "s/fluid1_viscosity .*/fluid1_viscosity $viscosity_r/g" ./multiphase_control.txt
    sed -i "s/fluid2_viscosity .*/fluid2_viscosity $viscosity_b/g" ./multiphase_control.txt
    sed -i "s/surface_tension .*/surface_tension $surface_tension/g" ./multiphase_control.txt
    sed -i "s/initial_interface_position .*/initial_interface_position $drop_radius/g" ./multiphase_control.txt
    sed -i "s/theta .*/theta $theta/g" ./multiphase_control.txt

    script_name=$folder_name
    mv run-LBM.sbatch $script_name   # easily identify the jobs
    sbatch $script_name

    cd ../../
    return 0
}

#code should be compiled with corresponding compiler and options consistent with the job script
#change of initial fluid distribution has to be done in the on Init_multiphase.F90 which requires recompile of the code
code_location=$(sed '2q;d' "path_info.txt")

sed -i "s|exec_location=.*|exec_location=\"$code_location/0.exec/MF_LBM.cpu\"|g" ./run-LBM.sbatch

mkdir 2d_drop_batch_example        
cd 2d_drop_batch_example  

#========================== batch case ==========================
viscosity_r=0.003  #red fluid
viscosity_b=0.3    #blue fluid
surface_tension=0.03
nx=200
ny=100
drop_radius=40
for theta in 30 60 90 120 150;    
do
    simulation_kernel $code_location $viscosity_r $viscosity_b $surface_tension $theta $nx $ny $drop_radius
done

viscosity_r=0.166667  #red fluid
viscosity_b=0.166667    #blue fluid
for theta in 30 60 90 120 150;    
do
    simulation_kernel  $code_location $viscosity_r $viscosity_b $surface_tension $theta $nx $ny $drop_radius
done