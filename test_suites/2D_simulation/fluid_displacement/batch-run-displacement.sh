function simulation_kernel () {
    geometry_name=$1    #geometry file name
    theta=$2
    ca=$3
    logM=$4
    vis_invade=$5

    folder_name=$geometry_name"_theta"$theta"ca"$ca"logM"$logM"vis"$vis_invade
    mkdir $folder_name
    #echo $folder_name
    cd $folder_name

    mkdir run_script
    cp ../../run-LBM.sbatch ./run_script  # job script should already be properly configured
    cp ../../multiphase_control.txt ./run_script
    cd run_script
    echo "new_simulation" > job_status.txt

    sed -i "s/geometry_name .*/geometry_name $geometry_name/g" ./multiphase_control.txt
    sed -i "s/fluid1_viscosity .*/fluid1_viscosity $vis_invade/g" ./multiphase_control.txt
    sed -i "s/fluid2_viscosity .*/fluid2_viscosity -1/g" ./multiphase_control.txt   #negative value so that viscosity_ratio_log is enabled
    sed -i "s/viscosity_ratio_log .*/viscosity_ratio_log $logM/g" ./multiphase_control.txt
    sed -i "s/theta .*/theta $theta/g" ./multiphase_control.txt
    sed -i "s/capillary_number .*/capillary_number "$ca"d-6/g" ./multiphase_control.txt
    
    script_name=$folder_name
    mv run-LBM.sbatch $script_name   # easily identify the jobs
    sbatch $script_name

    cd ../../
    return 0
}

#code should be compiled with corresponding compiler and options consistent with the job script
#change of initial fluid distribution has to be done in the on Init_multiphase.F90 which requires recompile of the code
code_location=$(sed '2q;d' "path_info.txt")
walldata_folder_location=$(sed '8q;d' "path_info.txt")  

sed -i "s|exec_location=.*|exec_location=\"$code_location/0.exec/MF_LBM.cpu\"|g" ./run-LBM.sbatch

mkdir 2d_displacement_batch_example    
cd 2d_displacement_batch_example   

mkdir walldata                   #default folder to store the geometry files
cp $walldata_folder_location/*.dat ./walldata

for geometry_name in "walls" "dog"
do
    #for theta in 10 40 70;
    for theta in 10 70;
    do
    #for ca in 01 05 25 125;    # multiply 10^-6
        for ca in 50;
        do
            #for M in -1.4 -0.7 0 0.7 1.4
            for M in -0.7;
            do
                #for vis_invade in 0.0017 0.0075 0.034;
                for vis_invade in 0.0017;
                do
                    simulation_kernel $geometry_name $theta $ca $M $vis_invade
                done
            done
        done
    done
done
