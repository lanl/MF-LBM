code_location=$(sed '2q;d' "path_info.txt")

sed -i "s|exec_location=.*|exec_location=\"$code_location/0.exec/MF_LBM.cpu\"|g" ./run-LBM.sbatch

mkdir 2d_drop_example    
cd 2d_drop_example 

mkdir single_case_simulation
cd single_case_simulation

mkdir run_script
cp ../../run-LBM.sbatch ./run_script  # job script should already be properly configured
cp ../../multiphase_control.txt ./run_script
cd run_script
echo "new_simulation" > job_status.txt

sbatch run-LBM.sbatch

cd ../../../
