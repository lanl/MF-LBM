# this script needs to modify the code to specify initial fluid distribution
# the original code will be copied here and recompiled

original_code_location=$(sed '2q;d' "path_info.txt")

mkdir 2d_layered_example    
cd 2d_layered_example

#copy the code to current location and modify the path
cp -r $original_code_location ./code
code_location="$PWD/code"
sed -i "s|exec_location=.*|exec_location=\"$code_location/0.exec/MF_LBM.cpu\"|g" ../run-LBM.sbatch

#modify the source code to change fluid configration
cd $code_location
sed -i "s/!fluid_interface_position_place_holder/if(y<=interface_x0.or.y>=(nyglobal-interface_x0+1))phi(i,j)=1d0/g" ./0.src/Init_multiphase.F90
make 
cd ../ 

mkdir single_case_simulation
cd single_case_simulation

mkdir run_script
cp ../../run-LBM.sbatch ./run_script  # job script should already be properly configured
cp ../../multiphase_control.txt ./run_script
cd run_script
echo "new_simulation" > job_status.txt

sbatch run-LBM.sbatch

cd ../../../
