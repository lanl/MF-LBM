#!/bin/bash
#*********************************************************************************#
#                   clean files/directories script
#*********************************************************************************#

echo "Are you sure you want to delete all simulation generated files/directories (except templates)? y/n"  
read user_choice

if [[ "$user_choice" == "y" ]] 
then
    rm -r out1.output
    rm -r out2.checkpoint
    rm -r out3.field_data
    mv irun.sh tobedel_irun.sh
    mv simulation_control.txt tobedel_simulation_control.txt
    mv job_status.txt tobedel_job_status.txt
    mv path_info.txt tobedel_path_info.txt
elif [[! "$user_choice" == "n" ]] 
then
    echo 'Wrong input! Please type y or n'
fi



