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
    rm irun.sh
    rm simulation_control.txt
    rm job_status.txt
    rm path_info.txt
elif [[! "$user_choice" == "n" ]] 
then
    echo 'Wrong input! Please type y or n'
fi



