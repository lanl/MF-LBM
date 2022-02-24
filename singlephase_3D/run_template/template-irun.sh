#!/bin/bash
#*********************************************************************************#
#                              run script
#*********************************************************************************#
# ------------ command line input ------------
job_status=$1   # new or old or blank; if not 'new' or 'old', the script will check job_status.txt file


# ------------ checking files/direcotries and job status ------------
status_file=job_status.txt
if [[ "$job_status" == "new" ]] 
then
    echo "Starting new simulation!" 
    echo "new_simulation" > job_status.txt
    if [ -d "out1.output" ];then
        rm -r out1.output
    fi
    if [ -d "out2.checkpoint" ];then
        rm -r out2.checkpoint
    fi
    if [ -d "out3.field_data" ];then
        rm -r out3.field_data
    fi
elif [[ "$job_status" == "old" ]] 
then
    echo "Restarting old simulation!" 
    echo "continue_simulation" > job_status.txt
else
    if [ -f "$status_file" ]
    then
        echo "Checking job_status.txt ..."
        line=$(head -n 1 job_status.txt)
        echo "$line"
        if [[ "$line" == "simulation_done" ]]
        then
            echo 'Previous simulation was already finished!'
            exit
        elif [[ "$line" == "simulation_reached_max_step" ]]
        then
            echo 'Previous simulation reached maximum step! Stop!'
            exit
        elif [[ "$line" == "simulation_failed" ]]
        then
            echo 'Previous simulation failed! Stop!'
            exit
        fi
    else
        echo "Error! $status_file does not exist. Stop!"
        exit
    fi
fi

# ------------ run command ------------
run_command


