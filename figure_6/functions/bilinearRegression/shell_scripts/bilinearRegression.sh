#!/bin/env bash
#
#SBATCH -p all                  # partition (queue)
#SBATCH -n 1                  # number of cores was24
#SBATCH -c 12
#SBATCH -t 350 #time (minutes)
#SBATCH -o logs/%N.%j.out       # STDOUT
#SBATCH -e logs/%N.%j.err       # STDERR
#SBATCH --mem=80000           # in MB was192
#SBATCH --mail-type=FAIL    # notifications for job fail
#SBATCH --mail-user=juliamc@princeton.edu # ADD YOUR EMAIL ADDRESS HERE TO GET EMAIL NOTIFICATIONS
#SBATCH --array=1-772
module load matlab/R2020b

matlab -nodisplay -r  " pause(randi(10)); regVer = 'v35'; savehere = 'dataset_single1'; rawFlag = 1; frameRate = 20;qFile='qLearn_session_all.mat'; try; bilinearEncodingModel($SLURM_ARRAY_TASK_ID,savehere,regVer,0,rawFlag,frameRate,qFile); 	catch me; fprintf('%s\n  %s\n %s\n %d\n',me.identifier,me.message,me.stack.file, me.stack.line); end; exit"

