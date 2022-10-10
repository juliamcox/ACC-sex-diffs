#!/bin/env bash
#
#SBATCH -p all                  # partition (queue)
#SBATCH -n 2                  # number of cores was24
#SBATCH -c 24
#SBATCH -t 700              # time (minutes)
#SBATCH -o logs/%N.%j.out       # STDOUT
#SBATCH -e logs/%N.%j.err       # STDERR
#SBATCH --mem=250000            # in MB was192
#SBATCH --mail-type=FAIL    # notifications for job done & fail

#SBATCH --mail-user=juliamc@princeton.edu # ADD YOUR EMAIL ADDRESS HERE TO GET EMAIL NOTIFICATIONS

echo "SBATCH job submitted for $1 $2 rawFlag $3 statMethod $4"

module load matlab/R2020b

matlab -nodisplay -r "try; imagingRegression_2020('$1','$2',1,$3,20,'$4',1000,1); 	catch me; fprintf('%s\n  %s\n %s\n %d\n',me.identifier,me.message,me.stack.file, me.stack.line); end; exit"

