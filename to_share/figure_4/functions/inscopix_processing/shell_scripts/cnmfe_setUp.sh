#!/bin/env bash
#
#SBATCH -p all                  # partition (queue)
#SBATCH -n 4                  # number of cores was24
#SBATCH -t 40                 # time (minutes)
#SBATCH -o logs/%N.%j.out       # STDOUT
#SBATCH -e logs/%N.%j.err       # STDERR
#SBATCH --mem=150000            # in MB was192
#SBATCH --mail-type=FAIL    # notifications for job done & fail

#SBATCH --mail-user=juliamc@princeton.edu # ADD YOUR EMAIL ADDRESS HERE TO GET EMAIL NOTIFICATIONS

echo "$1 $2"

module load matlab/R2018b

matlab -nojvm -nodisplay -r "try inscopixDffExtract_largeFile({'$1'}, {{'$2'}}, {'ACC_DMS_imaging'},1,0); catch me; fprintf('%s / %s\n',me.identifier,me.message); end; exit"
