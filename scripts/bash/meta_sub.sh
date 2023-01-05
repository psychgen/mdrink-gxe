#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=meta_sub_pgs
#
#Project:
#SBATCH --account=p471
#
#Wall clock limit (change based on resources required):
#SBATCH --time=00:15:00
#
#SBATCH --ntasks=1
#
#Output filename customization
#This specification is Jobname_User_JobID
#SBATCH --output=./%x_%u_%j.out
#
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=2G

## Set up job enviroment:
module init

mkdir -p ./reports

for FILE in ./*model.sh; do
echo ${FILE}
sbatch ${FILE}
sleep 1 # pause to be kind to the scheduler
done
