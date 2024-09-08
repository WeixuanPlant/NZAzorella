#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G 
#SBATCH --time=20-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --output="runsnaq_slurm%a.log"
#SBATCH --job-name="jointGeno_snq"
#SBATCH --array=0-10%2

## --array: to run multiple instances of this script,
##          one for each value in the array.
##          1 instance = 1 task
## -J job name
## -c number of cores (CPUs) per task

echo "slurm task ID = $SLURM_ARRAY_TASK_ID used as hmax"
echo "start of SNaQ parallel runs on $(hostname)"

module load curl/8.0.1-y76nj72
module load openblas/0.3.23-suvfvn7
module load git/2.40.0-dotnjuf
module load py-numpy/1.26.3-py310-gntgk2n
module load julia/1.8.3-py310-tgf4scr

# finally: launch the julia script, using Julia executable appropriate for slurm, with full paths:
julia --history-file=yes -- runSNaQ.jl $SLURM_ARRAY_TASK_ID 15 > net${SLURM_ARRAY_TASK_ID}_15runs.screenlog 2>&1
echo "end of SNaQ run ..."
