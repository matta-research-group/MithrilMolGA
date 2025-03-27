#!/bin/bash --login
#SBATCH -o test_run.out 
#SBATCH -e test_run.err 
#SBATCH --job-name=test_run 
#SBATCH -p cpu 
#SBATCH --ntasks=10
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=46:30:00 
 
module purge 
module load cuda/10.0.130-gcc-13.2.0 
python3 master_GA_script.py
