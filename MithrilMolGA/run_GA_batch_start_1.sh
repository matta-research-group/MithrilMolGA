#!/bin/bash --login
#SBATCH -o GA_manager_1.out 
#SBATCH -e GA_manager_1.err 
#SBATCH --job-name=GA_manager_achive_test
#SBATCH -p long_cpu
#SBATCH --mail-user=k2255489@kcl.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL    
#SBATCH --ntasks=10
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=99:00:00 
 
module purge 
module load cuda/10.0.130-gcc-13.2.0 
source /scratch/users/k2255489/miniconda3/etc/profile.d/conda.sh 
conda activate QCflow 
echo "Running on host: $(hostname)" 
echo "Python path: $(which python)" 
python --version 
python3 master_GA_script.py --run_start 1 --run_end 3
