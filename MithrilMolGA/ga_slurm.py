import subprocess
from QCflow.slurm import *

def write_ga_slurm(run_start, run_end, time=167, cpus=10):


    start_str = str(run_start)

    file_name = f'run_GA_batch_start_{start_str}.sh'
    
    calc_time = f'{time}:00:00'

    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o GA_manager_{start_str}.out \n')
        file.write(f'#SBATCH -e GA_manager_{start_str}.err \n')#
        file.write(f'#SBATCH --job-name=GA_manager_{start_str} \n')
        file.write(f'#SBATCH -p long_cpu \n')
        file.write(f'#SBATCH --mail-user=k2255489@kcl.ac.uk \n')
        file.write(f'#SBATCH --mail-type=BEGIN,END,FAIL \n')
        file.write(f'#SBATCH --ntasks={cpus}\n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --cpus-per-task=1 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time={calc_time} \n') #reduced to speed up queue time
        file.write(' \n')
        file.write(f'module purge \n')
        file.write(f'module load cuda/10.0.130-gcc-13.2.0 \n')
        file.write(f'source /scratch/users/k2255489/miniconda3/etc/profile.d/conda.sh \n')
        file.write(f'conda activate QCflow \n')
        file.write(f'echo "Running on host: $(hostname)" \n')
        file.write(f'echo "Python path: $(which python)" \n')
        file.write(f'python --version \n')
        file.write(f'python3 master_GA_script.py --run_start {run_start} --run_end {run_end} \n')

def submit_ga_slurm_job(run_start, max_retries=5, wait_seconds=30):
    """
    Submits a SLURM job using the specified job name and molecule name. Works on the KCL CREATE HPC.

    Parameters
    ----------
    run_start (float): The number of the run of the GA to submit at start of the batch

    max_retries (int): The maximum amount of times a job will attempt to submit. Deafult is 5.

    wait_seconds (int): How long python will go to sleep inbetween attempts to submit

    Returns
    -------
    bytes: The standard output from the SLURM job submission command.
    """

    start_str = str(run_start)

    string = f'sbatch run_GA_batch_start_{start_str}.sh'

    submission_name = f'run_GA_batch_start_{start_str}.sh'

    for attempt in range(1, max_retries + 1):
        if is_job_in_queue(submission_name):
            print(f"Job {submission_name} is already in the queue. Skipping submission.")
            return None

        try:
            process = subprocess.run(
                string,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                check=True
            )
            print(f"Successfully submitted GA_manager job")
            return process.stdout

        except subprocess.CalledProcessError as e:
            print(f"[Attempt {attempt}] Error submitting job GA_manager: {e.stderr.decode().strip()}")
            if attempt < max_retries:
                print(f"Waiting {wait_seconds} seconds before checking and retrying...")
                time.sleep(wait_seconds)
            else:
                print("Max retries reached. Moving on.")
                return None

