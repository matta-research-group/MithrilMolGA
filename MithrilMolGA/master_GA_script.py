import subprocess
import argparse
import os
import sys
from datetime import datetime
from MithrilMolGA.ga_slurm import *

# ==============================
# CLI options
# ==============================
options = {
    'run_start': {'default': 1},
    'run_end': {'default': 21},
}

parser = argparse.ArgumentParser()
for arg, opts in options.items():
    parser.add_argument(f'--{arg}', type=type(opts['default']), default=opts['default'])
args = parser.parse_args()

run_start = args.run_start
run_end = args.run_end

# ==============================
# Single GA run
# ==============================
def single_run(run_number, functional='b3lyp', basis_set='6-31g*', 
               time=6, number_of_cpus=10, EG_cutoff=3.2, 
               EG_rank_weight=1, planarity_rank_weight=1, SA_rank_weight=0.25,
               elite_value=25, anioinc_reorg_rank_weight=2, cationic_reorg_rank_weight=2, 
               eg_elite_value=2.5, planarity_elite_value=0.67, anioinc_reorg_elite_value=0.300, catioinc_reorg_elite_value=0.300,
               run_size=50):
    """
    Runs a single iteration of the GA workflow

    Parameters
    ----------
    run_number (int): The number of the run for the GA
    functional (str): The functional to use. Default is 'b3lyp'.
    basis_set (str): The basis_set to use. Default is '6-31g*'.
    time (int): How long each Psi4 calculation is given within the GA. Default is 6.
    number_of_cpus (int): How many cpus each Psi4 calculation is given within the GA. Default is 10.
    EG_cutoff (float): The energy gap cutoff for D-A matching. Default is 3.2 (eV)
    EG_rank_weight (float): How much ranking the energy gap holds, higher the value the more weight. Default is 1.
    planarity_rank_weight (float): How much ranking the planarity holds, higher the value the more weight. Default is 1.
    SA_rank_weight (float): How much ranking synthetic accessibility holds, higher the value the more weight. Default is 0.25
    elite_value (int): What percentage of each run's molecules are carried forward as elites for reorganisation calculations. Default is 25 (%)
    anioinc_reorg_rank_weight (float): How much ranking anionic reorganisation energy holds, higher the value the more weight. Default is 2.
    cationic_reorg_rank_weight (float): How much ranking cationic reorganisation energy holds, higher the value the more weight. Default is 2
    eg_elite_value (float): The cut-off for the elite df for energy gap. Default is 2.5 (eV).
    planarity_elite_value (float): The cut-off for the elite df for planarity. Default is 0.82.
    anioinc_reorg_elite_value (float): The cut-off for the elite df for anionic reorganisation energy. Default is 0.350 (eV).
    catioinc_reorg_elite_value (float): The cut-off for the elite df for cationic reorganisation energy. Default is 0.350 (eV).
    run_size (int): Minimum size of the next run of the GA. Default is 50.

    Returns
    -------
    None
    Executes one loop of the GA, submitting all calculations and updating dataframes.
    """

    base_path = os.path.join(os.getcwd(), "MithrilMolGA")  # absolute path to scripts

    # Monomer run
    subprocess.run([
        sys.executable, os.path.join(base_path, "monomer_run.py"),
        f"--run_num={run_number}",
        f"--functional={functional}",
        f"--basis_set={basis_set}",
        f"--time={time}",
        f"--cpus={number_of_cpus}"
    ], check=True)

    # Molecule run
    subprocess.run([
        sys.executable, os.path.join(base_path, "molecule_run.py"),
        f"--run_num={run_number}",
        f"--functional={functional}",
        f"--basis_set={basis_set}",
        f"--time={time}",
        f"--cpus={number_of_cpus}",
        f"--set_EG_value={EG_cutoff}"
    ], check=True)

    # GA data extraction
    subprocess.run([
        sys.executable, os.path.join(base_path, "GA_data_extraction_step.py"),
        f"--run_num={run_number}"
    ], check=True)

    # GA elite step
    subprocess.run([
        sys.executable, os.path.join(base_path, "GA_elite_step.py"),
        f"--run_num={run_number}",
        f"--functional={functional}",
        f"--basis_set={basis_set}",
        f"--time={time}",
        f"--cpus={number_of_cpus}",
        f"--EG_rank_weight={EG_rank_weight}",
        f"--planarity_rank_weight={planarity_rank_weight}",
        f"--SA_rank_weight={SA_rank_weight}",
        f"--elite_value={elite_value}",
        f"--anioinc_reorg_rank_weight={anioinc_reorg_rank_weight}",
        f"--cationic_reorg_rank_weight={cationic_reorg_rank_weight}",
        f"--eg_elite_value={eg_elite_value}",
        f"--planarity_elite_value={planarity_elite_value}",
        f"--anioinc_reorg_elite_value={anioinc_reorg_elite_value}",
        f"--catioinc_reorg_elite_value={catioinc_reorg_elite_value}"
    ], check=True)

    # Elite mutation
    subprocess.run([
        sys.executable, os.path.join(base_path, "elite_mutation.py"),
        f"--run_num={run_number}",
        f"--run_size={run_size}"
    ], check=True)

# ==============================
# Main GA loop
# ==============================
for i in range(run_start, run_end):
    single_run(
        i, functional='b3lyp', basis_set='6-31g*', time=6, number_of_cpus=6,
        EG_cutoff=3.2, EG_rank_weight=1, planarity_rank_weight=1, SA_rank_weight=0.25,
        elite_value=25, anioinc_reorg_rank_weight=2, cationic_reorg_rank_weight=2,
        eg_elite_value=2.5, planarity_elite_value=0.67,
        anioinc_reorg_elite_value=0.300, catioinc_reorg_elite_value=0.300,
        run_size=200
    )

    # Submit next batch if finished
    if i == (run_end - 1):
        progress_file_path = 'GA_status.txt'
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(progress_file_path, 'a') as file:
            file.write(f'GA complete batch, moving onto new batch at {current_time}.\n')
        
        new_start = run_end
        new_end = new_start + (run_end - run_start)

        write_ga_slurm(new_start, new_end, time=167, cpus=5)
        submit_ga_slurm_job(new_start, max_retries=5, wait_seconds=30)



