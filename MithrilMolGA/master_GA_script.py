import subprocess
import argparse

def single_run(run_number, functional='b3lyp', basis_set='6-31g*', 
               time=6, number_of_cpus=10, EG_cutoff=3.2, 
               EG_rank_weight=1, planarity_rank_weight=1, SA_rank_weight=0.25,
               elite_value=25, anioinc_reorg_rank_weight=2, cationic_reorg_rank_weight=2, 
               eg_elite_value=2.5, planarity_elite_value=0.82, anioinc_reorg_elite_value=0.350, catioinc_reorg_elite_value=0.350,
               run_size=50):
    """
    Runs a single interation of th GA workflow

    Parameters
    ----------
    run_number (int): The number of the run for the GA
    functional (str): The functional to use. Default is 'b3lyp'.
    basis_set (str): The basis_set to use. Default is '6-31g*'.
    time (int): How long each Psi4 calculation is given within the GA. Default is 6.
    number_of_cpus (int): How many cpus each Psi4 calculation is given within the GA. Default is 10.
    EG_cutoff (int): The energy gap cutoff for D-A matching. Default is 3.2 (eV)
    EG_rank_weight (int): How much ranking the energy gap holds, higher the value the more weight. Default is 1.
    planarity_rank_weight (int): How much ranking the planarity holds, higher the value the more weight. Default is 1.
    SA_rank_weight (int): How much ranking sysnethtic accessibility holds, higher the value the more weight. Default is 0.25
    elite_value (int): What percentage of each runs molecules are carried forward as elites for reorganisation calucltions. Default is 25 (%)
    anioinc_reorg_rank_weight (int): How much ranking anioinc reorganisation energy holds, higher the value the more weight. Defeault is 2.
    cationic_reorg_rank_weight (int): How much ranking catioinc reorganisation energy holds, higher the value the more weight. Defeault is 2
    eg_elite_value (int): The cut off for the elite df for energy gap. Default is 2.5 (eV).
    planarity_elite_value (int): The cut off for the elite df for planarity_elite_value. Default is 0.82.
    anioinc_reorg_elite_value (int): The cut off for the elite df for anioinc_reorg_elite_value. Default is 0.350 (eV).
    catioinc_reorg_elite_value (int): The cut off for the elite df for catioinc_reorg_elite_value. Default is 0.350 (eV).
    run_size (int): Minimum size of th next run of the GA. Default is 50.

    Returns
    -------
    Runs one loop of the GA and returns all the data into the specified folders.
    """
    subprocess.run(["python", "monomer_run.py", f"--run_num={run_number}", 
                    f"--functional={functional}", f"--basis_set={basis_set}", 
                    f"--time={time}", f"--cpus={number_of_cpus}"])  # runs the fragment systems
    
    subprocess.run(["python", "molecule_run.py", f"--run_num={run_number}", 
                    f"--functional={functional}", f"--basis_set={basis_set}", 
                    f"--time={time}", f"--cpus={number_of_cpus}", f"--set_EG_value={EG_cutoff}"])  # calculates D-A matching, submits molecules
    
    subprocess.run(["python", "GA_data_extraction_step.py", f"--run_num={run_number}"])  # extracts data and updates a df
    
    subprocess.run(["python", "GA_elite_step.py", f"--run_num={run_number}", 
                    f"--functional={functional}", f"--basis_set={basis_set}", 
                    f"--time={time}", f"--cpus={number_of_cpus}", f"--EG_rank_weight={EG_rank_weight}", 
                    f"--planarity_rank_weight={planarity_rank_weight}", f"--SA_rank_weight={SA_rank_weight}", 
                    f"--elite_value={elite_value}", f"--anioinc_reorg_rank_weight={anioinc_reorg_rank_weight}", 
                    f"--cationic_reorg_rank_weight={cationic_reorg_rank_weight}",
                    f"--eg_elite_value={eg_elite_value}", f"--planarity_elite_value={planarity_elite_value}", 
                    f"--anioinc_reorg_elite_value={anioinc_reorg_elite_value}", f"--catioinc_reorg_elite_value={catioinc_reorg_elite_value}"])  # does elite ordering, submits reorganisation calculation
    
    subprocess.run(["python", "elite_mutation.py", f"--run_num={run_number}", f"--run_size={run_size}"])  # extracts data and updates a df


#You may want to calculte how many total systems you can have and then have an elite that matches a certain percentage of that?

#run a test loop

for i in range(10, 31):
    single_run(i, functional='b3lyp', basis_set='6-31g*', 
               time=9, number_of_cpus=10, EG_cutoff=3.2, 
               EG_rank_weight=1, planarity_rank_weight=1, SA_rank_weight=0.25,
               elite_value=25, anioinc_reorg_rank_weight=2, cationic_reorg_rank_weight=2, 
               eg_elite_value=2.5, planarity_elite_value=0.82, anioinc_reorg_elite_value=0.300, catioinc_reorg_elite_value=0.300,
               run_size=100)
