import subprocess
import argparse

def single_run(run_number, functional='b3lyp', basis_set='6-31g*', 
               time=6, number_of_cpus=10, EG_cutoff=3.2, 
               EG_rank_weight=1, planarity_rank_weight=1, SA_rank_weight=4,
               elite_value=25, anioinc_reorg_rank_weight=0.5, cationic_reorg_rank_weight=0.5, elite_df_size=10):
    """
    Runs a single iteration of the GA workflow.
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
                    f"--elite_df_size={elite_df_size}"])  # does elite ordering, submits reorganisation calculation
    
    subprocess.run(["python", "elite_mutation.py", f"--run_num={run_number}"])  # extracts data and updates a df


#run a test loop

for i in range(4, 9):
    single_run(i, functional='b3lyp', basis_set='6-31g*', 
               time=5, number_of_cpus=10, EG_cutoff=3.2, 
               EG_rank_weight=1, planarity_rank_weight=1, SA_rank_weight=4,
               elite_value=50, anioinc_reorg_rank_weight=0.5, cationic_reorg_rank_weight=0.5, elite_df_size=20)