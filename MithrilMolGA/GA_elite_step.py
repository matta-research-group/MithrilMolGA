import os
import argparse
import pandas as pd
from datetime import datetime
from QCflow.load_gaussian import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *
from QCflow.write_psi4 import *
from QCflow.run_psi4 import *
from QCflow.energy_calculations import *
from MithrilMolGA.molecule_mutation import *
from MithrilMolGA.calculation_status import *

# ==============================
# Default variables
# ==============================
options = {
    'run_num': {'default': 0},
    'EG_rank_weight': {'default': 1.0},
    'planarity_rank_weight': {'default': 1.0},
    'SA_rank_weight': {'default': 0.25},
    'elite_value': {'default': 25},
    'anioinc_reorg_rank_weight': {'default': 2.0},
    'cationic_reorg_rank_weight': {'default': 2.0},
    'functional': {'default': 'b3lyp'},
    'basis_set': {'default': '6-31g*'},
    'time': {'default': 4},
    'cpus': {'default': 10},
    'eg_elite_value': {'default': 2.5},
    'planarity_elite_value': {'default': 0.82},
    'anioinc_reorg_elite_value': {'default': 0.350},
    'catioinc_reorg_elite_value': {'default': 0.350}
}

# ==============================
# Argument parser
# ==============================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for arg, opts in options.items():
        parser.add_argument(f'--{arg}', type=type(opts['default']), default=opts['default'])
    
    args, unknown = parser.parse_known_args()  # ignore unknown args passed from master script

    # ==============================
    # Store variables
    # ==============================
    run_num_str = str(args.run_num) if hasattr(args, 'run_num') and args.run_num else str(options['run_num']['default'])
    EG_rank_weight = args.EG_rank_weight
    planarity_rank_weight = args.planarity_rank_weight
    SA_rank_weight = args.SA_rank_weight
    elite_value = args.elite_value
    anioinc_reorg_rank_weight = args.anioinc_reorg_rank_weight
    cationic_reorg_rank_weight = args.cationic_reorg_rank_weight
    functional = args.functional
    basis_set = args.basis_set
    time = args.time
    cpus = args.cpus
    eg_elite_value = args.eg_elite_value
    planarity_elite_value = args.planarity_elite_value
    anioinc_reorg_elite_value = args.anioinc_reorg_elite_value
    catioinc_reorg_elite_value = args.catioinc_reorg_elite_value

    # ==============================
    # Load molecule dataframe
    # ==============================
    molecule_path = f'dataframes/run_{run_num_str}_data.csv'
    if os.path.exists(molecule_path):
        molecule_df = pd.read_csv(molecule_path)
    else:
        raise FileNotFoundError(f"{molecule_path} not found. Cannot proceed without molecule data.")

    molecule_df['Name'] = molecule_df['Name'].astype(str)

    # ==============================
    # Rank EG, Planarity, SA
    # ==============================
    molecule_df['EG Rank Order'] = molecule_df['EG /eV'].rank(ascending=True)
    molecule_df['Plan Rank Order'] = molecule_df['Planarity'].rank(ascending=False)
    molecule_df['SA Rank Order'] = molecule_df['SA Score'].rank(ascending=True)
    molecule_df['Rank Sum'] = (
        molecule_df['EG Rank Order'] * EG_rank_weight +
        molecule_df['Plan Rank Order'] * planarity_rank_weight +
        molecule_df['SA Rank Order'] * SA_rank_weight
    )

    # ==============================
    # Elite selection
    # ==============================
    sorted_molecule_df = molecule_df.sort_values(['Rank Sum'], ascending=True)
    elite_df = sorted_molecule_df.head(int(len(sorted_molecule_df)*(elite_value/100)))
    elite_df = elite_df[(elite_df['EG /eV'] <= eg_elite_value) & (elite_df['Planarity'] >= planarity_elite_value)]
    elite_smi_int = dict(zip(elite_df['Name'], elite_df['SMILES']))

    # ==============================
    # Archive handling
    # ==============================
    archive_path = 'archive_dataframes/all_ran_reorg.csv'
    if os.path.exists(archive_path):
        archive_df = pd.read_csv(archive_path)
    else:
        print(f"Warning: {archive_path} not found. Using empty DataFrame.")
        archive_df = pd.DataFrame()

    archive_elite_results = archive_df.iloc[0:0].copy()
    elite_smi = {}
    elite_archive_smi = {}

    for k, v in elite_smi_int.items():
        canon_smi = v
        if len(archive_df[archive_df['SMILES'] == canon_smi]) > 0:
            elite_archive_smi[k] = canon_smi
            ran_before = archive_df[archive_df['SMILES'] == canon_smi].copy()
            ran_before['Name'] = int(k)
            archive_elite_results = pd.concat([archive_elite_results, ran_before])
            print(k, 'Elite has been run before, data retrieved')
        else:
            elite_smi[str(k)] = v
            print(k, 'Elite has NOT been run before')

    save_dictionary(elite_archive_smi, f'run_dic/elite_archive_smi_{run_num_str}_molecules.json')

    if len(elite_smi) <= 0:
        elite_smi_int = dict(zip(sorted_molecule_df['Name'], sorted_molecule_df['SMILES']))
        length_of_elite = len(elite_df)
        subset_subset = dict(list(elite_smi_int.items())[(length_of_elite+1):length_of_elite+6])
        elite_smi = {str(int(k)): v for k, v in subset_subset.items()}
        elite_df = sorted_molecule_df

    print('Submit Calculations')

    # ==============================
    # Psi4 calculations
    # ==============================
    os.makedirs('data', exist_ok=True)
    os.chdir('data')
    for k, v in elite_smi.items():
        k = str(int(k))
        run_psi4('anion', k, v, time, cpus, functional, basis_set)
        run_psi4('cation', k, v, time, cpus, functional, basis_set)
        run_psi4('sp_c', k, v, time, cpus, functional, basis_set)
        run_psi4('sp_a', k, v, time, cpus, functional, basis_set)
    os.chdir('../')

    # ==============================
    # Update GA status
    # ==============================
    progress_file_path = 'GA_status.txt'
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(progress_file_path, 'a') as file:
        file.write(f'GA_elite_step complete for run {run_num_str} at {current_time}.\n')

    print('GA_elite_step complete')
