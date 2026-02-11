import rdkit
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdDepictor
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
import pandas as pd
import argparse
import QCflow
from QCflow.load_gaussian import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *
from QCflow.write_psi4 import *
from QCflow.run_psi4 import *
from QCflow.energy_calculations import *
from MithrilMolGA.molecule_mutation import *
from MithrilMolGA.calculation_status import *
import re
import itertools
import sys
from datetime import datetime
import os

# Ensure RDKit prefers CoordGen
rdDepictor.SetPreferCoordGen(True)

# Add RDKit contrib to path
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
from SA_Score import sascorer

# ==============================
# This file will first test if the calculations have been ran successfully
# Split into two dictionaries, one with the successful and another with unsuccessful
# The successful calculations will have their planarity, energy gap extracted and their synthetic accessibility score calculated
# The unsuccessful calculations will either be ran again or discarded (to be decided)
# Once all the data has been extracted they will be ordered into an elitism step
# With the top 25% undergo reorganisation calculations
# ==============================

options = {
    'run_num': {'default': 0},
}

# Create a parser for the arguments that can be changed by the user
parser = argparse.ArgumentParser()
for arg, opts in options.items():
    parser.add_argument(f'--{arg}', type=type(opts['default']), default=opts['default'])
args = parser.parse_args()

# Store run_num as a string for naming/retreaving files
run_num_str = str(getattr(args, 'run_num', options['run_num']['default']))

# Ensure directories exist
os.makedirs('run_dic', exist_ok=True)
os.makedirs('archive_dataframes', exist_ok=True)
os.makedirs('dataframes', exist_ok=True)
os.makedirs('failed_dic', exist_ok=True)

# ==============================
# Load dictionaries
# ==============================
ran_molecules_path = f'run_dic/ran_{run_num_str}_molecules.json'
archive_ran_molecules_path = f'run_dic/archive_ran_{run_num_str}_molecules.json'

ran_molecules = open_dictionary(ran_molecules_path) if os.path.exists(ran_molecules_path) else {}
archive_ran_molecules = open_dictionary(archive_ran_molecules_path) if os.path.exists(archive_ran_molecules_path) else {}


# Move into the data folder
os.makedirs('data', exist_ok=True)
os.chdir('data')

# ==============================
# Turn into a list of tasks for calculation_status
# ==============================
task_list = []
if ran_molecules:
    for k, v in ran_molecules.items():
        task = lambda k=k: is_file_present(f'{k}/{k}_opt_energy_and_gap.txt')
        task_name = f'{k}_opt'
        task_list.append((task_name, task()))

# ==============================
# Returns the failed and successful calculations
# ==============================
succesful_dict, failed_dict, attempts = calculations_status(task_list, sleep_time=5)

# ==============================
# Map failed molecules
# ==============================
failed_molecules = {}
if ran_molecules:
    for k in failed_dict:
        mol_key = k.split('_')[0]  # just the molecule number
        failed_molecules[mol_key] = ran_molecules[mol_key]

# ==============================
# Extract properties for successful molecules
# ==============================
SA_score_dict = {}
HOMO_dict = {}
LUMO_dict = {}
EG_dict = {}
planarity_dict = {}
succesful_dict_smi = {}

if ran_molecules:
    for k in succesful_dict:
        mol_key = str(int(k.split('_')[0]))
        mol_smi = ran_molecules[mol_key]
        succesful_dict_smi[mol_key] = mol_smi

        m = Chem.MolFromSmiles(mol_smi)
        # Run synthetic accessibility score
        SA_score_dict[mol_key] = sascorer.calculateScore(m)

        file_path = f'{mol_key}/{mol_key}_opt_energy_and_gap.txt'
        data = extract_data_from_txt(file_path) if os.path.exists(file_path) else {'homo': None, 'lumo': None, 'energy_gap': None}

        HOMO_dict[mol_key] = data['homo']
        LUMO_dict[mol_key] = data['lumo']
        EG_dict[mol_key] = data['energy_gap']

        linker_type = find_linker_type(m)
        # Planarity data
        planarity_dict[mol_key] = float(finding_planairty_psi4(mol_key, mol_smi, linker_type, 'opt'))

#leave data folder
os.chdir('../')

# ==============================
# Create DataFrame
# ==============================
run_x_df = pd.DataFrame({
    'Name': succesful_dict_smi.keys(),
    'SMILES': succesful_dict_smi.values(),
    'HOMO /eV': HOMO_dict.values(),
    'LUMO /eV': LUMO_dict.values(),
    'EG /eV': EG_dict.values(),
    'Planarity': planarity_dict.values(),
    'SA Score': SA_score_dict.values()
})

archive_results_path = f'archive_dataframes/archive_run_{run_num_str}_data.csv'
if os.path.exists(archive_results_path):
    archive_results = pd.read_csv(archive_results_path)
else:
    archive_results = pd.DataFrame()

run_x_final_df = pd.concat([run_x_df, archive_results], ignore_index=True)

# ==============================
# Save outputs
# ==============================
run_x_final_df.to_csv(f'dataframes/run_{run_num_str}_data.csv', index=False)
save_dictionary(failed_molecules, f'failed_dic/failed_molecules_run_{run_num_str}.json')

# ==============================
# Log progress
# ==============================
progress_file_path = 'GA_status.txt'
current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
with open(progress_file_path, 'a') as file:
    file.write(f'GA_data_extraction_step complete for run {run_num_str} at {current_time}.\n')

