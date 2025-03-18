import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
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
from molecule_mutation import *
from calculation_status import *
import re
import itertools
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
import sys
import os
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
from SA_Score import sascorer
import argparse

# This file will first test if the calculations have been ran successfully
# Split into two dictionaries, one with the successful and another with unsuccessful
# The successful calculations will have their planarity, energy gap extracted and their syenthic accessibility score calculated
# The unsuccessful calculations will either be ran again or discarded (to be decided)
# Once all the data has been extracted they will be ordered into an elitism step
# With the top 25% undergo reorganisation calculations

options = {
    'run_num': {'default': 0},
}


# Create a parser for the arguments that can be changed by the user
parser = argparse.ArgumentParser()
for arg, opts in options.items():
    parser.add_argument(f'--{arg}', type=type(opts['default']), default=opts['default'])
args = parser.parse_args()

# Store run_num as a string as it is mainly used for naming/retreaving files
if hasattr(args, 'run_num') and args.run_num:
    run_num_str = str(args.run_num)
else:
    run_num_str = str(options['run_num']['default'])


#dict of molecules and their SMILES
ran_molecules = open_dictionary(f'ran_{run_num_str}_molecules.json')

#turn into a list of tasks that calculation_status function can proccess
task_list = []
for k, v in ran_molecules.items():
    task = lambda: is_file_present(f'{k}/{k}_opt_energy_and_gap.txt')
    task_list.append((k, task()))

#returns the failed and successful calculations, keeps looping until all calculations are done
succesful_dict, failed_dict, attempts = calculations_status(task_list, sleep_time=5)

failed_molecules = {}
for k, v in failed_dict.items():
    failed_molecules[k] = ran_molecules[k]

SA_score_dict = {}
HOMO_dict = {}
LUMO_dict = {}
EG_dict = {}
planarity_dict = {}
succesful_dict_smi = {}
for k, v in succesful_dict.items():

    succesful_dict_smi[k] =  ran_molecules[k]

    m = Chem.MolFromSmiles(ran_molecules[k])
    #Run synethic accessibility score
    sa_score_val = sascorer.calculateScore(m)
    SA_score_dict[k] = sa_score_val

    file_path = f'{k}/{k}_opt_energy_and_gap.txt'

    data = extract_data_from_txt(file_path)
    #energy calcs
    mol_homo = data['homo']
    HOMO_dict[k] = mol_homo

    mol_lumo = data['lumo']
    LUMO_dict[k] = mol_lumo

    mol_EG = data['energy_gap']
    EG_dict[k] = mol_EG
    linker_type = find_linker_type(m)
    #planarity data
    mol_plan = finding_planairty_psi4(k, ran_molecules[k], linker_type, 'opt')
    planarity_dict[k] = mol_plan

run_x_df = pd.DataFrame()

run_x_df.insert(0, 'Name', succesful_dict_smi.keys())
run_x_df.insert(1, 'SMILES', succesful_dict_smi.values())
run_x_df.insert(2, 'HOMO /eV', HOMO_dict.values())
run_x_df.insert(3, 'LUMO /eV', LUMO_dict.values())
run_x_df.insert(4, 'EG /eV', EG_dict.values())
run_x_df.insert(5, 'Planarity', planarity_dict.values())
run_x_df.insert(6, 'SA Score', SA_score_dict.values())


run_x_df.to_csv(f'run_{run_num_str}_data.csv', index=False)
save_dictionary(failed_molecules, f'failed_molecules_run_{run_num_str}.json')

progress_file_path = 'GA_status.txt'

# Open the file in append mode and write some content
with open(progress_file_path, 'a') as file:
    file.write(f'GA_data_extraction_step complete for run {run_num_str}.\n')
