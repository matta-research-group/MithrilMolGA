import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
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

# This file will first test if the calculations have been ran successfully
# Split into two dictionaries, one with the successful and another with unsuccessful
# The successful calculations will have their planarity, energy gap extracted and their syenthic accessibility score calculated
# The unsuccessful calculations will either be ran again or discarded (to be decided)
# Once all the data has been extracted they will be ordered into an elitism step
# With the top 25% undergo reorganisation calculations

#load molecule df
molecule_df = pd.read_csv('fake_molecules_df.csv')

#dict of molecules and their SMILES
ran_molecules = dict(zip(molecule_df['Name'], molecule_df['SMILES']))

#turn into a list of tasks that calculation_status function can proccess
task_list = []
for k, v in ran_molecules.items():
    task = lambda: is_file_present(f'{k}_opt.txt')
    task_list.append((k, task))

#returns the failed and successful calculations, keeps looping until all calculations are done
succesful_dict, failed_dict, attempts = calculations_status(task_list, sleep_time=15)

failed_molecules = {}
for k, v in failed_dict.items():
    failed_molecules[k] = ran_molecules[k]

SA_score_dict = {}
HOMO_dict = {}
LUMO_dict = {}
EG_dict = {}
planarity_dict = {}
for k, v in succesful_dict.items():

    m = Chem.MolFromSmiles(molecule_study)
    #Run synethic accessibility score
    sa_score_val = sascorer.calculateScore(m)
    SA_score_dict[k] = sa_score_val

    data = extract_data_from_txt(file_path)
    #energy calcs
    mol_homo = data['homo']
    HOMO_dict[k] = mol_homo

    mol_lumo = data['lumo']
    LUMO_dict[k] = mol_lumo

    mol_EG = data['energy_gap']
    EG_dict[k] = mol_EG
    #planarity data
    mol_plan = finding_planairty_psi4(mol_name, molecule_study, linker_type, 'opt')
    planarity_dict[k] = mol_plan

run_x_df = pd.DatFrame()

run_x_df.insert(0, 'Name', succesful_dict.keys())
run_x_df.insert(1, 'SMILES', succesful_dict.values())
run_x_df.insert(2, 'HOMO /eV', HOMO_dict.values())
run_x_df.insert(3, 'LUMO /eV', LUMO_dict.values())
run_x_df.insert(4, 'EG /eV', EG_dict.values())
run_x_df.insert(5, 'Planarity', planarity_dict.values())
run_x_df.insert(6, 'SA Score', SA_score_dict.values())


run_x_df.to_csv(f'run_{X}_data.csv', index=False)
save_dictionary(f'failed_molecules, failed_molecules_run_{X}.json')