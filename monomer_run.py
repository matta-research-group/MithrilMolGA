import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from CombineMols.CombineMols import CombineMols
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
import argparse

# All monomers have to be CanonSmiles for retrieval from dataframes
# This file takes all the molecules that want to be run
# Reduces them to their monomers
# Checks if the monomers have been ran before
# If they have been ran before, then continue on
# If they haven't been ran before, then run them
# Wait for all the monomers to be ran and then update the monomer df
# Once all done, pass onto the next script

#FILE NAMES ARE ONLY PLACEHOLDERS AT THIS STAGE

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


#load molecules to run
potential_molecules = open_dictionary(f'molecules_to_run_{run_num_str}.json')

#load monomer df
monomer_df = pd.read_csv('monomer_df.csv')

#dict of monomers and their SMILES
monomer_already_run = dict(zip(monomer_df['Name'], monomer_df['SMILES']))

molecules_monomers = {}
for k, v in potential_molecules.items():
    #find the linker between the monomers
    linker_type = find_linker_type(Chem.MolFromSmiles(v))

    #fragment the molecule
    fragments = fragment_molecule(Chem.MolFromSmiles(v), linker_type)

    #One fragment system
    fragments_one_attach = [smi for smi in fragments if smi.count('I') == 1]
    # Remove the I from the fragment
    fragment_one = re.sub(r'\[I\]', '', fragments_one_attach[0])
    fragment_two = re.sub(r'\[I\]', '', fragments_one_attach[1])

    #adds the monomers to the dict, a and b to distinguish but naming doesn't matter
    molecules_monomers[f'{k}_a'] = fragment_one
    molecules_monomers[f'{k}_b'] = fragment_two

molecules_to_run = {}
for k1, v1 in molecules_monomers.items():
    match_found = False
    for k2, v2 in monomer_already_run.items():
        if Chem.CanonSmiles(v1) == Chem.CanonSmiles(v2):
            match_found = True
            break
    if match_found:
        print(f'{k1}: Match found')
    else:
        #gets the last key and value from the ran monomer dictionary
        last_key, last_value = list(monomer_already_run.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        make_num_str = str(make_num)
        #updates a new dictionary of monomers to run
        molecules_to_run[make_num_str] = v1
        #updates the ran dictionary so no overlap occures
        monomer_already_run[make_num_str] = v1

# if the data is missing, we need to run the psi4 calculations
if molecules_to_run is not None:
    for k, v in molecules_to_run.items():
        run_psi4('opt', k, v, time=4, cpus=10, functional, basis_set) #user set parameters

# Run the calculations and wait for the data to come back
#turn into a list of tasks that calculation_status function can proccess
task_list = []
for k, v in molecules_to_run.items():
    task = lambda: is_file_present(f'{k}_opt.txt')
    task_list.append((k, task))

#returns the failed and successful calculations, keeps looping until all calculations are done
succesful_dict, failed_dict, attempts = calculations_status(task_list, sleep_time=15)

#add failed monomers to a new dictionary
failed_monomers = {}
for k, v in failed_dict.items():
    failed_monomers[k] = molecules_to_run[k]

#extract data from the successful monomers
HOMO_dict = {}
LUMO_dict = {}
EG_dict = {}
for k, v in succesful_dict.items():

    data = extract_data_from_txt(file_path)
    #energy calcs
    mol_homo = data['homo']
    HOMO_dict[k] = mol_homo

    mol_lumo = data['lumo']
    LUMO_dict[k] = mol_lumo

    mol_EG = data['energy_gap']
    EG_dict[k] = mol_EG

succesful_dict_CanonSmiles = {}
for k, v in succesful_dict.items():
    succesful_dict_CanonSmiles[k] = Chem.CanonSmiles(v)
#make a new dataframe with the new monomers data
run_monomer_x_df = pd.DatFrame()

run_monomer_x_df.insert(0, 'Name', succesful_dict_CanonSmiles.keys())
run_monomer_x_df.insert(1, 'SMILES', succesful_dict_CanonSmiles.values())
run_monomer_x_df.insert(2, 'HOMO /eV', HOMO_dict.values())
run_monomer_x_df.insert(3, 'LUMO /eV', LUMO_dict.values())
run_monomer_x_df.insert(4, 'EG /eV', EG_dict.values())

#add the new monomer data to the dataframe of existing monomers
df_concat = pd.concat([monomer_df, run_monomer_x_df], ignore_index=True)

#override the old monomer dataframe with the new one with the data
df_concat.to_csv(f'monomer_df.csv', index=False)
save_dictionary(f'failed_monomers, failed_monomers_run_{run_num_str}.json')
    