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
from MithrilMolGA.molecule_mutation import *
from MithrilMolGA.calculation_status import *
import re
import itertools
import argparse
import os
from datetime import datetime
import shutil

# This script runs the potential molecules
# Its checks if the donor acceptor matching
# If the donor acceptor matching is good then it runs the molecule
# If the donor acceptor matching does not reach the threshold it is added to a failed list

options = {
    'run_num': {'default': 0},
    'set_EG_value': {'default': 3.2},
    'functional' : {'default': 'b3lyp'},
    'basis_set' : {'default': '6-31g*'},
    'time' : {'default': 6},
    'cpus' : {'default': 10}
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

# Store other variables and allow for them to be over written by the user
set_EG_value = args.set_EG_value if hasattr(args, 'set_EG_value') else options['set_EG_value']['default']
functional = args.functional if hasattr(args, 'functional') else options['functional']['default']
basis_set = args.basis_set if hasattr(args, 'basis_set') else options['basis_set']['default']
time = args.time if hasattr(args, 'time') else options['time']['default']
cpus = args.cpus if hasattr(args, 'cpus') else options['cpus']['default']


#load molecules to run
potential_molecules = open_dictionary(f'submission_dic/molecules_to_run_{run_num_str}.json')

total_molecules_ran = open_dictionary('ga_dic/total_molecules_ran.json')

total_molecules = total_molecules_ran | potential_molecules

save_dictionary(total_molecules, 'ga_dic/total_molecules_ran.json')

monomer_df = pd.read_csv('dataframes/monomer_df.csv') #df containing already ran monomers

monomer_smi = dict(zip(monomer_df['Name'], monomer_df['SMILES'])) #dict of monomers and their SMILES

#Load in the ran dataframe so can use the archived data
all_data = pd.read_csv('archive_dataframes/ran_all_data.csv')
#create a empty df with the same columns as all_data
archive_results = all_data.iloc[0:0].copy()

#change into data folder
os.chdir('data')

ran_molecules = {}
archive_ran_molecules = {}
failed_D_A_match = {}
all_molecules_D_A = {}
all_molecules_A_D = {}
donor_smiles = {}
acceptor_smiles = {}
for k, v in potential_molecules.items():
    k = str(int(k))  # Ensure k is a string for naming consistency
    #find the linker between the monomers
    linker_type = find_linker_type(Chem.MolFromSmiles(v))

    #fragment the molecule
    fragments = fragment_molecule(Chem.MolFromSmiles(v), linker_type)

    #One fragment system
    fragments_one_attach = [smi for smi in fragments if smi.count('I') == 1]
    # Remove the I from the fragment
    fragment_one = re.sub(r'\[I\]', '', fragments_one_attach[0])
    fragment_two = re.sub(r'\[I\]', '', fragments_one_attach[1])

    #monomer data for each fragment
    fragment_one_df = monomer_df[monomer_df['SMILES'] == Chem.CanonSmiles(fragment_one, useChiral=0)]
    fragment_two_df = monomer_df[monomer_df['SMILES'] == Chem.CanonSmiles(fragment_two, useChiral=0)]
    #get values for each fragment
    homo_fragment_one = pd.to_numeric(fragment_one_df['HOMO /eV']).values
    lumo_fragment_one = pd.to_numeric(fragment_one_df['LUMO /eV']).values

    homo_fragment_two = pd.to_numeric(fragment_two_df['HOMO /eV']).values
    lumo_fragment_two = pd.to_numeric(fragment_two_df['LUMO /eV']).values

    #calculate the predicted energy gap
    EG_D_A = np.abs((homo_fragment_one - lumo_fragment_two)[0])
    EG_A_D = np.abs((homo_fragment_two - lumo_fragment_one)[0])

    all_molecules_D_A[k] = EG_D_A
    all_molecules_A_D[k] = EG_A_D
    donor_smiles[k] = fragment_one
    acceptor_smiles[k] = fragment_two

    if (EG_D_A <= set_EG_value) or (EG_A_D <= set_EG_value):
        # run molecule_study

        #canon_smi = Chem.CanonSmiles(v, useChiral=0)
        canon_smi = v
        #check if its been ran before in the archive
        if len(all_data[all_data['SMILES'] == canon_smi]) > 0:
            archive_ran_molecules[k] = canon_smi
            #copy over from archive
            ran_before = all_data[all_data['SMILES'] == canon_smi].copy()
            #make sure the number is up to date

            name_value = all_data.loc[all_data['SMILES'] == canon_smi, 'Name'].values
            #get number
            name_str = str(name_value[0])

            if name_value[0] <= 357:
                #old data folder
                print('Located in run 9')
                src_dir = f'/scratch/prj/ch_mime/GA_pratice_runs/run_9/MithrilMolGA/MithrilMolGA/data/{name_str}'
            elif isinstance(name_value[0], str) and re.match(r'^\d+_\d+_$', name_value[0]):
                print('Located in funnel runs')
                src_dir = f'/scratch/prj/ch_mime/funnel_project_run_all/{name_str}'
            else:
                #old data folder
                print('Located in run 11')
                src_dir = f'/scratch/prj/ch_mime/GA_pratice_runs/run_11/MithrilMolGA/MithrilMolGA/data/{name_str}'
            #new data folder
            backup_dir = f'/scratch/prj/ch_mime/GA_new_bio_16_06_25/archive_testing_attempt_7/MithrilMolGA/MithrilMolGA/data/{k}'
            #change the name of the folder to the new number
            if os.path.exists(src_dir):
                shutil.copytree(src_dir, backup_dir)
                print(f'Copied {src_dir} to {backup_dir}')
                for filename in os.listdir(backup_dir):
                    if re.match(rf'^{name_str}\D', filename):
                        new_filename = re.sub(rf'^{name_str}', k, filename)
                        src = os.path.join(backup_dir, filename)
                        dst = os.path.join(backup_dir, new_filename)
                        os.rename(src, dst)
                        print(f'Renamed {filename} to {new_filename}')
            else:
                print(f'Source directory {src_dir} does not exist, skipping backup.')

            ran_before['Name'] = int(k)
            archive_results = pd.concat([archive_results, ran_before])
            print(k, 'Has been ran before, data retreived')
        else:
            #Run Psi4 calculations; planarity and energy gap
            ran_molecules[k] = v
            print(k, 'Has not been ran before, run the psi4 calculation') 
            run_psi4('opt', k, v, time, cpus, functional, basis_set) #user set parameters

    else:
        failed_D_A_match[k] = v

potential_molecules = {str(int(k)): v for k, v in potential_molecules.items()}

#back out of data folder
os.chdir('../')

d_a_matching_df = pd.DataFrame()

d_a_matching_df.insert(0, 'Name', potential_molecules.keys())
d_a_matching_df.insert(1, 'SMILES', potential_molecules.values())
d_a_matching_df.insert(2, 'Donor', donor_smiles.values())
d_a_matching_df.insert(3, 'Acceptor', acceptor_smiles.values())
d_a_matching_df.insert(4, 'D_A EG /eV', all_molecules_D_A.values())
d_a_matching_df.insert(5, 'A_D EG /eV', all_molecules_A_D.values())

d_a_df = pd.read_csv('dataframes/d_a_df.csv')

d_a_df_concat = pd.concat([d_a_df, d_a_matching_df], ignore_index=True)

d_a_df_concat.to_csv('dataframes/d_a_df.csv', index=False)

#save the ran molecules df
archive_results.to_csv(f'archive_dataframes/archive_run_{run_num_str}_data.csv', index=False)

save_dictionary(ran_molecules, f'run_dic/ran_{run_num_str}_molecules.json')
save_dictionary(archive_ran_molecules, f'run_dic/archive_ran_{run_num_str}_molecules.json')
save_dictionary(failed_D_A_match, f'failed_dic/failed_D_A_match_{run_num_str}_molecules.json')

progress_file_path = 'GA_status.txt'
current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

# Open the file in append mode and write some content
with open(progress_file_path, 'a') as file:
    file.write(f'molecule_run complete for run {run_num_str} at {current_time}.\n')
