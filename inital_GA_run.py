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
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
import sys
import os
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
from SA_Score import sascorer


#### Before this the pairs will have been checked to see if they already exist as a pair!!!!
#### This needs to be tested!!!!
#### This will be the part of the first box in the workflow (isolated fragment study)
#### This will be implanted within a input of a wider snakemake workflow

set_EG_value = 2.5 #eV of the cut off for donor-acceptor matching

mol_name = '0' #will be a number

molecule_study = 'COC1=C(C#CC2=CC(C(NC3=O)=O)=C(C4=C3C=CS4)S2)C5=C(C=C1OC)NC(C(O)=O)=C5'

monomer_df = pd.read_csv('fake_data_test/fake_data.csv') #df containing already ran monomers

monomer_smi = dict(zip(monomer_df['Monomer'], monomer_df['SMILES'])) #dict of monomers and their SMILES

# have these monomers been ran before so we can just extract that data?

#find the linker between the monomers
linker_type = find_linker_type(Chem.MolFromSmiles(molecule_study))

#fragment the molecule
fragmented = fragment_molecule(Chem.MolFromSmiles(molecule_study), linker_type)

#One fragment system
fragments_one_attach = [smi for smi in fragments if smi.count('I') == 1]
    # Remove the I from the fragment
fragment_one = re.sub(r'\[I\]', '', fragments_one_attach[0])
fragment_two = re.sub(r'\[I\]', '', fragments_one_attach[1])

data_exists = {}
data_missing = {}
updated_monomer_smi = {}
for k, v in monomer_smi.items():
    if Chem.CanonSmiles(v) == fragment_one:
        data_exists[k] = v
    else:
        last_key = int(list(monomer_smi.keys())[-1]) + 1 #adds to a missing dictionary and is 1 number greater than the last key
        data_missing[f'{last_key}'] = fragment_one
        updated_monomer_smi[f'{last_key}'] = fragment_one


for k, v in updated_monomer_smi.items():
    if Chem.CanonSmiles(v) == fragment_two:
        data_exists[k] = v
    else:
        last_key = int(list(updated_monomer_smi.keys())[-1]) + 1 #adds to a missing dictionary and is 1 number greater than the last key
        data_missing[f'{last_key}'] = fragment_two 
        updated_monomer_smi[f'{last_key}'] = fragment_two

# if the data is missing, we need to run the psi4 calculations

if data_missing is not None:
    for k, v in data_missing.items():
        run_psi4('opt', k, v, time=4, cpus=10, functional, basis_set) #user set parameters

# Run the calculations and wait for the data to come back
psi_homo = {}
psi_lumo = {}
psi_EG = {}
if data_missing is not None:
    for k, v in data_missing.items():
        file_path = f'{k}/{k}_opt_energy_and_gap.txt'
        wait_for_file(file_path, sleep_time=110, timeout=10)
        data = extract_data_from_txt(file_path)
        psi_homo = data['homo']
        psi_lumo = data['lumo']
        psi_EG = data['energy_gap']

if data_missing is not None:
    new_monomer_data = {
        'Monomer': list(data_missing.keys()),
        'SMILES': list(data_missing.values()),
        'HOMO': list(psi_homo.values()),
        'LUMO': list(psi_lumo.values()),
        'EG': list(psi_EG.values())
    }
    new_monomer_data_df = pd.DataFrame(new_monomer_data)
    monomer_df = pd.concat([monomer_df, new_monomer_data_df], ignore_index=True)
    monomer_df.to_csv('fake_data_test/fake_data.csv') #save the dataframe and overide the old one

data_finish = data_exists | data_missing

#Get data from the data frames
data_exists_homo = {}
data_exists_lumo = {}
for k, v in data_finish.items():
    specifc_monomer = monomer_df[monomer_df['Monomer'] == k]
    data_exists_homo[k] = float(specifc_monomer['HOMO'].values[0])
    data_exists_lumo[k] = float(specifc_monomer['LUMO'].values[0])



BG_D_A = {}
BG_A_D = {}
Donor = {}
Acceptor = {}
for (k1, v1), (k2, v2) in itertools.combinations(data_finish.items(), 2):
    if k1 != k2:
        homo_k1 = data_exists_homo[k1]
        homo_k2 = data_exists_homo[k2]
        lumo_k1 = data_exists_lumo[k1]
        lumo_k2 = data_exists_lumo[k2]
        
        BG_k1_k2 = homo_k1 - lumo_k2
        BG_k2_k1 = homo_k2 - lumo_k1

        BG_D_A[k1] = np.abs(BG_k1_k2)
        BG_A_D[k2] = np.abs(BG_k2_k1)

        Donor[k1] = v1
        Acceptor[k2] = v2

d_a_matching_GA = pd.read_csv('fake_data_test/d_a_matching_GA.csv')

# Create a new dataframe from the fake data dictionaries
new_energy_gap = {
    'Monomer_A': list(Donor.keys()),
    'Monomer_A_Smiles': list(Donor.values()),
    'Monomer_B': list(Acceptor.keys()),
    'Monomer_B_Smiles': list(Acceptor.values()),
    'D_A': list(BG_D_A.values()),
    'A_D': list(BG_A_D.values())
}

new_pair_df = pd.DataFrame(new_energy_gap)

# Concatenate the new dataframe with the existing df
d_a_matching_GA = pd.concat([d_a_matching_GA, new_pair_df], ignore_index=True)

# Save the updated dataframe
d_a_matching_GA.to_csv('fake_data_test/d_a_matching_GA.csv')

#set value is the value that the EG needs to be below to be a good pair
if (new_energy_gap['D_A'] <= set_EG_value) or (new_energy_gap['A_D'] <= set_EG_value):
    # run molecule_study
    m = Chem.MolFromSmiles(molecule_study)
    #Run synethic accessibility score
    sa_score_val = sascorer.calculateScore(m)
    #Run Psi4 calculations; planarity and energy gap
    run_psi4('opt', mol_name, molecule_study, time=4, cpus=10, functional, basis_set) #user set parameters
    #wait for file
    wait_for_file(file_path, sleep_time=110, timeout=10)
    #retreave data
    data = extract_data_from_txt(file_path)
    #energy calcs
    mol_homo = data['homo']
    mol_lumo = data['lumo']
    mol_EG = data['energy_gap']
    #planarity data
    mol_plan = finding_planairty_psi4(mol_name, molecule_study, linker_type, 'opt')
    #Adds all data to a df of ran systems
    #adds dtata to dataframe
else:
    #system is a bad match and is not ran but added to the overall dataset