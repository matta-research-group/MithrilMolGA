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
from workflow_initial import *
import re
import itertools


#### Before this the pairs will have been checked to see if they already exist as a pair!!!!
#### This needs to be tested!!!!
#### This will be the part of the first box in the workflow (isolated fragment study)
#### This will be implanted within a input of a wider snakemake workflow

bio_inspired_smi = 'COC1=CC2=C(C=C1OC)NC(C(O)=O)=C2' #DHICA

non_bio_inspired_smi = 'O=C(C1=C(C2=C3C=CS2)SC=C1)NC3=O' #None bio

bio_inspired_smi_CANON = Chem.CanonSmiles(bio_inspired_smi) #canonical SMILES

non_bio_inspired_molecule_CANON = Chem.CanonSmiles(non_bio_inspired_smi) #canonical SMILES

monomer_df = pd.read_csv('fake_data_test/fake_data.csv') #df containing already ran monomers

monomer_smi = dict(zip(monomer_df['Monomer'], monomer_df['SMILES'])) #dict of monomers and their SMILES

# have these monomers been ran before so we can just extract that data?
data_exists = {}
data_missing = {}
updated_monomer_smi = {}
for k, v in monomer_smi.items():
    if Chem.CanonSmiles(v) == bio_inspired_smi_CANON:
        data_exists[k] = v
    else:
        last_key = int(list(monomer_smi.keys())[-1])
        data_missing[f'{last_key}'] = bio_inspired_smi
        updated_monomer_smi[f'{last_key}'] = bio_inspired_smi


for k, v in updated_monomer_smi.items():
    if Chem.CanonSmiles(v) == non_bio_inspired_molecule_CANON:
        data_exists[k] = v
    else:
        last_key = int(list(updated_monomer_smi.keys())[-1])
        data_missing[f'{last_key}'] = non_bio_inspired_smi 
        updated_monomer_smi[f'{last_key}'] = non_bio_inspired_smi

# if the data is missing, we need to run the psi4 calculations

if data_missing is not None:
    for k, v in data_missing.items():
        run_psi4('opt', k, v, time=4, cpus=10, functional='b3lyp', basis_set='6-31g*')

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