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

# This script runs the potential molecules
# Its checks if the donor acceptor matching
# If the donor acceptor matching is good then it runs the molecule
# If the donor acceptor matching does not reach the threshold it is added to a failed list

#load molecules to run
potential_molecules = open_dictionary('molecules_to_run.json')

monomer_df = pd.read_csv('monomer_df.csv') #df containing already ran monomers

monomer_smi = dict(zip(monomer_df['Monomer'], monomer_df['SMILES'])) #dict of monomers and their SMILES

ran_molecules = {}
failed_D_A_match = {}
all_molecules_D_A = {}
all_molecules_A_D = {}
donor_smiles = {}
acceptor_smiles = {}
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

    #monomer data for each fragment
    fragment_one_df = monomer_df[monomer_df['SMILES'] == Chem.CanonSmiles(fragment_one)]
    fragment_two_df = monomer_df[monomer_df['SMILES'] == Chem.CanonSmiles(fragment_two)]
    #get values for each fragment
    homo_fragment_one = pd.to_numeric(fragment_one_df['HOMO']).values
    lumo_fragment_one = pd.to_numeric(fragment_one_df['LUMO']).values

    homo_fragment_two = pd.to_numeric(fragment_two_df['HOMO']).values
    lumo_fragment_two = pd.to_numeric(fragment_two_df['LUMO']).values

    #calculate the predicted energy gap
    EG_D_A = np.abs((homo_fragment_one - lumo_fragment_two)[0])
    EG_A_D = np.abs((homo_fragment_two - lumo_fragment_one)[0])

    all_molecules_D_A[k] = EG_D_A
    all_molecules_A_D[k] = EG_A_D
    donor_smiles[k] = fragment_one
    acceptor_smiles[k] = fragment_two

    if (EG_D_A <= set_EG_value) or (EG_A_D <= set_EG_value):
        # run molecule_study
        #Run Psi4 calculations; planarity and energy gap
        ran_molecules[k] = v
        run_psi4('opt', mol_name, molecule_study, time=4, cpus=10, functional, basis_set) #user set parameters
    else:
        failed_D_A_match[k] = v

d_a_matching_df = pd.DataFrame()

d_a_matching_df.insert(0, 'Name', potential_molecules.keys())
d_a_matching_df.insert(1, 'SMILES', potential_molecules.values())
d_a_matching_df.insert(2, 'Donor', donor_smiles.values())
d_a_matching_df.insert(3, 'Acceptor', acceptor_smiles.values())
d_a_matching_df.insert(4, 'D_A EG /eV', all_molecules_D_A.values())
d_a_matching_df.insert(5, 'A_D EG /eV', all_molecules_A_D.values())

d_a_df = pd.read_csv('d_a_df.csv')

d_a_df_concat = pd.concat([d_a_df, d_a_matching_df], ignore_index=True)

d_a_df_concat.to_csv('d_a_df.csv', index=False)

save_dictionary(ran_molecules, f'ran_{x}_molecules.json')
save_dictionary(failed_D_A_match, f'failed_D_A_match_{x}_molecules.json')