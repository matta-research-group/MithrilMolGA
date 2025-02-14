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
import random

#Mutate the elite molecules
#Includes new fragments, new linkers and even a whole new molecules
#User inputted balance for what type of mutations occure and how many new molecules are wanted
#Idea mutation weighting:
# 25% are new molecules, 25% have 1 biofragment fragment change, 25% non-biofragment change, 25% linker change

#load the elite df
elite_df = pd.read_csv(f'elite_25_run_{run_number}_df.csv')
#get th smi strings of the eilte_25
elite_smi = dict(zip(elite_25['Name'], elite_25['SMILES']))
#load all the molecules that have ever been ran in the GA
all_ran_df = pd.read_csv(f'ran_all_data.csv')
#get th nams and smiles so that the script does not produce duplicate molecules
all_ran_smi = dict(zip(all_ran_df['Name'], all_ran_df['SMILES']))

#Will have to determine which is the bioinspired fragment
#If more than one and the just randomly pick between the two to replace with a new bioinspired fragment
#Will have to determine which is the linker
#Have 4 possible choices and then pick at random between those choices
#Carry out the choice and then check if it has been done before, if so do again until new one found

new_study_molecules = {}
for k, v in elite_smi.items():
    #options to mutate
    choices = ['new_mol', 'new_bio', 'new_non_bio', 'new_linker']
    #which random one is chosen for this molecule
    selected_choice = random.choices(choices, k=1)[0]

    if selected_choice == 'new_mol':
        #do this
        #make a new molecule that does not have either fragment in it
    if selected_choice == 'new_bio':
        #do this
        #find the biofragment and replace it with a new biofragment
    if selected_choice == 'new_non_bio':
        #do this
        #find the non-biofragment and replace it with a new non-biofragment
    if selected_choice == 'new_linker':
        #do this
        #find the linker and replace it with a new linker that is not the same as old one