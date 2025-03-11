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

all_ran_smi_canon = {k: Chem.CanonSmiles(v) for k, v in all_ran_smi.items()}

#linker_dic
linker_dic = open_dictionary('linker_dic.json')

new_study_molecules = {}
for k, v in elite_smi.items():
    #options to mutate
    choices = ['new_mol', 'new_bio', 'new_non_bio', 'new_linker']
    #which random one is chosen for this molecule
    selected_choice = random.choices(choices, k=1)[0]

    if selected_choice == 'new_mol':
        #do this
        #make a new molecule that does not have either fragment in it

        #swap the bio fragment
        new_molecule_new_bio = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'bio')
        while Chem.CanonSmiles(new_molecule_new_bio) in all_ran_smi_canon.values():
            new_molecule = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'bio')

        #swap the non-bio fragment
        new_molecule_new_non_bio_and_bio = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'non_bio')
        while Chem.CanonSmiles(new_molecule_new_non_bio_and_bio) in all_ran_smi_canon.values():
            new_molecule_new_non_bio_and_bio = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'non_bio')

        #swap the linker
        linker_type = find_linker_type(Chem.MolFromSmiles(new_molecule_new_non_bio_and_bio))
        #fragment the molecule
        fragments = fragment_molecule(Chem.MolFromSmiles(new_molecule_new_non_bio_and_bio), linker_type)
        #replace the linker
        replaced_linker = replace_linker(fragments, linker_dic)

        last_key, last_value = list(all_ran_smi_canon.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        #creates the new name for the molecule
        make_num_str = str(make_num)
        #updates the new molecules list
        new_study_molecules[make_num_str] = replaced_linker
        #updates the ran dictionary so no overlap occures
        all_ran_smi[make_num_str] = replaced_linker

        
    if selected_choice == 'new_bio':
        #find the biofragment and replace it with a new biofragment
        new_molecule = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'bio')
        while Chem.CanonSmiles(new_molecule) in all_ran_smi_canon.values():
            new_molecule = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'bio')
        
        last_key, last_value = list(all_ran_smi_canon.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        #creates the new name for the molecule
        make_num_str = str(make_num)

        #updates the new molecules list
        new_study_molecules[make_num_str] = new_molecule
        #updates the ran dictionary so no overlap occures
        all_ran_smi[make_num_str] = new_molecule

    if selected_choice == 'new_non_bio':
        #do this
        #find the non-biofragment and replace it with a new non-biofragment
        new_molecule = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'non_bio')
        while Chem.CanonSmiles(new_molecule) in all_ran_smi_canon.values():
            new_molecule = swap_one_fragment(mol_smi, bio_dic, non_bio_dic, 'non_bio')
        
        last_key, last_value = list(all_ran_smi_canon.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        #creates the new name for the molecule
        make_num_str = str(make_num)
        #updates the new molecules list
        new_study_molecules[make_num_str] = new_molecule
        #updates the ran dictionary so no overlap occures
        all_ran_smi[make_num_str] = new_molecule

    if selected_choice == 'new_linker':
        #do this
        #find the linker and replace it with a new linker that is not the same as old one
        #swap the linker
        linker_type = find_linker_type(Chem.MolFromSmiles(v))
        #fragment the molecule
        fragments = fragment_molecule(Chem.MolFromSmiles(v), linker_type)
        #replace the linker
        replaced_linker = replace_linker(fragments, linker_dic)

        last_key, last_value = list(all_ran_smi_canon.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        #creates the new name for the molecule
        make_num_str = str(make_num)
        #updates the new molecules list
        new_study_molecules[make_num_str] = replaced_linker
        #updates the ran dictionary so no overlap occures
        all_ran_smi[make_num_str] = replaced_linker

save_dictionary(new_study_molecules, f'new_study_molecules_{run_number}.json')
