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
import argparse

# Define default variables
options = {
    'run_num': {'default': 0}
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



#Mutate the elite molecules
#Includes new fragments, new linkers and even a whole new molecules
#User inputted balance for what type of mutations occure and how many new molecules are wanted
#Idea mutation weighting:
# 25% are new molecules, 25% have 1 biofragment fragment change, 25% non-biofragment change, 25% linker change

#load the elite df
elite_df = pd.read_csv(f'elite_run_{run_num_str}_df.csv')
#get th smi strings of the eilte_25
elite_smi = dict(zip(elite_df['Name'], elite_df['SMILES']))
#load all the molecules that have ever been ran in the GA
all_ran_molecules_dic = open_dictionary('total_molecules_ran.json')

all_ran_smi_canon = {k: Chem.CanonSmiles(v) for k, v in all_ran_molecules_dic.items()}

#linker_dic
linker_dic = open_dictionary('linker_dic.json')
#bio_dic
bio_dic = open_dictionary('bio_dic.json')
#non_bio_dic
non_bio_dic = open_dictionary('non_bio_dic.json')

new_study_molecules = {}
for k, v in elite_smi.items():
    #options to mutate
    choices = ['new_bio', 'new_non_bio', 'new_linker']
    #which random one is chosen for this molecule
    selected_choice = random.choices(choices, k=1)[0]

    if selected_choice == 'new_bio':
        #find the biofragment and replace it with a new biofragment
        new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'bio')
        while Chem.CanonSmiles(new_molecule) in all_ran_smi_canon.values():
            new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'bio')
        
        last_key, last_value = list(all_ran_smi_canon.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        #creates the new name for the molecule
        make_num_str = str(make_num)

        #updates the new molecules list
        new_study_molecules[make_num_str] = new_molecule
        #updates the ran dictionary so no overlap occures
        all_ran_smi_canon[make_num_str] = new_molecule

    if selected_choice == 'new_non_bio':
        #do this
        #find the non-biofragment and replace it with a new non-biofragment
        new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
        while Chem.CanonSmiles(new_molecule) in all_ran_smi_canon.values():
            new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
        
        last_key, last_value = list(all_ran_smi_canon.items())[-1]
        #updates what the key will be by turning to int and then back to str
        make_num = int(last_key) + 1
        #creates the new name for the molecule
        make_num_str = str(make_num)
        #updates the new molecules list
        new_study_molecules[make_num_str] = new_molecule
        #updates the ran dictionary so no overlap occures
        all_ran_smi_canon[make_num_str] = new_molecule

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
        all_ran_smi_canon[make_num_str] = replaced_linker

#make a new list of molecules to run and make it same length as the molecule ran in list
current_run = open_dictionary(f'molecules_to_run_{run_num_str}.json')
length_of_run = len(current_run)
length_of_new_run = len(new_study_molecules)
#new_molecules_needed = length_of_run - length_of_new_run

#THIS IS JUST FOR TESTING PURPOSES
new_molecules_needed = 50 - length_of_new_run
#new molecules to make up the numbers lost via elite step
new_molecules = {}
for i in range(new_molecules_needed + 1):
    
    #pick a random linker
    random_linker_type = random.choice(list(linker_dic.keys()))
    random_linker = linker_dic[random_linker_type]

    #pick a random bio fragment
    random_bio = random.choice(list(bio_dic.keys()))
    fragment_1 = bio_dic[random_bio]

    #pick at random if the second fragment will be bio or non_bio
    weights = [0.25, 0.75]  # 70% chance for option1, 30% chance for option2
    choices = ['bio', 'non_bio']
    selected_choice = random.choices(choices, weights=weights, k=1)[0]
    if selected_choice == 'bio':
        random_bio_2 = random.choice(list(bio_dic.keys()))
        fragment_2 = bio_dic[random_bio_2]

    if selected_choice == 'non_bio':
        random_non_bio = random.choice(list(non_bio_dic.keys()))
        fragment_2 = non_bio_dic[random_non_bio]

    #connect them all together
    frag_1_and_linker = combine_structure(fragment_1, random_linker)
    final_mol = combine_structure(frag_1_and_linker, fragment_2)

    #make sure its a molecule that hasn't been ran before
    while Chem.CanonSmiles(final_mol) in all_ran_smi_canon.values():
        #pick a random linker
        random_linker_type = random.choice(list(linker_dic.keys()))
        random_linker = linker_dic[random_linker_type]

        #pick a random bio fragment
        random_bio = random.choice(list(bio_dic.keys()))
        fragment_1 = bio_dic[random_bio]

        #pick at random if the second fragment will be bio or non_bio
        weights = [0.25, 0.75]  # 25% chance for bio, 75% chance for non_bio
        choices = ['bio', 'non_bio']
        selected_choice = random.choices(choices, weights=weights, k=1)[0]
        if selected_choice == 'bio':
            random_bio_2 = random.choice(list(bio_dic.keys()))
            fragment_2 = bio_dic[random_bio_2]

        if selected_choice == 'non_bio':
            random_non_bio = random.choice(list(non_bio_dic.keys()))
            fragment_2 = non_bio_dic[random_non_bio]

        #connect them all together
        frag_1_and_linker = combine_structure(fragment_1, random_linker)
        final_mol = combine_structure(frag_1_and_linker, fragment_2)

    new_molecules[str(int(i))] = final_mol


# Extract the last key from the reference dictionary
last_key = list(all_ran_molecules_dic.keys())[-1]

# Determine the base number from the last key
base_number = int(last_key)

# Renumber the new dictionary
new_molecules_renumbered = {f"{base_number + 1 + i}": v for i, (k, v) in enumerate(new_molecules.items())}

if len(new_molecules_renumbered) == 0:
    molecules_to_run = new_study_molecules

else:

    # Extract the last key from the reference dictionary
    last_key = list(new_molecules_renumbered.keys())[-1]

    # Determine the base number from the last key
    base_number = int(last_key)

    # Renumber the new dictionary
    renumbered_dict_mutation = {f"{base_number + 1 + i}": v for i, (k, v) in enumerate(new_study_molecules.items())}

    molecules_to_run = new_molecules_renumbered | renumbered_dict_mutation

run_num_int = int(run_num_str)
new_run_num = run_num_int + 1
new_run_num_str = str(new_run_num)

save_dictionary(molecules_to_run, f'molecules_to_run_{new_run_num_str}.json')

progress_file_path = 'GA_status.txt'

# Open the file in append mode and write some content
with open(progress_file_path, 'a') as file:
    file.write(f'elite_mutation complete for run {run_num_str}.\n')
