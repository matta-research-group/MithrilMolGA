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
from MithrilMolGA.molecule_mutation import *
from MithrilMolGA.calculation_status import *
from datetime import datetime
import random
import argparse
import subprocess
import os

# ==============================
# Defaults at import time
# ==============================
run_num_str = '0'
run_size = 50
elite_df = None  # placeholder, safe at import

# ==============================
# Command-line argument parser
# ==============================
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_num", type=int, default=0)
    parser.add_argument("--run_size", type=int, default=50)
    return parser.parse_args()

# ==============================
# Safe function to load the dataframe
# ==============================
def load_elite_df(run_num):
    """Load the elite dataframe for the given run number."""
    df_path = os.path.join("dataframes", f"elite_run_{run_num}_df.csv")
    if os.path.exists(df_path):
        return pd.read_csv(df_path)
    else:
        print(f"No dataframe found at {df_path}, skipping load.")
        return None

# ==============================
# Main GA function
# ==============================
def run_elite_mutation(run_num_str, run_size):
    global elite_df

    # ==============================
    # Mutate the elite molecules
    # Includes new fragments, new linkers and even a whole new molecules
    # User inputted balance for what type of mutations occure and how many new molecules are wanted
    # Idea mutation weighting:
    #  25% are new molecules, 25% have 1 biofragment fragment change, 25% non-biofragment change, 25% linker change
    # ==============================

    # Load the elite df safely
    elite_df = load_elite_df(run_num_str)
    if elite_df is None:
        print("No elite dataframe loaded, skipping GA mutation step.")
        return

    # get the smi strings of the elite_25
    elite_smi = dict(zip(elite_df['Name'], elite_df['SMILES']))
    # load all the molecules that have ever been ran in the GA
    all_ran_molecules_dic = open_dictionary('ga_dic/total_molecules_ran.json')
    all_ran_smi_canon = {k: Chem.CanonSmiles(v, useChiral=0) for k, v in all_ran_molecules_dic.items()}

    # linker_dic
    linker_dic = open_dictionary('ga_dic/linker_dic.json')
    # bio_dic
    bio_dic = open_dictionary('ga_dic/bio_dic.json')
    # non_bio_dic
    non_bio_dic = open_dictionary('ga_dic/non_bio.json')

    new_study_molecules = {}

    for k, v in elite_smi.items():
        # options to mutate
        choices = ['new_bio', 'new_non_bio', 'new_linker']
        # which random one is chosen for this molecule
        selected_choice = random.choices(choices, k=1)[0]

        if selected_choice == 'new_bio':
            # find the biofragment and replace it with a new biofragment
            new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'bio')

            # there are less bio fragments so could run out of them
            # this is a fail safe to make sure it doesn't get stuck
            attempts = 0
            while Chem.CanonSmiles(new_molecule, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'bio')
                attempts += 1
                print(f'New Bio attempt {attempts}')

            if attempts == 15:
                new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
                print(f'New Bio attempt maxed out, using non bio fragment')

                attempts = 0
                while Chem.CanonSmiles(new_molecule, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                    attempts += 1
                    print(f'New non bio attempt {attempts}')
                if attempts == 15:
                    continue
                    print(f'New Bio attempt maxed out, moving on')

            last_key, last_value = list(all_ran_smi_canon.items())[-1]
            make_num = int(last_key) + 1
            make_num_str = str(make_num)
            new_study_molecules[make_num_str] = new_molecule
            all_ran_smi_canon[make_num_str] = new_molecule

        if selected_choice == 'new_non_bio':
            # find the non-biofragment and replace it with a new non-biofragment
            new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
            attempts = 0
            while Chem.CanonSmiles(new_molecule, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
                attempts += 1
                print(f'New non Bio attempt {attempts}')

            if attempts == 15:
                new_molecule = swap_one_fragment(v, bio_dic, non_bio_dic, 'bio')
                print(f'New non Bio attempt maxed out, using bio fragment')

                attempts = 0
                while Chem.CanonSmiles(new_molecule, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                    attempts += 1
                    print(f'New bio attempt {attempts}')
                if attempts == 15:
                    continue
                    print(f'New non Bio attempt maxed out, moving on')

            last_key, last_value = list(all_ran_smi_canon.items())[-1]
            make_num = int(last_key) + 1
            make_num_str = str(make_num)
            new_study_molecules[make_num_str] = new_molecule
            all_ran_smi_canon[make_num_str] = new_molecule

        if selected_choice == 'new_linker':
            # find the linker and replace it with a new linker that is not the same as old one
            linker_type = find_linker_type(Chem.MolFromSmiles(v))
            fragments = fragment_molecule(Chem.MolFromSmiles(v), linker_type)
            replaced_linker = replace_linker(fragments, linker_dic)

            attempts = 0
            while Chem.CanonSmiles(replaced_linker, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                replaced_linker = replace_linker(fragments, linker_dic)
                attempts += 1
                print(f'New linker attempt {attempts}')

            if attempts == 15:
                replaced_linker = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
                print(f'New linker attempt maxed out, using non bio fragment')

                attempts = 0
                while Chem.CanonSmiles(replaced_linker, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                    replaced_linker = swap_one_fragment(v, bio_dic, non_bio_dic, 'non_bio')
                    attempts += 1
                    print(f'New bio attempt {attempts}')
                if attempts == 15:
                    continue
                    print(f'New non Bio attempt maxed out for linker, moving on')

            last_key, last_value = list(all_ran_smi_canon.items())[-1]
            make_num = int(last_key) + 1
            make_num_str = str(make_num)
            new_study_molecules[make_num_str] = replaced_linker
            all_ran_smi_canon[make_num_str] = replaced_linker

    # ==============================
    # Generate molecules to maintain run size
    # ==============================
    current_run = open_dictionary(f'submission_dic/molecules_to_run_{run_num_str}.json')
    length_of_run = len(current_run)
    length_of_new_run = len(new_study_molecules)
    new_molecules_needed = run_size - length_of_new_run

    new_molecules = {}
    if new_molecules_needed > 0:
        for i in range(new_molecules_needed + 1):
            random_linker_type = random.choice(list(linker_dic.keys()))
            random_linker = linker_dic[random_linker_type]

            random_bio = random.choice(list(bio_dic.keys()))
            fragment_1 = bio_dic[random_bio]

            weights = [0.25, 0.75]
            choices = ['bio', 'non_bio']
            selected_choice = random.choices(choices, weights=weights, k=1)[0]

            if selected_choice == 'bio':
                random_bio_2 = random.choice(list(bio_dic.keys()))
                fragment_2 = bio_dic[random_bio_2]
            else:
                random_non_bio = random.choice(list(non_bio_dic.keys()))
                fragment_2 = non_bio_dic[random_non_bio]

            frag_1_and_linker = combine_structure(fragment_1, random_linker)
            final_mol = combine_structure(frag_1_and_linker, fragment_2)

            attempts = 0
            while Chem.CanonSmiles(final_mol, useChiral=0) in all_ran_smi_canon.values() and attempts < 15:
                attempts += 1
                print(f'New molecule {attempts}')
                random_linker_type = random.choice(list(linker_dic.keys()))
                random_linker = linker_dic[random_linker_type]

                random_bio = random.choice(list(bio_dic.keys()))
                fragment_1 = bio_dic[random_bio]

                selected_choice = random.choices(choices, weights=weights, k=1)[0]

                if selected_choice == 'bio':
                    random_bio_2 = random.choice(list(bio_dic.keys()))
                    fragment_2 = bio_dic[random_bio_2]
                else:
                    random_non_bio = random.choice(list(non_bio_dic.keys()))
                    fragment_2 = non_bio_dic[random_non_bio]

                frag_1_and_linker = combine_structure(fragment_1, random_linker)
                final_mol = combine_structure(frag_1_and_linker, fragment_2)

            if attempts == 15:
                print(f'New molecule maxed out')
                continue

            new_molecules[str(int(i))] = final_mol

    # ==============================
    # Renumber dictionaries and save
    # ==============================
    last_key = list(all_ran_molecules_dic.keys())[-1]
    base_number = int(last_key)
    new_molecules_renumbered = {f"{base_number + 1 + i}": v for i, (k, v) in enumerate(new_molecules.items())}

    if len(new_molecules_renumbered) == 0:
        molecules_to_run = new_study_molecules
        run_num_int = int(run_num_str)
        new_run_num_str = str(run_num_int + 1)
        save_dictionary(molecules_to_run, f'submission_dic/molecules_to_run_{new_run_num_str}.json')
    else:
        last_key = list(new_molecules_renumbered.keys())[-1]
        base_number = int(last_key)
        renumbered_dict_mutation = {f"{base_number + 1 + i}": v for i, (k, v) in enumerate(new_study_molecules.items())}
        molecules_to_run = new_molecules_renumbered | renumbered_dict_mutation
        run_num_int = int(run_num_str)
        new_run_num_str = str(run_num_int + 1)
        save_dictionary(molecules_to_run, f'submission_dic/molecules_to_run_{new_run_num_str}.json')

    # ==============================
    # Log GA progress
    # ==============================
    progress_file_path = 'GA_status.txt'
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    with open(progress_file_path, 'a') as file:
        file.write(f'elite_mutation complete for run {run_num_str} at {current_time}.\n')

    with open('talk_to_GA.txt', 'r') as f:
        content = f.read().strip()
        if content == 'end':
            file.write("The file says 'end'")
            subprocess.run(['scancel', '-p', 'long_cpu', '--me'])
        else:
            print("The file does not say 'end' continue the GA")

# ==============================
# Only run GA if executed directly
# ==============================
if __name__ == "__main__":
    args = parse_args()
    run_num_str = str(args.run_num)
    run_size = args.run_size
    run_elite_mutation(run_num_str, run_size)

