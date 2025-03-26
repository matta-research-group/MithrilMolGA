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
import argparse

#Open the data df
#Sort top EG, top Plan, top SA
#Submit top 25% for reorganisation calculation (need to do some Psi4 stuff)
#Can definitely be all done within a few psi4 calculations rather than the gaussian splitting across multiple files
#Retreave reorganisation data
#Make a new df with the reorganisation data and then combine with old elite 25% from previous step and then sort again

#The user inputted weighting of the different parameters calculated

# Define default variables
options = {
    'run_num': {'default': 0},
    'EG_rank_weight': {'default': 1},
    'planarity_rank_weight': {'default': 1},
    'SA_rank_weight': {'default': 4},
    'elite_value': {'default': 25},
    'anioinc_reorg_rank_weight': {'default': 0.5},
    'cationic_reorg_rank_weight': {'default': 0.5},
    'elite_df_size': {'default': 10},
    'functional' : {'default': 'b3lyp'},
    'basis_set' : {'default': '6-31g*'},
    'time' : {'default': 4},
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
EG_rank_weight = args.EG_rank_weight if hasattr(args, 'EG_rank_weight') else options['EG_rank_weight']['default']
planarity_rank_weight = args.planarity_rank_weight if hasattr(args, 'planarity_rank_weight') else options['planarity_rank_weight']['default']
SA_rank_weight = args.SA_rank_weight if hasattr(args, 'SA_rank_weight') else options['SA_rank_weight']['default']
elite_value = args.elite_value if hasattr(args, 'elite_value') else options['elite_value']['default']
anioinc_reorg_rank_weight = args.anioinc_reorg_rank_weight if hasattr(args, 'anioinc_reorg_rank_weight') else options['anioinc_reorg_rank_weight']['default']
cationic_reorg_rank_weight = args.cationic_reorg_rank_weight if hasattr(args, 'cationic_reorg_rank_weight') else options['cationic_reorg_rank_weight']['default']
elite_df_size = args.elite_df_size if hasattr(args, 'elite_df_size') else options['elite_df_size']['default']
functional = args.functional if hasattr(args, 'functional') else options['functional']['default']
basis_set = args.basis_set if hasattr(args, 'basis_set') else options['basis_set']['default']
time = args.time if hasattr(args, 'time') else options['time']['default']
cpus = args.cpus if hasattr(args, 'cpus') else options['cpus']['default']

#load molecule df
molecule_df = pd.read_csv(f'run_{run_num_str}_data.csv')


molecule_df['EG Rank Order'] = molecule_df['EG /eV'].rank(ascending=True)
molecule_df['Plan Rank Order'] = molecule_df['Planarity'].rank(ascending=False)
molecule_df['SA Rank Order'] = molecule_df['SA Score'].rank(ascending=True)

#having SA score not have as much weight
molecule_df['Rank Sum'] = (molecule_df['EG Rank Order']/EG_rank_weight) + (molecule_df['Plan Rank Order']/planarity_rank_weight) + (molecule_df['SA Rank Order']/SA_rank_weight)

#rank the combined planarity, energy gap and SA score rankings to produce the best balanced molecule
sorted_molecule_df = molecule_df.sort_values(['Rank Sum'], ascending=True)

#retreave the elite 25% of the molecules
elite_df = sorted_molecule_df.head(int(len(sorted_molecule_df)*(elite_value/100)))

#dict of elite 25 monomers and their SMILES
elite_smi_int = dict(zip(elite_df['Name'], elite_df['SMILES']))

elite_smi = {}
for k, v in elite_smi_int.items():
    elite_smi[str(k)] = v

#make sure the new potential elites is long enough
if len(elite_smi) <= elite_df_size:
    elite_smi_int = dict(zip(sorted_molecule_df['Name'], sorted_molecule_df['SMILES']))

    elite_smi = {}
    for k, v in elite_smi_int.items():
        elite_smi[str(k)] = v

    elite_df = sorted_molecule_df

print('Submit Calculations')
# run reorganisation calucltions for elite 25%
# THE psi4 scripts have not been written for these yet
# The opt_c and opt_a will also hav to contain the n_c_geo and n_a_geo as those calucltions rely off the coordinates of the optimised geometry
for k, v in elite_smi.items():
    run_psi4('anion', str(k), v, time, cpus, functional, basis_set) #user set parameters
    run_psi4('cation', str(k), v, time, cpus, functional, basis_set) 
    run_psi4('sp_c', str(k), v, time, cpus, functional, basis_set) 
    run_psi4('sp_a', str(k), v, time, cpus, functional, basis_set)

print('Reorgansaition calculations submited')

#turn into a list of tasks that calculation_status function can proccess
task_list = []
for k, v in elite_smi.items():
    k = str(k)
    task_opt_c = lambda: is_file_present(f'{k}/{k}_opt_c_energy_and_gap.txt', 'opt_c')
    task_opt_a = lambda: is_file_present(f'{k}/{k}_opt_a_energy_and_gap.txt', 'opt_a')
    task_sp_c = lambda: is_file_present(f'{k}/{k}_sp_c_energy_and_gap.txt')
    task_sp_a = lambda: is_file_present(f'{k}/{k}_sp_a_energy_and_gap.txt')
    task_n_c_geo = lambda: is_file_present(f'{k}/{k}_n_c_geo_energy_and_gap.txt', 'n_c_geo')
    task_n_a_geo = lambda: is_file_present(f'{k}/{k}_n_a_geo_energy_and_gap.txt', 'n_a_geo')
    #add all the tasks to the task list
    task_list.append((f'{k}_opt_c', task_opt_c()))
    task_list.append((f'{k}_opt_a', task_opt_a()))
    task_list.append((f'{k}_sp_c', task_sp_c()))
    task_list.append((f'{k}_sp_a', task_sp_a()))
    task_list.append((f'{k}_n_c_geo', task_n_c_geo()))
    task_list.append((f'{k}_n_a_geo', task_n_a_geo()))

print('Testing calculations status')

#returns the failed and successful calculations, keeps looping until all calculations are done
succesful_dict, failed_dict, attempts = calculations_status(task_list, sleep_time=5)

print('Calculations status test finished')

failed_molecules = {}
for k, v in failed_dict.items():
    k = k.split('_')[0]
    failed_molecules[k] = elite_smi[k]

succesful_molecules = {k: v for k, v in elite_smi.items() if k not in failed_molecules}

#How this data is extracted needs to be determined by psi4 input but this is sudo code as follows
reorganisation_anionic = {}
reorganisation_cationic = {}
for k, v in succesful_molecules.items():
    #extract data from the successful monomers
    data_opt = extract_data_from_txt(f'{k}/{k}_opt_energy_and_gap.txt')
    data_opt_c = extract_data_from_txt(f'{k}/{k}_opt_c_energy_and_gap.txt')
    data_opt_a = extract_data_from_txt(f'{k}/{k}_opt_a_energy_and_gap.txt')
    data_sp_c = extract_data_from_txt(f'{k}/{k}_sp_c_energy_and_gap.txt')
    data_sp_a = extract_data_from_txt(f'{k}/{k}_sp_a_energy_and_gap.txt')
    data_n_c_geo = extract_data_from_txt(f'{k}/{k}_n_c_geo_energy_and_gap.txt')
    data_n_a_geo = extract_data_from_txt(f'{k}/{k}_n_a_geo_energy_and_gap.txt')
    #energy calcs
    cation_reorg = cal_reorg(data_opt, data_sp_c, data_opt_c, data_n_c_geo, calculation_software='Psi4')
    anion_reorg = cal_reorg(data_opt, data_sp_a, data_opt_a, data_n_a_geo, calculation_software='Psi4')

    reorganisation_anionic[k] = anion_reorg
    reorganisation_cationic[k] = cation_reorg

print('Calculating reorg complete')

#Making sure that if just one of the reorganisation energies have failed then the molecules still progresses
if len(failed_dict.items()) != 0:
    for k, v in elite_smi.items():
        if all(f'{k}_{job_name}' not in failed_dict for job_name in ['n_c_geo', 'opt_c', 'sp_c']) and any(f'{k}_{job_name}' in failed_dict for job_name in ['n_a_geo', 'opt_a', 'sp_a']):
            print(f'Cation reorganisation did not fail for {k}')
            data_opt = extract_data_from_txt(f'{k}/{k}_opt_energy_and_gap.txt')
            data_opt_c = extract_data_from_txt(f'{k}/{k}_opt_c_energy_and_gap.txt')
            data_sp_c = extract_data_from_txt(f'{k}/{k}_sp_c_energy_and_gap.txt')
            data_n_c_geo = extract_data_from_txt(f'{k}/{k}_n_c_geo_energy_and_gap.txt')

            cation_reorg = cal_reorg(data_opt, data_sp_c, data_opt_c, data_n_c_geo, calculation_software='Psi4')
            reorganisation_cationic[k] = cation_reorg
            reorganisation_anionic[k] = 100000 #stupidly large number so to still allow for it to be present but not affect data
        
        elif all(f'{k}_{job_name}' not in failed_dict for job_name in ['n_a_geo', 'opt_a', 'sp_a']) and any(f'{k}_{job_name}' in failed_dict for job_name in ['n_c_geo', 'opt_c', 'sp_c']):
            print(f'Anioinc reorganisation did not fail for {k}')
            data_opt = extract_data_from_txt(f'{k}/{k}_opt_energy_and_gap.txt')
            data_opt_a = extract_data_from_txt(f'{k}/{k}_opt_a_energy_and_gap.txt')
            data_sp_a = extract_data_from_txt(f'{k}/{k}_sp_a_energy_and_gap.txt')
            data_n_a_geo = extract_data_from_txt(f'{k}/{k}_n_a_geo_energy_and_gap.txt')

            anion_reorg = cal_reorg(data_opt, data_sp_a, data_opt_a, data_n_a_geo, calculation_software='Psi4')
            reorganisation_anionic[k] = anion_reorg
            reorganisation_cationic[k] = 100000 #stupidly large number so to still allow for it to be present but not affect data
        
        elif any(f'{k}_{job_name}' in failed_dict for job_name in ['n_c_geo', 'opt_c', 'sp_c']) and any(f'{k}_{job_name}' in failed_dict for job_name in ['n_a_geo', 'opt_a', 'sp_a']):
             reorganisation_cationic[k] = 100000
             reorganisation_anionic[k] = 100000

#reorder to make sure it matches up for df
reorganisation_cationic = reorder_dict(elite_smi, reorganisation_cationic)
reorganisation_anionic = reorder_dict(elite_smi, reorganisation_anionic)


print('semi failed reorg molecules calculated complete')

#create dataframe of elite 25% with reorganmsaition energy
elite_df = elite_df.drop(['EG Rank Order', 'Plan Rank Order', 'SA Rank Order', 'Rank Sum'], axis=1)
elite_df.insert(7, 'Anionic Reorganisation Energy /eV', reorganisation_anionic.values())
elite_df.insert(8, 'Cationic Reorganisation Energy /eV', reorganisation_cationic.values())

#combine old elite 25% with new elite 25% and then sort again
run_num_float = float(run_num_str)
previous_run_num = str(int(run_num_float - 1))
old_elite_df = pd.read_csv(f'elite_run_{previous_run_num}_df.csv') #previous elite 25% dataframe

#combine the two dataframes
combined_elite_df = pd.concat([old_elite_df, elite_df])

#sort the combined dataframe
combined_elite_df['EG Rank Order'] = combined_elite_df['EG /eV'].rank(ascending=True)
combined_elite_df['Plan Rank Order'] = combined_elite_df['Planarity'].rank(ascending=False)
combined_elite_df['SA Rank Order'] = combined_elite_df['SA Score'].rank(ascending=True)
combined_elite_df['Anionic Reorg Rank Order'] = combined_elite_df['Anionic Reorganisation Energy /eV'].rank(ascending=True)
combined_elite_df['Cationic Reorg Rank Order'] = combined_elite_df['Cationic Reorganisation Energy /eV'].rank(ascending=True)

#good at anioinc reorganisation energy
combined_elite_df['Rank Sum Anionic'] = (combined_elite_df['EG Rank Order']/EG_rank_weight) + (combined_elite_df['Plan Rank Order']/planarity_rank_weight) + (combined_elite_df['SA Rank Order']/SA_rank_weight) + (combined_elite_df['Anionic Reorg Rank Order']/anioinc_reorg_rank_weight)
#good at cationic reorganisation energy
combined_elite_df['Rank Sum Cationic'] = (combined_elite_df['EG Rank Order']/EG_rank_weight) + (combined_elite_df['Plan Rank Order']/planarity_rank_weight) + (combined_elite_df['SA Rank Order']/SA_rank_weight) + (combined_elite_df['Cationic Reorg Rank Order']/cationic_reorg_rank_weight)

#rank the anionic reorganisation energy molecules
sorted_anionic_reorg_df = combined_elite_df.sort_values(['Rank Sum Anionic'], ascending=True)
#rank the cationic reorganisation energy molecules
sorted_cationic_reorg_df = combined_elite_df.sort_values(['Rank Sum Cationic'], ascending=True)

#retreave the elite 25% of the molecules
elite_anionic = sorted_anionic_reorg_df.head(int(len(sorted_anionic_reorg_df)*(elite_value/100)))
elite_cationic = sorted_cationic_reorg_df.head(int(len(sorted_cationic_reorg_df)*(elite_value/100)))

new_elite_df = pd.concat([elite_anionic, elite_cationic])
new_elite_df = new_elite_df.drop_duplicates()

#make sure the elite doesn't get too small and can be set by the user to be at least a certain size
if len(new_elite_df) <= elite_df_size:
    new_elite_df = pd.concat([sorted_anionic_reorg_df, sorted_cationic_reorg_df]).drop_duplicates()

else:
    new_elite_df = new_elite_df

#similarity comparison
old_elite_smi = dict(zip(old_elite_df['Name'], old_elite_df['SMILES']))
new_elite_smi = dict(zip(new_elite_df['Name'], new_elite_df['SMILES']))

same_as_old = {}
for k1, v1 in old_elite_smi.items():
    for k2, v2 in new_elite_smi.items():
        if Chem.CanonSmiles(v1) == Chem.CanonSmiles(v2):
            same_as_old[k2] = v2

length_of_old = len(old_elite_smi)
length_of_same = len(same_as_old)

#how similar in percentage are the two elite 25% lists
if length_of_old == 0 and length_of_same == 0:
    similarity_percentage = 0
else:
    similarity_percentage = (length_of_same/length_of_old)*100

similarity_percentage_dic = {f'Run {run_num_str}': similarity_percentage}

current_run_df = pd.DataFrame()
current_run_df.insert(0, 'Name', similarity_percentage_dic.keys())
current_run_df.insert(1, 'Similarity Percentage', similarity_percentage_dic.values())

run_df = pd.read_csv(f'run_df.csv') #this dataframe tracks the progress of the GA

run_df_combined = pd.concat([run_df, current_run_df])
run_df_combined.to_csv(f'run_df.csv', index=False)

#save new elite dataframe without the scoring
new_elite_df = new_elite_df.drop(['EG Rank Order', 'Plan Rank Order', 'SA Rank Order', 'Anionic Reorg Rank Order', 'Cationic Reorg Rank Order', 'Rank Sum Anionic', 'Rank Sum Cationic'], axis=1)
new_elite_df.to_csv(f'elite_run_{run_num_str}_df.csv', index=False)

#combine all the runs into one big dataframe
all_ran_molecules = pd.read_csv(f'ran_all_data.csv')

clean_molecule_df = molecule_df.drop(['EG Rank Order', 'Plan Rank Order', 'SA Rank Order', 'Rank Sum'], axis=1)

# Create a boolean mask
mask = clean_molecule_df['Name'].isin(elite_df['Name'])

# Apply the mask to filter to remove the duplicates
filtered_clean_molecule_df = clean_molecule_df[~mask]

# Print the filtered DataFrame

adding_new_runs = pd.concat([all_ran_molecules, filtered_clean_molecule_df, elite_df])
adding_new_runs_no_dup = adding_new_runs.drop_duplicates()
adding_new_runs_no_dup.to_csv(f'ran_all_data.csv', index=False)


progress_file_path = 'GA_status.txt'

# Open the file in append mode and write some content
with open(progress_file_path, 'a') as file:
    file.write(f'GA_elite_step complete for run {run_num_str}.\n')
