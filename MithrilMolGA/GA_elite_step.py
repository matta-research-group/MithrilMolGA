import os
import argparse
import pandas as pd
from datetime import datetime
from QCflow.load_gaussian import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *
from QCflow.write_psi4 import *
from QCflow.run_psi4 import *
from QCflow.energy_calculations import *
from MithrilMolGA.molecule_mutation import *
from MithrilMolGA.calculation_status import *

# ==============================
# Default variables
# ==============================
options = {
    'run_num': {'default': 0},
    'EG_rank_weight': {'default': 1.0},
    'planarity_rank_weight': {'default': 1.0},
    'SA_rank_weight': {'default': 0.25},
    'elite_value': {'default': 25},
    'anioinc_reorg_rank_weight': {'default': 2.0},
    'cationic_reorg_rank_weight': {'default': 2.0},
    'functional': {'default': 'b3lyp'},
    'basis_set': {'default': '6-31g*'},
    'time': {'default': 4},
    'cpus': {'default': 10},
    'eg_elite_value': {'default': 2.5},
    'planarity_elite_value': {'default': 0.82},
    'anioinc_reorg_elite_value': {'default': 0.350},
    'catioinc_reorg_elite_value': {'default': 0.350}
}

# ==============================
# Argument parser
# ==============================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for arg, opts in options.items():
        parser.add_argument(f'--{arg}', type=type(opts['default']), default=opts['default'])
    
    args, unknown = parser.parse_known_args()  # ignore unknown args passed from master script

    # ==============================
    # Store variables
    # ==============================
    run_num_str = str(args.run_num) if hasattr(args, 'run_num') and args.run_num else str(options['run_num']['default'])
    EG_rank_weight = args.EG_rank_weight
    planarity_rank_weight = args.planarity_rank_weight
    SA_rank_weight = args.SA_rank_weight
    elite_value = args.elite_value
    anioinc_reorg_rank_weight = args.anioinc_reorg_rank_weight
    cationic_reorg_rank_weight = args.cationic_reorg_rank_weight
    functional = args.functional
    basis_set = args.basis_set
    time = args.time
    cpus = args.cpus
    eg_elite_value = args.eg_elite_value
    planarity_elite_value = args.planarity_elite_value
    anioinc_reorg_elite_value = args.anioinc_reorg_elite_value
    catioinc_reorg_elite_value = args.catioinc_reorg_elite_value

    # ==============================
    # Load molecule dataframe
    # ==============================
    molecule_path = f'dataframes/run_{run_num_str}_data.csv'
    if os.path.exists(molecule_path):
        molecule_df = pd.read_csv(molecule_path)
    else:
        raise FileNotFoundError(f"{molecule_path} not found. Cannot proceed without molecule data.")

    molecule_df['Name'] = molecule_df['Name'].astype(str)

    # ==============================
    # Rank EG, Planarity, SA
    # ==============================
    molecule_df['EG Rank Order'] = molecule_df['EG /eV'].rank(ascending=True)
    molecule_df['Plan Rank Order'] = molecule_df['Planarity'].rank(ascending=False)
    molecule_df['SA Rank Order'] = molecule_df['SA Score'].rank(ascending=True)
    molecule_df['Rank Sum'] = (
        molecule_df['EG Rank Order'] * EG_rank_weight +
        molecule_df['Plan Rank Order'] * planarity_rank_weight +
        molecule_df['SA Rank Order'] * SA_rank_weight
    )

    # ==============================
    # Elite selection
    # ==============================
    sorted_molecule_df = molecule_df.sort_values(['Rank Sum'], ascending=True)
    elite_df = sorted_molecule_df.head(int(len(sorted_molecule_df)*(elite_value/100)))
    elite_df = elite_df[(elite_df['EG /eV'] <= eg_elite_value) & (elite_df['Planarity'] >= planarity_elite_value)]
    elite_smi_int = dict(zip(elite_df['Name'], elite_df['SMILES']))

    # ==============================
    # Archive handling
    # ==============================
    archive_path = 'archive_dataframes/all_ran_reorg.csv'
    if os.path.exists(archive_path):
        archive_df = pd.read_csv(archive_path)
    else:
        print(f"Warning: {archive_path} not found. Using empty DataFrame.")
        archive_df = pd.DataFrame()

    archive_elite_results = archive_df.iloc[0:0].copy()
    elite_smi = {}
    elite_archive_smi = {}

    for k, v in elite_smi_int.items():
        canon_smi = v
        if len(archive_df[archive_df['SMILES'] == canon_smi]) > 0:
            elite_archive_smi[k] = canon_smi
            ran_before = archive_df[archive_df['SMILES'] == canon_smi].copy()
            ran_before['Name'] = int(k)
            archive_elite_results = pd.concat([archive_elite_results, ran_before])
            print(k, 'Elite has been run before, data retrieved')
        else:
            elite_smi[str(k)] = v
            print(k, 'Elite has NOT been run before')

    save_dictionary(elite_archive_smi, f'run_dic/elite_archive_smi_{run_num_str}_molecules.json')

    if len(elite_smi) <= 0:
        elite_smi_int = dict(zip(sorted_molecule_df['Name'], sorted_molecule_df['SMILES']))
        length_of_elite = len(elite_df)
        subset_subset = dict(list(elite_smi_int.items())[(length_of_elite+1):length_of_elite+6])
        elite_smi = {str(int(k)): v for k, v in subset_subset.items()}
        elite_df = sorted_molecule_df

    print('Submit Calculations')

    # ==============================
    # Psi4 calculations
    # ==============================
    os.makedirs('data', exist_ok=True)
    os.chdir('data')
    for k, v in elite_smi.items():
        k = str(int(k))
        run_psi4('anion', k, v, time, cpus, functional, basis_set)
        run_psi4('cation', k, v, time, cpus, functional, basis_set)
        run_psi4('sp_c', k, v, time, cpus, functional, basis_set)
        run_psi4('sp_a', k, v, time, cpus, functional, basis_set)
    os.chdir('../')

task_list = []
for k, v in elite_smi.items():
    k = str(int(k))
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
succesful_dict, failed_dict, attempts = calculations_status(task_list, sleep_time=10)

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
    k = str(int(k))
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
        k = str(int(k))
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

#leave data folder
os.chdir('../')

#create dataframe of elite 25% with reorganmsaition energy
elite_df = elite_df.drop(['EG Rank Order', 'Plan Rank Order', 'SA Rank Order', 'Rank Sum'], axis=1)
elite_df.insert(7, 'Anionic Reorganisation Energy /eV', reorganisation_anionic.values())
elite_df.insert(8, 'Cationic Reorganisation Energy /eV', reorganisation_cationic.values())

#combine old elite 25% with new elite 25% and then sort again
run_num_float = float(run_num_str)
previous_run_num = str(int(run_num_float - 1))
old_elite_df = pd.read_csv(f'dataframes/elite_run_{previous_run_num}_df.csv') #previous elite 25% dataframe

#combine the three dataframes (old elite, the new elite that has been calculated and the archive elite run results)
combined_elite_df = pd.concat([old_elite_df, elite_df, archive_elite_run_results])

#sort the combined dataframe
combined_elite_df['EG Rank Order'] = combined_elite_df['EG /eV'].rank(ascending=True)
combined_elite_df['Plan Rank Order'] = combined_elite_df['Planarity'].rank(ascending=False)
combined_elite_df['SA Rank Order'] = combined_elite_df['SA Score'].rank(ascending=True)
combined_elite_df['Anionic Reorg Rank Order'] = combined_elite_df['Anionic Reorganisation Energy /eV'].rank(ascending=True)
combined_elite_df['Cationic Reorg Rank Order'] = combined_elite_df['Cationic Reorganisation Energy /eV'].rank(ascending=True)

#good at anioinc reorganisation energy
combined_elite_df['Rank Sum Anionic'] = (combined_elite_df['EG Rank Order']*EG_rank_weight) + (combined_elite_df['Plan Rank Order']*planarity_rank_weight) + (combined_elite_df['SA Rank Order']*SA_rank_weight) + (combined_elite_df['Anionic Reorg Rank Order']*anioinc_reorg_rank_weight)
#good at cationic reorganisation energy
combined_elite_df['Rank Sum Cationic'] = (combined_elite_df['EG Rank Order']*EG_rank_weight) + (combined_elite_df['Plan Rank Order']*planarity_rank_weight) + (combined_elite_df['SA Rank Order']*SA_rank_weight) + (combined_elite_df['Cationic Reorg Rank Order']*cationic_reorg_rank_weight)

#rank the anionic reorganisation energy molecules
sorted_anionic_reorg_df = combined_elite_df.sort_values(['Rank Sum Anionic'], ascending=True)
#rank the cationic reorganisation energy molecules
sorted_cationic_reorg_df = combined_elite_df.sort_values(['Rank Sum Cationic'], ascending=True)

#retreave the elite 25% of the molecules
#elite_anionic = sorted_anionic_reorg_df.head(int(len(sorted_anionic_reorg_df)*(elite_value/100)))
#elite_cationic = sorted_cationic_reorg_df.head(int(len(sorted_cationic_reorg_df)*(elite_value/100)))

elite_anionic = sorted_anionic_reorg_df
elite_cationic = sorted_cationic_reorg_df

new_elite_df = pd.concat([elite_anionic, elite_cationic, archive_elite_run_results]) #add the archive elite run results to the new elite df
new_elite_df = new_elite_df.drop_duplicates()

#make sure the elite doesn't get too small and can be set by the user to be at least a certain size
if len(new_elite_df) <= 0:
    new_elite_df = pd.concat([sorted_anionic_reorg_df, sorted_cationic_reorg_df]).drop_duplicates()

else:
    new_elite_df = new_elite_df

#similarity comparison
old_elite_smi = dict(zip(old_elite_df['Name'], old_elite_df['SMILES']))
new_elite_smi = dict(zip(new_elite_df['Name'], new_elite_df['SMILES']))

same_as_old = {}
for k1, v1 in old_elite_smi.items():
    for k2, v2 in new_elite_smi.items():
        if Chem.CanonSmiles(v1, useChiral=0) == Chem.CanonSmiles(v2, useChiral=0):
            same_as_old[k2] = v2

length_of_old = len(old_elite_smi)
length_of_same = len(same_as_old)
length_of_new = len(new_elite_smi)

#how similar in percentage are the two elite 25% lists
if length_of_old == 0 or length_of_same == 0:
    similarity_percentage = 0
else:
    similarity_percentage = (length_of_same/length_of_new)*100

similarity_percentage_dic = {f'Run {run_num_str}': similarity_percentage}

current_run_df = pd.DataFrame()
current_run_df.insert(0, 'Name', similarity_percentage_dic.keys())
current_run_df.insert(1, 'Similarity Percentage', similarity_percentage_dic.values())

run_df = pd.read_csv(f'dataframes/run_df.csv') #this dataframe tracks the progress of the GA

run_df_combined = pd.concat([run_df, current_run_df])
run_df_combined.to_csv(f'dataframes/run_df.csv', index=False)

#save new elite dataframe without the scoring
new_elite_df = new_elite_df.drop(['EG Rank Order', 'Plan Rank Order', 'SA Rank Order', 'Anionic Reorg Rank Order', 'Cationic Reorg Rank Order', 'Rank Sum Anionic', 'Rank Sum Cationic'], axis=1)

#removing poor elites that don't match the threshold values
elite_df_plan_eg_filter = new_elite_df[(new_elite_df['EG /eV'] <= eg_elite_value) & (new_elite_df['Planarity'] >= planarity_elite_value )]

elite_df_filter_final = elite_df_plan_eg_filter[(elite_df_plan_eg_filter['Anionic Reorganisation Energy /eV'] <= anioinc_reorg_elite_value) | (elite_df_plan_eg_filter['Cationic Reorganisation Energy /eV'] <= catioinc_reorg_elite_value)]

elite_df_filter_final.to_csv(f'dataframes/elite_run_{run_num_str}_df.csv', index=False)

#combine all the runs into one big dataframe
all_ran_molecules = pd.read_csv(f'dataframes/ran_all_data.csv')

clean_molecule_df = molecule_df.drop(['EG Rank Order', 'Plan Rank Order', 'SA Rank Order', 'Rank Sum'], axis=1)

# Create a boolean mask
mask = clean_molecule_df['Name'].isin(elite_df['Name'])

# Apply the mask to filter to remove the duplicates
filtered_clean_molecule_df = clean_molecule_df[~mask]

# Print the filtered DataFrame

adding_new_runs = pd.concat([all_ran_molecules, filtered_clean_molecule_df, elite_df])
adding_new_runs_no_dup = adding_new_runs.drop_duplicates()
adding_new_runs_no_dup.to_csv(f'dataframes/ran_all_data.csv', index=False)


# Open the file in append mode and write some content
with open(progress_file_path, 'a') as file:
    file.write(f'GA_elite_step complete for run {run_num_str} at {current_time}.\n')
    # ==============================
    # Update GA status
    # ==============================
    progress_file_path = 'GA_status.txt'
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(progress_file_path, 'a') as file:
        file.write(f'GA_elite_step complete for run {run_num_str} at {current_time}.\n')

    print('GA_elite_step complete')
