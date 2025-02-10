# GA_TSJ_project
Code base of the genetic algorithm project

## Notes
This contains the scripts that will be used within the GA. Alot of the file names are currently just place holders.

## Dataframes
run_df.csv is the dataframe containing the similarity results of each run. This is used to calculate the average similarity of each run.

elite_25_run_{run_number}_df.csv contains the yop 25% for a specific run. 

ran_all_data.csv contains the data for all the runs. This is to keep track of every molecule that has been studied to avoid duplictes.

run_{X}_data.csv just contains the data for molecules in the current run.

## Scripts
'monomer_run.py' this takes all the potential molecules, checks if the monomers have been run and if they haven't it runs them

'molecule_run.py' takes the monomer data and does D-A matching, if they are a good match it runs the geometry optimisation to get HOMO, LUMO, Energy Gap and Planarity

'GA_data_extraction_step.py' extracts the data from the molecule run step, this runs in th background and waits for all the molecules to finish

'GA_elite_step.py' ranks all the molecules in the current run, takes th top 25% runs reorganisation energy calculations on them. Combines the top 25% with the old 25% and then reranks and records the similarity to the old one.

'mutation_of_elite.py' pending

## Function scripts

'calculation_status.py' - monitors the calculation

'molecule_mutaion.py' - functions to mutate molecules