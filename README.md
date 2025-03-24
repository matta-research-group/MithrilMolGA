# GA_TSJ_project
Code base of the genetic algorithm project

## QCflow

This GA uses a package created by Tristan Stephens-Jones as part of the Matta Research Group. It uses QCflow and more specifically the beta version of the Psi4 version (https://github.com/matta-research-group/QCflow/tree/qcflow-psi4). All the calculations submitted using this GA use Psi4 and are not adapted for Gaussian16 as Psi4 is more accessible and was found to be faster and produce highly similar results compared to G16.

Tristan is the principle developer and any questions should be directed to him either via the issues feature of GitHub or via email.

## Notes
This contains the scripts that will be used within the GA. Alot of the file names are currently just place holders.

## Scripts
`monomer_run.py` this takes all the potential molecules, checks if the monomers have been run and if they haven't it runs them

`molecule_run.py` takes the monomer data and does D-A matching, if they are a good match it runs the geometry optimisation to get HOMO, LUMO, Energy Gap and Planarity

`GA_data_extraction_step.py` extracts the data from the molecule run step, this runs in th background and waits for all the molecules to finish

`GA_elite_step.py` ranks all the molecules in the current run, takes the top X% runs reorganisation energy calculations on them. Combines the top X% with the old X% and then reranks and records the similarity to the old one.

`mutation_of_elite.py` This step mutates the elite molecules be either changing the bio, non-bio or the linker. It then creates new molecules so that the amount of molecules for the next run matches the amount used in the previous run

`master_GA_script.py` This script overseas and runs all the previous scripts and is what is submitted to CREATE. This houses all the user defined parameters and is what should be adjusted by the user if they want certain filter parameters.

## Function Scripts

`calculation_status.py` - monitors the calculation

`molecule_mutaion.py` - functions to mutate molecules

## DataFrames and Dictionaries

`run_df.csv` - contains the similarity results of each run

`d_a_df.csv` - contains the D-A matching results

`monomer_df.csv` - contains data for all the monomers ran

`ran_all_data.csv` - contains data for all the molecules ran

`elite_run_0_df.csv` - empty df to start the GA

`bio_dic.json` - dictionary of all the bio monomers with attachment points denoted by `I`

`non_bio_dic.json` - dictionary of all the non-bio monomers with attachment points denoted by `I`

`linker_dic.json` - dictionary of all the linkers with attachment points denoted by `I`

`molecule_to_run_1.json` - dictionary of all molecules to run for the first run of the GA

## Progress Checking

When the GA is submitted a `GA_status.txt` file is created. This file tracks when each script has been successful completed and shows for which run it is completed for. This is helpful for the user to know where the GA is currently at.

## User Determined Values & Their Defaults

- functional = `b3lyp`
- basis_set = `6-31g*`
- time = `6` (hours) for calculation submissions within the GA
- number_of_cpus = `10` for calculation submissions within the GA
- EG_cutoff = `3.2` (eV) for energy gap cutoff 
- EG_rank_weight = `1` weighting for energy gap in the ranking
- planarity_rank_weight = `1` weighting for planarity in the ranking
- SA_rank_weight = `4` weighting for synthetic score (1/4 that of Eg and planairty)
- anioinc_reorg_rank_weight = `0.5` (worth double that of Eg and planarity and 8 times of SA score)
- cationic_reorg_rank_weight = `0.5` (worth double that of Eg and planarity and 8 times of SA score)
- elite_value = `25` (percent) percentage of molecules to be considered elite
- elite_df_size = `10` makes sure the elite df has a minimum size.


## Future Work

- Add a function to check similarity and to finish the GA if there is continued similarity
- Add the time to the `GA_status.txt` file so the user can see if it is stuck in a loop
- Adapt so to work with a workflow manager as need to be able to run for more than 2 days
- Change from pandas to polars for dataframes
- Add documentation