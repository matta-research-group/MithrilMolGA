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
import re
import itertools
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
import sys
import os
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
from SA_Score import sascorer

# This file will first test if the calculations have been ran successfully
# Split into two dictionaries, one with the successful and another with unsuccessful
# The successful calculations will have their planarity, energy gap extracted and their syenthic accessibility score calculated
# The unsuccessful calculations will either be ran again or discarded (to be decided)
# Once all the data has been extracted they will be ordered into an elitism step
# With the top 25% undergo reorganisation calculations


data = extract_data_from_txt(file_path)
    #energy calcs
mol_homo = data['homo']
mol_lumo = data['lumo']
mol_EG = data['energy_gap']
    #planarity data
mol_plan = finding_planairty_psi4(mol_name, molecule_study, linker_type, 'opt')