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

#Open the data df
#Sort top EG, top Plan, top SA
#Submit top 25% for reorganisation calculation (need to do some Psi4 stuff)
#Can definitely be all done within a few psi4 calculations rather than the gaussian splitting across multiple files
#Retreave reorganisation data
#Make a new df with the reorganisation data and then combine with old elite 25% from previous step and then sort again