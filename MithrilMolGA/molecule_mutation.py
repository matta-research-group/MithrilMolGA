import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from CombineMols.CombineMols import CombineMols
import QCflow
from QCflow.load_gaussian import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *
import re
import random
from datetime import datetime

def find_linker_type(mol):
    """
    Determines the type of linker present in a given molecule.

    Parameters
    ----------
    mol (rdkit.Chem.Mol): The molecule to analyze.

    Returns
    -------
    str: The type of linker found in the molecule. Possible values are:
        - 'single': Single bond linker
        - 'double': Double bond linker
        - 'imine': Imine linker
        - 'thio': Thiophene linker
        - 'triple': Triple bond linker

    Note
    ----
    The function uses SMARTS patterns to identify the linker types.
    """
    # old thio = [R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]

    linker_type_smarts = {
    'single' : '[R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]',
    'double' : '[R!$(*#*)&!D1]C=C-!@[R!$(*#*)&!D1]',
    'imine' : '[R!$(*#*)&!D1]N=C-!@[R!$(*#*)&!D1]',
    'thio' : '[#6]-[#6]1:[#6]:[#6](-[#8]-[#6]):[#6](-[#6]):[#16]:1',
    'triple' : '[*R1!$(*#*)!D1]C#C-!@[*R1!$(*#*)!D1]'}

    for k, v in linker_type_smarts.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(v)):
            linker_type =  k

    if k == 'single':
        thio_test = mol.HasSubstructMatch(Chem.MolFromSmarts(linker_type_smarts['thio']))
        if thio_test == True:
            linker_type = 'thio'

    return linker_type


def find_fragment_type(mol_smi, bio_dic):
    """
    Determines the type of a molecular fragment based on its SMILES representation.
    Parameters
    ----------
    mol_smi (str): The SMILES string of the molecule to be classified.
    bio_dic (dict): A dictionary where keys are identifiers and values are SMILES strings of known bio fragments.

    Returns
    -------
    str: 'bio' if the molecule matches any SMILES in bio_dic, otherwise 'non_bio'.
    """

    fragment_one_type = []
    for k1, v1 in bio_dic.items():
        match_found = False
        if Chem.CanonSmiles(v1, useChiral=0) == Chem.CanonSmiles(mol_smi, useChiral=0):
            match_found = True
            fragment_one_type.append('bio')
            break
    
    if not fragment_one_type:
        fragment_one_type.append('non_bio')

    return fragment_one_type[0]

def fragment_molecule(mol, linker):
    """
    Fragments a molecule at specified linker bonds and returns the SMILES representation of the fragments.
    
    Parameters
    ----------
    mol (rdkit.Chem.Mol): The molecule to be fragmented.
    linker (str): The type of linker bond to fragment. Supported values are 'double', 'imine', 'triple', 'thio', and 'single'.

    Returns
    -------
    list: A list of SMILES strings representing the fragments with attachment points replaced by a generic atom '[I]'.

    Notes
    -----
    - For 'double', 'imine', 'triple', and 'thio' linkers, the molecule is fragmented at two bonds.
    - For 'single' linkers, the molecule is fragmented at one bond.
    - The function uses regular expressions to replace attachment points with a generic atom '[I]'.
    """

    linker_bonds = getBondLinkers(mol, linker) #finding atoms involved in bond

    fragmented_decorated_smi = []

    atoms_in_bond = []
    if (linker == 'double') or (linker == 'imine') or (linker == 'triple'):
        bond = linker_bonds[0]
        atoms_in_bond.append(bond)

    if (linker == 'thio'):
        bond = linker_bonds[0] + linker_bonds[1]
        atoms_in_bond.append(bond)
    
    if (linker == 'single'):
        bond = linker_bonds[0]
        atoms_in_bond.append(bond)

    atoms_in_bond_cut = atoms_in_bond[0]

    if (linker == 'double') or (linker == 'imine') or (linker == 'triple') or (linker == 'thio'): 

        bond_1_idx = mol.GetBondBetweenAtoms(atoms_in_bond_cut[0], atoms_in_bond_cut[1]).GetIdx() #index of bond 1

        fragment_1 = Chem.FragmentOnBonds(mol, (bond_1_idx,)) #fragment with bond 1 cut and linker still attached

        bond_2_idx = fragment_1.GetBondBetweenAtoms(atoms_in_bond_cut[2], atoms_in_bond_cut[3]).GetIdx() #index of bond 2

        fragment_2 = Chem.FragmentOnBonds(fragment_1, (bond_2_idx,)) #final fragment

        # getting smi of the fragments
        fragmented_smi = Chem.MolToSmiles(fragment_2)

        splitting = fragmented_smi.split('.') #splitting up the fragments

        #could remove linker via --> [smi for smi in splitting if smi.count('*') == 1]

        for smi in splitting:
            modified_smi = re.sub(r'\[\d+\*\]', '[I]', smi) #replacing the attachment point with a generic atom
            fragmented_decorated_smi.append(modified_smi)

    if (linker == 'single'):

        bond_1_idx = mol.GetBondBetweenAtoms(atoms_in_bond_cut[0], atoms_in_bond_cut[1]).GetIdx()

        fragment_1 = Chem.FragmentOnBonds(mol, (bond_1_idx,))

        # getting smi of the fragments
        fragmented_smi = Chem.MolToSmiles(fragment_1)

        splitting = fragmented_smi.split('.') #splitting up the fragments

        #could remove linker via --> [smi for smi in splitting if smi.count('*') == 1]

        for smi in splitting:
            modified_smi = re.sub(r'\[\d+\*\]', '[I]', smi) #replacing the attachment point with a generic atom
            fragmented_decorated_smi.append(modified_smi)

    return fragmented_decorated_smi

def replace_linker(fragments, linker_dic):
    """
    Replaces the linker in a fragmented molecule with a new linker from the provided dictionary.

    Parameters
    ----------
    fragments (list): A list of SMILES strings representing the fragments of the molecule.
    linker_dic (dict): A dictionary of possible linkers to replace the old linker.

    Returns
    -------
    str: A SMILES string representing the new molecule with the replaced linker.

    Notes
    -----
    - The function randomly selects a new linker from the dictionary that is not the same as the old linker.
    - The function combines the fragments with the new linker to form the new molecule.
    """

    # List of the fragments with attachment points
    molecule_fragments = [smi for smi in fragments if smi.count('I') == 1]
    # List of the old linker
    old_linker = [smi for smi in fragments if smi.count('I') == 2]
    # RDKit mol object of old linker
    if not old_linker:
        old_linker_canon = Chem.CanonSmiles('II')
    else:
        old_linker_canon = Chem.CanonSmiles(old_linker[0])
    # Randomly selecting a linker that is not the old linker
    random_linker = random.choice([v for k, v in linker_dic.items() if Chem.CanonSmiles(v, useChiral=0) != old_linker_canon])
    # Combining the fragments with the random linker
    frag_linker = combine_structure(molecule_fragments[0], random_linker)
    # Combine the remaining fragment with the old linker that is attached to the other fragment
    frag_linker_frag = combine_structure(molecule_fragments[1], frag_linker)
    # The new linker is added to a dictionary but doesn't need to be like this
    return frag_linker_frag


def find_replacement_fragment(fragments, bio_dic, non_bio_dic, fragment_replace):
    """
     Chooses which fragment to replace.

    Parameters
    ----------
    fragments (list): The list of fragmented molecule SMILES strings.
    bio_dic (dict): The dictionary of bio-inspired fragments.
    non_bio_dic (dict): The dictionary of non-bio-inspired fragments.
    fragment_replace (str): The type of fragment to replace. Possible values are 'bio' or 'non_bio'

    Returns
    -------
    dict: A dictionary containing the old fragment, the new fragment, and the type of fragment.
    """
    # List of the fragments with one attachment point
    fragments_one_attach = [smi for smi in fragments if smi.count('I') == 1]
    # List of the old linker
    linkage = [smi for smi in fragments if smi.count('I') == 2]
    # Remove the I from the fragment
    #fragment_one = re.sub(r'\[I\]', '', fragments_one_attach[0])
    #fragment_two = re.sub(r'\[I\]', '', fragments_one_attach[1])
    fragment_one = re.sub(fragments_one_attach[0])
    fragment_two = re.sub(fragments_one_attach[1])

    # Compare fragments to bio_dictionary and non-bio_dictionary
    fragment_types = {}
    fragment_smi = {}
    if find_fragment_type(fragment_one, bio_dic) == 'True':
        fragment_types['fragment_one'] = 'bio'
        fragment_smi['fragment_one'] = fragment_one
    else:
        fragment_types['fragment_one'] = 'non_bio'
        fragment_smi['fragment_one'] = fragment_one
    
    if find_fragment_type(fragment_two, bio_dic) == 'True':
        fragment_types['fragment_two'] = 'bio'
        fragment_smi['fragment_two'] = fragment_two
    else:
        fragment_types['fragment_two'] = 'non_bio'
        fragment_smi['fragment_two'] = fragment_two

    swap_fragment = {}
    # Randomly decide which fragment to replace
    random_replace = random.choice(list(fragment_types.keys()))
    # Store the old fragment
    swap_fragment[f'old_{fragment_types[random_replace]}'] = fragment_smi[random_replace]

    if fragment_types[random_replace] == 'bio':
        old_fragment_connon = Chem.CanonSmiles(fragment_smi[random_replace], useChiral=0)
        #pick a random new fragment but it cannot be the same as the old one
        new_fragment_name = random.choice([k for k, v in bio_dic.items() if Chem.CanonSmiles(v, useChiral=0) != old_fragment_connon])
        swap_fragment['new_bio'] = bio_dic[new_fragment_name]
    else:
        old_fragment_connon = Chem.CanonSmiles(fragment_smi[random_replace])
        #pick a random new fragment but it cannot be the same as the old one
        new_fragment_name = random.choice([k for k, v in non_bio_dic.items() if Chem.CanonSmiles(v, useChiral=0) != old_fragment_connon])
        swap_fragment['new_non_bio'] = non_bio_dic[new_fragment_name]

    return swap_fragment

def replace_fragment(new_fragment, old_fragment, fragments):
    """
    Replaces a fragment in a list of fragments with a new fragment.

    Parameters
    ----------
    new_fragment (str): The new fragment to replace the old fragment.
    old_fragment (str): The old fragment to be replaced.
    fragments (list): The list of fragments to be modified.

    Returns
    -------
    list: A list of SMILES strings representing the fragments with the old fragment replaced by the new fragment.

    Notes
    -----
    - The function identifies the fragment to be replaced and the fragment to be kept.
    - The new fragment is modified to include an attachment point.
    - The function combines the constant fragment and the new fragment with the linker to form new molecules.
    """

    # List of the fragments with one attachment point
    fragments_one_attach = [smi for smi in fragments if smi.count('I') == 1]
    # List of the old linker
    linkage = [smi for smi in fragments if smi.count('I') == 2]
    # Remove the I from the fragment
    fragment_one = re.sub(r'\[I\]', '', fragments_one_attach[0])
    fragment_two = re.sub(r'\[I\]', '', fragments_one_attach[1])

    # Keeping fragment that is not being replaced
    if Chem.CanonSmiles(fragment_one, useChiral=0) != Chem.CanonSmiles(old_fragment, useChiral=0):
        constant_fragment = fragments_one_attach[0]
    else:
        constant_fragment = fragments_one_attach[1]

    # Adding attachment point to new fragment
    new_fragment_w_attach = adding_attach(new_fragment, find='[cH;^2]', get_rid='C([I])')

    new_molecules = []

    constant_frag_and_linker = combine_structure(constant_fragment, linkage[0])

    for fragment in new_fragment_w_attach:
        constant_frag_linker_new_frag = combine_structure(constant_frag_and_linker, fragment)
        new_molecules.append(constant_frag_linker_new_frag)

    return new_molecules


def swap_one_fragment(mol_smi, bio_dic, non_bio_dic, fragment_replace):
    """
    Replaces one of the fragments depending on the users input.

    Parameters
    ----------
    fragments (list): The list of fragmented molecule SMILES strings.
    bio_dic (dict): The dictionary of bio-inspired fragments.
    non_bio_dic (dict): The dictionary of non-bio-inspired fragments.
    fragment_replace (str): The type of fragment to replace. Possible values are 'bio' or 'non_bio'

    Returns
    -------
    str: smi str of the new molecule.
    """
    #find the linker type
    linker_type = find_linker_type(Chem.MolFromSmiles(mol_smi))
    #fragment the molecule
    fragments = fragment_molecule(Chem.MolFromSmiles(mol_smi), linker_type)
    # List of the fragments with one attachment point
    fragments_one_attach = [smi for smi in fragments if smi.count('I') == 1]
    # List of the old linker
    linkage = [smi for smi in fragments if smi.count('I') == 2]

    #if the linker is a single bond it is unable to separate it so have to account for it
    if len(linkage) == 0:
        linkage = ['II']
    # Remove the I from the fragment
    fragment_one = re.sub(r'\[I\]', '', fragments_one_attach[0])
    fragment_two = re.sub(r'\[I\]', '', fragments_one_attach[1])

    # Compare fragments to bio_dictionary and non-bio_dictionary
    fragment_one_type = find_fragment_type(fragments_one_attach[0], bio_dic)
    fragment_two_type = find_fragment_type(fragments_one_attach[1], bio_dic)
    
    replaced_smi = []
    if fragment_replace == 'bio':
        if fragment_one_type == 'bio':
            replaced_smi.append(fragments_one_attach[0])
            random_key = random.choice(list(bio_dic.keys()))
            random_value = bio_dic[random_key]
            while Chem.CanonSmiles(random_value, useChiral=0) == Chem.CanonSmiles(fragments_one_attach[0], useChiral=0): #cannot be the same fragmeant over and over again
                random_key = random.choice(list(bio_dic.keys()))
                random_value = bio_dic[random_key]
        else:
            replaced_smi.append(fragments_one_attach[1])
            random_key = random.choice(list(bio_dic.keys()))
            random_value = bio_dic[random_key]
            while Chem.CanonSmiles(random_value, useChiral=0) == Chem.CanonSmiles(fragments_one_attach[1], useChiral=0): #cannot be the same fragmeant over and over again
                random_key = random.choice(list(bio_dic.keys()))
                random_value = bio_dic[random_key]

    if fragment_replace == 'non_bio':
        if fragment_one_type == 'non_bio':
            replaced_smi.append(fragments_one_attach[0])
            random_key = random.choice(list(non_bio_dic.keys()))
            random_value = non_bio_dic[random_key]
            while Chem.CanonSmiles(random_value, useChiral=0) == Chem.CanonSmiles(fragments_one_attach[0], useChiral=0): #cannot be the same fragmeant over and over again
                random_key = random.choice(list(non_bio_dic.keys()))
                random_value = non_bio_dic[random_key]
        else:
            replaced_smi.append(fragments_one_attach[1])
            random_key = random.choice(list(non_bio_dic.keys()))
            random_value = non_bio_dic[random_key]
            while Chem.CanonSmiles(random_value, useChiral=0) == Chem.CanonSmiles(fragments_one_attach[1], useChiral=0): #cannot be the same fragmeant over and over again
                random_key = random.choice(list(non_bio_dic.keys()))
                random_value = non_bio_dic[random_key]

    if Chem.CanonSmiles(replaced_smi[0], useChiral=0) == Chem.CanonSmiles(fragments_one_attach[0], useChiral=0):
        new_fragment = combine_structure(linkage[0], fragments_one_attach[1])
        new_molecule = combine_structure(new_fragment, random_value)
        #new_fragment = fragments_one_attach[1]
        #new_fragment_2 = random_value
        #new_linker = linkage
    else:
        new_fragment = combine_structure(linkage[0], fragments_one_attach[0])
        new_molecule = combine_structure(new_fragment, random_value)
        #new_fragment = fragments_one_attach[0]
        #new_fragment_2 = random_value
        #new_linker = linkage

    return new_molecule
