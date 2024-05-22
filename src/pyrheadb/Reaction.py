__author__ = "Anastasia Sveshnikova"
__email__ = "anastasia.sveshnikova@sib.swiss"
__status__ = "Prototype"

import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

from collections import Counter

RDLogger.DisableLog('rdApp.*')

class Reaction:
	
	def __init__(self):
		pass
	
	def read_rxnfile(self, rxnfile):
		"""
		Use rdkit to read .rxn file provided on Rhea FTP as a rdkit reaction (rdChemReaction)
		Args:
			rxnfile: .rxn file
		Returns:
			rdkit rxn object
		"""
		pattern = r'\n*\s*\$MOL'
		with open(rxnfile) as f:
			content = f.read()
			content = re.sub(pattern, '\n$MOL', content)
		with  open(rxnfile, 'w') as w:
			w.write(content)
		rxn = AllChem.ReactionFromRxnFile(rxnfile)
		try:
			rxn.Initialize()
		except:
			# Temporal hack - might generate recursive behaviour
			# TO FIX
			with open(rxnfile) as f:
				content = f.read()
				content = content.replace('R#','R')
			with open(rxnfile, 'w') as w:
				w.write(content)
			return self.read_rxnfile(rxnfile)
		return rxn
	
	def stoichiometric_equation(self, identifiers_equation):
		reactants = identifiers_equation.split('>>')[0].split('.')
		products = identifiers_equation.split('>>')[1].split('.')
		
		reactants = Counter(reactants)
		products = Counter(products)
		reactants = [str(reactants[cid]) + ' ' + cid for cid in reactants]
		products = [str(products[cid]) + ' ' + cid for cid in products]
		return ' + '.join(reactants) + ' <=> ' + ' + '.join(products)
	
	def check_reaction_balance(self, reaction_string, format='smiles'):
		"""
		Checks if the reaction represented by the reaction string is balanced.
		Supports both SMILES and InChI formats.

		Parameters:
		- reaction_string (str): The reaction string in either SMILES or InChI format.
		- format (str): The format of the reaction string ('smiles' or 'inchi').

		Returns:
		- bool: True if the reaction is balanced, False otherwise.
		"""
		if format == 'smiles':
			reactants, products = reaction_string.split(">>")
			reactant_split, product_split = '.', '.'
			mol_func = Chem.MolFromSmiles
		elif format == 'inchi':
			reactants, products = reaction_string.split(">>")
			reactant_split, product_split = ' + ', ' + '
			mol_func = Chem.MolFromInchi
		
		reactants = [reactant.strip() for reactant in reactants.split(reactant_split)]
		products = [product.strip() for product in products.split(product_split)]
		
		reactant_mols = [mol_func(reactant) for reactant in reactants]
		product_mols = [mol_func(product) for product in products]
		
		reactant_counts_total, product_counts_total = self.calculate_atom_counts(
			reactant_mols), self.calculate_atom_counts(product_mols)
		
		return reactant_counts_total == product_counts_total
	
	def calculate_atom_counts(self, mols):
		"""
		Calculates the total count of each atomic number in the given list of molecule objects.

		Parameters:
		- mols (list): List of molecule objects.

		Returns:
		- dict: Dictionary of atom counts keyed by atomic number.
		"""
		atom_counts_total = dict()
		for mol in mols:
			if mol is not None:  # Handle the case where molecule conversion fails
				atom_counts = {atom.GetAtomicNum(): mol.GetNumAtoms(atom.GetAtomicNum())
							   for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1}
				atom_counts_total = self.sum_dicts(atom_counts_total, atom_counts)
		return atom_counts_total
	
	def sum_dicts(self, dict1, dict2):
		result_dict = {}
		
		for key, value in dict1.items():
			result_dict[key] = result_dict.get(key, 0) + value
		
		for key, value in dict2.items():
			result_dict[key] = result_dict.get(key, 0) + value
		
		return result_dict
	# def check_reaction_balance_smiles(self, reaction_smiles):
	# 	reactants, products = reaction_smiles.split(">>")
	# 	reactants = [reactant.strip() for reactant in reactants.split('.')]
	# 	products = [product.strip() for product in products.split('.')]
	#
	# 	reactant_mols = [Chem.MolFromSmiles(reactant) for reactant in reactants]
	# 	product_mols = [Chem.MolFromSmiles(product) for product in products]
	#
	# 	reactant_counts_total = dict()
	# 	for reactant_mol in reactant_mols:
	# 		reactant_atoms = [atom.GetAtomicNum() for atom in reactant_mol.GetAtoms() if atom.GetAtomicNum() != 1]
	# 		reactant_counts = dict(zip(set(reactant_atoms), [reactant_atoms.count(atom) for atom in set(reactant_atoms)]))
	# 		reactant_counts_total = self.sum_dicts(reactant_counts_total, reactant_counts)
	#
	# 	product_counts_total = dict()
	# 	for product_mol in product_mols:
	# 		product_atoms = [atom.GetAtomicNum() for atom in product_mol.GetAtoms() if atom.GetAtomicNum() != 1]
	# 		product_counts = dict(zip(set(product_atoms), [product_atoms.count(atom) for atom in set(product_atoms)]))
	# 		product_counts_total = self.sum_dicts(product_counts_total, product_counts)
	# 	return reactant_counts == product_counts
	#
	# def check_reaction_balance(self, reaction_inchi):
	# 	reactants, products = reaction_inchi.split(">>")
	# 	reactants = [reactant.strip() for reactant in reactants.split(' + ')]
	# 	products = [product.strip() for product in products.split(' + ')]
	#
	# 	reactant_mols = [Chem.MolFromInchi(reactant) for reactant in reactants]
	# 	product_mols = [Chem.MolFromInchi(product) for product in products]
	#
	# 	reactant_counts_total = dict()
	# 	for reactant_mol in reactant_mols:
	# 		reactant_atoms = [atom.GetAtomicNum() for atom in reactant_mol.GetAtoms() if atom.GetAtomicNum() != 1]
	# 		reactant_counts = dict(
	# 			zip(set(reactant_atoms), [reactant_atoms.count(atom) for atom in set(reactant_atoms)]))
	# 		reactant_counts_total = self.sum_dicts(reactant_counts_total, reactant_counts)
	#
	# 	product_counts_total = dict()
	# 	for product_mol in product_mols:
	# 		product_atoms = [atom.GetAtomicNum() for atom in product_mol.GetAtoms() if atom.GetAtomicNum() != 1]
	# 		product_counts = dict(zip(set(product_atoms), [product_atoms.count(atom) for atom in set(product_atoms)]))
	# 		product_counts_total = self.sum_dicts(product_counts_total, product_counts)
	# 	return reactant_counts == product_counts



