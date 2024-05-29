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
	
	def check_reaction_balance(self, reaction_string, format='smiles', account_H=False, debug_mode=False):
		"""
		Checks if the reaction represented by the reaction string is balanced.
		Supports both SMILES and InChI formats.

		Parameters:
		- reaction_string (str): The reaction string in either SMILES or InChI format.
		- format (str): The format of the reaction string ('smiles' or 'inchi').
		- account_H (bool): if H is counted in the balance
		- debug_mode (bool): if the current run is to debug and print intermediate values

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
		
		if debug_mode==True:
			print('reactants', reactants)
			print('products', products)
		reactant_mols = [mol_func(reactant) for reactant in reactants]
		product_mols = [mol_func(product) for product in products]
		
		reactant_counts_total, product_counts_total = self.calculate_atom_counts(
			reactant_mols, debug_mode=debug_mode), self.calculate_atom_counts(product_mols, debug_mode=debug_mode)
		
		# Remove hydrogens from count - option
		if account_H==False and 1 in reactant_counts_total.keys():
			del reactant_counts_total[1]
		if account_H==False and 1 in product_counts_total.keys():
			del product_counts_total[1]
		return reactant_counts_total == product_counts_total
	
	def calculate_atom_counts(self, mols, debug_mode=False):
		"""
		Calculates the total count of each atomic number in the given list of molecule objects.

		Parameters:
		- mols (list): List of molecule objects.

		Returns:
		- dict: Dictionary of atom counts keyed by atomic number.
		"""
		atom_counts_total = dict()
		for mol in mols:
			if debug_mode == True:
				print(Chem.MolToSmiles(mol))
			if mol is not None:  # Handle the case where molecule conversion fails
				atom_counts = dict()
				for atom in mol.GetAtoms():
					if atom.GetAtomicNum() not in atom_counts:
						atom_counts[atom.GetAtomicNum()]=0
					atom_counts[atom.GetAtomicNum()]+=1
				if debug_mode == True:
					print(atom_counts)
				atom_counts_total = self.sum_dicts(atom_counts_total, atom_counts)
		return atom_counts_total
	
	def sum_dicts(self, dict1, dict2):
		result_dict = {}
		
		for key, value in dict1.items():
			result_dict[key] = result_dict.get(key, 0) + value
		
		for key, value in dict2.items():
			result_dict[key] = result_dict.get(key, 0) + value
		
		return result_dict
	
	def clean_rxnsmiles_of_redox(self, rxnsmiles):
		substrates = rxnsmiles.split('>>')[0].split('.')
		products = rxnsmiles.split('>>')[1].split('.')
		self.remove_molecules(substrates, products, '[H]*[H]', min(substrates.count('[H]*[H]'), products.count('*')))
		self.remove_molecules(products, substrates, '[H]*[H]', min(products.count('[H]*[H]'), substrates.count('*')))
		return '.'.join(substrates) + '>>' + '.'.join(products)
	
	def remove_molecules(self, substrates, products, molecule, count):
		for _ in range(count):
			substrates.remove(molecule)
			products.remove('*')
	
	def remove_A_AH_pattern(self, rxnsmiles):
		
		oxidant_replacement = {
			r"\.\*\.": '.',
			r"^\*\.": '',
			r">>\*\.": '>>',
			r"\.\*>>": '>>',
			r"\.\*$": '',
		}
		# Oxidant patterns
		for pattern, replacement in oxidant_replacement.items():
			found_pattern = re.findall(pattern, rxnsmiles)
			if found_pattern:
				for pattern in found_pattern:
					rxnsmiles = rxnsmiles.replace(pattern, replacement)
		
		# Reducing patterns
		
		reducing_replacement = {
			r"\.\[H\]\*\[H\]\.": '.',
			r"^\[H\]\*\[H\]\.": '',
			r">>\[H\]\*\[H\]\.": '>>',
			r"\.\[H\]\*\[H\]>>": '>>',
			r"\.\[H\]\*\[H\]$": ''
		}
		for pattern, replacement in reducing_replacement.items():
			found_pattern = re.findall(pattern, rxnsmiles)
			if found_pattern:
				for pattern in found_pattern:
					rxnsmiles = rxnsmiles.replace(pattern, replacement)
		
		rxnsmiles = self.remove_extra_dots_smiles(rxnsmiles)
		return rxnsmiles
	
	def remove_extra_dots_smiles(self, rxnsmiles):
		reactants = [r for r in rxnsmiles.split('>>')[0].split('.') if r]
		products = [p for p in rxnsmiles.split('>>')[1].split('.') if p]
		return '.'.join(reactants) + '>>' + '.'.join(products)
	

	def star_to_isotopic_label(self, rxn_smiles):
		"""
		Modifies the reaction SMILES to replace placeholders with isotopically labelled carbon.
		This is necessary to be able to handle reaction SMILES with RXNMapper.
		The order of replacements is important! '*': '[13C]' should be the last.
		:param rxn_smiles: The reaction SMILES string
		:return: Modified SMILES string
		"""
		replacements = {
			'[1*:0]': '[13C]', '[1*]': '[13C]', '[2*]': '[13C]', '[3*]': '[13C]', '[4*]': '[13C]', '[5*]': '[13C]',
			'[6*]': '[13C]', '[7*]': '[13C]', '[8*]': '[13C]', '[9*]': '[13C]', '[*-]': '[13C-]',
			'[*:0]': '[13C]', '*': '[13C]'
		}
		for old, new in replacements.items():
			rxn_smiles = rxn_smiles.replace(old, new)
		return rxn_smiles


