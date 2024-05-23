import json
import os
import re
from tqdm import tqdm

from .AtomMapper import AtomMapper
from .Reaction import Reaction

from rdkit.Chem import rdChemReactions

tqdm.pandas()

class ReactionSmartsConverter():
	
	def __init__(self, rhea_db):
		"""
		Initializes the ReactionSmarts class which handles conversion of reaction SMILES to SMARTS.
		ReactionSmarts class handles everything related to reaction SMARTS
		- atom mapped version of reactions that can act as reaction patterns.
		The purpose of this class is to provide interface between Rhea reactions
		as structures (originally stored and .sdf and converted to SMILES by pyrheadb)
		and reaction patters as SMARTS that can be used by rdkit reaction class to generate
		new reactions.
		
		"""
		self.RDBv_loc = rhea_db.RDBv_loc
		self.rxn_mapper = AtomMapper()
		self.df_smiles = rhea_db.df_smiles_master_id.copy()
		
	def convert_all_rhea_smiles_to_smarts(self):
		"""
		Main function of the class:
		Converts SMILES to SMARTS and outputs the result into the text files
		:return: output into text files
		"""
		self.generate_rhea_smarts_for_all_generic_reactions()
		self.write_smarts_files()
	
	def generate_rhea_smarts_for_all_generic_reactions(self):
		"""
		Add SMARTS column to the copy of the pandas DataFrame with SMILES and Rhea MASTER_ID
		:param df_smiles: pandas df with MASTER_ID (rhea id no direction) and rxn (reaction smiles) columns
		"""
		
		# Cleaning the DataFrame of reactions not fitting into the SMARTS generation
		self.df_smiles['rxn_clean_redox'] = self.df_smiles.apply(self.clean_rxn_smiles_of_redox, axis=1)
		self.df_smiles['star_count'] = self.df_smiles['rxn_clean_redox'].apply(lambda x: x.count('*'))
		
		self.df_smiles = self.df_smiles[self.df_smiles['star_count'] > 0]
		self.df_smiles = self.df_smiles[
			~self.df_smiles['star_count'].isin([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 39, 41])]
		
		print('Converting reaction SMILES into SMARTS')
		self.df_smiles['SMARTS'] = self.df_smiles.progress_apply(self.generate_smarts, axis=1)
		
	def write_smarts_files(self):
		"""
		:param df_smiles: pandas DataFrame with SMARTS generated for each possible Rhea reaction
		:return: outputs files (I/O function)
		"""
		smarts_directory = os.path.join(self.RDBv_loc, 'smarts')
		os.makedirs(smarts_directory, exist_ok=True)
		self.df_smiles.to_csv(os.path.join(smarts_directory, 'rheaSmarts.tsv'), sep='\t', index=False,
						 columns=['MASTER_ID', 'rheaid', 'template_type', 'chebi_equation', 'SMARTS', 'rxnsmiles'])
		rhea_smarts_file = os.path.join(smarts_directory, 'rheaSmarts.json')
		rhea_smarts = dict(zip(self.df_smiles['MASTER_ID'], self.df_smiles['SMARTS']))
		
		with open(rhea_smarts_file, "w") as file:
			json.dump(rhea_smarts, file)
	
	def remove_molecules(self, substrates, products, molecule, count):
		for _ in range(count):
			substrates.remove(molecule)
			products.remove('*')
	
	def clean_rxn_smiles_of_redox(self, row):
		substrates = row['rxn'].split('>>')[0].split('.')
		products = row['rxn'].split('>>')[1].split('.')
		self.remove_molecules(substrates, products, '[H]*[H]', min(substrates.count('[H]*[H]'), products.count('*')))
		self.remove_molecules(products, substrates, '[H]*[H]', min(products.count('[H]*[H]'), substrates.count('*')))
		return '.'.join(substrates) + '>>' + '.'.join(products)
	
	def remove_free_hydrogen_pattern(self, input_string):
		# Define the regular expression pattern to match "[F:*]"
		pattern = r"\.\[H\+:\d{1,3}\]$"
		
		# Use the re.sub() function to replace the pattern with an empty string
		result_string = re.sub(pattern, "", input_string)
		result_string = result_string.replace('.[H+]', '')
		return result_string
	
	def remove_A_AH_pattern(self, result_string):
		
		oxidant_replacement = {
			r"\.\*\.": '.',
			r"^\*\.": '',
			r">>\*\.": '>>',
			r"\.\*>>": '>>',
			r"\.\*$": '',
		}
		# Oxidant patterns
		for pattern, replacement in oxidant_replacement.items():
			found_pattern = re.findall(pattern, result_string)
			if found_pattern:
				for pattern in found_pattern:
					result_string = result_string.replace(pattern, replacement)
		
		# Reducing patterns
		
		reducing_replacement = {
			r"\.\[H\]\*\[H\]\.": '.',
			r"^\[H\]\*\[H\]\.": '',
			r">>\[H\]\*\[H\]\.": '>>',
			r"\.\[H\]\*\[H\]>>": '>>',
			r"\.\[H\]\*\[H\]$": ''
		}
		for pattern, replacement in reducing_replacement.items():
			found_pattern = re.findall(pattern, result_string)
			if found_pattern:
				for pattern in found_pattern:
					result_string = result_string.replace(pattern, replacement)
		
		result_string = self.remove_extra_dots_smiles(result_string)
		return result_string
	
	def remove_extra_dots_smiles(self, input_string):
		reactants = [r for r in input_string.split('>>')[0].split('.') if r]
		products = [p for p in input_string.split('>>')[1].split('.') if p]
		return '.'.join(reactants) + '>>' + '.'.join(products)
	
	def generate_smarts(self, row):
		"""
		:param row: Row of dataframe of rhea reaction smiles
		:return: SMARTS of Rhea reaction with Atom Mapping numbers that can be in the next step converted to rdkit reactions
		"""
		# Logic from another impementation of this
		# substrates = row['rxn_clean_redox'].split('>>')[0].split('.')
		# products = row['rxn_clean_redox'].split('>>')[1].split('.')
		# substrates = '.'.join([i for i in substrates if '*' in i])
		# products = '.'.join([i for i in products if '*' in i])
		#
		# rxns = [substrates + '>>' + products]
		#
		# if len(rxns[0]) < 512:
		#
		# 	if not all(['*' in i for i in row['rxn'].split('>>')]):
		# 		# print('unbalanced R', row['MASTER_ID'])
		# 		return 'unbalanced R'
		
		reaction = Reaction()
		
		balance = reaction.check_reaction_balance(row['rxn'])
		if not balance:
			return 'unbalanced reaction'
		
		rxn_smiles = self.remove_A_AH_pattern(row['rxn'])
		
		if len(row['rxn']) > 512:
			return 'reaction too big for the atom mapper', None
		
		#if '*' in row['rxn'] and row['MASTER_ID'] not in (63388):
		
		# Rhea reaction 63388 : A +  Fe(II)-heme o + H2O = AH2 + Fe(II)-heme a is not processable using RXNmapper
		# A - AH2 - too generic, Fe(II)-heme o - too big molecule
		
		rxns = [self.modify_rxn_smiles(row['rxn'])]
		atommap_results = self.rxn_mapper.map_many_reactions(rxns)
		smarts = atommap_results['mapped_rxn']
		smarts = self.remove_pattern(smarts)
		smarts = self.remove_free_hydrogen_pattern(smarts)
		return smarts

	
	def modify_rxn_smiles(self, rxn_smiles):
		"""
		Modifies the reaction SMILES to replace placeholders with isotopically labelled carbon.
		This is necessary to be able to handle reaction SMILES with RXNMapper
		:param rxn_smiles: The reaction SMILES string
		:return: Modified SMILES string
		"""
		replacements = {
			'[1*]': '[13C]', '[2*]': '[13C]', '[3*]': '[13C]', '[4*]': '[13C]', '[5*]': '[13C]',
			'[6*]': '[13C]', '[7*]': '[13C]', '[8*]': '[13C]', '[9*]': '[13C]', '[*-]': '[13C-]',
			'[*:0]': '[13C]', '*': '[13C]'
		}
		for old, new in replacements.items():
			rxn_smiles = rxn_smiles.replace(old, new)
		return rxn_smiles
	
	def remove_pattern(self, input_string):
		"""
		Remove the pattern that was previously introduced to overcome inability of RXNmapper to handle dummy atoms (*)
		:param input_string: SMARTS string
		:return:
		"""
		result_string = input_string.replace('13C', '*')
		return result_string
	
	def smarts_to_smiles(self, rxn_smarts):

		rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
		rdChemReactions.SanitizeRxn(rxn)
		rxn_smiles = rdChemReactions.ReactionToSmiles(rxn)
		return rxn_smiles
