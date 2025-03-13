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
		self.rhea_db_version_location = rhea_db.rhea_db_version_location
		self.rxn_mapper = AtomMapper()
		self.reactionobj = Reaction()
		self.df_smiles = rhea_db.df_reactions.copy()
		
	def convert_all_rhea_smiles_to_smarts(self, debug=False):
		"""
		Main function of the class:
		Converts SMILES to SMARTS and outputs the result into the text files
		:return: output into text files
		"""
		self.generate_rhea_smarts_for_all_generic_reactions(debug_mode=debug)
		self.write_smarts_files()
	
	def generate_rhea_smarts_for_all_generic_reactions(self, debug_mode=False):
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
		self.df_smiles['SMARTS_with_halogens'] = self.df_smiles.progress_apply(self.generate_smarts, axis=1, debug_mode=debug_mode)
		self.df_smiles['SMARTS'] = self.df_smiles.apply(self.clean_rxn_smarts_of_halogen_star, axis=1)
		self.df_smiles = self.df_smiles[~self.df_smiles['SMARTS'].isna()]
		
	def write_smarts_files(self):
		"""
		:param df_smiles: pandas DataFrame with SMARTS generated for each possible Rhea reaction
		:return: outputs files (I/O function)
		"""
		smarts_directory = os.path.join(self.rhea_db_version_location, 'smarts')
		os.makedirs(smarts_directory, exist_ok=True)
		
		# Remove unbalanced reactions
		self.df_smiles=self.df_smiles[self.df_smiles['SMARTS']!='unbalanced reaction']
		# Remove macromolecular reactions and pathway-reactions (lumped reactions)
		self.df_smiles=self.df_smiles[self.df_smiles['SMARTS']!='macromolecular reaction']
		# Remove not pattern reactions - reactions in which A / AH were the only compound with star
		self.df_smiles=self.df_smiles[self.df_smiles['SMARTS']!='Not pattern after removing A / AH']
		
		self.df_smiles.to_csv(os.path.join(smarts_directory, 'rheaSmarts.tsv'), sep='\t', index=False,
						 columns=['MASTER_ID', 'chebi_equation', 'SMARTS', 'rxnsmiles'])
		rhea_smarts_file = os.path.join(smarts_directory, 'rheaSmarts.json')
		rhea_smarts = dict(zip(self.df_smiles['MASTER_ID'], self.df_smiles['SMARTS']))
		
		with open(rhea_smarts_file, "w") as file:
			json.dump(rhea_smarts, file)
	
	def clean_rxn_smiles_of_redox(self, row):
		rxnsmiles=row['rxnsmiles']
		return self.reactionobj.clean_rxnsmiles_of_redox(rxnsmiles)

	def clean_rxn_smarts_of_halogen_star(self, row):
		"""
		CHEBI:16042 - halide anion
		"""
		if 'CHEBI:16042' in row['chebi_equation']:
			return self.replace_isolated_star_with_halogens(row['SMARTS_with_halogens'])
		return row['SMARTS_with_halogens']

	def replace_isolated_star_with_halogens(self, equation):
		# Identify all isolated atoms from both sides of the equation
		isolated_atoms = self.find_isolated_atoms(equation)

		# Replace isolated atoms in the equation based on identified atom numbers
		for number in isolated_atoms:
			equation = equation.replace(f"[*:{number}]", f"[F,Cl,Br,I:{number}]")
			equation = equation.replace(f"[*-:{number}]", f"[F,Cl,Br,I-:{number}]")
		return equation
	
	def find_isolated_atoms(self, equation):
		# Split the equation into left and right sides based on ">>"
		left_side, right_side = equation.split(">>")
		
		# Identify isolated atoms on both sides using regex
		isolated_left = set(re.findall(r'\[\*\-\:(\d+)\]', left_side))
		isolated_right = set(re.findall(r'\[\*\-\:(\d+)\]', right_side))
		
		# Combine all unique isolated atoms from both sides
		all_isolated_atoms = isolated_left.union(isolated_right)
		return all_isolated_atoms

	def remove_free_hydrogen_pattern(self, smarts):
		# Define the regular expression pattern to match "[F:*]"
		pattern = r"\.\[H\+:\d{1,3}\]$"
		
		# Use the re.sub() function to replace the pattern with an empty string
		smarts = re.sub(pattern, "", smarts)
		smarts = smarts.replace('.[H+]', '')
		return smarts
	
	def generate_smarts(self, row, debug_mode=False):
		"""
		:param row: Row of dataframe of rhea reaction smiles
		:return: SMARTS of Rhea reaction with Atom Mapping numbers that can be in the next step converted to rdkit reactions
		"""
		
		balance = self.reactionobj.check_reaction_balance(row['rxnsmiles'], format='smiles', account_H=False, debug_mode=debug_mode)
		
		if not balance:
			return 'unbalanced reaction'
		
		rxn_smiles = self.reactionobj.remove_A_AH_pattern(row['rxnsmiles'])
		
		if not '*' in rxn_smiles:
			return 'Not pattern after removing A / AH'
		
		rxn_smiles = self.reactionobj.star_to_isotopic_label(rxn_smiles)
		
		atommap_results = self.rxn_mapper.map_one_reaction(rxn_smiles)
		
		# Exclude reactions that are not likely to be sensible patterns -
		# star in those reactions likely defines a protein or other macromolecule - not small molecule part
		# Besides, such reactions cannot be processed with atom mapper
		if atommap_results[0] == 'tokenlength_error':
			return 'macromolecular reaction'
		
		smarts = atommap_results[0]
		if type(smarts)!=str:
			print(atommap_results)
			return None
		smarts = self.isotope_to_star_pattern(smarts)
		smarts = self.remove_free_hydrogen_pattern(smarts)
		return smarts

	
	def isotope_to_star_pattern(self, rxnsmiles):
		"""
		Remove the pattern that was previously introduced to overcome inability of RXNmapper to handle dummy atoms (*)
		:param input_string: SMARTS string
		:return:
		"""
		return rxnsmiles.replace('13C', '*')
	
	def smarts_to_smiles(self, rxn_smarts):

		rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
		rdChemReactions.SanitizeRxn(rxn)
		rxn_smiles = rdChemReactions.ReactionToSmiles(rxn)
		return rxn_smiles

