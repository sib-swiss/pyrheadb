import os.path

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
import io
import json
from tqdm import tqdm
from PIL import Image

tqdm.pandas()

from .AtomMapper import AtomMapper
from .RheaDB import RheaDB

class ReactionSmarts(RheaDB):
	
	def __init__(self):
		"""
		ReactionSmarts class handles everything related to reaction SMARTS
		- atom mapped version of reactions that can act as reaction patterns
		"""
		pass
	
	def generateRheaSmartsForAllGenericReactions(self, df):
		"""
		
		:param df: pandas df with MASTER_ID (rhea id no direction) and rxn (reaction smiles) columns
		:return:
		"""
		if not os.path.exists(self.RDBv_loc+'/smarts'):
			os.mkdir(self.RDBv_loc+'/smarts')
		rheaSmartsFile = self.RDBv_loc+'/smarts/rheaSmarts.json'
		self.rxn_mapper = AtomMapper()
		
		print('Converting reaction SMILES into SMARTS')
		df['SMARTS'] = df.progress_apply(self.generateSmarts, axis=1)
		
		rheaSmarts = dict(zip(df['MASTER_ID'], df['SMARTS']))
		
		with open(rheaSmartsFile, "w") as file:
			json.dump(rheaSmarts, file)
		
		return rheaSmarts
	
	def generateSmarts(self, row):
		"""
		:param row: row of dataframe of rhea reaction smiles
		:return: Smarts of Rhea reaction with Atom Mapping numbers that can be in the next step converted to rdkit reactions
		"""
		if '*' in row['rxn'] and len(row['rxn']) < 512 and row['MASTER_ID'] not in (63388):
			# Rhea reaction 63388 : A +  Fe(II)-heme o + H2O = AH2 + Fe(II)-heme a is not processable using RXNmapper
			# A - AH2 - too generic, Fe(II)-heme o - too big molecule
			
			rxns = [row['rxn']]
			
			rxns = [i.replace('[1*]', '[13C]').replace('[2*]', '[13C]').replace('[3*]', '[13C]')
					.replace('[4*]', '[13C]').replace('[5*]', '[13C]').replace('[6*]', '[13C]')
					.replace('[7*]', '[13C]').replace('[8*]', '[13C]').replace('[9*]', '[13C]')
					.replace('[*-]', '[13C-]').replace('[*:0]', '[13C]')
					.replace('*', '[13C]') for i in rxns]
			atommap_results = self.rxn_mapper.map_many_reactions(rxns)
			smarts = atommap_results['mapped_rxn']
			smarts = self.remove_pattern(smarts)
		else:
			smarts = None
		return smarts

	
	def remove_pattern(self, input_string):
		"""
		Remove the pattern that was previously introduced to overcome inability of RXNmapper to handle dummy atoms (*)
		:param input_string: SMARTS string
		:return:
		"""
		result_string = input_string.replace('13C', '*')
		return result_string
	
	def setImageOutputFolder(self, image_folder):
		"""
		:param image_folder: folder where the generated images are going to be stored
		:return:
		"""
		self.image_folder = image_folder
	
	def drawReactionSmarts(self, rxn_smarts, rxn_name, bw=False):
		"""
		:param rxn_smarts: SMARTS of the reaction
		:param rxn_name: name of the reaction (e.g. Rhea ID)
		:param bw: black and white, True / False
		:return: image is being saved to the defined folder
		"""
		smarts = rxn_smarts.replace('][', ']-[').replace(')[', ')-[')
		rxn = rdChemReactions.ReactionFromSmarts(smarts)
		drawer = rdMolDraw2D.MolDraw2DCairo(700, 300)
		if bw:
			drawer.drawOptions().useBWAtomPalette()
		drawer.DrawReaction(rxn)
		drawer.FinishDrawing()
		bio = io.BytesIO(drawer.GetDrawingText())
		img = Image.open(bio)
		img.save(f'{self.image_folder}/rxn_smarts_{rxn_name}.png')
	
	def drawCompoundSmarts(self, compounds_smarts, name, bw=False):
		"""
		
		:param compounds_smarts: SMARTS of compounds separated by dot (.)
		:param name: name of the file that will be found in reaction_png
		:param bw: black and white, True / False
		:return: image is being saved to the defined folder
		"""
		compounds_smarts = compounds_smarts.split('.')
		ms = [Chem.MolFromSmiles(x) for x in compounds_smarts]
		drawer = rdMolDraw2D.MolDraw2DCairo(500, 180, 200, 180)
		if bw:
			drawer.drawOptions().useBWAtomPalette()
		drawer.DrawMolecules(ms)
		drawer.FinishDrawing()
		bio = io.BytesIO(drawer.GetDrawingText())
		img = Image.open(bio)
		img.save(f'{self.image_folder}/cmps_smarts_{name}.png')