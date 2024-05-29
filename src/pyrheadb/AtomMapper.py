from rxnmapper import RXNMapper
import pandas as pd

class AtomMapper:
	def __init__(self):
		"""
		It is a wrapper for RXNMapper, created by Schwaller et al.
		"""
		self.rxn_mapper = RXNMapper()
	
	def map_one_reaction(self, rxnsmiles):
		"""
		Use function from RXNMapper to get atom map for one reaction.
		Atom mapper has 512 token length limitation. This function handles the exception resulting from this limitaiton.
		:param rxnsmiles: reaction smiles
		:return: atom map for one reaction
		"""
		try:
			res = self.rxn_mapper.get_attention_guided_atom_maps([rxnsmiles])[0]
			return pd.Series([res['mapped_rxn'], res['confidence']])
		except:
			return pd.Series(['tokenlength_error', 'err'])
	
	def map_many_reactions(self, rxnsmileslist):
		"""
		
		:param rxnsmileslist: list of reaction smiles
		:return: list of atom maps
		"""
		return [self.map_one_reaction(rxn) for rxn in rxnsmileslist]

# Conditions to consider:
# 1.  if '*' in row['rxn'] and row['MASTER_ID'] not in (63388):
# Rhea reaction 63388 : A +  Fe(II)-heme o + H2O = AH2 + Fe(II)-heme a is not processable using RXNmapper
# A - AH2 - too generic, Fe(II)-heme o - too big molecule
#
#
