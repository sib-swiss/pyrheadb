from rxnmapper import RXNMapper

class AtomMapper:
	def __init__(self):
		"""
		This class is a wrapper for RXNMapper, created by Schwaller et al.
		"""
		self.rxn_mapper = RXNMapper()
	
	def map_one_reaction(self, rxnsmiles):
		return self.rxn_mapper.get_attention_guided_atom_maps([rxnsmiles])[0]
	
	def map_many_reactions(self, rxnsmileslist):
		return [self.map_one_reaction(rxn) for rxn in rxnsmileslist]
