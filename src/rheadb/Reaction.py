__author__ = "Anastasia Sveshnikova"
__email__ = "anastasia.sveshnikova@sib.swiss"
__status__ = "Prototype"

import re
from rdkit.Chem import AllChem
from rdkit import RDLogger

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


