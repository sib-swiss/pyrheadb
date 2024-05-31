import os.path
import pickle
import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem
import pandas as pd
from tqdm import tqdm

from .AtomMapper import AtomMapper
from .Reaction import Reaction

tqdm.pandas()

class AtomTransitionNetwork():
	
	def __init__(self, rheadata=None):
		"""
		This class is implementation of atom mapping of the reaction network of RheaDB to transform it into atom transition network.
		:param rheadata:
		"""
		self.rheadata = rheadata  # RheaDB object
		self.rxn_mapper = AtomMapper()
		self.reactionobj = Reaction()
		os.makedirs(os.path.join(rheadata.rhea_db_version_location, 'atom_transition_network'), exist_ok=True)
		
		self.filter_rhea_data()
		self.calculate_dict_rheaid_comp_id()
		
	def filter_rhea_data(self):
		"""
		Atom transition network works based on the fully defined compounds.
		Certain definitions used in the Rhea database to not fit the purpose of atom transition network.
		1. R/star(*)-reactions : do not fit since they represent a class of reactions and not a particular reaction.
		2. generic reducing agent - generic oxidation agent : do not fit since they do not transfer any non-hydrogen atom.
		:return:
		"""
		self.reactiondata = self.rheadata.df_smiles_master_id.copy()
		print('Calculating reaction balance')
		self.reactiondata['balance']=self.reactiondata['rxnsmiles'].progress_apply(self.reactionobj.check_reaction_balance, format='smiles', account_H=False)
		# drop unbalanced reactions
		self.reactiondata=self.reactiondata[self.reactiondata['balance']==True]
		# drop reactions with star compounds
		self.reactiondata['rxnsmilesnoredox']=self.reactiondata['rxnsmiles'].progress_apply(self.reactionobj.clean_rxnsmiles_of_redox)
		self.reactiondata.drop(columns=['rxnsmiles'], inplace=True)
		self.reactiondata.rename(columns={'rxnsmilesnoredox':'rxnsmiles'}, inplace=True)
		self.reactiondata['star_no_redox'] = self.reactiondata['rxnsmiles'].apply(lambda x: "*" in x)
		self.reactiondata = self.reactiondata[self.reactiondata['star_no_redox'] == False]
		# drop reactions with polymers
		self.reactiondata['polymer_reaction'] = self.reactiondata['chebi_equation'].apply(lambda x: 'POLYMER' in x)
		self.reactiondata = self.reactiondata[self.reactiondata['polymer_reaction '] == False]
		
	def get_dict_compid_to_smiles(self):
		"""
		OBSOLETE GET DIRECTLY FROM RHEA DATA
		:param df_rhea_rheactions_smiles_chebi:
		:return:
		"""
		chebiids = self.rheadata.df_chebi_smiles['chebiid'].to_list()
		smiles = self.rheadata.df_chebi_smiles['smiles'].to_list()
		dict_compid_to_smiles = dict(zip(chebiids, smiles))
		# for index, row in self.reactiondata.iterrows():
		# 	smiles = row['rxnsmiles'].replace('>>', '.').split('.')
		# 	chebi_equation = row['chebi_equation'].replace('>>', '.').split('.')
		# 	assert len(smiles) == len(chebi_equation)
		# 	for i in range(len(smiles)):
		# 		dict_compid_to_smiles[chebi_equation[i]] = smiles[i]
		
		self.test_all_compounds_can_be_molobj(dict_compid_to_smiles)
		
		return dict_compid_to_smiles
	
	def test_all_compounds_can_be_molobj(self, dict_compid_to_smiles):
		# check generate mol objects for all chebi ids
		for compid in dict_compid_to_smiles:
			try:
				Chem.MolFromSmiles(dict_compid_to_smiles[compid])
			except Exception as e:
				print(e, compid, dict_compid_to_smiles[compid])
	def number_atoms_in_dict(self, dict_compid_to_smiles):
		"""
		This function numbers the atoms in ChEBI compounds (independently of reactions)
		:param dict_compid_to_smiles: dictionary containing ChEBI compound ID with corresponding SMILES
		:return:
		"""
		dict_numbered = dict()
		for compid, smiles in dict_compid_to_smiles.items():
			dict_numbered[compid] = self.__assign_atom_number(smiles)
		with open(os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network', 'chebi_smiles_numbered.pkl'),
				  'wb') as fp:
			pickle.dump(dict_numbered, fp)
		return dict_numbered
	
	def __assign_atom_number(self, smiles):
		"""
		Called by number_atoms_in_dict() to assign numbers to atoms in ChEBI reactions
		:param smiles: smiles with atoms not numbered
		:return: mol with atoms numbered
		"""
		mol = Chem.MolFromSmiles(smiles)
		mol = self.__mol_with_atom_index(mol)
		return Chem.MolToSmiles(mol)
	
	def __mol_with_atom_index(self, mol):
		"""
		Called by number_atoms_in_dict() -> assign_atom_number() to assign numbers to atoms in ChEBI reactions
		:param mol: rdkit Mol object
		:return: mol object with atoms numbered
		"""
		for atom in mol.GetAtoms():
			atom.SetAtomMapNum(atom.GetIdx() + 1)  # necessary to add one since 0 cannot be atom map number
		return mol
	
	def match_mapped_rxn_to_comp_ids(self, row):
		"""
		Since rxnmapper atom mapping changes the order of compounds in the smiles, to be able to link them to CHEBI ID,
		we have to map the compounds back to the original reaction order based on structure
		:param row:
		:return:
		"""
		mapped_smiles = row['mapped_rxn']
		dict_compid_to_mol = row['relevant_comp_dict']
		reactants = [Chem.MolFromSmiles(sm) for sm in mapped_smiles.split('>>')[0].split('.')]
		products = [Chem.MolFromSmiles(sm) for sm in mapped_smiles.split('>>')[1].split('.')]
		reactant_comp_ids = [self.find_compid_in_compid_dict(mol, dict_compid_to_mol) for mol in reactants]
		product_comp_ids = [self.find_compid_in_compid_dict(mol, dict_compid_to_mol) for mol in products]
		return '.'.join(reactant_comp_ids) + '>>' + '.'.join(product_comp_ids)
	
	def find_compid_in_compid_dict(self, check_mol, dict_compid_to_mol):
		"""
		find compound based on its structure
		:param check_mol:
		:param dict_compid_to_mol:
		:return:
		"""
		for comp_id, mol in dict_compid_to_mol.items():
			if self.same_mol(mol, check_mol):
				return comp_id
		return 'not_found'
	
	def same_mol(self, mol1, mol2):
		"""
		Trick to identify if the two molecules are the same: first molecule is a substructure of the second, and
		second molecule is a substructure of the first
		:param mol1: rdkit Mol object
		:param mol2: rdkit Mol object
		:return:
		"""
		return mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1)
	
	def match_compound_atoms(self, row, dict_compid_to_numbered_smiles):
		"""
		Function that helps to match the atoms of the reactants and products
		with their universal (not dependent on reaction) numbering
		:param equation_compounds:
		:param equation_mapped:
		:param dict_compid_to_numbered_smiles:
		:return:
		"""
		equation_compounds = row['mapped_comp_ids']
		equation_mapped = row['mapped_rxn']
		
		reactants_comp_id = equation_compounds.split('>>')[0].split('.')
		products_comp_id = equation_compounds.split('>>')[1].split('.')
		
		reactants_mapped = equation_mapped.split('>>')[0].split('.')
		products_mapped = equation_mapped.split('>>')[1].split('.')
		
		reactants_atoms = dict()
		for reactant_index in range(len(reactants_comp_id)):
			react_comp_id = reactants_comp_id[reactant_index]
			react_mapped = reactants_mapped[reactant_index]
			react_comp_univ_mapping = dict_compid_to_numbered_smiles[reactants_comp_id[reactant_index]]
			# Get structure match between universal atom number of compound and atom map number specific to reaction
			structure_match = self.get_structure_match(react_comp_univ_mapping, react_mapped, react_comp_id)
			if structure_match is None:
				print('No structure match in Rhea ID:', row['rheaid'])
				return None
			# Add this atom number
			for atom in structure_match:
				reactants_atoms[atom] = react_comp_id + '-' + structure_match[atom]
		
		reaction_atom_edges = list()
		for product_index in range(len(products_comp_id)):
			prod_comp_id = products_comp_id[product_index]
			prod_mapped = products_mapped[product_index]
			prod_comp_univ_mapping = dict_compid_to_numbered_smiles[products_comp_id[product_index]]
			structure_match = self.get_structure_match(prod_comp_univ_mapping, prod_mapped, prod_comp_id)
			if structure_match is None:
				print('No structure match in Rhea ID:', row['rheaid'])
				return None
			for atom in structure_match:
				if atom.split(':')[0] != 'H':
					try:
						reaction_atom_edges.append((reactants_atoms[atom], prod_comp_id + '-' + structure_match[atom]))
					except Exception as e:
						return e
		return reaction_atom_edges
	
	def get_structure_match(self, comp_univ_mapping, comp_mapped, comp_id):
		try:
			structure_match = self.structure_match(comp_univ_mapping, comp_mapped)
		except Exception as e:
			print(e)
			return None
		# print('structure_match reactant', structure_match)
		if structure_match is None or structure_match == 'num atoms match error':
			print("if structure_match is None or structure_match == 'num atoms match error':  -->")
			print(comp_univ_mapping)
			print(comp_id)
			print(comp_mapped)
			print(" <-- if structure_match is None or structure_match == 'num atoms match error':")
			return None
		return structure_match
	
	def structure_match(self, smile1, smile2):
		
		qmol = Chem.MolFromSmiles(smile1)
		mol = Chem.MolFromSmiles(smile2)
		
		if not qmol.GetNumAtoms() == mol.GetNumAtoms():
			return 'num atoms match error'
		
		ind_map_qmol = {}
		map_qmol = {}
		for atom in qmol.GetAtoms():
			map_num = atom.GetAtomMapNum()
			if map_num and atom.GetSymbol() != 'H':
				ind_map_qmol[atom.GetSymbol() + ':' + str(map_num)] = atom.GetIdx()
				map_qmol[atom.GetIdx()] = atom.GetSymbol() + ':' + str(map_num)
		
		map_list_qmol = [ind_map_qmol[x] for x in sorted(ind_map_qmol)]
		
		map_mol = {}
		
		for atom in mol.GetAtoms():
			map_num = atom.GetAtomMapNum()
			if map_num and atom.GetSymbol() != 'H':
				# reverse atom index and atom map number for the second compound
				map_mol[atom.GetIdx()] = atom.GetSymbol() + ':' + str(map_num)
		
		for match in mol.GetSubstructMatches(qmol):
			match_map = {}
			mas = [match[x] for x in map_list_qmol]
			for i, m in enumerate(mas):
				match_map[map_mol[i]] = map_qmol[m]
			
			return match_map
	
	def number_atoms_in_chebi_molecules(self, dict_compid_to_smiles):
		file_chebi_numbered = os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network',
										   'chebi_smiles_numbered.pkl')
		if not os.path.exists(file_chebi_numbered):
			# number atoms of all chebi compounds independently of reactions
			dict_compid_to_numbered_smiles = self.number_atoms_in_dict(dict_compid_to_smiles)
		else:
			with open(file_chebi_numbered, 'rb') as pickle_file:
				dict_compid_to_numbered_smiles = pickle.load(pickle_file)
				
		return dict_compid_to_numbered_smiles
	def get_rhea_participants_dict(self, rheadata):
		dict_rhea_chebi = dict()
		for index, row in rheadata.df_smiles_master_id.iterrows():
			dict_rhea_chebi[row['MASTER_ID']] = row['chebi_equation'].replace('>>', '.').split('.')
		return dict_rhea_chebi
	
	def get_only_relevant_comps_dict(self, rheaid, dict_all_comp_id=dict()):
		forbidden_chebis = ['CHEBI:4194', 'CHEBI:10545', 'CHEBI:24431', 'CHEBI:59720', 'CHEBI:84139', 'CHEBI:131859',
							'CHEBI:140159', 'CHEBI:137328', 'CHEBI:30212']
		return {key: Chem.MolFromSmiles(dict_all_comp_id[key]) for key in self.dict_rheaid_comp_id[rheaid] if
				key not in forbidden_chebis}  # photon, electron, etc

	def calculate_dict_rheaid_comp_id(self):
		df = self.rheadata.rhea_reaction_long_format_smiles_chebi.copy()
		self.dict_rheaid_comp_id = df.groupby("MASTER_ID")["chebiid"].apply(list).to_dict()
		
	def create_atom_transition_network(self):
		
		dict_rhea_id_comp_id = self.get_rhea_participants_dict(self.rheadata)
		
		list_all_participants = []
		for participants in dict_rhea_id_comp_id.values():
			list_all_participants.extend(participants)
		
		# get dict of chebi ids to smiles
		dict_compid_to_smiles = self.get_dict_compid_to_smiles()
		
		# number atoms to stabilise a number for each atom in compound
		dict_compid_to_numbered_smiles = self.number_atoms_in_chebi_molecules(dict_compid_to_smiles)
		
		atom_mapped_df_file = os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network', 'rheadf_atom_mapped.tsv')
		if not os.path.exists(atom_mapped_df_file):
			self.calculate_atom_mapped_reaction_data(atom_mapped_df_file)
		
		df = pd.read_csv(atom_mapped_df_file, sep='\t')
		
		# for each reaction get relevant dict
		df['relevant_comp_dict'] = \
			df['MASTER_ID'].progress_apply(self.get_only_relevant_comps_dict, dict_all_comp_id=dict_compid_to_numbered_smiles)
		
		# match mapped_rxns_to_comp_ids
		df['mapped_comp_ids'] = df.progress_apply(self.match_mapped_rxn_to_comp_ids, axis=1)
		
		# Save intermediate step as pickle
		df.to_pickle(os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network', 'rheadf_atom_mapped.pkl'))
		
		# Identify chebiids that could not be matched - should be empty list -> review the list and fix the problem!
		df['not_found_chebi_id'] = df.apply(lambda row: 'not_found' in row['mapped_comp_ids'], axis=1)
		df[df['not_found_chebi_id'] == True]. \
			to_csv(os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network', 'not_found_chebi_ids.tsv'),
				   sep='\t', index=False)
		
		# proceed only with reactions for which every chebiid was found
		df = df[df['not_found_chebi_id'] == False]
		
		# add edges data
		df['atom_edges_reaction'] = df.progress_apply(self.match_compound_atoms, axis=1,
															   args=[dict_compid_to_numbered_smiles, ])
		# Save as pickle
		df.to_pickle(
			os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network', 'rheadf_atom_mapped_with_edges.pkl'))
		
		self.load_networkx_graph()
	
	# edges_reaction = match_compound_atoms(comp_id_match_mapped_reaction, equation_mapped, dict_compid_to_numbered_smiles)
	
	def calculate_atom_mapped_reaction_data(self, atom_mapped_df_file):
		
		df = self.reactiondata.copy()
		
		df['rxnsmiles_no_h'] = df['rxnsmiles'].apply(
				lambda x: x.replace('.[H+]', '').replace('[H+].', '').replace('[H][H].', '').replace('.[H][H]', ''))
		
		# atom map rhea reactions
		df[['mapped_rxn', 'confidence']] = df['rxnsmiles_no_h'].progress_apply(self.rxn_mapper.map_one_reaction)
		
		# remove the reactions with no atom mapping
		df = df[df['mapped_rxn'] != 'tokenlength_error']
		
		df.to_csv(atom_mapped_df_file, sep='\t', index=False)
	
	def load_networkx_graph(self):
		df = pd.read_pickle(
			os.path.join(self.rheadata.rhea_db_version_location, 'atom_transition_network',
						 'rheadf_atom_mapped_with_edges.pkl'))
		
		list_edges = []
		
		for edges in df['atom_edges_reaction']:
			if type(edges) == list:
				list_edges.extend(edges)
		
		self.G = nx.from_edgelist(list_edges)
		
	def plot_subgraph_size_vs_count(self, plotid='default'):
		"""
		This function saves the plot of disconnected subgraph sizes (number of compound-atoms included) vs
		number of subgraphs of this size
		:return:
		"""
		
		print('Number of subgraphs:', nx.number_connected_components(self.G))
		
		sizes = []
		for comp in nx.connected_components(self.G):
			sizes.append(len(comp))
		print('Number of compound-assigned atoms in the biggest subgraph:', max(sizes))
		
		len_comp_small = [i for i in sizes if i <= 10]
		len_comp_medium = [i for i in sizes if i < 10000 and i > 10]
		
		print('Atoms disconnected from biggest subnetwork:', sum(len_comp_small) + sum(len_comp_medium))
		
		# Define custom bins
		bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
				17, 18, 19, 20, 50, 100, max(sizes)]
		bin_labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
					  '13', '14', '15', '16', '17', '18', '19', '20-49', '50-99',
					  f'100-{max(sizes)}']
		
		# Count the sizes in each bin range
		counts = [0] * (len(bins) - 1)
		for size in sizes:
			for i in range(len(bins) - 1):
				if bins[i] <= size < bins[i + 1]:
					counts[i] += 1
					break
		
		# Add the count for the max bin
		counts[-1] += sizes.count(max(sizes))
		
		# Plotting the bar chart with custom bins
		fig, ax = plt.subplots(figsize=(12, 8))
		ax.bar(range(len(counts)), counts, tick_label=bin_labels, edgecolor='black', color='black')
		
		# Customizing ticks and labels
		ax.tick_params(axis='x', rotation=45, labelsize=15)
		ax.tick_params(axis='y', labelsize=15)
		
		ax.set_yscale('log')
		
		ax.set_xlabel('Number of compound-atoms (e.g. CHEBI:57681-C:2) in subgraph', fontsize=20)
		ax.set_ylabel('Number of subgraphs', fontsize=20)
		ax.set_title('Connectivity of the atom transition network', fontsize=30)
		
		ax.grid(axis='y', linestyle='--', alpha=0.7)
		fig.tight_layout()
		
		os.makedirs('plots', exist_ok=True)
		plt.savefig(os.path.join('plots', f'atom_transition_network_subgraphs_size_vs_count_{plotid}.png'))
		
		plt.close()


# import igraph as ig
# from itertools import product
# reactions_file = 'reactions_enumeration.tsv'
# import matplotlib.pyplot as plt
#
# def produce_graph(reactions_file):
# 	with open(reactions_file) as f, open('graph.ncol', 'w') as g:
# 		header = f.readline().split('\t')
# 		equ_index = header.index('stoichiometric_equation')
# 		for line in f:
# 			eq = line.split('\t')[equ_index]
# 			reacts = [i.split()[1] for i in eq.split(' <=> ')[0].split(' + ')]
# 			reacts = [i for i in reacts if i.startswith('SLM:')]
# 			products = [i.split()[1] for i in eq.split(' <=> ')[1].split(' + ')]
# 			products = [i for i in products if i.startswith('SLM:')]
# 			for combination in product(reacts, products):
# 				g.write(f'{combination[0]} {combination[1]}\n')
#
# #produce_graph(reactions_file)
# g = ig.Graph.Read_Ncol('graph.ncol')
# sum = 0
# comp_sizes=[]
# for component in g.components():
# 	sum+=len(component)
# 	comp_sizes.append(len(component))
# print(sum)
