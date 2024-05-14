from rxnmapper import RXNMapper
import os.path
import pickle
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem

tqdm.pandas()

class AtomMapper:
	def __init__(self, rheadata=None):
		"""
		This class is implementation of atom mapping of the reaction network of RheaDB.
		It is based on RXNMapper, created by Schwaller et al.
		"""
		self.rxn_mapper = RXNMapper()
		self.rheadata = rheadata # RheaDB object
	
	def map_one_reaction(self, rxnsmiles):
		"""
		Use function from RXNMapper to get atom map for one reaction
		:param rxnsmiles: reaction smiles
		:return: atom map for one reaction
		"""
		return self.rxn_mapper.get_attention_guided_atom_maps([rxnsmiles])[0]
	
	def map_many_reactions(self, rxnsmileslist):
		"""
		
		:param rxnsmileslist: list of reaction smiles
		:return: list of atom maps
		"""
		return [self.map_one_reaction(rxn) for rxn in rxnsmileslist]

	def get_dict_compid_to_smiles(self):
		"""
		OBSOLETE GET DIRECTLY FROM RHEA DATA
		:param df_rhea_rheactions_smiles_chebi:
		:return:
		"""
		dict_compid_to_smiles = dict()
		for index, row in self.rheadata.df_rhea_reactions_smiles_chebi.iterrows():
			smiles = row['rxnsmiles'].replace('>>', '.').split('.')
			chebi_equation = row['chebi_equation'].replace('>>', '.').split('.')
			assert len(smiles) == len(chebi_equation)
			for i in range(len(smiles)):
				dict_compid_to_smiles[chebi_equation[i]] = smiles[i]
		return dict_compid_to_smiles
	
	def number_atoms_in_dict(self, dict_compid_to_smiles):
		"""
		This function numbers the atoms in ChEBI compounds (independently from reactions)
		:param dict_compid_to_smiles: dictionary containing ChEBI compound ID with corresponding SMILES
		:return:
		"""
		dict_numbered = dict()
		for compid, smiles in dict_compid_to_smiles.items():
			dict_numbered[compid] = self.assign_atom_number(smiles)
		with open(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','chebi_smiles_numbered.pkl'), 'wb') as fp:
			pickle.dump(dict_numbered, fp)
		return dict_numbered
	
	def assign_atom_number(self, smiles):
		"""
		Called by number_atoms_in_dict() to assign numbers to atoms in ChEBI reactions
		:param smiles: smiles with atoms not numbered
		:return: mol with atoms numbered
		"""
		mol = Chem.MolFromSmiles(smiles)
		mol = self.mol_with_atom_index(mol)
		return Chem.MolToSmiles(mol)
	
	def mol_with_atom_index(self, mol):
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
			try:
				structure_match = self.structure_match(react_comp_univ_mapping, react_mapped)
			except Exception as e:
				print(e)
				return 'failed structure match'
			# print('structure_match reactant', structure_match)
			if structure_match is None or structure_match == 'num atoms match error':
				print("if structure_match is None or structure_match == 'num atoms match error':  -->")
				print(react_comp_univ_mapping)
				print(row['rheaid'])
				print(reactants_comp_id)
				print(reactants_mapped)
				print(" <-- if structure_match is None or structure_match == 'num atoms match error':")
				return None
			for atom in structure_match:
				reactants_atoms[atom] = react_comp_id + '-' + structure_match[atom]
		# print('reactants_atoms', reactants_atoms)
		
		reaction_atom_edges = list()
		for product_index in range(len(products_comp_id)):
			prod_comp_id = products_comp_id[product_index]
			prod_mapped = products_mapped[product_index]
			prod_comp_univ_mapping = dict_compid_to_numbered_smiles[products_comp_id[product_index]]
			try:
				structure_match = self.structure_match(prod_comp_univ_mapping, prod_mapped)
			except:
				return 'failed structure match'
			
			if structure_match is None or structure_match == 'num atoms match error':
				print(react_comp_univ_mapping)
				print(row['rheaid'])
				print(reactants_comp_id)
				print(reactants_mapped)
				return None
			
			for atom in structure_match:
				if atom.split(':')[0] != 'H':
					try:
						reaction_atom_edges.append((reactants_atoms[atom], prod_comp_id + '-' + structure_match[atom]))
					except Exception as e:
						return e
		return reaction_atom_edges
	
	def structure_match(self, smile1, smile2):
		
		qmol = Chem.MolFromSmarts(smile1)
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
	
	def debug_atommapper(self, rxnsmiles):
		"""
		Atom mapper has 512 token length limitation. This function overcomes the limitaiton
		:param rxnsmiles: input reaction smiles
		:return:
		"""
		try:
			res = self.rxn_mapper.map_one_reaction(rxnsmiles)
			return pd.Series([res['mapped_rxn'], res['confidence']])
		except:
			return pd.Series(['tokenlength_error', 'err'])
	
	def get_rhea_participants_dict(self, rheadata):
		dict_rhea_chebi = dict()
		for index, row in rheadata.df_smiles_master_id.iterrows():
			dict_rhea_chebi[row['MASTER_ID']] = row['chebi_equation'].replace('>>', '.').split('.')
		return dict_rhea_chebi
	
	def get_only_relevant_comps_dict(self, rheaid, dict_all_comp_id, dict_rheaid_comp_id):
		forbidden_chebis = ['CHEBI:4194', 'CHEBI:10545', 'CHEBI:24431', 'CHEBI:59720', 'CHEBI:84139', 'CHEBI:131859',
							'CHEBI:140159', 'CHEBI:137328', 'CHEBI:30212']
		return {key: dict_all_comp_id[key] for key in dict_rheaid_comp_id[rheaid] if
				key not in forbidden_chebis}  # photon, electron, etc
	
	def calculate_atommapped_reaction_data(self):
		# Initialize rhea data
		
		dict_rhea_id_comp_id = self.get_rhea_participants_dict(self.rheadata)
		
		list_all_participants = []
		for participants in dict_rhea_id_comp_id.values():
			list_all_participants.extend(participants)
		
		# get dict of chebi ids to smiles
		dict_compid_to_smiles = self.get_dict_compid_to_smiles(self.rheadata.df_smiles_chebi_equation)
		
		# generate mol objects for all chebi ids
		for compid in dict_compid_to_smiles:
			try:
				Chem.MolFromSmiles(dict_compid_to_smiles[compid])
			except Exception as e:
				print(e, compid, dict_compid_to_smiles[compid])
				
		dict_compid_to_mol = dict([(compid, Chem.MolFromSmiles(dict_compid_to_smiles[compid]))
								   for compid in dict_compid_to_smiles])
		
		file_chebi_numbered = os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','chebi_smiles_numbered.pkl')
		if not os.path.exists(file_chebi_numbered):
			# number atoms of all chebi compounds independently of reactions
			dict_compid_to_numbered_smiles = self.number_atoms_in_dict(dict_compid_to_smiles)
		else:
			with open(file_chebi_numbered, 'rb') as pickle_file:
				dict_compid_to_numbered_smiles = pickle.load(pickle_file)
		
		if not os.path.exists(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','rheadf_atom_mapped.tsv')):
			
			self.rheadata.df_one_dir_only_defined_rhea_smiles['rxnsmiles_no_h'] = \
			self.rheadata.df_one_dir_only_defined_rhea_smiles['rxnsmiles']. \
				apply(
				lambda x: x.replace('.[H+]', '').replace('[H+].', '').replace('[H][H].', '').replace('.[H][H]', ''))
			
			# atom map rhea reactions
			self.rheadata.df_one_dir_only_defined_rhea_smiles[['mapped_rxn', 'confidence']] = \
				self.rheadata.df_one_dir_only_defined_rhea_smiles['rxnsmiles_no_h'].progress_apply(self.debug_atommapper)
			self.rheadata.df_one_dir_only_defined_rhea_smiles.to_csv(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','rheadf_atom_mapped.tsv'), sep='\t', index=False)
		
		# remove the reactions with no atom mapping
		self.rheadata.df_one_dir_only_defined_rhea_smiles = self.rheadata.df_one_dir_only_defined_rhea_smiles[
			self.rheadata.df_one_dir_only_defined_rhea_smiles['mapped_rxn'] != 'tokenlength_error'
			]
		
		# for each reaction get relevant dict
		self.rheadata.df_one_dir_only_defined_rhea_smiles['relevant_comp_dict'] = \
			self.rheadata.df_one_dir_only_defined_rhea_smiles['MASTER_ID'].progress_apply(self.get_only_relevant_comps_dict)

		# match mapped_rxns_to_comp_ids
		self.rheadata.df_one_dir_only_defined_rhea_smiles['mapped_comp_ids'] = \
			self.rheadata.df_one_dir_only_defined_rhea_smiles.progress_apply(self.match_mapped_rxn_to_comp_ids, axis=1)
		
		self.rheadata.df_one_dir_only_defined_rhea_smiles.to_csv(
			os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','rheadf_atom_mapped.tsv'), sep='\t', index=False)
		
		self.rheadata.df_one_dir_only_defined_rhea_smiles[
			'not_found_chebi_id'] = self.rheadata.df_one_dir_only_defined_rhea_smiles.apply(
			lambda row: 'not_found' in row['mapped_comp_ids'], axis=1
		)
		self.rheadata.df_one_dir_only_defined_rhea_smiles[
			self.rheadata.df_one_dir_only_defined_rhea_smiles['not_found_chebi_id'] == True]. \
			to_csv(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','not_found_chebi_ids.tsv'), sep='\t', index=False)
		
		self.rheadata.df_one_dir_only_defined_rhea_smiles = self.rheadata.df_one_dir_only_defined_rhea_smiles[
			self.rheadata.df_one_dir_only_defined_rhea_smiles['not_found_chebi_id'] == False]
		
		# add edges data
		self.rheadata.df_one_dir_only_defined_rhea_smiles['atom_edges_reaction'] = \
			self.rheadata.df_one_dir_only_defined_rhea_smiles.progress_apply(self.match_compound_atoms, axis=1,
																		args=[dict_compid_to_numbered_smiles, ])
		self.rheadata.df_one_dir_only_defined_rhea_smiles.to_pickle(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','rheadf_atom_mapped_with_edges.pkl'))
	
	# edges_reaction = match_compound_atoms(comp_id_match_mapped_reaction, equation_mapped, dict_compid_to_numbered_smiles)

	def visualise_graph(self):
		df = pd.read_pickle(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','rheadf_atom_mapped_with_edges.pkl'))
		G = nx.Graph()
		for index, row in df.iterrows():
			if type(row['atom_edges_reaction']) == list:
				for edge in row['atom_edges_reaction']:
					G.add_edge(edge[0], edge[1])
		
		print(nx.number_connected_components(G))
		
		len_components = list()
		for component in nx.connected_components(G):
			len_components.append(len(component))
		len_comp_small = [i for i in len_components if i <= 10]
		len_comp_medium = [i for i in len_components if i < 10000 and i > 10]
		len_comp_large = [i for i in len_components if i >= 10000]
		print('Atoms connected:', len_comp_large)
		print('Atoms disconnected:', sum(len_comp_small) + sum(len_comp_medium))
		plt.hist(len_comp_small)
		plt.show()
		
		plt.hist(len_comp_medium)
		plt.show()

	def visualise_atom_mapped_graph(self):
		
		df = pd.read_pickle(os.path.join(self.rheadata.rhea_db_version_location,'atom_mapping','rheadf_atom_mapped_with_edges.pkl'))
		
		list_edges = []
		
		for edges in df['atom_edges_reaction']:
			if type(edges) == list:
				list_edges.extend(edges)
		
		G = nx.from_edgelist(list_edges)
		print('graph loaded')
		
		print(nx.number_connected_components(G))
		
		sizes = []
		for comp in nx.connected_components(G):
			sizes.append(len(comp))
		print(max(sizes))
		
		plt.hist(sizes)
		plt.xlabel('number of pairs in the graph component')
		plt.ylabel('number of components')
		plt.show()
		
		sizes = [i for i in sizes if i < 10]
		sizes.sort(reverse=True)
		
		print(len(sizes))
		plt.hist(sizes)
		plt.xlabel('number of pairs in the graph component')
		plt.ylabel('number of components')
		plt.show()
	
