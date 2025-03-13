import networkx as nx
import pandas as pd
import os

class ReactionNetwork():
	def __init__(self, rdb, type='chebiid', hub_compounds_from='from_count'):
		"""
		: param type : describes if network should be built from chebi ids, ichikeys or inchikey14L as nodes.
		processes long_format_reaction_participant_table: table that contains information about which participant is in which reaction
		"""
		self.long_format_reaction_participant_table = rdb.rhea_reaction_long_format_smiles_chebi
		if type == 'chebiid':
			self.build_chebiid_network_from_long_format_table(hub_compounds_from=hub_compounds_from)
		elif type == 'inchikey':
			self.build_inchikey_network_from_long_format_table(hub_compounds_from=hub_compounds_from)
		elif type == 'inchikey14L':
			self.build_inchikey14L_network_from_long_format_table(hub_compounds_from=hub_compounds_from)

	def build_chebiid_network_from_long_format_table(self, hub_compounds_from='from_count'):
		
		self.filter_hub_compounds_from_long_format_table( hub_compounds_from=hub_compounds_from)
		
		self.df_long_format['chebipair'] = self.df_long_format.apply(lambda row: '.'.join(sorted([row['chebiid_x'], row['chebiid_y']])), axis=1)
		grouped = self.df_long_format.groupby('chebipair', as_index=False, dropna=False).agg(lambda x: set(list(x)))
		grouped['chebiid1']=grouped['chebipair'].apply(lambda x: x.split('.')[0])
		grouped['chebiid2']=grouped['chebipair'].apply(lambda x: x.split('.')[1])
		self.participants_graph = nx.from_pandas_edgelist(grouped, 'chebiid1', 'chebiid2', ['MASTER_ID'])
	
	def build_inchikey_network_from_long_format_table(self,  hub_compounds_from='from_count'):
		
		self.filter_hub_compounds_from_long_format_table(hub_compounds_from=hub_compounds_from)
		
		self.df_long_format = self.df_long_format[self.df_long_format['inchikey_x']!='NA-NA-NA']
		self.df_long_format = self.df_long_format[self.df_long_format['inchikey_y']!='NA-NA-NA']
		
		self.df_long_format['inchikeypair'] = self.df_long_format.apply(lambda row: '.'.join(sorted([row['inchikey_x'], row['inchikey_y']])), axis=1)
		grouped = self.df_long_format.groupby('inchikeypair', as_index=False, dropna=False).agg(lambda x: set(list(x)))
		grouped['inchikey1']=grouped['inchikeypair'].apply(lambda x: x.split('.')[0])
		grouped['inchikey2']=grouped['inchikeypair'].apply(lambda x: x.split('.')[1])
		self.participants_graph = nx.from_pandas_edgelist(grouped, 'inchikey1', 'inchikey2', ['MASTER_ID'])	

	def build_inchikey14L_network_from_long_format_table(self,  hub_compounds_from='from_count'):
		
		self.filter_hub_compounds_from_long_format_table( hub_compounds_from=hub_compounds_from)
				
		self.df_long_format = self.df_long_format[self.df_long_format['inchikey_x']!='NA-NA-NA']
		self.df_long_format = self.df_long_format[self.df_long_format['inchikey_y']!='NA-NA-NA']
		
		self.df_long_format['inchikeypair'] = self.df_long_format.apply(lambda row: '.'.join(sorted([row['inchikey14L_x'], row['inchikey14L_y']])), axis=1)
		grouped = self.df_long_format.groupby('inchikeypair', as_index=False, dropna=False).agg(lambda x: set(list(x)))
		grouped['inchikey14L1']=grouped['inchikeypair'].apply(lambda x: x.split('.')[0])
		grouped['inchikey14L2']=grouped['inchikeypair'].apply(lambda x: x.split('.')[1])
		self.participants_graph = nx.from_pandas_edgelist(grouped, 'inchikey14L1', 'inchikey14L2', ['MASTER_ID'])	
		
	def filter_hub_compounds_from_long_format_table(self, hub_compounds_from='from_count'):
		"""
		Removes hub compounds and creates a networkx graph
		:param hub_compounds_from: 'from_count' or 'from_list' or 'all_compounds'
		:return: adds reaction graph to self
		"""
		counts_chebiids = self.long_format_reaction_participant_table['chebiid'].value_counts().rename_axis('chebiid').reset_index(name='counts')
		
		# identify non hub compounds
		# by count
		if hub_compounds_from == 'from_count':
			df_non_hubs = pd.DataFrame(counts_chebiids[counts_chebiids['counts']<100]['chebiid'])
		elif hub_compounds_from == 'from_list':
			# by manual list of hub compounds
			filepath = os.path.dirname(__file__)
			with open(os.path.join(filepath, 'biochemical_assumptions', 'rhea_chebi_hub.tsv')) as f:
				hub_compounds = ['CHEBI:'+i.strip() for i in f.readlines() if i]
			df_non_hubs = pd.DataFrame(counts_chebiids[~counts_chebiids['chebiid'].isin(hub_compounds)]['chebiid'])
		elif hub_compounds_from == 'all_compounds':
			df_non_hubs = counts_chebiids['chebiid']
		# Remove hub compounds from the network
		self.long_format_reaction_participant_table = self.long_format_reaction_participant_table.merge(df_non_hubs, on='chebiid', how='inner')
		# Remove polymer compounds from the network
		self.long_format_reaction_participant_table['POLYMER'] = self.long_format_reaction_participant_table['chebiid'].apply(lambda x: x.startswith('POLYMER'))
		
		self.long_format_reaction_participant_table = self.long_format_reaction_participant_table[self.long_format_reaction_participant_table['POLYMER']==False]
		# create chebiid to chebiid table
		df_filtered = self.long_format_reaction_participant_table.merge(self.long_format_reaction_participant_table, on='MASTER_ID')
		# Drop same reaction side connections
		df_filtered = df_filtered[df_filtered['reaction_side_x'] != df_filtered['reaction_side_y']]
		#print(df_filtered[df_filtered['MASTER_ID'].isin([69340, 69348, 69352])])
		df_filtered.drop(columns=['reaction_side_x',  'smiles_x', 'inchi_x',
		                   'reaction_side_y',  'smiles_y', 'inchi_y'], inplace=True)
		self.df_long_format = df_filtered
	
	def get_disconnected_compounds(self):
		components = [comp for comp in nx.connected_components(self.participants_graph)]
		component_sizes = [len(comp) for comp in components]

		compounds = []

		for comp in components:
			if len(comp) < max(component_sizes):
				compounds.extend(comp)
		return compounds

	def analyze_basic_properties(self):
		"""
		1. Basic Network Properties. Analyzes basic properties of the participants graph.
		This method would provide general statistics about the network,
		including the number of nodes, number of edges, average degree, density of the network,
		and whether the network is connected.
		"""
		properties = {
			'Number of Nodes': len(self.participants_graph.nodes()),
			'Number of Edges': len(self.participants_graph.edges()),
			'Average Degree': sum(dict(self.participants_graph.degree()).values()) / len(
				self.participants_graph.nodes()),
			'Density': nx.density(self.participants_graph),
			'Is Connected': nx.is_connected(self.participants_graph)
		}
		return properties

	def analyze_components(self):
		"""
		2. Component Analysis. Analyzes the connected components of the participants graph.
		Identify and analyze connected components within the network.
		This is particularly useful for understanding isolated subnetworks within your reaction network.
		
		"""
		components = [comp for comp in nx.connected_components(self.participants_graph)]
		component_sizes = [len(comp) for comp in components]
		largest_component = max(components, key=len)
		return {
			'Number of Components': len(components),
			'Sizes of Components': component_sizes,
			'Largest Component Size': len(largest_component)
		}

	def analyze_centrality(self):
		"""
		3. Centrality Measures. Computes centrality measures for nodes in the participants graph.
		Compute various centrality measures to find the most influential nodes (compounds) in your network.
		Measures like degree centrality, betweenness centrality, and closeness centrality can be useful.
		"""
		centrality = {
			'Degree Centrality': nx.degree_centrality(self.participants_graph),
			'Betweenness Centrality': nx.betweenness_centrality(self.participants_graph),
			'Closeness Centrality': nx.closeness_centrality(self.participants_graph)
		}
		return centrality

	def analyze_clustering(self):
		"""
		4. Clustering Coefficient. Analyzes the clustering coefficient of the participants graph.
		This analysis provides insights into the clustering behavior of nodes, indicating how nodes tend to cluster together.
		"""
		clustering = nx.clustering(self.participants_graph)
		clustering_inverse = [(value, key) for key, value in clustering.items()]
		clustering_inverse.sort()
		average_clustering = nx.average_clustering(self.participants_graph)
		return {'Clustering Coefficients': clustering, 'Average Clustering': average_clustering, 'Max 20 Clustering':clustering_inverse[-20:]}

	def analyze_shortest_paths(self):
		"""
		5. Path Analysis. Analyzes shortest paths in the participants graph.
		Examine shortest paths within the network to understand the connectivity and distance between different compounds.
		"""
		all_pairs_shortest_path = dict(nx.all_pairs_shortest_path_length(self.participants_graph))
		return all_pairs_shortest_path
