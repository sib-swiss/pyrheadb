import networkx as nx
import pandas as pd
import os

class ReactionNetwork():
	def __init__(self):
		pass
	
	def build_network_from_long_format_table(self, long_format_reaction_participant_table, hub_compounds_from='from_count'):
		"""
		Removes hub compounds and creates a networkx graph
		:param long_format_reaction_participant_table: table that contains information about which participant is in which reaction
		:param hub_compounds_from: 'from_count' or 'from_list'
		:return: adds reaction graph to self
		"""
		counts_chebiids = long_format_reaction_participant_table['chebiid'].value_counts().rename_axis('chebiid').reset_index(name='counts')
		
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
		long_format_reaction_participant_table = long_format_reaction_participant_table.merge(df_non_hubs, on='chebiid', how='inner')
		# create chebiid to chebiid table
		df_m = long_format_reaction_participant_table.merge(long_format_reaction_participant_table, on='MASTER_ID')

		# Drop transports
		df_m = df_m[df_m['reaction_side_x'] != df_m['reaction_side_y']]
		print(df_m[df_m['MASTER_ID'].isin([69340, 69348, 69352])])
		df_m['chebipair'] = df_m.apply(lambda row: '.'.join(sorted([row['chebiid_x'], row['chebiid_y']])), axis=1)
		df_m.drop(columns=['reaction_side_x', 'chebiid_x', 'smiles_x',
		'reaction_side_y', 'chebiid_y', 'smiles_y'], inplace=True)
		grouped = df_m.groupby('chebipair', as_index=False, dropna=False).agg(lambda x: set(list(x)))
		grouped['chebiid1']=grouped['chebipair'].apply(lambda x: x.split('.')[0])
		grouped['chebiid2']=grouped['chebipair'].apply(lambda x: x.split('.')[1])
		self.participants_graph = nx.from_pandas_edgelist(grouped, 'chebiid1', 'chebiid2', ['MASTER_ID'])
	
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
