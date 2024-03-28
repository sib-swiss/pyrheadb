import networkx as nx
import pandas as pd
import os

class ReactionNetwork():
	def __init__(self):
		pass
	
	def buildNetworkFromLongFormatTable(self, long_format_reaction_participant_table, hub_compounds_from='from_count'):
		"""
		Removes hub compounds and gretes a networkx graph
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
			with open(filepath+'/rhea_chebi_hub.tsv') as f:
				hub_compounds = ['CHEBI:'+i.strip() for i in f.readlines() if i]
			df_non_hubs = pd.DataFrame(counts_chebiids[~counts_chebiids['chebiid'].isin(hub_compounds)]['chebiid'])
		
		# Remove hub compounds from the network
		long_format_reaction_participant_table = long_format_reaction_participant_table.merge(df_non_hubs, on='chebiid', how='inner')
		# create chebiid to chebiid table
		df_m = long_format_reaction_participant_table.merge(long_format_reaction_participant_table, on='MASTER_ID')

		# Drop transports
		df_m = df_m[df_m['reaction_side_x'] != df_m['reaction_side_y']]
		df_m['chebipair'] = df_m.apply(lambda row: '.'.join(sorted([row['chebiid_x'], row['chebiid_y']])), axis=1)
		df_m.drop(columns=['reaction_side_x', 'chebiid_x', 'smiles_x',
		'reaction_side_y', 'chebiid_y', 'smiles_y'], inplace=True)
		grouped = df_m.groupby('chebipair', as_index=False, dropna=False).agg(lambda x: set(list(x)))
		grouped['chebiid1']=grouped['chebipair'].apply(lambda x: x.split('.')[0])
		grouped['chebiid2']=grouped['chebipair'].apply(lambda x: x.split('.')[1])
		self.participants_graph = nx.from_pandas_edgelist(grouped, 'chebiid1', 'chebiid2', ['MASTER_ID'])
		
		