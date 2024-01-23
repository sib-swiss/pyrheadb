import networkx as nx
import pandas as pd


class ReactionNetwork():
	def __init__(self):
		pass
	
	def buildNetworkFromLongFormatTable(self, long_format_reaction_participant_table):
		"""
		Removes hub compounds and gretes a networkx graph
		:param long_format_reaction_participant_table: table that contains information about which participant is in which reaction
		:return: adds reaction graph to self
		"""
		counts_chebiids = long_format_reaction_participant_table['chebiid'].value_counts().rename_axis('chebiid').reset_index(name='counts')
		
		# identify non hub compounds
		non_hub_compounds = counts_chebiids[counts_chebiids['counts']<100]['chebiid']
		df_non_hubs = pd.DataFrame(non_hub_compounds)
		
		# Remove hub compounds from the network
		long_format_reaction_participant_table = long_format_reaction_participant_table.merge(df_non_hubs, on='chebiid', how='inner')
		
		# create chebiid to chebiid table
		df_m = long_format_reaction_participant_table.merge(long_format_reaction_participant_table, on='MASTER_ID')
		
		# Drop transports
		df_m = df_m[df_m['reaction_side_x'] != df_m['reaction_side_y']]
		
		self.participants_graph = nx.from_pandas_edgelist(df_m, 'chebiid_x', 'chebiid_y')
		
		