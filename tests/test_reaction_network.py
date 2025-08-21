import unittest

from pyrheadb.ReactionNetwork import ReactionNetwork

class TestReactionNetwork(unittest.TestCase):
    def setUp(self):
        self.reaction_net = ReactionNetwork()
	
    def test_build_network_from_reaction_smiles_list(self):
        reactions = ['CCO.C[CH]C>>C[CH]CC', 'CCCC>>CCC.C', 'CCO.O>>CCO[O-]']
        self.reaction_net.build_network_from_reaction_smiles_list(reactions)
        print("Nodes:", self.reaction_net.participant_graph.nodes())
        print("Edges:", self.reaction_net.participant_graph.edges())
        #self.assertTrue(result)
		
    def test_get_disconnected_compounds(self):
        reactions = ['CCO.C[CH]C>>C[CH]CC', 'CCCC>>CCC.C', 'CCO.O>>CCO[O-]']
        self.reaction_net.build_network_from_reaction_smiles_list(reactions)
        disconnected_compounds = self.reaction_net.get_disconnected_compounds()
        print(disconnected_compounds)