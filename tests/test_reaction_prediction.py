# This is the test of the ReactionPrediction class

from rdkit import Chem
import unittest
from unittest.mock import patch, MagicMock, mock_open

from pyrheadb.ReactionPrediction import ReactionPrediction

class TestReactionPrediction(unittest.TestCase):
    def setUp(self):
        self.prediction = ReactionPrediction()

    @patch('builtins.open', new_callable=unittest.mock.mock_open, read_data='{"SMARTS": ["C>>CC"]}')
    def test_load_smarts_data(self):
        self.prediction.load_smarts_data()
        self.assertIn("SMARTS", self.prediction.smarts_data)
        self.assertEqual(self.prediction.smarts_data["SMARTS"][0], "C>>CC")

    def test_set_substrates_incorrect_input(self):
        with patch('rdkit.Chem.AllChem.MolFromSmiles', side_effect=[None]):
            self.prediction.set_substrates("InvalidSMILES")
            self.assertEqual(self.prediction.substrates, [])

    def test_predict_products_no_substrates_set(self):
        with self.assertRaises(ValueError):
            self.prediction.predict_products()

    @patch('rdkit.Chem.AllChem.ReactionFromSmarts')
    def test_predict_products_successful_prediction(self, mock_rxn_from):
        # Using actual RDKit functionality to convert SMILES to Mol
        self.prediction.smarts_data = {'SMARTS': ["[*:49][CH2:50][OH:51]>>[*:49][CH:50]=[O:51]"]}
        self.prediction.substrates = ['CC(O)CO']

        # Set up a mock reaction that simulates RDKit reaction behavior
        mock_rxn = MagicMock()
        product_mol = Chem.MolFromSmiles('CC(O)C=O')
        mock_rxn.RunReactants.return_value = [(product_mol,)]
        mock_rxn_from.return_value = mock_rxn

        self.prediction.predict_products()

        # Assert that reactions are populated
        self.assertGreater(len(self.prediction.reactions), 0)
        # Verify that reaction results as expected
        expected_product_smiles = Chem.MolToSmiles(product_mol)
        actual_product_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(self.prediction.reactions[0]['products'][0]))
        self.assertEqual(actual_product_smiles, expected_product_smiles)

    def test_check_balance_balanced(self):
        self.prediction.reactions = [{'substrate': 'C', 'products': ['CC'], 'smarts': 'C>>CC'}]
        with patch('pyrheadb.Reaction.Reaction.check_reaction_balance', return_value=True):
            balanced = self.prediction.get_balanced_reactions()
            self.assertEqual(len(balanced), 1)

    @patch('builtins.open', mock_open())
    def test_output_results(self):
        # Assuming the structure of each reaction in the list
        self.prediction.reactions = [
            {'substrate': 'C', 'products': ['CC'], 'smarts': 'C>>CC'}
        ]
    
        self.prediction.output_results('output.txt')
    
        # Verify that file was opened in write mode
        open.assert_called_once_with('output.txt', 'w')
    
        # Get the mock file handle from the open function
        file_handle = open()
    
        # Since there's one reaction, the write method should be called once for each line of output
        expected_write_calls = [((reaction['substrate'] + '>>' + str(reaction['products']) + '\n',),) for reaction in
                                self.prediction.reactions]
        file_handle.write.assert_has_calls(expected_write_calls)

if __name__ == '__main__':
    unittest.main()