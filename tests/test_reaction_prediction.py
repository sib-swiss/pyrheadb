# This is the test of the ReactionPrediction class

from rdkit import Chem
import unittest
from unittest.mock import patch, MagicMock, mock_open
import os
from pyrheadb.ReactionPrediction import ReactionPrediction
from pyrheadb.RheaDB import RheaDB

current_dir = os.path.dirname(__file__)  # Gets the directory where this test file resides
test_rhea_db_path = os.path.join(current_dir, 'test_rhea_db')

rhea_db_version = 130
rdb = RheaDB(rhea_versions_folder_location=test_rhea_db_path, rhea_version=rhea_db_version)

class TestReactionPrediction(unittest.TestCase):
    def setUp(self):
        self.prediction = ReactionPrediction(rdb)

    def test_predict_products_no_substrates_set(self):
        with self.assertRaises(TypeError):
            self.prediction.predict_products()

    def test_predict_products_successful_prediction(self):
        # Using actual RDKit functionality to convert SMILES to Mol
        self.prediction.smarts_data = {'SMARTS': "[*:1][C:2][O:3]>>[*:1][C:2]=[O:3]"}
        self.prediction.parse_all_smarts_into_rdkit_rxn()
        # Set up a mock reaction that simulates RDKit reaction behavior
        # mock_rxn = MagicMock()
        product_mol = Chem.MolFromSmiles('CC=O')
        # mock_rxn.RunReactants.return_value = [(product_mol,)]
        # mock_rxn_from.return_value = mock_rxn

        predicted_products = self.prediction.predict_products('CCO')
        
        # # Assert that reactions are populated - incorrect, not the right structure, will return [('SMARTS', None)] even if
        # # no reactions are generated
        # self.assertGreater(len(predicted_products), 0)
        
        # Verify that reaction results as expected
        expected_product_smiles = Chem.MolToSmiles(product_mol)
        actual_product_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(predicted_products[0][1][0][0]))
        self.assertEqual(actual_product_smiles, expected_product_smiles)

    def test_add_explicit_zero_charge(self):
        # Example usage
        smarts_strings = [
            ('[CH2:1][O-:2][O-:3]>>[CH2:1][O:2][O:3]','[CH2+0:1][O-:2][O-:3]>>[CH2+0:1][O+0:2][O+0:3]'),
            ('[C:1][O-:2]>>[C:1][O:2]','[C+0:1][O-:2]>>[C+0:1][O+0:2]'),
            ('[N+1:1][Cl-1:2]>>[N:1][Cl:2]','[N+1:1][Cl-1:2]>>[N+0:1][Cl+0:2]'),
            ('[S+2:1][O-1:2]>>[S:1][O:2]','[S+2:1][O-1:2]>>[S+0:1][O+0:2]'),
        ]
    
        # Apply the function to multiple SMARTS strings
        for smarts_original, smart_updated in smarts_strings:
            calc_updated_smarts = self.prediction.add_explicit_zero_charge(smarts_original)
            self.assertEqual(smart_updated, calc_updated_smarts)

if __name__ == '__main__':
    unittest.main()