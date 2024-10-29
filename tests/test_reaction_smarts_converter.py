# This is the test of the ReactionSmartsConverter class

import unittest
from unittest.mock import MagicMock, patch
import pandas as pd

from pyrheadb.ReactionSmartsConverter import ReactionSmartsConverter

class TestReactionSmartsConverter(unittest.TestCase):

    def setUp(self):
        # Set up a test instance and other recurring objects
        self.rhea_db_mock = MagicMock()
        self.rhea_db_mock.rhea_db_version_location = 'mock/location'
        self.rhea_db_mock.df_reactions = pd.DataFrame({
            'MASTER_ID': [1, 2, 3],
            'rxn': ['C>>C', 'C.C>>CC', 'CC>>C.C'],
        })

        self.converter = ReactionSmartsConverter(self.rhea_db_mock)

    def test_init(self):
        # Test initialization and attribute setup
        self.assertEqual(self.converter.rhea_db_version_location, 'mock/location')
        self.assertIsNotNone(self.converter.rxn_mapper)
        self.assertTrue(isinstance(self.converter.df_smiles, pd.DataFrame))

    @patch('os.makedirs')
    @patch('pandas.DataFrame.to_csv')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    @patch('json.dump')
    
    def test_write_smarts_files(self, mock_json_dump, mock_open, mock_to_csv, mock_makedirs):
        # Prepare a DataFrame to test file writing
        self.converter.df_smiles = pd.DataFrame({
            'MASTER_ID': [1],
            'SMARTS': ['C>>C']
        })
        self.converter.write_smarts_files()
        mock_makedirs.assert_called_once_with('mock/location/smarts', exist_ok=True)
        mock_to_csv.assert_called_once()
        mock_open.assert_called_once()
        mock_json_dump.assert_called_once()

    def test_clean_rxn_smiles_of_redox(self):
        # Example to test the redox cleaning process
        input_row = {'rxnsmiles': 'C.C.[H]*[H]>>CC.*'}
        expected_output = 'C.C>>CC'
        result = self.converter.clean_rxn_smiles_of_redox(input_row)
        self.assertEqual(result, expected_output)

    def test_generate_smarts(self):
        # Mocking the Reaction class and its methods
        with patch('pyrheadb.Reaction', autospec=True) as mock_reaction_cls:
            mock_reaction = mock_reaction_cls.return_value
            mock_reaction.check_reaction_balance.return_value = True

            # Test 1
            self.util_test_one_smiles_case('C>>C', 'Not pattern after removing A / AH')
            # Test 2
            self.util_test_one_smiles_case('C>>CC', 'unbalanced reaction')
            # Test 3
            self.util_test_one_smiles_case("*CNC(C)=O.[H]O[H]>>CC(=O)[O-].*C[NH3+]",
                                           '[CH3:1][C:2](=[O:3])[NH:7][CH2:6][*:5].[OH2:4]>>[CH3:1][C:2](=[O:3])[O-:4].[*:5][CH2:6][NH3+:7]')
            # Test 4
            self.util_test_one_smiles_case("*C(=O)OC[C@H](O)CO.[H]O[H]>>*C(=O)[O-].OCC(O)CO.[H+]",
                                           '[OH2:10].[O:1]([CH2:2][C@H:3]([OH:4])[CH2:5][OH:6])[C:8]([*:7])=[O:9]>>[OH:1][CH2:2][CH:3]([OH:4])[CH2:5][OH:6].[*:7][C:8](=[O:9])[O-:10]')
    
    def util_test_one_smiles_case(self, smiles, expected_smarts):
        input_row = {'rxnsmiles': smiles, 'MASTER_ID': 1}
        expected_smarts = expected_smarts
        self.converter.rxn_mapper.map_many_reactions = MagicMock(return_value={'mapped_rxn': expected_smarts})
        result = self.converter.generate_smarts(input_row)
        self.assertEqual(result, expected_smarts)
        

    def test_isotope_to_star_pattern(self):
        # Test the removal of specific patterns
        input_string = '13C>>13C'
        expected_output = '*>>*'
        result = self.converter.isotope_to_star_pattern(input_string)
        self.assertEqual(result, expected_output)

# Add more tests as necessary for each method.

if __name__ == '__main__':
    unittest.main()