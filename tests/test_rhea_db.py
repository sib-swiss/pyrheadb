import unittest
from unittest.mock import patch, MagicMock
import os

import pandas as pd
from pyrheadb.RheaDB import RheaDB  # Adjust this import according to your actual module structure.

class TestRheaDB(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures, if any."""
        self.env_var = 'RHEADB_LOC'
        self.current_dir = os.path.dirname(__file__)  # Gets the directory where this test file resides
        test_rhea_db_path = os.path.join(self.current_dir, 'test_rhea_db')
        os.environ[self.env_var] = test_rhea_db_path
        self.rhea_db_version = 130
        self.rhea_db = RheaDB(rhea_versions_folder_location=test_rhea_db_path, rhea_version=self.rhea_db_version)

    def test_init_sets_properties_correctly(self):
        """Test that the __init__ method sets properties correctly."""
        self.assertEqual(self.rhea_db.rhea_versions_folder_location, os.getenv(self.env_var))
        self.assertEqual(self.rhea_db.rhea_db_version, self.rhea_db_version)


    @patch('urllib.request.urlretrieve')
    def test_get_current_rhea_version(self, mock_urlretrieve):
        """Test fetching the current Rhea version from the properties file."""
        mock_urlretrieve.return_value = None
        version_data = "rhea-version=132"  # String mimicking file content
        with patch('builtins.open', unittest.mock.mock_open(read_data=version_data)) as mock_file:
            version = self.rhea_db._RheaDB__get_current_rhea_version()
            mock_file.assert_called_with('rhea-release.properties')
            self.assertEqual(version, 132)

    # @patch('tarfile.open')
    # def test_extract_rhea_tar(self, mock_tar_open):
    #     """Test the extraction of a Rhea tar file."""
    #     # Create a mock for the tar file object
    #     mock_tar = MagicMock()
    #     mock_tar.__enter__ = MagicMock(return_value=mock_tar)
    #     mock_tar.__exit__ = MagicMock(return_value=None)
    #     mock_tar.extractall = MagicMock()
    #     mock_tar_open.return_value = mock_tar
    
    #     # Call the method under test
    #     self.rhea_db._RheaDB__extract_rhea_tar()
    
    #     # Check that tarfile.open was called with the correct parameters
    #     mock_tar_open.assert_called_once_with(f'{self.rhea_db.rhea_db_version}.tar.bz2', 'r:bz2')
    
    #     # Ensure that the extractall method was called correctly
    #     mock_tar.extractall.assert_called_once_with(f'{self.rhea_db.rhea_versions_folder_location}', 'rhea-versions')

    @patch('pandas.read_csv')
    def test_read_rhea_files(self, mock_read_csv):
        """Test loading Rhea files into dataframes."""
        self.rhea_db._RheaDB__read_rhea_files()
        self.assertTrue(mock_read_csv.called)
        self.assertEqual(mock_read_csv.call_count, 4)  # Should match the number of files being loaded

    @patch('pandas.DataFrame.to_csv')
    def test_generate_reaction_compound_names_file(self,  mock_to_csv):
        """Test generating the file for reaction compound names."""
        self.rhea_db.df_smiles_master_id=pd.read_csv(
            os.path.join(self.current_dir, 'test_rhea_db', 'rhea-versions', f'{self.rhea_db_version}', 'tsv', 'rhea-reaction-smiles-master-id.tsv'), sep='\t')
        self.rhea_db.df_chebi_cmpname = pd.read_csv(
            os.path.join(self.current_dir, 'test_rhea_db', 'rhea-versions', f'{self.rhea_db_version}', 'tsv', 'chebi_cmpname.tsv'), sep='\t',
                               names=['chebiid','cmpname'])
        with patch('os.path.exists') as mock_exists:
            mock_exists.return_value = False
            self.rhea_db._RheaDB__generate_reaction_compound_names_file()
            mock_to_csv.assert_called_once()

    # More tests can be added for the remaining methods following similar patterns.

if __name__ == '__main__':
    unittest.main()