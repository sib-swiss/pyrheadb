import unittest
from unittest.mock import MagicMock, patch
from pyrheadb.AtomTransitionNetwork import AtomTransitionNetwork
from pyrheadb.RheaDB import RheaDB
from rdkit import Chem

class TestAtomTransitionNetwork(unittest.TestCase):
    def setUp(self):
        rhea_db = RheaDB(rhea_versions_folder_location='/scratch/rheadata/', rhea_version=130)
        self.atn = AtomTransitionNetwork(rheadata=rhea_db)
        self.smiles_test = 'CCCCC(N)=O.[H]O[H]>>CCCCC(=O)[O-].[H][N+]([H])([H])[H]'
        self.dict_smiles_to_compid = {
            'CCCCC(N)=O': 'comp:1', '[H]O[H]': 'comp:2',
            'CCCCC(=O)[O-]': 'comp:3', '[H][N+]([H])([H])[H]': 'comp:4'
        }
        self.dict_compid_to_smiles = {
            'comp:1': 'CCCCC(N)=O', 'comp:2': '[H]O[H]',
            'comp:3': 'CCCCC(=O)[O-]', 'comp:4': '[H][N+]([H])([H])[H]'
        }
        self.dict_compid_to_mol = {
            compid: Chem.MolFromSmiles(smile)
            for compid, smile in self.dict_compid_to_smiles.items()
        }
        self.atn.map_one_reaction = MagicMock(return_value={
            'confidence': 0.9965346908959984,
            'mapped_rxn': '[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[NH2:8].[OH2:7]>>[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[O-:7].[NH4+:8]'
        })
    
    def test_number_atoms_in_dict(self):
        result = self.atn.number_atoms_in_dict(self.dict_compid_to_smiles)
        expected_numbered_smiles = {
            'comp:1': '[CH3:1][CH2:2][CH2:3][CH2:4][C:5]([NH2:6])=[O:7]',
            'comp:2': '[OH2:1]',
            'comp:3': '[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[O-:7]',
            'comp:4': '[NH4+:1]'
        }
        self.assertEqual(result, expected_numbered_smiles)
    
    def test_match_mapped_rxn_to_comp_ids(self):
        mapped_rxn = self.atn.map_one_reaction(self.smiles_test)['mapped_rxn']
        mapped_data = {'mapped_rxn':mapped_rxn, 'relevant_comp_dict':self.dict_compid_to_mol}
        result = self.atn.match_mapped_rxn_to_comp_ids(mapped_data)
        expected_comp_ids = 'comp:1.comp:2>>comp:3.comp:4'
        self.assertEqual(result, expected_comp_ids)
    
    def test_match_compound_atoms(self):
        equation_mapped = self.atn.map_one_reaction(self.smiles_test)['mapped_rxn']
        mapped_data = {'mapped_rxn':equation_mapped, 'relevant_comp_dict':self.dict_compid_to_mol}
        comp_ids = self.atn.match_mapped_rxn_to_comp_ids(mapped_data)
        dict_compid_to_numbered_smiles = self.atn.number_atoms_in_dict(self.dict_compid_to_smiles)
        mapped_data={'mapped_comp_ids':comp_ids,'mapped_rxn':equation_mapped}
        result = self.atn.match_compound_atoms(mapped_data, dict_compid_to_numbered_smiles)
        self.assertTrue(isinstance(result, list))  # assuming list of atom edges
    
    # Add more tests for other methods

if __name__ == '__main__':
    unittest.main()
