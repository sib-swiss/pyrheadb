import unittest
import os

from rdkit import Chem

from pyrheadb.Reaction import Reaction

class TestReaction(unittest.TestCase):
	def setUp(self):
		self.reaction = Reaction()
	
	def test_check_reaction_balance_smiles_balanced(self):
		reaction_string = "C(C)C>>C(C)C"
		result = self.reaction.check_reaction_balance(reaction_string, format='smiles')
		self.assertTrue(result)
	
	def test_check_reaction_balance_smiles_unbalanced(self):
		reaction_string = "C(C)C>>CC(C)C"
		result = self.reaction.check_reaction_balance(reaction_string, format='smiles')
		self.assertFalse(result)
	
	def test_check_reaction_balance_inchi_balanced(self):
		reaction_string = "InChI=1S/C2H6/c1-2/h1-2H3>>InChI=1S/C2H6/c1-2/h1-2H3"
		result = self.reaction.check_reaction_balance(reaction_string, format='inchi')
		self.assertTrue(result)
	
	def test_check_reaction_balance_inchi_balanced_2(self):
		reaction_string = "InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)/p-3/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1" \
						  " + " \
						  "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3 " \
						  ">> " \
						  "InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)/p-3/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1" \
						  " + " \
						  "InChI=1S/C3H6O2/c1-3(5)2-4/h4H,2H2,1H3"
		result = self.reaction.check_reaction_balance(reaction_string, format='inchi')
		self.assertTrue(result)
		
	def test_check_reaction_balance_inchi_unbalanced(self):
		reaction_string = "InChI=1S/C2H6/c1-2/h1-2H3>>InChI=1S/C3H8/c1-3-2/h3H,1-2H3"
		result = self.reaction.check_reaction_balance(reaction_string, format='inchi')
		self.assertFalse(result)
	
	def test_calculate_atom_counts(self):
		# Testing with a simple molecule: methane (CH4)
		methane = Chem.MolFromSmiles('C')
		counts = self.reaction.calculate_atom_counts([methane])
		expected_counts = {6: 1}  # Carbon atom count
		self.assertEqual(counts, expected_counts)
	
	def test_read_rxnfile(self):
		# Construct the path to the .rxn file
		current_dir = os.path.dirname(__file__)  # Gets the directory where this test file resides
		rxnfile_path = os.path.join(current_dir, 'test_files', '10001.rxn')  # Builds the path to the .rxn file
	
		rxn = self.reaction.read_rxnfile(rxnfile_path)
		self.assertIsNotNone(rxn)
		self.assertTrue(rxn.IsInitialized())

	def test_star_to_isotopic_label(self):
		# Test the reaction SMILES modification method
		input_smiles = '[1*]C[2*]>>[3*]CC[4*]'
		expected_output = '[13C]C[13C]>>[13C]CC[13C]'
		result = self.reaction.star_to_isotopic_label(input_smiles)
		self.assertEqual(result, expected_output)

if __name__ == '__main__':
	unittest.main()