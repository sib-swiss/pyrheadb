import unittest

from pyrheadb.RInChI import RInChI

class TestRInChI(unittest.TestCase):
    def setUp(self):
        self.rinchi = RInChI()
	
    def test_inchis_from_rinchi(self):
        
        rinchi_nostructs = "RInChI=1.00.1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3/t3-,4?/m0/s1<>C4H9BrO/c1-3(5)4(2)6/h3-4,6H,1-2H3/t3-,4+/m1/s1!Na.H2O/h;1H2/q+1;/p-1/d-"
        inchi_components = self.rinchi.inchis_from_rinchi(rinchi_nostructs, "")
        #Expected result:
        # {'Direction': '-', 
        # 'No-Structures': [0, 0, 0], 
        # 'Reactants': [('InChI=1S/C4H9BrO/c1-3(5)4(2)6/h3-4,6H,1-2H3/t3-,4+/m1/s1', ''), ('InChI=1S/Na.H2O/h;1H2/q+1;/p-1', '')], 
        # 'Products': [('InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3/t3-,4?/m0/s1', '')], 
        # 'Agents': []}
        self.assertEqual(inchi_components['Direction'], '-')
        self.assertEqual(inchi_components['No-Structures'], [0, 0, 0])
        self.assertEqual(inchi_components['Reactants'], [('InChI=1S/C4H9BrO/c1-3(5)4(2)6/h3-4,6H,1-2H3/t3-,4+/m1/s1', ''), ('InChI=1S/Na.H2O/h;1H2/q+1;/p-1', '')])
        self.assertEqual(inchi_components['Products'], [('InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3/t3-,4?/m0/s1', '')])
        self.assertEqual(inchi_components['Agents'], [])
