from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
import io
import os
from PIL import Image

class Visualization():
    
    def __init__(self, image_folder=None):
        """
        Initializes the Visualization instance with an optional image folder path.
        :param image_folder: Optional. Folder where the generated images will be stored. Default is None.
        """
        self.image_folder = image_folder if image_folder else "default_image_folder"
        os.makedirs(self.image_folder, exist_ok=True)
		
    def set_image_output_folder(self, image_folder):
        """
        :param image_folder: folder where the generated images are going to be stored
        :return:
        """
        self.image_folder = image_folder
        os.makedirs(self.image_folder, exist_ok=True)

    def draw_reaction_smarts(self, rxn_smarts, rxn_name, bw=False):
        """
        :param rxn_smarts: SMARTS of the reaction
        :param rxn_name: name of the reaction (e.g. Rhea ID)
        :param bw: black and white, True / False
        :return: image is being saved to the defined folder
        """
        smarts = rxn_smarts.replace('][', ']-[').replace(')[', ')-[')
        rxn = rdChemReactions.ReactionFromSmarts(smarts)
        drawer = rdMolDraw2D.MolDraw2DCairo(700, 300)
        if bw:
            drawer.drawOptions().useBWAtomPalette()
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        bio = io.BytesIO(drawer.GetDrawingText())
        img = Image.open(bio)
        img.save(os.path.join(f'{self.image_folder}',f'rxn_smarts_{rxn_name}.png'))
        
    def draw_compound_smarts(self, compounds_smarts, name, bw=False):
        """
    
        :param compounds_smarts: SMARTS of compounds separated by dot (.)
        :param name: name of the file that will be found in reaction_png
        :param bw: black and white, True / False
        :return: image is being saved to the defined folder
        """
        compounds_smarts = compounds_smarts.split('.')
        ms = [Chem.MolFromSmiles(x) for x in compounds_smarts]
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 180, 200, 180)
        if bw:
            drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecules(ms)
        drawer.FinishDrawing()
        bio = io.BytesIO(drawer.GetDrawingText())
        img = Image.open(bio)
        img.save(os.path.join(f'{self.image_folder}',f'cmps_smarts_{name}.png'))

    def example(self):
        rxn_smarts = "[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[NH2:8].[OH2:7]\
        >>[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[O-:7].[NH4+:8]"
        self.draw_reaction_smarts(rxn_smarts, 'test_rxn')
        
        self.draw_compound_smarts("[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[NH2:7].[OH2:1]", 'test_reactants')
        self.draw_compound_smarts("[CH3:1][CH2:2][CH2:3][CH2:4][C:5](=[O:6])[O-:7].[NH4+:1]", 'test_products')