from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
import io
import os
from PIL import Image
import matplotlib.pyplot as plt

class Visualization():
    
    def __init__(self, image_folder=None):
        """
        Initializes the Visualization instance with an optional image folder path.
        :param image_folder: Optional. Folder where the generated images will be stored. Default is None.
        """
        self.set_image_output_folder(image_folder)
        
    def set_image_output_folder(self, image_folder):
        """
        :param image_folder: folder where the generated images are going to be stored
        :return:
        """
        self.image_folder =  image_folder if image_folder else "default_image_folder"
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

    def create_and_save_histogram(df, column_name, bins=10, output_image_path='histogram.png'):
        """
        Creates a histogram of the specified column in the DataFrame and saves it as an image.
        
        Parameters:
            df (pd.DataFrame): The DataFrame containing the data.
            column_name (str): The name of the column to create the histogram for.
            bins (int): Number of bins for the histogram. Default is 10.
            output_image_path (str): The path where the histogram image will be saved. Default is 'histogram.png'.
            
        Returns:
            None

        Raises:
            ValueError: If the specified column does not exist in the DataFrame.
        """
        # Ensure the column exists in the DataFrame
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' does not exist in the DataFrame.")
        
        # Filter values between 0 and 1
        data = df[column_name].dropna()  # Remove NaN values
        data = data[(data >= 0) & (data <= 1)]
        
        # Create the histogram
        plt.figure(figsize=(8, 6))
        plt.hist(data, bins=bins, edgecolor='black', color='blue')
        plt.title(f'Histogram of {column_name}')
        plt.xlabel(column_name)
        plt.ylabel('Frequency')
        
        # Save the histogram as an image
        plt.savefig(output_image_path)
        plt.close()  # Close the figure to free memory