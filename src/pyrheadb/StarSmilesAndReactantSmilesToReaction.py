import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

from .AtomMapper import AtomMapper
from .Reaction import Reaction
from .Reaction import Reaction
from .RInChI import RInChI

# from rdkit import RDLogger
# RDLogger.DisableLog('rdApp.*')

#####################################################################################################################
#   ReactionPrediction class that takes a compound or a set of compounds as potential substrates and tests          #
#   the resulting reaction SMARTS on these compounds, gets products and checks balance for the resulting reactions. #
#####################################################################################################################

class StarSmilesAndReactantSmilesToReaction():
    
    def __init__(self):
        """
        Initializes the ReactionSmarts class which handles conversion of reaction SMILES to SMARTS.
        ReactionSmarts class handles everything related to reaction SMARTS
        - atom mapped version of reactions that can act as reaction patterns.
        The purpose of this class is to provide interface between Rhea reactions
        as structures (originally stored and .sdf and converted to SMILES by pyrheadb)
        and reaction patters as SMARTS that can be used by rdkit reaction class to generate
        new reactions.
        
        """
        #self.rhea_db_version_location = rhea_db.rhea_db_version_location
        self.rxn_mapper = AtomMapper()
        self.reactionobj = Reaction()
        self.rinchiobj = RInChI()
        #self.df_smiles = rhea_db.df_reactions.copy()

    def generate_smarts(self, rxn_smiles, debug_mode=False):
        """
        :param row: Row of dataframe of rhea reaction smiles
        :return: SMARTS of Rhea reaction with Atom Mapping numbers that can be in the next step converted to rdkit reactions
        """
        rxn_smiles = self.reactionobj.star_to_isotopic_label(rxn_smiles)
        atommap_results = self.rxn_mapper.map_one_reaction(rxn_smiles)

        if atommap_results[0] == 'tokenlength_error':
            return 'macromolecular reaction'
        
        smarts = atommap_results[0]
        if type(smarts)!=str:
            print(atommap_results)
            return None
        smarts = self.isotope_to_star_pattern(smarts)
        smarts, defined_reactant_templates, defined_product_templates = self.rework_smarts(smarts)
        return smarts, defined_reactant_templates, defined_product_templates
    
    def isotope_to_star_pattern(self, rxnsmiles):
        """
        Remove the pattern that was previously introduced to overcome inability of RXNmapper to handle dummy atoms (*)
        :param input_string: SMARTS string
        :return:
        """
        return rxnsmiles.replace('13C', '*')
    
    def smarts_to_smiles(self, rxn_smarts):
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
        rdChemReactions.SanitizeRxn(rxn)

        reactant_smiles = []
        product_smiles = []

        # Clean reactants
        for mol in rxn.GetReactants():
            reactant_smiles.append(self.clean_mol(mol))

        # Clean products
        for mol in rxn.GetProducts():
            product_smiles.append(self.clean_mol(mol))

        return '.'.join(reactant_smiles) + '>>' + '.'.join(product_smiles)
    
    # Helper to clean each molecule and turn from smarts format into smiles
    def clean_mol(self, mol):
        mol = Chem.RemoveHs(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.SanitizeMol(mol)
        smiles = Chem.MolToSmiles(mol, canonical=True)
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    
    def parse_one_smarts_into_rdkit_rxn(self, smarts):
        """
        This function transforms smarts into rdkit rxn
        https://rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html#rdkit.Chem.rdChemReactions.ReactionFromSmarts
        :param smarts: reaction SMARTS
        :return:
        """
        rxn = AllChem.ReactionFromSmarts(smarts)

        return rxn
    
    def count_substrate_product_templates_per_rxn(self, rxn):
        num_reactants = rxn.GetNumReactantTemplates()
        num_products = rxn.GetNumProductTemplates()
        return num_reactants, num_products
    
    def rework_smarts(self, smarts):
        """
        Prepare original SMARTS from rhea in Rhea order to be imported into rxn with params
        :param smarts:
        :return:
        """
        undefined_reactant_templates, \
            defined_reactant_templates, \
            undefined_product_templates, \
            defined_product_templates = self.__split_reaction_smarts(smarts)
        
        smarts = '.'.join(undefined_reactant_templates) + '>>' + \
                     '.'.join(undefined_product_templates)
        
        smarts = self.__add_explicit_zero_charge(smarts)
        return smarts, defined_reactant_templates, defined_product_templates
    
    def __split_reaction_smarts(self, smarts):
        reactant_templates = smarts.split('>>')[0].split('.')
        product_templates = smarts.split('>>')[1].split('.')
    
        undefined_reactant_templates = [t for t in reactant_templates if '*' in t]
        defined_reactant_templates = [t for t in reactant_templates if '*' not in t]
        undefined_product_templates = [t for t in product_templates if '*' in t]
        defined_product_templates = [t for t in product_templates if '*' not in t]
        
        defined_reactant_smiles = [self.clean_mol(Chem.MolFromSmiles(i)) for i in defined_reactant_templates]
        defined_product_smiles = [self.clean_mol(Chem.MolFromSmiles(i)) for i in defined_product_templates]

        return undefined_reactant_templates, \
               defined_reactant_smiles, \
               undefined_product_templates, \
               defined_product_smiles
    
    def __add_explicit_zero_charge(self, smarts):
        # This regex pattern captures elements possibly followed by hydrogens and digits, optional charges, and indices
        pattern = r'(\[)([A-Z][a-z]?H?\d*)([+-]?\d*)?(:\d+)?(\])'
    
        def replace(match):
            start_bracket = match.group(1)
            element_hydrogen = match.group(2)  # Includes element and any hydrogen notation
            charge = match.group(3)
            index = match.group(4) if match.group(4) else ''
            end_bracket = match.group(5)
        
            # Check if a charge is already specified, if not, append '+0'
            if not charge:
                charge = '+0'
        
            return f'{start_bracket}{element_hydrogen}{charge}{index}{end_bracket}'
    
        # Split reactants and products using '>>'
        reactants, products = smarts.split('>>')
        # Apply regex replacement to both reactants and products
        updated_reactants = re.sub(pattern, replace, reactants)
        updated_products = re.sub(pattern, replace, products)
    
        return '>>'.join([updated_reactants, updated_products])
    
    def load_defined_cofactor_smarts(self):
        for rxnid, smarts in self.smarts_data.items():
            self.defined_cofactor_smarts[rxnid] = self.get_defined_substrates_and_products(smarts)
            
    def get_defined_substrates_and_products(self, smarts):
        _, defined_reactant_templates, _, defined_product_templates = self.__split_reaction_smarts(smarts)
        reactant_smiles = [self.smiles_from_smarts(sm) for sm in defined_reactant_templates]
        product_smiles = [self.smiles_from_smarts(sm) for sm in defined_product_templates]
        return (reactant_smiles, product_smiles)
    
    def smiles_from_smarts(self, smarts):
        mol = Chem.MolFromSmarts(smarts)
        # Solution to remove atom map numbers as decribed at https://github.com/rdkit/rdkit/discussions/3837
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(mol)
    
    def check_substrate(self, substrate_smiles):
        """
        Set one substrate for reaction prediction.
        :param substrate_smiles: String or list of SMILES representing substrates.
        """
        if not Chem.MolFromSmiles(substrate_smiles):
            print(f'Incorrect SMILES: impossible to import as rdkit mol object: {substrate_smiles}')
            return None
        return substrate_smiles


    ##########################
    #  Prediction functions  #
    ##########################
        
    def predict_products_one_reactant_one_reaction(self, substrate_smiles, rxn):
        """
        Individual block of reaction prediction: one set of substrates + one reaction
        :param input_substrates: SMILES
        :param rxn: rdkit ChemicalReaction object
        :return: list of SMILES of products
        """
        if not self.check_substrate(substrate_smiles):
            print(f'Subtrate quality check not passed for {substrate_smiles}')
            return None
        
        if rxn.GetNumReactantTemplates() != 1:
            return None
        
        ps = rxn.RunReactants([Chem.MolFromSmiles(substrate_smiles)])
        products_all = [self.mol_products_to_smiles(product_set) for product_set in ps]
        return products_all
    
    def mol_products_to_smiles(self, product_set):
        # The molecules have not been sanitized, so it’s a good idea to at least update the valences before continuing:
        for p in product_set:
            p.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(p, Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        # Convert to SMILES
        products = [Chem.MolToSmiles(p) for p in product_set]

        return products
    
    def predict_reactions(self, input_substrates, smarts, defined_cosubstrate_smiles=[], defined_coproduct_smiles=[]):
        """
        Predict reactions for input substrates
        :param input_substrates:
        :return: df_rheaid_rxnsmiles - pandas.DataFrame with source rheaid and generated reactions per id
        """
        rxn = self.parse_one_smarts_into_rdkit_rxn(smarts)

        products = self.predict_products_one_reactant_one_reaction(input_substrates, rxn)
        
        # Handle the situation when no products are predicted
        products = [i for i in products if i is not None]
        if not products:
            print('No reactions possible for input substrates ', input_substrates)
            return None
        
        # Harmonise handling of input substrates for SMILES as str and list of SMILES for several substrates
        if isinstance(input_substrates, str):
            input_substrates=[input_substrates]
        
        reactions = []
        # Export all product sets as reactions
        for product_set in products:
            reaction_smiles = self.reaction_smiles_from_reactants_and_products(substrates=input_substrates, products=product_set, defined_cosubstrate_smiles=defined_cosubstrate_smiles, defined_coproduct_smiles=defined_coproduct_smiles)
            reactions.append((reaction_smiles))
        df_rxnsmiles = pd.DataFrame(reactions,columns=['rxnsmiles'])
        df_rxnsmiles['balance']=df_rxnsmiles['rxnsmiles'].apply(self.check_balance)
        print('QC: number of unbalanced reactions generated:', len(df_rxnsmiles[df_rxnsmiles['balance']==False]))
        return df_rxnsmiles
        
    def reaction_smiles_from_reactants_and_products(self, substrates=[], products=[], defined_cosubstrate_smiles=[], defined_coproduct_smiles=[]):
        """
        convert substrates and products lists into reaction SMILES
        :param substrates: SMILES of compounds
        :param products: SMILES of compounds
        :return: reaction SMILES
        """

        substrates=substrates+defined_cosubstrate_smiles
        products=products+defined_coproduct_smiles
        return '.'.join(substrates)+'>>'+'.'.join(products)
    
    def check_balance(self, smiles):
        """
        Use check_reaction_balance function from
        :param smiles:
        :return:
        """
        return self.reactionobj.check_reaction_balance(smiles, account_H=True)
    
    def group_predicted_reactions_based_on_rinchikey(self, df_predicted_reactions):
        """
        Careful - the direction of the reaction will be lost since RInChI is used (non-directional)
        :param df_predicted_reactions: pandas.DataFrame of predicted reactions
        :return: filtered df
        """

        print('Calculating Reaction InChiKeys')
        df_predicted_reactions[['RInChI','Web-RInChIKey']] = \
             df_predicted_reactions.apply(lambda row: self.rinchiobj.error_handle_wrap_rinchi(row['rxnsmiles']), axis=1, result_type='expand')
        result = df_predicted_reactions.groupby('Web-RInChIKey').agg({'rxnsmiles': 'first'}).reset_index()
        return result
        