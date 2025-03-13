import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import json
import os
import re

# from rdkit import RDLogger
# RDLogger.DisableLog('rdApp.*')

from .Reaction import Reaction
from .RInChI import RInChI

#####################################################################################################################
#   ReactionPrediction class that takes a compound or a set of compounds as potential substrates and tests          #
#   the resulting reaction SMARTS on these compounds, gets products and checks balance for the resulting reactions. #
#####################################################################################################################
class ReactionPrediction:
    def __init__(self, rhea_db):
        """
        Initializes the ReactionPrediction class.
        """
        self.rdb = rhea_db
        self.rhea_db_version_location = rhea_db.rhea_db_version_location 
        self.rxnclass = Reaction()
        self.account_H_in_balance = False
        self.star_smarts_only = True # use only template/changeable part (compound with *) or full smarts (with defined cofactors, set this to False)
        self.rdkit_stereo_rxn_data = dict()
        self.rdkit_flat_rxn_data = dict()  # not stereo - remove all @ from reactions
        self.defined_cofactor_smarts=dict()
        self.load_smarts_data()
    
    def set_account_H_in_balance(self, new_account_H_in_balance_option):
        """
        Reset the defaut option for handling hydrogen in calculating reaction balance
        :param account_H_in_balance_option: default False
        :return:
        """
        self.account_H_in_balance = new_account_H_in_balance_option
    
    def set_star_smarts_only(self, new_star_smarts_only_option):
        self.star_smarts_only = new_star_smarts_only_option
        
    def load_smarts_data(self):
        """
        Correct the location of SMARTS json
        self.smarts_data: DataFrame or dictionary containing reaction SMARTS patterns.
        """
        smarts_directory = os.path.join(self.rhea_db_version_location, 'smarts')
        with open(os.path.join(smarts_directory,'rheaSmarts.json')) as j:
            self.smarts_data = json.load(j)
        self.parse_all_smarts_into_rdkit_rxn()
        if self.star_smarts_only:
            self.load_defined_cofactor_smarts()
    
    def parse_all_smarts_into_rdkit_rxn(self):
        """
        Parse all SMARTS extracted from Rhea reactions and create rdkit reactions
        class rdkit.Chem.rdChemReactions.ChemicalReaction((object)self)
        see https://rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html for more info about the class
        :return:
        """
        
        for rxnid, smarts in self.smarts_data.items():
            try:
                self.rdkit_stereo_rxn_data[rxnid] = self.parse_one_smarts_into_rdkit_rxn(smarts)
                self.rdkit_flat_rxn_data[rxnid] = self.parse_one_smarts_into_rdkit_rxn(smarts.replace('@', '')) # could be not enough, need to see use cases when stereo info is still left
            except:
                print(rxnid, smarts)
            
    def parse_one_smarts_into_rdkit_rxn(self, smarts):
        """
        This function transforms smarts into rdkit rxn
        https://rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html#rdkit.Chem.rdChemReactions.ReactionFromSmarts
        :param smarts: reaction SMARTS
        :return:
        """
        rxn = AllChem.ReactionFromSmarts(
            self.rework_smarts(smarts)
        )
        # Chem.rdChemReactions.SanitizeRxn(rxn) # commented since does not help
        # The operations carried out by default are:
        #         fixRGroups(): sets R group labels on mapped dummy atoms when possible
        #         fixAtomMaps(): attempts to set atom maps on unmapped R groups
        #         adjustTemplate(): calls adjustQueryProperties() on all reactant templates
        #         fixHs(): merges explicit Hs in the reactant templates that don’t map to heavy atoms
        # None of these operations ensures that the charge is interpreted correctly


        return rxn
    
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
        
        if self.star_smarts_only == False:
            # full_rxn_smarts
            smarts = '.'.join(undefined_reactant_templates + defined_reactant_templates) + '>>' + \
                     '.'.join(undefined_product_templates + defined_product_templates)
        if self.star_smarts_only == True:
            # only_changeable_part_rxn_smarts
            smarts = '.'.join(undefined_reactant_templates) + '>>' + \
                     '.'.join(undefined_product_templates)
        
        smarts = self.add_explicit_zero_charge(smarts)
        return smarts
    
    def __split_reaction_smarts(self, smarts):
        reactant_templates = smarts.split('>>')[0].split('.')
        product_templates = smarts.split('>>')[1].split('.')
    
        undefined_reactant_templates = [t for t in reactant_templates if '*' in t]
        defined_reactant_templates = [t for t in reactant_templates if '*' not in t]
        undefined_product_templates = [t for t in product_templates if '*' in t]
        defined_product_templates = [t for t in product_templates if '*' not in t]
        
        return undefined_reactant_templates, \
               defined_reactant_templates, \
               undefined_product_templates, \
                defined_product_templates
    

    def add_explicit_zero_charge(self, smarts):
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

    def predict_products(self, input_substrates):
        """
        Predicts products using defined and undefined templates.
        input_substrates can be list (of SMILES) or str (SMILES).
        If list, the total number should not be more than 3! (biochemically impossible)
        For list, RunReactants rdkit function is used.
        For str, Run Reactant rdkit function is used.
        return: reactions_to_products - list of rheaids with product versions
        """
        if input_substrates is None:
            raise ValueError("Substrates not set. Please set substrates before predicting products.")
        
        # Stereochemistry Handling
        # The handling of stereochemistry (@ symbols in SMILES) based on whether the substrate has stereochemistry
        # is a specific but important detail, especially in pharmaceutical applications.
        # This should also be integrated to ensure the correct chemical context is maintained.
        if any(['@' in i for i in input_substrates]): rxn_used = self.rdkit_stereo_rxn_data
        else: rxn_used = self.rdkit_flat_rxn_data

        # Produce list of reaction ids with product options as
        # [(rheaid1:[products1-1, products1-2]), (rheaid2:[products2-1, products2-2]), ...]
        if isinstance(input_substrates, str):
            reactions_to_products = self.predict_products_one_reactant(input_substrates, rxn_used)
        if isinstance(input_substrates, list):
            reactions_to_products = self.predict_products_several_reactants(input_substrates, rxn_used)
        return reactions_to_products

    def predict_products_one_reactant(self, substrate_smiles, rxn_used):
        """
        predict potential products of singe reactants
        :param input_substrates: substrates of the reaction (SMILES)
        :param rxn_used: rdkit ChemicalReaction objects - selected is stereo or flat
        :return: list of products
        """
        if not self.check_substrate(substrate_smiles):
            print(f'Subtrate quality check not passed for {substrate_smiles}')
            return None
        rheaid_to_products = []
        for rheaid, rxn in rxn_used.items():
            try:
                products = self.predict_products_one_reactant_one_reaction(substrate_smiles, rxn)
            except Exception as e:
                print('Problematic rheaid', rheaid)
                print(e)
            rheaid_to_products.append((rheaid, products))
        return rheaid_to_products
    
    def predict_products_one_reactant_one_reaction(self, substrate_smiles, rxn):
        """
        Individual block of reaction prediction: one set of substrates + one reaction
        :param input_substrates: SMILES
        :param rxn: rdkit ChemicalReaction object
        :return: list of SMILES of products
        """
        if rxn.GetNumReactantTemplates() != 1:
            return None
        ps = rxn.RunReactants([Chem.MolFromSmiles(substrate_smiles)])
        products_all = [self.mol_products_to_smiles(product_set) for product_set in ps]
        return products_all
    
    def predict_products_several_reactants(self, input_substrates, rxn_used):
        """
        predict potential products of mixing several reactants (in one reaction, so theoretically the number should not be more than 2)
        :param input_substrates: substrates of the reaction (list of SMILES)
        :param rxn_used: rdkit ChemicalReaction objects - selected is stereo or flat
        :return: list of products
        """
        if not all([self.check_substrate(substrate_smiles) for substrate_smiles in input_substrates]):
            print(f'Subtrate quality check not passed for {input_substrates}')
            return None
        rheaid_to_products = []
        for rheaid, rxn in rxn_used.items():
            products = self.predict_products_several_reactants_one_reaction(input_substrates, rxn)
            rheaid_to_products.append((rheaid, products))
        return rheaid_to_products
    
    def predict_products_several_reactants_one_reaction(self, input_substrates, rxn):
        """
        Individual block of reaction prediction: one set of substrates + one reaction
        :param input_substrates: list of SMILES
        :param rxn: rdkit ChemicalReaction object
        :return: list of SMILES of products
        """
        if rxn.GetNumReactantTemplates() != len(input_substrates):
            return None
        ps = rxn.RunReactants([Chem.MolFromSmiles(substrate) for substrate in input_substrates])
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
    
    def predict_reactions(self, input_substrates):
        """
        Predict reactions for input substrates
        :param input_substrates:
        :return: df_rheaid_rxnsmiles - pandas.DataFrame with source rheaid and generated reactions per id
        """
        rheaid2products = self.predict_products(input_substrates)
        
        # Handle the situation when no products are predicted
        rheaid2products = [i for i in rheaid2products if i[1] is not None]
        if not rheaid2products:
            print('No reactions possible for input substrates ', input_substrates)
            return None
        
        # Harmonise handling of input substrates for SMILES as str and list of SMILES for several substrates
        if isinstance(input_substrates, str):
            input_substrates=[input_substrates]
        
        rheaid2reaction = []
        # Export all product sets as reactions
        for rheaid, products in rheaid2products:
            for product_set in products:
                reaction_smiles = self.reaction_smiles_from_reactants_and_products(rheaid=rheaid, substrates=input_substrates, products=product_set)
                rheaid2reaction.append((rheaid, reaction_smiles))
        df_rheaid_rxnsmiles = pd.DataFrame(rheaid2reaction,columns=['rheaid', 'rxnsmiles'])
        df_rheaid_rxnsmiles['balance']=df_rheaid_rxnsmiles['rxnsmiles'].apply(self.check_balance)
        print('QC: number of unbalanced reactions generated:', len(df_rheaid_rxnsmiles[df_rheaid_rxnsmiles['balance']==False]))
        return df_rheaid_rxnsmiles
        
    def reaction_smiles_from_reactants_and_products(self, rheaid='rheaid', substrates=[], products=[]):
        """
        convert substrates and products lists into reaction SMILES
        :param substrates: SMILES of compounds
        :param products: SMILES of compounds
        :return: reaction SMILES
        """
        if self.star_smarts_only==True:
            substrates=substrates+self.defined_cofactor_smarts[rheaid][0]
            products=products+self.defined_cofactor_smarts[rheaid][1]
        return '.'.join(substrates)+'>>'+'.'.join(products)
    
    def check_balance(self, smiles):
        """
        Use check_reaction_balance function from
        :param smiles:
        :return:
        """
        return self.rxnclass.check_reaction_balance(smiles, account_H=self.account_H_in_balance)
    
    def group_predicted_reactions_based_on_rinchikey(self, df_predicted_reactions):
        """
        Careful - the direction of the reaction will be lost since RInChI is used (non-directional)
        :param df_predicted_reactions: pandas.DataFrame of predicted reactions
        :return: filtered df
        """
        rinchiobj = RInChI()

        print('Calculating Reaction InChiKeys')
        df_predicted_reactions[['RInChI','Web-RInChIKey']] = \
             df_predicted_reactions.apply(lambda row: rinchiobj.error_handle_wrap_rinchi(row['rxnsmiles']), axis=1, result_type='expand')
        result = df_predicted_reactions.groupby('Web-RInChIKey').agg({
            'rheaid': lambda x: ', '.join(x),
            'rxnsmiles': 'first'
        }).reset_index()
        return result
        
    def count_substrate_templates(self):
        counts_reactants = []
        counts_products = []
        counts_together = []
        for rheaid, rxn in self.rdkit_stereo_rxn_data.items():
            num_reactants, num_products = self.count_substrate_product_templates_per_rxn(rxn)
            if num_reactants>3 or num_products>3:
                print('4+ substrate templates in', rheaid)
            counts_products.append(num_products)
            counts_reactants.append(num_reactants)
            counts_together.extend([num_reactants, num_products])
        return counts_reactants, counts_products, counts_together
        
    def count_substrate_product_templates_per_rxn(self, rxn):
        num_reactants = rxn.GetNumReactantTemplates()
        num_products = rxn.GetNumProductTemplates()
        return num_reactants, num_products