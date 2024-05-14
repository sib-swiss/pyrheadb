from rdkit import Chem
from rdkit.Chem import AllChem
import json

# from rdkit import RDLogger
# RDLogger.DisableLog('rdApp.*')

from .Reaction import Reaction

#####################################################################################################################
# ReactionPrediction class that takes a compound or a set of compounds as potential substrates and tests
# the resulting reaction SMARTS on these compounds, gets products and checks balance for the resulting reactions.
#####################################################################################################################
class ReactionPrediction:
    def __init__(self):
        """
        Initializes the ReactionPrediction class.
        """
        self.substrates = None
        self.reactions = []

    def load_smarts_data(self):
        """
        Correct the location of SMARTS json
        self.smarts_data: DataFrame or dictionary containing reaction SMARTS patterns.
        """
        with open('data/rheaSmarts.json') as j:
            self.smarts_data = json.load(j)
            
    def set_substrates(self, substrate_smiles):
        """
        Set the substrates for reaction prediction.
        :param substrate_smiles: String or list of SMILES representing substrates.
        """
        if isinstance(substrate_smiles, str):
            substrate_smiles = [substrate_smiles]
        correct_substrate_smiles = []
        for smiles in substrate_smiles:
            if not Chem.MolFromSmiles(smiles):
                print(f'Incorrect SMILES: impossible to import as rdkit mol object{smiles}')
            else:
                correct_substrate_smiles.append(smiles)
        self.substrates = correct_substrate_smiles

    def predict_products(self):
        """
        Predicts products using defined and undefined templates.
        """
        if self.substrates is None:
            raise ValueError("Substrates not set. Please set substrates before predicting products.")
    
        for smarts in self.smarts_data['SMARTS']:
            reactant_templates = smarts.split('>>')[0].split('.')
            product_templates = smarts.split('>>')[1].split('.')
        
            undefined_reactant_templates = [t for t in reactant_templates if '*' in t]
            defined_reactant_templates = [t for t in reactant_templates if '*' not in t]
            undefined_product_templates = [t for t in product_templates if '*' in t]
            defined_product_templates = [t for t in product_templates if '*' not in t]
        
            rxn_smarts = '.'.join(undefined_reactant_templates + defined_reactant_templates) + '>>' + \
                         '.'.join(undefined_product_templates + defined_product_templates)
                
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            for substrate in self.substrates:
                # Stereochemistry Handling
                # The handling of stereochemistry (@ symbols in SMILES) based on whether the substrate has stereochemistry
                # is a specific but important detail, especially in pharmaceutical applications.
                # This should also be integrated to ensure the correct chemical context is maintained.
                if not '@' in substrate:
                    rxn = AllChem.ReactionFromSmarts(rxn_smarts.replace('@', ''))
                ps = rxn.RunReactants((Chem.MolFromSmiles(substrate),))
                for p_set in ps:
                    products = [Chem.MolToSmiles(p) for p in p_set]
                    self.reactions.append(
                        {'substrate': substrate, 'products': products, 'smarts': rxn_smarts})

    def check_balance(self):
        """
        Check the balance of the reactions based on atom counts.
        """
        balanced_reactions = []
        for reaction_info in self.reactions:
            substrate_smiles = reaction_info['substrate']
            for product_smiles in reaction_info['products']:
                reaction_smiles = substrate_smiles + '>>' + product_smiles
                reaction_obj = Reaction()
                balance = reaction_obj.check_reaction_balance(reaction_smiles)
                if balance:
                    balanced_reactions.append({
                        'substrate': substrate_smiles,
                        'products': product_smiles,
                        'smarts': reaction_info['smarts']
                    })
        return balanced_reactions

    def get_reactions(self):
        """
        Get a list of reactions with substrates, products, and the used SMARTS.
        """
        return self.reactions
    
    def get_balanced_reactions(self):
        return self.check_balance()
        
    def output_results(self, filepath):
        """
        Outputs the results to a file.
        """
        with open(filepath, 'w') as file:
            for reaction in self.reactions:
                file.write(f"{reaction['substrate']}>>{reaction['products']}\n")


###########################
# Apply SMARTS
############################


# import json
#
#
#
# def enumerate_one_reaction_template(rhea_smarts, undefined_reactant_smiles, direction='forward'):
#
#     # print('rhea_smarts',rhea_smarts)
#     reactant_templates = rhea_smarts.split('>>')[0].split('.')
#     defined_reactant_templates = [templ for templ in reactant_templates if not '*' in templ]
#     undefined_reactant_templates = [templ for templ in reactant_templates if '*' in templ]
#
#     # print("defined_reactant_templates, undefined_reactant_templates")
#     # print(defined_reactant_templates, undefined_reactant_templates)
#
#     product_templates = rhea_smarts.split('>>')[1].split('.')
#     defined_product_templates = [templ for templ in product_templates if not '*' in templ]
#     undefined_product_templates = [templ for templ in product_templates if '*' in templ]
#
#     if direction == 'forward':
#         rxn_smarts = '.'.join(undefined_reactant_templates + defined_reactant_templates) + '>>' + \
#                      '.'.join(undefined_product_templates + defined_product_templates)
#     # print("rxn_smarts direction")
#     # print(rxn_smarts)
#
#     if not '@' in undefined_reactant_smiles:
#         stereo = False
#
#     if stereo == False:
#         rxn_smarts = rxn_smarts.replace('@', '')
#
#     rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
#
#     # Assembling the input for forward enumeration of the library
#     # print("undefined_reactant_smiles")
#     # print(undefined_reactant_smiles)
#     # print("defined_reactant_templates")
#     # print(defined_reactant_templates)
#     undefined_reactant_mols = [Chem.MolFromSmiles(undefined_reactant_smiles)]
#     defined_reactant_mols = [Chem.MolFromSmiles(m) for m in defined_reactant_templates]
#     reactant_mols = undefined_reactant_mols + defined_reactant_mols
#     # print(lib)
#
#     # Enumerating forward library
#     results_total = rxn.RunReactants(reactant_mols)
#     options_products_smiles = []
#     for result_option in results_total:
#         try:
#             product_smiles = [Chem.MolToSmiles(m) for m in result_option]
#             product_smiles.sort()
#             options_products_smiles.append('.'.join(product_smiles))
#         except:
#             pass
#     options_products_smiles = set(options_products_smiles)
#     reactions = []
#     if options_products_smiles:
#         reactant_smiles = [Chem.MolToSmiles(m) for m in reactant_mols]
#         reactant_smiles.sort()
#         reactant_smiles = '.'.join(reactant_smiles)
#         for option in options_products_smiles:
#             reaction_smiles = reactant_smiles + '>>' + option
#             reaction = Reaction()
#             balance = reaction.check_reaction_balance(reaction_smiles)
#             if balance:
#                 reactions.append(reaction)
#     return reactions
#
# def enumerate_one_compound(compound):
#     with open(f'individual_compound_results/{compound}.tsv', 'w') as w:
#         w.write('rheaid\treaction_option\n')
#         # count = 0
#         # count_products = 0
#         for rheaid, rhea_smarts in data_rxn_smarts.items():
#             # print(rheaid)
#             if rhea_smarts:
#                 reactants = rhea_smarts.split('>>')[0].split('.')
#                 star_reactants = [r for r in reactants if '*' in r]
#                 if len(star_reactants) == 1:
#                     # print(star_reactants)
#                     reactions = enumerate_one_reaction_template(rhea_smarts, compound, direction='forward')
#                     # count +=1
#                     for reaction in reactions:
#                         w.write(f"{rheaid}\t{reaction}\n")
#                         # count_products += 1
#         # print(count)
#         # print(count_products)
#
# undefined_reactant_smiles = ['CC(O)CO', 'CCC(O)CO', 'NCC(O)CO']
#
# enumerate_one_compound(undefined_reactant_smiles[0])
