import pandas as pd
import os

from cobra import Model, Reaction, Metabolite

from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import MolFromSmiles

from tqdm import tqdm
tqdm.pandas()

class RheaCobraModel():
	
	def __init__(self, model_name, rdb, solver=None):
		self.model_name = model_name
		
		# preprocessing dfs
		self.rdb = rdb
		self.df_smiles_chebi_equation_master_id = self.rdb.df_reactions.copy()
		self.assign_subsystems_to_reactions()
		self.assign_names_to_reactions()
		self.set_bounds()
		
		# creating the model
		self.model = Model(model_name)
		if solver:
			self.model.solver = solver
		self.add_all_rhea_reactions_to_model()
		self.add_boundaries()
	
	def assign_subsystems_to_reactions(self):
		"""
		Add subsystems
		TO DO:
		The logic has to be improved based on the SPARQL query to UniProt to identify compartments
		:return:
		"""
		self.df_smiles_chebi_equation_master_id['reaction_subsystem']='dummy_subsystem'
		
	def set_bounds(self):
		"""
		Add bounds to reactions
		TO DO:
		values of lower and upper bounds should be decided in the further modeling -
		at the moment most permissive bounds are introduced
		:return:
		"""
		self.df_smiles_chebi_equation_master_id['lower_bound'] = -1000.
		self.df_smiles_chebi_equation_master_id['upper_bound'] = 1000.
		
	def assign_names_to_reactions(self):
		"""
		Reactions in cobra models require names
		Example name: '3 oxoacyl acyl carrier protein synthase n C140 '
		Since we do not have names for reactions, I use Web-RInChIKeys
		:return:
		"""
		len_before = len(self.df_smiles_chebi_equation_master_id)
		print('Total number of reactions in the rhea reactions set:', len_before)
		self.df_smiles_chebi_equation_master_id.dropna(subset=['Web-RInChIKey'], inplace=True)
		len_after = len(self.df_smiles_chebi_equation_master_id)
		print('Number of reactions after excluding those for which it was impossible to calculate Web-RInChIKey:', len_after)
		print(len_before-len_after, ' reactions were excluded and not added to the model since their structures were not defined')
		self.df_smiles_chebi_equation_master_id.drop_duplicates(subset=['Web-RInChIKey'], inplace=True)
		print('Remained:', len(self.df_smiles_chebi_equation_master_id), 'non-duplicate reactions')
		self.df_smiles_chebi_equation_master_id['webinchi_name'] = self.df_smiles_chebi_equation_master_id['Web-RInChIKey'].apply(
			lambda x: x.replace('Web-RInChIKey=', '').replace('-', '_'))
		
	def add_all_rhea_reactions_to_model(self):
		"""
		Add all Rhea reactions with fully defined structure to the model one by one
		:return:
		"""
		print('Adding reactions to model')
		self.df_smiles_chebi_equation_master_id.progress_apply(self.add_reaction_to_model, axis=1)
		
	def add_reaction_to_model(self, row):
		"""
		Example input reactions:
		10001	[H]O[H].CCCCC(N)=O>>[H][N+]([H])([H])[H].CCCCC(=O)[O-]	CHEBI:15377.CHEBI:16459>>CHEBI:28938.CHEBI:31011
		10009	*N[C@@H](CS)C(*)=O.*N[C@@H](CS)C(*)=O.*OO>>*N[C@@H](CSSC[C@H](N*)C(*)=O)C(*)=O.*O.[H]O[H]	CHEBI:29950.CHEBI:29950.CHEBI:35924>>CHEBI:50058.CHEBI:30879.CHEBI:15377
		10069	CC(=O)N[C@H]1[C@@H](OP(=O)([O-])OP(=O)([O-])OC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(\C)CC/C=C(\C)CCC=C(C)C)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](OP(=O)([O-])OC[C@H](O)COP(=O)([O-])OC[C@@H](O)[C@@H](O)[C@@H](O)COP(=O)([O-])OC[C@@H](O)[C@@H](O)[C@@H](O)CO)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H]2O)C(=O)N1>>CC(=O)N[C@H]1[C@@H](OP(=O)([O-])OP(=O)([O-])OC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(\C)CC/C=C(\C)CCC=C(C)C)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](OP(=O)([O-])OC[C@H](O)COP(=O)([O-])OC[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])OC[C@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H](O)CO)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O.[H+].[H+].[H+].[H+].[H+].O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1.O=C1C=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1	POLYMER:12833.CHEBI:58885.CHEBI:58885.CHEBI:58885.CHEBI:58885.CHEBI:58885>>POLYMER:14260.CHEBI:15378.CHEBI:15378.CHEBI:15378.CHEBI:15378.CHEBI:15378.CHEBI:58223.CHEBI:58223.CHEBI:58223.CHEBI:58223.CHEBI:58223
		
		Example reaction addition from Cobra:
		reaction = Reaction('R_3OAS140')
		reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
		reaction.subsystem = 'Cell Envelope Biosynthesis'
		reaction.lower_bound = 0.  # This is the default
		reaction.upper_bound = 1000.  # This is the default
		--> define metabolites like ACP_c
		reaction.add_metabolites({
			malACP_c: -1.0,
			h_c: -1.0,
			ddcaACP_c: -1.0,
			co2_c: 1.0,
			ACP_c: 1.0,
			omrsACP_c: 1.0
		})
		:param reaction_as_chebi_ids_equation:
		:return:
		"""
		reaction_as_chebi_ids_equation = row['chebi_equation']
		reaction_smiles = row['rxnsmiles']
		reaction_id = f"RHEA_{row['MASTER_ID']}"
		reaction_subsystem = row['reaction_subsystem']
		reaction_name = row['webinchi_name']
		
		reaction = Reaction(reaction_id)
		reaction.name = reaction_name
		reaction.subsystem = reaction_subsystem
		reaction.lower_bound = row['lower_bound'] # permissive bounds since rhea reactions to not have direction
		reaction.upper_bound = row['upper_bound'] # permissive bounds since rhea reactions to not have direction
		
		# Now we have to define metabolites and add stoichiometric coefficient to each metabolite
		reactant_chebis = reaction_as_chebi_ids_equation.split('>>')[0].split('.')
		reactant_smiles = reaction_smiles.split('>>')[0].split('.')
		reactants = self.get_stoich_metabolite_tuple(reactant_chebis, reactant_smiles)
		product_chebis = reaction_as_chebi_ids_equation.split('>>')[1].split('.')
		product_smiles = reaction_smiles.split('>>')[1].split('.')
		products = self.get_stoich_metabolite_tuple(product_chebis, product_smiles)
		products = [(met,-1*coef) for met,coef in products] # change stoich coef to negative for products
		reactants.extend(products)
		reaction.add_metabolites(dict(reactants))
		
		# Add reaction to model
		self.model.add_reactions([reaction])

	def get_stoich_metabolite_tuple(self, side_chebis, side_smiles):
		"""
		translate [CHEBI:29950,CHEBI:29950,CHEBI:35924] to [(Metabolite(CHEBI:29950), 2), (Metabolite(CHEBI:35924), 1)]
		:return:
		"""
		unique_metabolite_ids = set(side_chebis)
		metabolites = []
		for metabolite_id in unique_metabolite_ids:
			stoich_coef = side_chebis.count(metabolite_id)
			metabolite_obj = self.define_metabolite(metabolite_id=metabolite_id, metabolite_SMILES=side_smiles[side_chebis.index(metabolite_id)])
			metabolites.append((metabolite_obj, stoich_coef))
		return metabolites
		
	def define_metabolite(self, metabolite_id, metabolite_SMILES=None, compartment='e'):
		"""
		Example from Cobra:
		ACP_c = Metabolite(
		'ACP_c',
		formula='C11H21N2O7PRS',
		name='acyl-carrier-protein',
		compartment='c')
		"""
		# Lookup metabolite name from CHEBI ID to the name
		metabolite_name_df = self.rdb.df_chebi_cmpname[self.rdb.df_chebi_cmpname['chebiid']==metabolite_id]['cmpname']
		if len(metabolite_name_df)>0:
			metabolite_name = metabolite_name_df.iloc[0]
		else:
			print('No name found for metabolite, report the case:', metabolite_id)
			metabolite_name = 'noname'
		if not metabolite_SMILES:
			metabolite_SMILES = self.rdb.df_chebi_smiles[self.rdb.df_chebi_smiles['chebiid']==metabolite_id]['smiles'].iloc[0]
		calculated_formula = rdMolDescriptors.CalcMolFormula(MolFromSmiles(metabolite_SMILES))
		# Change : to _ since cobra does not take strings with : as ids
		metabolite_id = metabolite_id.replace(':','_')
		return Metabolite(metabolite_id,
						  formula=calculated_formula,
						  name=metabolite_name,
						  compartment=compartment)
	
	def inspect_model(self):
		"""
		Provide basic stats of the model
		:return:
		"""
		print(f'{len(self.model.reactions)} reactions')
		print(f'{len(self.model.metabolites)} metabolites')
		print(f'{len(self.model.variables)} variables')
	

	def add_boundaries(self):
		"""
		TO REVIEW
		Add boundary metabolites to the model.
		From Cobra doc:
		"There are three different types of pre-defined boundary reactions: exchange, demand, and sink reactions.
		All of them are unbalanced pseudo reactions, that means they fulfill a function for modeling by adding
		to or removing metabolites from the model system but are not based on real biology.
		An exchange reaction is a reversible reaction that adds to or removes an extracellular metabolite from the extracellular compartment.
		A demand reaction is an irreversible reaction that consumes an intracellular metabolite.
		A sink is similar to an exchange but specifically for intracellular metabolites, i.e.,
		a reversible reaction that adds or removes an intracellular metabolite.
		...
		Boundary reactions are defined on metabolites."
		:return:
		"""
		self.load_boundary_metabolites()
		for boundary in self.boundary_metabolites:
			try:
				self.model.add_boundary(self.model.metabolites.get_by_id(boundary.replace(':','_')), type="exchange")
			except KeyError:
				self.model.add_metabolites([self.define_metabolite(boundary, compartment='e')])
				self.model.add_boundary(self.model.metabolites.get_by_id(boundary.replace(':','_')), type="exchange")
		print(f'Added {len(self.model.exchanges)} exchanges:')
		
	def load_boundary_metabolites(self, boundary_metabolite_file=None):
		"""
		This function loads metabolites that should be available as source to start modeling all reactions
		:param boundary_metabolite_file: file in which boundary metabolites are defined
		boundary_metabolite_file has 3 columns:
		chebiid	smiles	cmpname
		:return:
		"""
		if boundary_metabolite_file:
			self.df_boundary_metabolites = pd.read_csv(boundary_metabolite_file, sep='\t')
		else:
			# read the default boundaries - set of inorganics
			filepath = os.path.dirname(__file__)
			self.df_boundary_metabolites = pd.read_csv(
				os.path.join(filepath, 'biochemical_assumptions',
							 'default_metabolic_precursor_molecules_with_manual_check.tsv'),
				sep='\t')
			
		self.boundary_metabolites = self.df_boundary_metabolites['chebiid']
		# replace : with _ since optlang in cobra does not take comps ids with :
		print('Total boundary metabolites:', len(self.boundary_metabolites))

