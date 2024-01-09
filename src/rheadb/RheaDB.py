__author__ = "Anastasia Sveshnikova"
__email__ = "anastasia.sveshnikova@sib.swiss"
__version__ = "0.0.1"
__status__ = "Prototype"

import pandas as pd
import bioversions
from pathlib import Path
import os
import urllib
import tarfile

from rdkit import Chem
from rdkit.Chem.rdChemReactions import PreprocessReaction

from .Reaction import Reaction

################## DATA #####################

class RheaDB:
    """
    Class for handling the data stored in the RheaDB
    """
    
    def __init__(self):
        """
        Initiate Rhea DB by setting the current version
        """
        pass
    
    def set_rhea_data_location(self, rhea_versions_folder_location):
        self.RDBv_loc = rhea_versions_folder_location
        
    def set_rhea_version(self):
        """
        Use bioversions to store your rhea version and to not mix it up
        rhea_versions_folder_location : location at which your rhea data is stored
        :return:
        """
        
        version = self.get_current_rhea_version()
        if not os.path.exists(f'{self.RDBv_loc}/rhea-versions'):
            os.mkdir(f'{self.RDBv_loc}/rhea-versions')
        rhea_versions = os.listdir(f'{self.RDBv_loc}/rhea-versions')
        if version not in rhea_versions:
            if not os.path.exists(version):
                os.mkdir(f'{self.RDBv_loc}/rhea-versions/{version}')
                os.mkdir(f'{self.RDBv_loc}/rhea-versions/{version}/tsv')
        print(f'Your Rhea DB version is {version}')
        self.rhea_db_version = version

    def get_current_rhea_version(self):
        """
        Get current version of the Rhea DB
        :return:
        """
        version = bioversions.get_version("rhea")
        return version

    def loadRhea(self):
        """
        Load Rhea DB - download files and generate tsv from .rxn
        :return:
        """
        self.download_rhea_structure()
        self.generateSmilesChebiReactionEquationFile()
        self.df_smiles_chebi_equation = pd.read_csv(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-smiles-chebi.tsv', sep='\t')
        
        self.download_rhea_files()
        self.add_master_id_to_hierarchy()
        self.add_master_id_to_rxnsmiles()
        
    def download_rhea_files(self):
        """
        Download all the relevant Rhea table from the FTP
        :return:
        """
        self.df_hierarchy = self.load_df('rhea-relationships.tsv')
        self.df_directions = self.load_df('rhea-directions.tsv')
        # self.df_smiles = self.load_df('rhea-reaction-smiles.tsv', columnsnames=['rheaid', 'rxnsmiles'])
        # self.df_chebi_smiles = self.load_df('rhea-chebi-smiles.tsv', columnsnames=['chebiid', 'smiles'])
    
    def download_rhea_structure(self):
        """
        Download from Rhea FTP files that store the rhea reactions as .rxn files
        - with ChEBI IDs and .mol structures
        :return:
        """
        if not os.path.exists(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/rxn'):
            urllib.request.urlretrieve('https://ftp.expasy.org/databases/rhea/ctfiles/rhea-rxn.tar.gz',
                                       'rhea-rxn.tar.gz')
            # open file
            file = tarfile.open('rhea-rxn.tar.gz')
            # extracting file
            file.extractall(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}')
            
            if os.path.exists('rhea-rxn.tar.gz'):
                os.remove('rhea-rxn.tar.gz')
        else:
            print('Using previously downloaded Rhea version')
            
    def load_df(self, filename, columnsnames=[]):
        """
        Load the .tsv dataframes from Rhea FTP
        :param filename:
        :param columnsnames:
        :return:
        """
        rhea_tsv_file = f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/{filename}'
        if not os.path.exists(rhea_tsv_file):
            if columnsnames:
                df = pd.read_csv(f'https://ftp.expasy.org/databases/rhea/tsv/{filename}', sep='\t',
                                   names=columnsnames)
            else:
                df = pd.read_csv(f'https://ftp.expasy.org/databases/rhea/tsv/{filename}', sep='\t')
            df.to_csv(rhea_tsv_file, index=False, sep='\t')
            return df
        else:
            return pd.read_csv(rhea_tsv_file, sep='\t')
            
    def add_master_id_to_hierarchy(self):
        """
        Load the data for hierarchy.
        
        Returns:
           pd.DataFrame: The loaded data as a pandas DataFrame.
        """
        if not Path(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-relationships-master-id.tsv').exists():
            self.df_hierarcy_master_id = self.df_hierarchy.copy()
            self.df_hierarcy_master_id[['FROM_REACTION_ID_MASTER_ID', 'DIR_FROM']] =\
                self.df_hierarcy_master_id.apply(self.findMasterID, axis=1,
                                                args=['FROM_REACTION_ID', ], result_type='expand')
            self.df_hierarcy_master_id[['TO_REACTION_ID_MASTER_ID', 'DIR_TO']] = \
                self.df_hierarcy_master_id.apply(self.findMasterID, axis=1,
                                                args=['TO_REACTION_ID', ], result_type='expand')
            self.df_hierarcy_master_id = self.df_hierarcy_master_id[~self.df_hierarcy_master_id['FROM_REACTION_ID_MASTER_ID'].isna()]
            self.df_hierarcy_master_id = self.df_hierarcy_master_id[~self.df_hierarcy_master_id['TO_REACTION_ID_MASTER_ID'].isna()]
            self.df_hierarcy_master_id = self.df_hierarcy_master_id.astype({'FROM_REACTION_ID_MASTER_ID': 'int32', 'TO_REACTION_ID_MASTER_ID': 'int32'})
            self.df_hierarcy_master_id.to_csv(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-relationships-master-id.tsv', sep='\t', index=False)
        else:
            self.df_hierarcy_master_id = pd.read_csv(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-relationships-master-id.tsv', sep='\t')

    def add_master_id_to_rxnsmiles(self):
        """
		Load Rhea reaction smiles
		:return: pandas df of rxn smiles
		"""
        if not Path(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-smiles-master-id.tsv').exists():
            self.df_smiles_master_id = self.df_smiles_chebi_equation.copy()
            self.df_smiles_master_id[['MASTER_ID', 'DIR']] = self.df_smiles_master_id.apply(self.findMasterID, axis=1, args=['rheaid'],
                                                              result_type='expand')
            self.df_smiles_master_id = self.df_smiles_master_id[~self.df_smiles_master_id['MASTER_ID'].isna()]
            self.df_smiles_master_id.drop_duplicates(subset=['MASTER_ID'], inplace=True)
            self.df_smiles_master_id.to_csv(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-smiles-master-id.tsv', sep='\t', index=False)
        else:
            self.df_smiles_master_id = pd.read_csv(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-smiles-master-id.tsv', sep='\t')
    
    def findMasterID(self, row, id):
        """
        Assign the rhea Master ID (undirectional) to avoid LR - RL duplication
        :param row: row of smiles dataframe
        :param df_dir: dataframe of directions in rhea (Mapping LRand RL to undirectional)
        :param id: column name used as rhea id
        :return:
        """
        if self.df_directions.loc[self.df_directions['RHEA_ID_MASTER'] == row[id]]['RHEA_ID_MASTER'].to_list():
         return row[id], 'MASTER'
        if self.df_directions.loc[self.df_directions['RHEA_ID_LR'] == row[id]]['RHEA_ID_MASTER'].to_list():
         return self.df_directions.loc[self.df_directions['RHEA_ID_LR'] == row[id]]['RHEA_ID_MASTER'].to_list()[0], 'LR'
        if self.df_directions.loc[self.df_directions['RHEA_ID_RL'] == row[id]]['RHEA_ID_MASTER'].to_list():
         return self.df_directions.loc[self.df_directions['RHEA_ID_RL'] == row[id]]['RHEA_ID_MASTER'].to_list()[0], 'RL'
        if self.df_directions.loc[self.df_directions['RHEA_ID_BI'] == row[id]]['RHEA_ID_MASTER'].to_list():
         return self.df_directions.loc[self.df_directions['RHEA_ID_BI'] == row[id]]['RHEA_ID_MASTER'].to_list()[0], 'BI'
        # None is for obsolete reactions
        return None
    
    def generateSmilesChebiReactionEquationFile(self):
        """
        Download the reactions as .rxn folder from Rhea ftp
        Transform them into reaction smiles and corresponding ChEBI ids in the matching order
        Returns:
        Generates the .tsv with smiles and ChEBI equation in .tsv folder
        """
        rgen = Reaction() # generic reaction class object to access the processing method
        
        file_rhea_db = f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-smiles-chebi.tsv'
        if not os.path.exists(file_rhea_db):
            rxns = os.listdir(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/rxn')
            with open(file_rhea_db, 'w') as w:
                w.write('rheaid\trxnsmiles\tchebi_equation\n')
                for rxnid in rxns:
                
                    with open(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/rxn/{rxnid}') as f:
                        chebis = []
                        for line in f:
                            if line.startswith('CHEBI') or line.startswith('POLYMER'):
                                chebis.append(line.strip())
                
                    rxn = rgen.read_rxnfile(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/rxn/{rxnid}')
                    nWarn, nError, nReacts, nProds, reactantLabels = PreprocessReaction(rxn)
                
                    if not nReacts + nProds == len(chebis):
                        print(f'{rxn} did not match the number of chebi ids and number of compounds:')
                        print(nReacts, nProds, chebis)
                
                    reactant_smiles = []
                    reactant_chebis = []
                    product_smiles = []
                    product_chebis = []
                
                    for i in range(nReacts):
                        react = rxn.GetReactantTemplate(i)
                        reactant_smiles.append(Chem.MolToSmiles(react))
                        reactant_chebis.append(chebis[i])
                
                    for i in range(nProds):
                        p = rxn.GetProductTemplate(i)
                        product_smiles.append(Chem.MolToSmiles(p))
                        product_chebis.append(chebis[nReacts + i])
                
                    smiles_equation = '.'.join(reactant_smiles) + '>>' + '.'.join(product_smiles)
                    chebi_equation = '.'.join(reactant_chebis) + '>>' + '.'.join(product_chebis)
                
                    w.write(f'{rxnid.split(".")[0]}\t{smiles_equation}\t{chebi_equation}\n')
                    