__author__ = "Anastasia Sveshnikova"
__email__ = "anastasia.sveshnikova@sib.swiss"
__status__ = "Prototype"

import pandas as pd
import bioversions
from pathlib import Path
import os
import urllib
import tarfile
import gzip
import shutil

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
        
    def set_rhea_version(self, rheaversion = False):
        """
        rheaversion : optional, the version of Rhea db to use
        Use bioversions to store your rhea version and to not mix it up
        rhea_versions_folder_location : location at which your rhea data is stored
        :return:
        """
        if not rheaversion:
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
        self.download_rhea_compound_sdf()
        self.parseSDF()
        self.chebiId_name = pd.read_csv(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/compound_names.tsv', sep='\t')
        self.generateReactionCompoundNamesFile()
        
    def download_rhea_files(self):
        """
        Download all the relevant Rhea table from the FTP
        :return:
        """
        self.df_hierarchy = self.load_df('rhea-relationships.tsv')
        self.df_directions = self.load_df('rhea-directions.tsv')
        self.df_smiles = self.load_df('rhea-reaction-smiles.tsv', columnsnames=['rheaid', 'rxnsmiles'])
        self.df_chebi_smiles = self.load_df('rhea-chebi-smiles.tsv', columnsnames=['chebiid', 'smiles'])

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
            print('Using previously downloaded Rhea .rxn file version')

    def download_rhea_compound_sdf(self):
        """
        Download from Rhea FTP files that store the rhea reactions as .rxn files
        - with ChEBI IDs and .mol structures
        :return:
        """
        if not os.path.exists(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/sdf'):
            os.mkdir(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/sdf')
            urllib.request.urlretrieve('https://ftp.expasy.org/databases/rhea/ctfiles/rhea.sdf.gz',
                                       'rhea.sdf.gz')
            # open file
            with gzip.open('rhea.sdf.gz', 'rb') as file_in:
                with open(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/sdf/rhea.sdf', 'wb') as file_out:
                    shutil.copyfileobj(file_in, file_out)
        
            if os.path.exists('rhea.sdf.gz'):
                os.remove('rhea.sdf.gz')
        else:
            print('Using previously downloaded Rhea compounds .sdf version')
            
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
    
    def setUndefinedCompoundStarFlag(self):
        """
        Mark the reactions with compounds that do not have fully defined structures in the dataframe
        :return:
        """
        self.df_smiles_master_id['star'] = self.df_smiles_master_id['rxnsmiles'].apply(lambda x: "*" in x)

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
                    
    def generateReactionCompoundNamesFile(self):
        """
        :return: create self.df_reaction_participants_names - table that stores the reactions as compounds names
        in format of:
        reactant_name1 + reactant_name2 = reactant_name3
        """
        filename_df_reaction_participants_names = f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-participant-names.tsv'
        if not os.path.exists(filename_df_reaction_participants_names):
            df_temp = self.df_smiles_master_id[['rheaid','chebi_equation', 'MASTER_ID']].copy()
            chebi_dict = dict(zip(self.chebiId_name['chebiid'].to_list(), self.chebiId_name['cmpname'].to_list()))
            df_temp['reaction_participant_names'] = df_temp.apply(self.get_reaction_in_names, axis=1, args = [chebi_dict,])
            df_temp.to_csv(filename_df_reaction_participants_names, sep='\t', index=False, columns=['rheaid', 'MASTER_ID', 'reaction_participant_names'])
        self.df_reaction_participants_names = pd.read_csv(filename_df_reaction_participants_names, sep='\t')
    
    def get_reaction_in_names(self, row, chebi_dict):
        """
        subfunction of generateReactionCompoundNamesFile(self)
        :param row: row of self.df_smiles_chebi_equation copy dataframe
        :return:
        """
        chebi_equation = row['chebi_equation']
        chebi_reactants = chebi_equation.split('>>')[0].split('.')
        chebi_products =  chebi_equation.split('>>')[1].split('.')
        
        # Test I was running to check that no reactant names were lost
        # extra_reactants = set(chebi_reactants).union(set(chebi_products))-set(chebi_dict.keys())
        # if len(extra_reactants)>0:
        #     print('Extra', extra_reactants)
        #     return None
        
        reactant_names = [chebi_dict[chebiid] for chebiid in chebi_reactants]
        product_names = [chebi_dict[chebiid] for chebiid in chebi_products]
        return ' + '.join(reactant_names) + ' = ' + ' + '.join(product_names)
        
    def loadLongTableReactionParticipats(self):
        """
        This function transforms the wide format of storing reaction-compounds relationship
        (wide format example: RHEA:10008    CHEBI:29950.CHEBI:29950.CHEBI:35924>>CHEBI:50058.CHEBI:30879.CHEBI:15377)
        into long format of storing reaction-compounds relationship:
        RHEA:10008    CHEBI:29950
        RHEA:10008    CHEBI:35924
        RHEA:10008    CHEBI:50058
        RHEA:10008    CHEBI:30879
        RHEA:10008    CHEBI:15377
        
        Uses self.df_smiles_master_id table
        :return: generates self.rhea_reaction_long_format_smiles_chebi pabdas dataframe
        and self.rhea_reaction_long_format_smiles_chebi.tsv file
        """
        filename_rhea_reaction_smiles_chebi = f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-long-format-smiles-chebi.tsv'
        if not os.path.exists(filename_rhea_reaction_smiles_chebi):
            with open(filename_rhea_reaction_smiles_chebi, 'w') as w:
                w.write('MASTER_ID\treaction_side\tchebiid\tsmiles\n')
                for index, row in self.df_smiles_master_id.iterrows():
                    chebi_equation = row['chebi_equation']
                    rxnsmiles = row['rxnsmiles']
                    rheaid = row['MASTER_ID']
                    
                    reactant_chebis = chebi_equation.split('>>')[0].split('.')
                    reactant_smiles = rxnsmiles.split('>>')[0].split('.')
                    assert len(reactant_chebis) == len(reactant_smiles)
                    for i in range(len(reactant_chebis)):
                        w.write(f'{rheaid}\t{rheaid}_L\t{reactant_chebis[i]}\t{reactant_smiles[i]}\n')

                    product_chebis = chebi_equation.split('>>')[1].split('.')
                    product_smiles = rxnsmiles.split('>>')[1].split('.')
                    assert len(product_chebis) == len(product_smiles)
                    for i in range(len(product_chebis)):
                        w.write(f'{rheaid}\t{rheaid}_R\t{product_chebis[i]}\t{product_smiles[i]}\n')
                        
        self.rhea_reaction_long_format_smiles_chebi = pd.read_csv(
            f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/rhea-reaction-long-format-smiles-chebi.tsv', sep='\t')
        self.rhea_reaction_long_format_smiles_chebi.drop_duplicates(inplace=True)

    def parseSDF(self):
        """
        Parse the downloaded .sdf file that contains the structures and names for all compounds
        Extract only names per id as this ids are used in the .rxn files and not all chebiids are
        present in .tsv on the website
        :return:
        """
        compound_id_name_dict = dict()
        new_compound = True
        name_attention = False
        with open(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/sdf/rhea.sdf') as f:
            for line in f:
                if new_compound == True:
                    compound_id = line.strip()
                    new_compound = False
                if line.startswith('$$$$'):
                    new_compound = True
                if name_attention == True:
                    compound_id_name_dict[compound_id] = line.strip()
                    name_attention = False
                if line.startswith('> <Rhea_ascii_name>'):
                    name_attention = True
        with open(f'{self.RDBv_loc}/rhea-versions/{self.rhea_db_version}/tsv/compound_names.tsv', 'w') as w:
            w.write('chebiid\tcmpname\n')
            for compound_id, compound_name in compound_id_name_dict.items():
                w.write(f'{compound_id}\t{compound_name}\n')

