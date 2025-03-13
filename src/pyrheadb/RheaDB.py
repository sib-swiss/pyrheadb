__author__ = "Anastasia Sveshnikova"
__email__ = "anastasia.sveshnikova@sib.swiss"
__status__ = "Prototype"

import pandas as pd
from pathlib import Path
import os
import urllib.request as request
import tarfile
import gzip
import shutil

from rdkit import Chem
from rdkit.Chem.rdChemReactions import PreprocessReaction

from .Reaction import Reaction
from .RInChI import RInChI

################## DATA #####################

class RheaDB:
    """
    Class for handling the data stored in the RheaDB
    """
    
    def __init__(self, rhea_versions_folder_location=None, rhea_version=None):
        """
        Initiate Rhea DB by setting the current version
        """
        if not rhea_versions_folder_location:
            self.rhea_versions_folder_location = os.getenv('RHEADB_LOC')
        else:
            self.rhea_versions_folder_location = rhea_versions_folder_location
        self.__set_rhea_version(rhea_version=rhea_version)
        self.__load_rhea()
        
    def __set_rhea_version(self, rhea_version=None):
        """
        rheaversion : optional, the version of Rhea db to use
        Use bioversions to store your rhea version and to not mix it up
        rhea_versions_folder_location : location at which your rhea data is stored
        :return:
        """
        if not rhea_version:
            version = self.__handle_ulr_error()
        else:
            version = rhea_version
        
        if not version:
            print('Impossible to access any Rhea version.')
            print('Connect to internet and try again.')
            exit()
        os.makedirs(os.path.join(f'{self.rhea_versions_folder_location}','rhea-versions'), exist_ok=True)

        print(f'Your Rhea DB version is {version}')
        self.rhea_db_version = version
        self.rhea_db_version_location = os.path.join(f'{self.rhea_versions_folder_location}','rhea-versions', f'{version}')

    def __handle_ulr_error(self):
        try:
            return self.__get_current_rhea_version()
        except Exception as e:
            print(e)
            print('It appears you are not connected to internet and we cannot download the database.')
            print('We will check what is the latest version on your computer.')
            versions = os.listdir(os.path.join(f'{self.rhea_versions_folder_location}','rhea-versions'))
            versions = [int(i) for i in versions if all([c.isdigit() for c in i])]
            if versions:
                return (max(versions))
            
    def __get_current_rhea_version(self):
        """
        Get current version of the Rhea DB
        :return:
        """
        request.urlretrieve('https://ftp.expasy.org/databases/rhea/rhea-release.properties',
                                   'rhea-release.properties')
        with open('rhea-release.properties') as f:
            version = int(f.readline().split('=')[1].strip())
        Path('rhea-release.properties').unlink(missing_ok=True)
        return version
    
    def __load_rhea(self):
        """
        Load Rhea DB - download files and generate tsv from .rxn
        :return:
        """
        self.__download_rhea_version_full()
        self.__extract_rhea_structure()
        self.__extract_rhea_compound_sdf()
        self.__read_rhea_files()
        self.__add_master_id_to_hierarchy()
        
        pyrheadb_reactions_file = os.path.join(f'{self.rhea_db_version_location}','tsv','pyrheadb_reactions.tsv')
        if not os.path.exists(pyrheadb_reactions_file):
            self.__generate_smiles_chebi_reaction_equation_file()
            self.df_smiles_chebi_equation = self.__read_tsv_to_pandas('rhea-reaction-smiles-chebi.tsv')
            self.__add_master_id_to_rxnsmiles()
            self.__generate_reaction_compound_names_file()
            self.__merge_all_intermediate_reaction_files_into_one()
            self.__separate_protein_residue_from_unstructured_compound_reactions()
            self.__set_class_reaction_flag()
            self.__set_polymer_reaction_flag()
            self.__add_rinchikey()
            self.df_reactions.to_csv(pyrheadb_reactions_file, sep='\t', index=False)
        else:
            self.df_reactions = pd.read_csv(pyrheadb_reactions_file, sep='\t')

        self.__load_long_table_reaction_participants()
    
    def __read_tsv_to_pandas(self, tsv_filename, columnsnames=[]):
        """
        Load the .tsv dataframes from Rhea FTP
        :param filename:
        :param columnsnames:
        :return:
        """
        rhea_tsv_file = os.path.join(f'{self.rhea_db_version_location}','tsv',f'{tsv_filename}')
        if columnsnames:
            df = pd.read_csv(rhea_tsv_file, sep='\t',
                               names=columnsnames)
        else:
            df = pd.read_csv(rhea_tsv_file, sep='\t')
        return df
    
    def __download_rhea_version_full(self):
        """
        Download all Rhea version from expasy FTP
        :return:
        """
        if not os.path.exists(f'{self.rhea_db_version_location}'):
            request.urlretrieve(f'https://ftp.expasy.org/databases/rhea/old_releases/{self.rhea_db_version}.tar.bz2',
                                       f'{self.rhea_db_version}.tar.bz2')
            self.__extract_rhea_tar()
            Path(f'{self.rhea_db_version}.tar.bz2').unlink(missing_ok=True)
        else:
            print('Using previously downloaded Rhea version')
    
    def __extract_rhea_tar(self):
        # open file
        file = tarfile.open(f'{self.rhea_db_version}.tar.bz2', 'r:bz2')
        # extracting file
        file.extractall(os.path.join(f'{self.rhea_versions_folder_location}', 'rhea-versions'))
        
    def __read_rhea_files(self):
        """
        Load all the relevant Rhea tables from the downloaded FTP folder
        :return:
        """
        self.df_hierarchy = self.__read_tsv_to_pandas('rhea-relationships.tsv')
        self.df_directions = self.__read_tsv_to_pandas('rhea-directions.tsv')
        self.df_smiles = self.__read_tsv_to_pandas('rhea-reaction-smiles.tsv', columnsnames=['rheaid', 'rxnsmiles'])

    def __extract_rhea_structure(self):
        """
        Download from Rhea FTP files that store the rhea reactions as .rxn files
        - with ChEBI IDs and .mol structures
        :return:
        """
        if not os.path.exists(os.path.join(f'{self.rhea_db_version_location}','rxn')):
            # open file
            file = tarfile.open(os.path.join(f'{self.rhea_db_version_location}','ctfiles','rhea-rxn.tar.gz'))
            # extracting file
            file.extractall(os.path.join(f'{self.rhea_db_version_location}'))

    def __extract_rhea_compound_sdf(self):
        """
        Download from Rhea FTP files that store the rhea reactions as .rxn files
        - with ChEBI IDs and .mol structures
        :return:
        """
        if not os.path.exists(os.path.join(f'{self.rhea_db_version_location}','sdf')):
            os.mkdir(os.path.join(f'{self.rhea_db_version_location}','sdf'))
            # open file
            with gzip.open(os.path.join(f'{self.rhea_db_version_location}','ctfiles','rhea.sdf.gz'), 'rb') as file_in:
                with open(os.path.join(f'{self.rhea_db_version_location}','sdf','rhea.sdf'), 'wb') as file_out:
                    shutil.copyfileobj(file_in, file_out)
        if not os.path.exists(os.path.join(f'{self.rhea_db_version_location}','tsv','chebi_in_rhea.tsv')):
            self.__parse_sdf()
        
        self.df_chebi_cmpname = self.__read_tsv_to_pandas('chebi_in_rhea.tsv')

        
    def __add_master_id_to_hierarchy(self):
        """
        Load the data for hierarchy.
        
        Returns:
           pd.DataFrame: The loaded data as a pandas DataFrame.
        """
        rhea_hierarchy_file = os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-relationships-reaction-hierarchy.tsv')
        if not Path(rhea_hierarchy_file).exists():
            self.df_hierarchy_master_id = self.df_hierarchy.copy()
            self.df_hierarchy_master_id[['FROM_REACTION_ID_MASTER_ID', 'DIR_FROM']] =\
                self.df_hierarchy_master_id.apply(self.__find_master_id, axis=1,
                                                args=['FROM_REACTION_ID', ], result_type='expand')
            self.df_hierarchy_master_id[['TO_REACTION_ID_MASTER_ID', 'DIR_TO']] = \
                self.df_hierarchy_master_id.apply(self.__find_master_id, axis=1,
                                                args=['TO_REACTION_ID', ], result_type='expand')
            self.df_hierarchy_master_id = self.df_hierarchy_master_id[~self.df_hierarchy_master_id['FROM_REACTION_ID_MASTER_ID'].isna()]
            self.df_hierarchy_master_id = self.df_hierarchy_master_id[~self.df_hierarchy_master_id['TO_REACTION_ID_MASTER_ID'].isna()]
            self.df_hierarchy_master_id = self.df_hierarchy_master_id.astype({'FROM_REACTION_ID_MASTER_ID': 'int32', 'TO_REACTION_ID_MASTER_ID': 'int32'})
            self.df_hierarchy_master_id.to_csv(rhea_hierarchy_file, sep='\t', index=False)
        else:
            self.df_hierarchy_master_id = pd.read_csv(rhea_hierarchy_file, sep='\t')

    def __add_master_id_to_rxnsmiles(self):
        """
		Load Rhea reaction smiles
		:return: pandas df of rxn smiles
		"""
        rhea_reaction_smiles_master_id_file = os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-smiles-master-id.tsv')
        if not Path(rhea_reaction_smiles_master_id_file).exists():
            self.df_smiles_master_id = self.df_smiles_chebi_equation.copy()
            self.df_smiles_master_id[['MASTER_ID', 'DIR']] = self.df_smiles_master_id.apply(self.__find_master_id, axis=1, args=['rheaid'],
                                                              result_type='expand')
            self.df_smiles_master_id = self.df_smiles_master_id[~self.df_smiles_master_id['MASTER_ID'].isna()]
            self.df_smiles_master_id.drop_duplicates(subset=['MASTER_ID'], inplace=True)
            self.df_smiles_master_id.to_csv(rhea_reaction_smiles_master_id_file, sep='\t', index=False)
        else:
            self.df_smiles_master_id = pd.read_csv(rhea_reaction_smiles_master_id_file, sep='\t')
            
    def __find_master_id(self, row, id):
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
    
    def __set_class_reaction_flag(self):
        """
        Mark the reactions with compounds that do not have fully defined structures in the dataframe
        :return:
        """
        reactionobj = Reaction()
        self.df_reactions['rxn_smiles_no_A_AH']= self.df_reactions['rxnsmiles'].apply(reactionobj.remove_A_AH_pattern)
        self.df_reactions['class_reaction_flag'] = self.df_reactions['rxn_smiles_no_A_AH'].apply(lambda x: "*" in x)

    def __set_polymer_reaction_flag(self):
        """
        Mark the reactions with compounds that do not have fully defined structures in the dataframe
        :return:
        """
        self.df_reactions['polymer_reaction_flag'] = self.df_reactions['chebi_equation'].apply(lambda x: "POLYMER" in x)
    #
    # def rhea_atommap_reaction_reader(self):
    #     atommapped_file = 'data/rheadf_atom_mapped.tsv'
    #     if not os.path.exists(atommapped_file):
    #         self.df_smiles_master_id['star'] = self.df_smiles_master_id['rxnsmiles'].apply(lambda x: "*" in x)
    #         print('Reactions with * :', len(self.df_smiles_master_id[self.df_smiles_master_id['star'] == True]))
    #         print('Reactions with compounds defined :',
    #               len(self.df_smiles_master_id[self.df_smiles_master_id['star'] == False]))
    #         self.df_one_dir_only_defined_rhea_smiles = self.df_smiles_master_id[
    #             self.df_smiles_master_id['star'] == False].copy()
    #         self.df_one_dir_only_defined_rhea_smiles.drop(columns=['star', 'DIR'], inplace=True)
    #         self.df_smiles_master_id.drop(columns=['star'], inplace=True)
    #     else:
    #         self.df_one_dir_only_defined_rhea_smiles = pd.read_csv(atommapped_file, sep='\t')

        
    def __generate_smiles_chebi_reaction_equation_file(self):
        """
        Download the reactions as .rxn folder from Rhea ftp
        Transform them into reaction smiles and corresponding ChEBI ids in the matching order
        Returns:
        Generates the .tsv with smiles and ChEBI equation in .tsv folder
        """
       
        rhea_reaction_smiles_chebi_file = os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-smiles-chebi.tsv')
        if not os.path.exists(rhea_reaction_smiles_chebi_file):
            rgen = Reaction()  # generic reaction class object to access the processing method
    
            rxns = os.listdir(os.path.join(f'{self.rhea_db_version_location}','rxn'))
            with open(rhea_reaction_smiles_chebi_file, 'w') as w:
                w.write('rheaid\trxnsmiles\tchebi_equation\n')
                for rxnid in rxns:
                
                    with open(os.path.join(f'{self.rhea_db_version_location}','rxn',f'{rxnid}')) as f:
                        chebis = []
                        for line in f:
                            if line.startswith('CHEBI') or line.startswith('POLYMER'):
                                chebis.append(line.strip())
                
                    rxn = rgen.read_rxnfile(os.path.join(f'{self.rhea_db_version_location}','rxn',f'{rxnid}'))
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
                    
    def __generate_reaction_compound_names_file(self):
        """
        :return: create self.df_reaction_participants_names - table that stores the reactions as compounds names
        in format of:
        reactant_name1 + reactant_name2 = reactant_name3
        """
        filename_df_reaction_participants_names = os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-participant-names.tsv')
        if not os.path.exists(filename_df_reaction_participants_names):
            df_temp = self.df_smiles_master_id[['rheaid','chebi_equation', 'MASTER_ID']].copy()
            chebi_dict = dict(zip(self.df_chebi_cmpname['chebiid'].to_list(), self.df_chebi_cmpname['cmpname'].to_list()))
            df_temp['reaction_participant_names'] = df_temp.apply(self.__get_reaction_in_names, axis=1, args = [chebi_dict,])
            df_temp.to_csv(filename_df_reaction_participants_names, sep='\t', index=False, columns=['rheaid', 'MASTER_ID', 'reaction_participant_names'])
        self.df_reaction_participants_names = pd.read_csv(filename_df_reaction_participants_names, sep='\t')
    
    def __get_reaction_in_names(self, row, chebi_dict):
        """
        subfunction of generateReactionCompoundNamesFile(self)
        :param row: row of self.df_smiles_chebi_equation copy dataframe
        :return:
        """
        chebi_equation = row['chebi_equation']
        chebi_reactants = chebi_equation.split('>>')[0].split('.')
        chebi_products =  chebi_equation.split('>>')[1].split('.')
        
        # Test I was running to check that no reactant names were lost
        extra_reactants = set(chebi_reactants).union(set(chebi_products))-set(chebi_dict.keys())
        if len(extra_reactants)>0:
            for reactant in extra_reactants:
                chebi_dict[reactant]='noname'
        
        reactant_names = [chebi_dict[chebi] for chebi in chebi_reactants]
        product_names = [chebi_dict[chebi] for chebi in chebi_products]
        return ' + '.join(reactant_names) + ' = ' + ' + '.join(product_names)
    
    def __merge_all_intermediate_reaction_files_into_one(self):
        """
        This function takes all the intermediate reaction files that were generated and merges them into one "df_reactions" dataframe.
        This function also deletes all intermediate files, and only leaves the non-redundant dataframe files
        """
        
        self.df_reaction_participants_names.drop(columns=['rheaid'],inplace=True)
        self.df_smiles_master_id.drop(columns=['DIR'],inplace=True)
        self.df_reactions = self.df_reaction_participants_names.merge(self.df_smiles_master_id, on='MASTER_ID')
        os.remove(os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-smiles-master-id.tsv'))
        os.remove(os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-participant-names.tsv'))
        os.remove(os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-smiles-chebi.tsv'))
    
    def __separate_protein_residue_from_unstructured_compound_reactions(self):
        """
        Some Rhea reactions have * because they have protein, DNA etc molecule attached.
        These reactions cannot be used as templates and they have to be substituted by dummy atom (here At)
        """
        self.df_reactions.rename(columns={'rxnsmiles':'rxnsmiles_no_residue_correction'}, inplace=True)
        self.df_reactions[['residue_rxn_flag', 'rxnsmiles']] = self.df_reactions.apply(self.__identify_and_process_residue_row, axis=1, result_type='expand')
    
    def __identify_and_process_residue_row(self, row):
        if not 'residue' in str(row['reaction_participant_names']):
            return False, row['rxnsmiles_no_residue_correction']
        
        reactant_names = row['reaction_participant_names'].split(' = ')[0].split(' + ')
        product_names = row['reaction_participant_names'].split(' = ')[1].split(' + ')

        reactant_smiles = row['rxnsmiles_no_residue_correction'].split('>>')[0].split('.')
        product_smiles = row['rxnsmiles_no_residue_correction'].split('>>')[1].split('.')

        new_smiles_reactants = self.__process_side_identify_residues(reactant_names, reactant_smiles)
        new_smiles_products = self.__process_side_identify_residues(product_names, product_smiles)

        return True, '.'.join(new_smiles_reactants)+'>>'+'.'.join(new_smiles_products)

    def __process_side_identify_residues(self, names, smiles):
        assert len(names) == len(smiles)
        new_smiles_reactants = []
        for i in range(len(names)):
            if 'residue' in names[i]:
                new_smiles_reactants.append(smiles[i].replace('*','At'))
            else:
                new_smiles_reactants.append(smiles[i])
        return new_smiles_reactants

    def __add_rinchikey(self):
        """
        add_rinchikey_column
        :param df_rxnsmiles: pandas dataframe with rxnsmiles columns
        :return: dataframe with RInChI column added
        """
        rinchiobj = RInChI()
        print('Calculating Reaction InChiKeys')
        self.df_reactions[['RInChI','Web-RInChIKey']] = \
            self.df_reactions.apply(lambda row: rinchiobj.error_handle_wrap_rinchi(row['rxnsmiles']), axis=1, result_type='expand')
          
    def __load_long_table_reaction_participants(self):
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
        and rhea_reaction_long_format_smiles_chebi.tsv file
        """
        filename_rhea_reaction_smiles_chebi = os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-long-format-smiles-chebi.tsv')
        if not os.path.exists(filename_rhea_reaction_smiles_chebi):
            with open(filename_rhea_reaction_smiles_chebi, 'w') as w:
                w.write('MASTER_ID\treaction_side\tchebiid\tsmiles\tinchi\tinchikey\tinchikey14L\n')
                for index, row in self.df_smiles_master_id.iterrows():
                    chebi_equation = row['chebi_equation']
                    rxnsmiles = row['rxnsmiles']
                    rheaid = row['MASTER_ID']
                    
                    reactant_chebis = chebi_equation.split('>>')[0].split('.')
                    reactant_smiles = rxnsmiles.split('>>')[0].split('.')
                    assert len(reactant_chebis) == len(reactant_smiles)
                    for i in range(len(reactant_chebis)):
                        reactant_inchi, reactant_inchikey = self.__get_inchi_inchikey(reactant_smiles[i])
                        w.write(f"{rheaid}\t{rheaid}_L\t{reactant_chebis[i]}\t{reactant_smiles[i]}\t{reactant_inchi}\t{reactant_inchikey}\t{reactant_inchikey.split('-')[0]}\n")

                    product_chebis = chebi_equation.split('>>')[1].split('.')
                    product_smiles = rxnsmiles.split('>>')[1].split('.')
                    assert len(product_chebis) == len(product_smiles)
                    for i in range(len(product_chebis)):
                        product_inchi, product_inchikey = self.__get_inchi_inchikey(product_smiles[i])
                        w.write(f"{rheaid}\t{rheaid}_R\t{product_chebis[i]}\t{product_smiles[i]}\t{product_inchi}\t{product_inchikey}\t{product_inchikey.split('-')[0]}\n")
                        
        self.rhea_reaction_long_format_smiles_chebi = pd.read_csv(
            os.path.join(f'{self.rhea_db_version_location}','tsv','rhea-reaction-long-format-smiles-chebi.tsv'), sep='\t')
        self.rhea_reaction_long_format_smiles_chebi.drop_duplicates(inplace=True)

    def __get_inchi_inchikey(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)
        if not inchi:
            return 'NA', 'NA-NA-NA'
        return inchi, inchikey

    def __parse_sdf(self):
        """
        Parse the downloaded .sdf file that contains the structures and names for all compounds
        Extract only names per id as this ids are used in the .rxn files and not all chebiids are
        present in .tsv on the website
        :return:
        """
        compound_id_name_dict = dict()
        new_compound = True
        name_attention = False
        compound_sdf = ''
        with open(os.path.join(f'{self.rhea_db_version_location}','sdf','rhea.sdf')) as f, \
            open(os.path.join(f'{self.rhea_db_version_location}','tsv','chebi_in_rhea.tsv'), 'w') as w:
            w.write('chebiid\tcmpname\tsmiles\tinchi\tinchikey\n')
            for line in f:
                compound_sdf += line
                if new_compound == True:
                    compound_id = line.strip()
                    new_compound = False
                if line.startswith('$$$$'):
                    m = Chem.MolFromMolBlock(compound_sdf)
                    if not m:
                        inchi = None
                        inchikey=None
                        smiles=None
                    else:
                        inchi = Chem.MolToInchi(m)
                        inchikey = Chem.MolToInchiKey(m)
                        smiles = Chem.MolToSmiles(m)
                    w.write(f'{compound_id}\t{compound_name}\t{smiles}\t{inchi}\t{inchikey}\n')
                    compound_sdf = ''
                    new_compound = True
                if name_attention == True:
                    compound_name = line.strip()
                    name_attention = False
                if line.startswith('> <Rhea_ascii_name>'):
                    name_attention = True
                

    def data_overview(self):
        """
        Function for summary of the RheaDB:
        prints all the dataframes that can be used for other applications
        :return:
        """
        print('1.')
        print("self.df_hierarchy_master_id")
        print(self.df_hierarchy_master_id.columns)
        print(self.df_hierarchy_master_id.shape)
        print()
        print('2.')
        print("self.df_reactions")
        print(self.df_reactions.columns)
        print(self.df_reactions.shape)
        print()
        print('3.')
        print("self.rhea_reaction_long_format_smiles_chebi")
        print(self.rhea_reaction_long_format_smiles_chebi.columns)
        print(self.rhea_reaction_long_format_smiles_chebi.shape)
        print()
        print(4.)
        print("self.df_chebi_in_rhea")
        print(self.df_chebi_cmpname.columns)
        print(self.df_chebi_cmpname.shape)
