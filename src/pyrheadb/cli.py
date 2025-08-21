#!/usr/bin/env python3

import argparse
from .ReactionSmartsConverter import ReactionSmartsConverter
from .ReactionPrediction import ReactionPrediction
from .ReactionNetwork import ReactionNetwork
from .RheaDB import RheaDB
import pandas as pd

rhea_db = RheaDB()

def setup_smarts_converter(args):
    # Uses the default Rhea database path
    converter = ReactionSmartsConverter(rhea_db)
    converter.convert_all_rhea_smiles_to_smarts()
    print('Rhea templates converted to SMARTS')

def setup_reaction_prediction(args):
    # Uses the default Rhea database path
    predictor = ReactionPrediction(rhea_db)
    substrates = args.substrates.split(',')
    df_predicted = predictor.predict_reactions(substrates)
    if type(df_predicted)==pd.DataFrame and len(df_predicted)>0:
        df_predicted.to_csv(f'predicted_reactions_{args.substrates}.tsv', sep='\t', index=False)

def setup_disconnected_compound_identification(args):
    #
    rxnnet = ReactionNetwork()
    rxnnet.build_network_from_rhea_tsv(rhea_db, type=args.network_type)
    comps = rxnnet.get_disconnected_compounds()
    print('Total disconnected compounds identified:', len(comps))
    return comps

def main():
    parser = argparse.ArgumentParser(description="CLI tool for RheaDB reaction utilities")
    parser.add_argument('--smarts', action='store_true', help='Convert SMILES to SMARTS')
    parser.add_argument('--predict', action='store_true', help='Predict reaction products from substrates')
    parser.add_argument('--dxd_compounds', action='store_true', help='Identify disconnected compound IDs')
    parser.add_argument('--substrates', type=str, help='Comma-separated list of substrate SMILES for product prediction')
    parser.add_argument('--network_type', type=str, help='chebiid OR inchikey OR inchikey14L : what will be used as a node in the network analysis')

    args = parser.parse_args()

    if args.smarts:
        setup_smarts_converter(args)
    elif args.predict:
        setup_reaction_prediction(args)
    elif args.dxd_compounds:
        setup_disconnected_compound_identification(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
