# From SWAT4HCLS course
# https://github.com/sib-swiss/sparql-training/blob/master/rhea/SWAT4HCLS_2019/rhea_tutorial_SWAT4HCLS_2019.ipynb
from SPARQLWrapper import SPARQLWrapper, JSON
import pandas as pd


#SPARQL endpoints
sparql_uniprot_url = "https://sparql.uniprot.org/sparql/"
sparql_rhea_url = "https://sparql.rhea-db.org/sparql"

def sparql2pandas(sparql_query, sparql_service_url):
    """
    Use SPARQLWrapper to query a SPARQL endpoint (sparql_service_url) with a given SPARQL query string
    (sparql_query) that uses SELECT query form.
    Return the SPARQL query result as a pandas Dataframe.
    """
    #
    # run the SPARQL query
    #
    sparql = SPARQLWrapper(sparql_service_url)
    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)
    res = sparql.query().convert()
    
    # convert the JSON result (res) in a pandas dataframe
    res_sparql_df = pd.json_normalize(res["results"]["bindings"])
    
    # distinguish .type and .value
    col_type = [c for c in res_sparql_df.columns.tolist() if ".type" in c]
    col_value = [c for c in res_sparql_df.columns.tolist() if ".value" in c]
    col_datatype = [c for c in res_sparql_df.columns.tolist() if ".datatype" in c]
    
    # Remove .type columns
    res_sparql_df.drop(col_type, axis=1, inplace=True)
    # Remove .datatype columns
    res_sparql_df.drop(col_datatype, axis=1, inplace=True)
    # Remove ".value" from column names
    res_sparql_df = res_sparql_df.rename(columns=lambda col: col.replace(".value", ""))
    
    return res_sparql_df

def get_all_human_reactions():
    """
    Use uniprot taxid to get all associated rhea reactions
    :return:
    """
    # Homo sapiens (human), species
    taxid = '9606'
        
    # Display the 5 first rows
    return get_all_taxid_reactions(taxid)

def get_all_taxid_reactions(taxid):
    """
    Use uniprot taxid to get all associated rhea reactions
    :return:
    """

    sparql_Q6 = """
    #endpoint:https://sparql.rhea-db.org/sparql

    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX ch:<http://purl.obolibrary.org/obo/>

    SELECT DISTINCT ?protein
                    ?proteinId
                    ?proteinName
                    ?reaction
    WHERE {
      ?reaction rdfs:subClassOf rh:Reaction .

      SERVICE <https://sparql.uniprot.org/sparql/> {
        ?protein up:mnemonic ?proteinId .
        ?protein up:recommendedName ?rn .
        ?rn up:fullName ?proteinName .
        ?protein up:organism taxon:""" + taxid + """ .
        ?protein up:annotation ?a1 .
        ?a1 a up:Catalytic_Activity_Annotation .
        ?a1 up:catalyticActivity ?ca .
        ?ca up:catalyzedReaction ?reaction .
      }
    }
    """
    
    try:
        df = sparql2pandas(sparql_Q6, sparql_rhea_url)
    except Exception as e:
        df = None
        print('ERROR in SPARQL query')
        print(e)
    
    if 'reaction' in df.columns:
        df['MASTER_ID'] = df['reaction'].apply(lambda x: x.replace('http://rdf.rhea-db.org/', ''))
    
    return df
