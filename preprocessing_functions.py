import math
import networkx as nx
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from collections import defaultdict
from scipy.stats import norm

# Get DrugBank DDI dataframe
# INPUT:
#   None
# OUTPUT:
#   drugbank_ddi_df: (DataFrame) the DrugBank drug-drug interaction data
def get_drugbank_ddi():
    drugbank_ddi_fp = 'data/DrugBank/u_Brown_DrugBank_CSV/structured_drug_interactions.csv'
    drugbank_ddi_df = pd.read_csv(drugbank_ddi_fp)

    # Answer a few questions about the dataset to begin with?
    print('How many interactions are in each severity category in the drugbank database?')
    print(drugbank_ddi_df['severity'].value_counts())

    # make all drug names lowercase
    drugbank_ddi_df['subject_drug_name'] = drugbank_ddi_df['subject_drug_name'].str.lower()
    drugbank_ddi_df['affected_drug_name'] = drugbank_ddi_df['affected_drug_name'].str.lower()

    # create a dataframe mapping 
    subject_drug_df = drugbank_ddi_df[['subject_drug_id', 'subject_drug_drugbank_id', 'subject_drug_name']]
    subject_drug_df = subject_drug_df.drop_duplicates()
    subject_drug_df = subject_drug_df.reset_index(drop=True)
    subject_drug_df.columns = ['drug_id', 'drugbank_id', 'drug_name']
    affected_drug_df = drugbank_ddi_df[['affected_drug_id', 'affected_drug_drugbank_id', 'affected_drug_name']]
    affected_drug_df = affected_drug_df.drop_duplicates()
    affected_drug_df = affected_drug_df.reset_index(drop=True)
    affected_drug_df.columns = ['drug_id', 'drugbank_id', 'drug_name']
    drugbank_drug_ids_to_names = pd.concat([subject_drug_df, affected_drug_df], axis=0)
    drugbank_drug_ids_to_names = drugbank_drug_ids_to_names.drop_duplicates()
    drugbank_drug_ids_to_names = drugbank_drug_ids_to_names.reset_index(drop=True)
    drugbank_drug_ids_to_names.to_csv('data_processed/drugbank_drug_ids_to_names.csv', index=False)

    # drop rows where subject_drug_name or affected_drug_name is NaN
    drugbank_ddi_df = drugbank_ddi_df.dropna(subset=['subject_drug_name', 'affected_drug_name'])
    # drop rows where subject_drug_name or affected_drug_name is empty
    drugbank_ddi_df = drugbank_ddi_df[drugbank_ddi_df['subject_drug_name'] != '']
    drugbank_ddi_df = drugbank_ddi_df[drugbank_ddi_df['affected_drug_name'] != '']
    
    # print the shape of the dataframe
    print('Shape of drugbank ddi dataframe: ', drugbank_ddi_df.shape)

    drugbank_ddi_df.to_csv('data_processed/drugbank_ddi.csv', index=False)

    return drugbank_ddi_df, drugbank_drug_ids_to_names

# Get the drug comb data
# INPUT:
#   bliss: (bool) whether to drop NA values for bliss modality
#   hsa: (bool) whether to drop NA values for hsa modality
#   loewe: (bool) whether to drop NA values for loewe modality
#   s_mean: (bool) whether to drop NA values for s_mean modality
#   s_max: (bool) whether to drop NA values for s_max modality
#   s_sum: (bool) whether to drop NA values for s_sum modality
#   zip: (bool) whether to drop NA values for zip modality
# OUTPUT:
#   drug_comb_data: (DataFrame) the drug comb data
def get_drug_comb_data(bliss=False, hsa=False, loewe=False, s_mean=False, s_max=False, s_sum=False, zip=False):
    # Read the drug comb data
    drugcomb_df = pd.read_csv('data/DrugComb/drugcomb_summary_v_1_5.csv', sep=',', index_col=False)
    print("Original shape of drugcomb data: ", drugcomb_df.shape)

    # Filter out NaN values in the drug_row and drug_col and cell line columns
    drugcomb_df = drugcomb_df.dropna(subset=['drug_row', 'drug_col', 'cell_line_name'])

    # Drop values where malaria tissue is present or was studied in SARS-CoV-2 or Ebola data
    malaria_filter = drugcomb_df['tissue_name'] != 'malaria'
    MOTT_filter = drugcomb_df['study_name'] != 'MOTT'
    ELLINGER_filter = drugcomb_df['study_name'] != 'ELLINGER'
    GORDON_filter = drugcomb_df['study_name'] != 'GORDON'
    TOURET_filter = drugcomb_df['study_name'] != 'TOURET'
    NCATS_SARS_COV_2DPI_filter = drugcomb_df['study_name'] != 'NCATS_SARS-COV-2DPI'
    BOBROWSKI_filter = drugcomb_df['study_name'] != 'BOBROWSKI'
    DYALL_filter = drugcomb_df['study_name'] != 'DYALL'

    drugcomb_df = drugcomb_df[malaria_filter & MOTT_filter & ELLINGER_filter & GORDON_filter \
                            & TOURET_filter & NCATS_SARS_COV_2DPI_filter & BOBROWSKI_filter \
                            & DYALL_filter]
    #print("Shape after filtering out malaria and SARS-CoV-2 and Ebola data: ", drugcomb_df.shape)
    
    # Do not filter out antagonistic values, only filter out NaN values
    if bliss:
        # drop any rows where synergy_bliss cannot be converted to a number
        drugcomb_df['synergy_bliss'] = pd.to_numeric(drugcomb_df['synergy_bliss'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_bliss'])
    if hsa:
        drugcomb_df['synergy_hsa'] = pd.to_numeric(drugcomb_df['synergy_hsa'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_hsa'])
    if loewe:
        drugcomb_df['synergy_loewe'] = pd.to_numeric(drugcomb_df['synergy_loewe'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_loewe'])
    if s_mean:
        drugcomb_df['S_mean'] = pd.to_numeric(drugcomb_df['S_mean'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['S_mean'])
    if s_max:
        drugcomb_df['S_max'] = pd.to_numeric(drugcomb_df['S_max'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['S_max'])
    if s_sum:
        drugcomb_df['S_sum'] = pd.to_numeric(drugcomb_df['S_sum'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['S_sum'])
    if zip:
        drugcomb_df['synergy_zip'] = pd.to_numeric(drugcomb_df['synergy_zip'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_zip'])

    # Convert all values in drug_row and drug_col to lowercase
    drugcomb_df['drug_row'] = drugcomb_df['drug_row'].str.lower()
    drugcomb_df['drug_col'] = drugcomb_df['drug_col'].str.lower()
        
    print("Final shape of filtered drugcomb data: ", drugcomb_df.shape)

    drugcomb_df.to_csv('data_processed/drugcomb_data.csv', index=False)

    return drugcomb_df


# Get the DDInter data
# INPUT:
#   None
# OUTPUT:
# ddinter_df: (DataFrame) the DDInter severity data
def get_ddinter_data():
    # Get the DDInter data
    ddinter_files = [
        'data/DDInter/ddinter_downloads_code_A.csv',
        'data/DDInter/ddinter_downloads_code_B.csv',
        'data/DDInter/ddinter_downloads_code_D.csv',
        'data/DDInter/ddinter_downloads_code_H.csv',
        'data/DDInter/ddinter_downloads_code_L.csv',
        'data/DDInter/ddinter_downloads_code_P.csv',
        'data/DDInter/ddinter_downloads_code_R.csv',
        'data/DDInter/ddinter_downloads_code_V.csv',
    ]
    ddinter_df = None
    for file in ddinter_files:
        if ddinter_df is None:
            ddinter_df = pd.read_csv(file, sep=',')
        else:
            ddinter_df = pd.concat([ddinter_df, pd.read_csv(file, sep=',')], ignore_index=True)

    # Convert all values in Drug_A and Drug_B columns to lowercase
    ddinter_df['Drug_A'] = ddinter_df['Drug_A'].str.lower()
    ddinter_df['Drug_B'] = ddinter_df['Drug_B'].str.lower()

    # Remove duplicates
    ddinter_df = ddinter_df.drop_duplicates()

    # Drop rows where Drug_A or Drug_B is NaN
    ddinter_df = ddinter_df.dropna(subset=['Drug_A', 'Drug_B'])
    # Drop rows where Drug_A or Drug_B is empty
    ddinter_df = ddinter_df[ddinter_df['Drug_A'] != '']
    ddinter_df = ddinter_df[ddinter_df['Drug_B'] != '']
    # Drop rows where Drug_A or Drug_B is 'nan'
    ddinter_df = ddinter_df[ddinter_df['Drug_A'] != 'nan']
    ddinter_df = ddinter_df[ddinter_df['Drug_B'] != 'nan']

    print("DDInter Shape without NA values: " + str(ddinter_df.shape))
    ddinter_df.to_csv('data_processed/ddinter_data.csv', index=False)

    return ddinter_df


# Get the DrugBank Target Data
# INPUT:
#   xml_file: (str) the path to the DrugBank XML file
# OUTPUT:
#   drugbank_target_df: (DataFrame) the DrugBank target data
def parse_drugbank_xml(xml_file='data/DrugBank/full database.xml'):
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Define the namespace used in DrugBank XML
    ns = {'db': 'http://www.drugbank.ca'}
    
    # List to store all drug-target pairs
    drug_target_pairs = []

    # Dictionary to store SMILES for each drug
    drug_smiles_dict = {}
    
    # Iterate through all drugs
    for drug in root.findall('db:drug', ns):
        drug_name = drug.find('db:name', ns).text
        
        # Get all targets for the drug
        targets = drug.findall('.//db:target', ns)
        
        for target in targets:
            target_data = {
                'drug_name': drug_name,
                'target_name': None,
                'target_DrugBank_ID': None,
                'GenBank_Protein_ID': None,
                'GenBank_Gene_ID': None,
                'UniProtKB_ID': None,
                'GenAtlas_ID': None,
                'HGNC_ID': None,
                'SMILES': None,
            }
            
            # Get target name
            polypeptide = target.find('.//db:polypeptide', ns)
            if polypeptide is not None:
                target_name = polypeptide.find('db:name', ns)
                if target_name is not None:
                    target_data['target_name'] = target_name.text
            
            # Get target DrugBank ID
            target_id = target.find('.//db:id', ns)
            if target_id is not None:
                target_data['target_DrugBank_ID'] = target_id.text
            
            # Get external identifiers
            if polypeptide is not None:
                external_ids = polypeptide.findall('.//db:external-identifier', ns)
                for ext_id in external_ids:
                    resource = ext_id.find('db:resource', ns).text
                    identifier = ext_id.find('db:identifier', ns).text
                    
                    if resource == 'GenBank Protein Database':
                        target_data['GenBank_Protein_ID'] = identifier
                    elif resource == 'GenBank Gene Database':
                        target_data['GenBank_Gene_ID'] = identifier
                    elif resource == 'UniProtKB':
                        target_data['UniProtKB_ID'] = identifier
                    elif resource == 'GenAtlas':
                        target_data['GenAtlas_ID'] = identifier
                    elif resource == 'HUGO Gene Nomenclature Committee (HGNC)':
                        target_data['HGNC_ID'] = identifier
            
                # Get SMILES
                drug_smiles = drug.find('.//db:calculated-properties/db:property[db:kind="SMILES"]/db:value', ns)
                if drug_smiles is not None:
                    smiles_text = drug_smiles.text
                    if smiles_text != '' and smiles_text is not None:
                        target_data['SMILES'] = smiles_text
                        if drug_name not in drug_smiles:
                            drug_smiles_dict[drug_name] = smiles_text
                        else:
                            if drug_smiles_dict[drug_name] != smiles_text:
                                print(f"SMILES mismatch for {drug_name}")
            
            drug_target_pairs.append(target_data)
    
    # Create DataFrame
    df = pd.DataFrame(drug_target_pairs)

    df['drug_name'] = df['drug_name'].str.lower() # lower case

    df.to_csv('data_processed/drugbank_drug_targets.csv', index=False)

    return df


# Get the Reactome Data
# INPUT:
#   None
# OUTPUT:
#   lowest_pathways_df: (DataFrame) the Reactome data with the lowest pathways
#   all_pahtways_df: (DataFrame) the Reactome data with all pathways
def get_reactome_data():
    lowest_pathway_reactome_file_path = 'data/Reactome/UniProt2Reactome.tsv'
    all_pathways_reactome_file_path = 'data/Reactome/UniProt2Reactome_All_Levels.tsv'

    lowest_pathways_df = pd.read_csv(lowest_pathway_reactome_file_path, sep='\t', names=['UniProtKB_ID', 'Reactome_ID', 'Reactome_URL', 'Pathway_Name', 'Evidence_Code', 'Species'])
    print("Original lowest pathways shape: " + str(lowest_pathways_df.shape))
    lowest_pathways_df = lowest_pathways_df[lowest_pathways_df['Species'] == 'Homo sapiens']

    all_pathways_df = pd.read_csv(all_pathways_reactome_file_path, sep='\t', names=['UniProtKB_ID', 'Reactome_ID', 'Reactome_URL', 'Pathway_Name', 'Evidence_Code', 'Species'])
    print("Original all pathways shape: " + str(all_pathways_df.shape))
    all_pathways_df = all_pathways_df[all_pathways_df['Species'] == 'Homo sapiens']

    # Drop the Reactome URL, Evidence Code, and Species columns (unnecessary)
    lowest_pathways_df.drop(columns=['Reactome_URL', 'Evidence_Code', 'Species'], inplace=True)
    all_pathways_df.drop(columns=['Reactome_URL', 'Evidence_Code', 'Species'], inplace=True)

    # Save to CSV file
    lowest_pathways_df.to_csv('data_processed/reactome_lowest_pathways_homo_sapiens.csv', index=None)
    all_pathways_df.to_csv('data_processed/reactome_all_pathways_homo_sapiens.csv', index=None)

    print("Filtered lowest pathways shape restricting to Homo sapiens and removing unnecessary columns: " + str(lowest_pathways_df.shape))
    print("Filtered all pathways shape restricting to Homo sapiens and removing unnecessary columns: " + str(all_pathways_df.shape))

    return lowest_pathways_df, all_pathways_df


# Get the STRING graph
# INPUT:
#   string_fp: (str) - File path to STRING edge list
# OUTPUT:
#   STRING_G: (NetworkX graph) - graph with the STRING PPIN network
def get_STRING_graph(string_fp='data/STRING/9606.protein.physical.links.detailed.v12.0.txt'):
    string_edge_list_df = pd.read_csv(string_fp, sep=' ')
    print('Original shape of STRING edge list, physical detailed: ' + str(string_edge_list_df.shape))
    string_edge_list_df = string_edge_list_df[string_edge_list_df['experimental'] != 0]
    STRING_G = nx.from_pandas_edgelist(string_edge_list_df, 'protein1', 'protein2')
    STRING_G = STRING_G.to_undirected()
    return STRING_G


# Filter the drugcomb data -- remove the drugs that are not in the ddinter data
# INPUT:
#   drugcomb_df: (DataFrame) the drug comb data
#   ddinter_df: (DataFrame) the CID to drug name mapping
# OUTPUT:
#   drug_syntox_df: (DataFrame) the filtered drug comb data
#   intersection_major_pairs: (set) - drugcomb pairs with major severity in ddinter
#   intersection_moderate_pairs: (set) - drugcomb pairs with moderate severity in ddinter 
#   intersection_minor_pairs: (set) - drugcomb pairs with minor severity in ddinter
#   intersection_unknown_pairs: (set) - drugcomb pairs with unknown severity in ddinter
def find_drugcomb_ddinter_intersect(drugcomb_df, ddinter_df):
    # Get the set of drugs in the drugcomb data and the ddinter data
    drugcomb_drugs = set(drugcomb_df['drug_row']).union(set(drugcomb_df['drug_col']))
    ddinter_drugs = set(ddinter_df['Drug_A']).union(set(ddinter_df['Drug_B']))
    common_drugs = ddinter_drugs.intersection(drugcomb_drugs)

    # Create filtered drugcomb data that includes only drugs that are in the intersection and the synergy scores
    drug_syntox_df = pd.DataFrame(columns=[
        'drug_row',
        'drug_col',
        'cell_line_name',
        'synergy_bliss',
        'synergy_hsa',
        'synergy_loewe',
        'synergy_zip',
        'S_max',
        'S_mean',
        'S_sum',
        'toxicity_category',
    ])

    # Print how many drugs are common between drugcomb and ddinter
    print("Number of drugs in common between drugcomb and ddinter [lowercase enforced]: ", len(common_drugs))

    # Want to find the common drug pairs between drugcomb and ddinter
    ddinter_major_pairs_both_ways = set()
    ddinter_moderate_pairs_both_ways = set()
    ddinter_minor_pairs_both_ways = set()
    ddinter_unknown_pairs_both_ways = set()
    for index, row in ddinter_df.iterrows():
        drug_a = row['Drug_A']
        drug_b = row['Drug_B']
        if row['Level'] == 'Unknown':
            ddinter_unknown_pairs_both_ways.add((drug_a, drug_b))
            ddinter_unknown_pairs_both_ways.add((drug_b, drug_a))
        elif row['Level'] == 'Minor':
            ddinter_minor_pairs_both_ways.add((drug_a, drug_b))
            ddinter_minor_pairs_both_ways.add((drug_b, drug_a))
        elif row['Level'] == 'Moderate':
            ddinter_moderate_pairs_both_ways.add((drug_a, drug_b))
            ddinter_moderate_pairs_both_ways.add((drug_b, drug_a))
        elif row['Level'] == 'Major':
            ddinter_major_pairs_both_ways.add((drug_a, drug_b))
            ddinter_major_pairs_both_ways.add((drug_b, drug_a))
        else:
            print("Level is: " + str(row['Level']))
    
    intersection_major_pairs = set()
    intersection_moderate_pairs = set()
    intersection_minor_pairs = set()
    intersection_unknown_pairs = set()
    for index, row in drugcomb_df.iterrows():
        drug_pair = (row['drug_row'], row['drug_col'])
        reverse_drug_pair = (row['drug_col'], row['drug_row'])
        partial_row = [
                row['drug_row'],
                row['drug_col'],
                row['cell_line_name'],
                row['synergy_bliss'],
                row['synergy_hsa'],
                row['synergy_loewe'],
                row['synergy_zip'],
                row['S_max'],
                row['S_mean'],
                row['S_sum'],
        ]
        if drug_pair in ddinter_major_pairs_both_ways and reverse_drug_pair not in intersection_major_pairs:
            intersection_major_pairs.add(drug_pair)
            row_to_add = partial_row + ['Major']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        elif drug_pair in ddinter_moderate_pairs_both_ways and reverse_drug_pair not in intersection_moderate_pairs:
            intersection_moderate_pairs.add(drug_pair)
            row_to_add = partial_row + ['Moderate']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        elif drug_pair in ddinter_minor_pairs_both_ways and reverse_drug_pair not in intersection_minor_pairs:
            intersection_minor_pairs.add(drug_pair)
            row_to_add = partial_row + ['Minor']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        elif drug_pair in ddinter_unknown_pairs_both_ways and reverse_drug_pair not in intersection_unknown_pairs:
            intersection_unknown_pairs.add(drug_pair)
            row_to_add = partial_row + ['Unknown']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        else:
            continue
    
    print("Major pairs in both DrugComb and in DDInter: ", len(intersection_major_pairs))
    print("Moderate pairs in both DrugComb and in DDInter: ", len(intersection_moderate_pairs))
    print("Minor pairs in both DrugComb and in DDInter: ", len(intersection_minor_pairs))
    print("Unknown toxicity pairs in both DrugComb and in DDInter: ", len(intersection_unknown_pairs))

    total_num_common_pairs = len(intersection_unknown_pairs) + len(intersection_major_pairs) + \
          len(intersection_moderate_pairs) + len(intersection_minor_pairs)

    print("Total common pairs: ", total_num_common_pairs)
    print("Total known pairs: ", total_num_common_pairs - len(intersection_unknown_pairs))

    # Save the known toxicity data to a csv file
    known_syntox_df = drug_syntox_df[drug_syntox_df['toxicity_category'] != 'Unknown']
    known_syntox_df.to_csv('data_processed/ddinter_syntox_known.csv', index=False)

    return known_syntox_df, intersection_major_pairs, intersection_moderate_pairs, intersection_minor_pairs, intersection_unknown_pairs


# Filter the drugcomb data -- remove the drugs that are not in the DrugBank data
# INPUT:
#   drugcomb_df: (DataFrame) the DrugComb data
#   drugbank_ddi_df: (DataFrame) the DrugBank drug-drug interaction data
# OUTPUT:
#   drugbank_syntox_df: (DataFrame) the filtered drug comb data
#   intersection_major_pairs: (set) - drugcomb pairs with major severity in DrugBank
#   intersection_moderate_pairs: (set) - drugcomb pairs with moderate severity in DrugBank 
#   intersection_minor_pairs: (set) - drugcomb pairs with minor severity in DrugBank
#   intersection_unknown_pairs: (set) - drugcomb pairs with unknown severity in DrugBank
def find_drugcomb_drugbankddi_intersect(drugcomb_df, drugbank_ddi_df):
    # Get the set of drugs in the drugcomb data and the drugbank data
    drugcomb_drugs = set(drugcomb_df['drug_row']).union(set(drugcomb_df['drug_col']))
    drugbank_drugs = set(drugbank_ddi_df['subject_drug_name']).union(set(drugbank_ddi_df['affected_drug_name']))
    common_drugs = drugbank_drugs.intersection(drugcomb_drugs)

    # Create filtered drugcomb data that includes only drugs that are in the intersection and the synergy scores
    drugbank_syntox_df = pd.DataFrame(columns=[
        'drug_row',
        'drug_col',
        'cell_line_name',
        'synergy_bliss',
        'synergy_hsa',
        'synergy_loewe',
        'synergy_zip',
        'S_max',
        'S_mean',
        'S_sum',
        'toxicity_category',
    ])

    # Print how many drugs are common between drugcomb and drugbank
    print("Number of drugs in common between drugcomb and drugbank [lowercase enforced]: ", len(common_drugs))

    # Want to find the common drug pairs between drugcomb and drugbank
    drugbank_major_pairs_both_ways = set()
    drugbank_moderate_pairs_both_ways = set()
    drugbank_minor_pairs_both_ways = set()
    for index, row in drugbank_ddi_df.iterrows():
        drug_a = row['subject_drug_name']
        drug_b = row['affected_drug_name']
        if row['severity'] == 0:
            drugbank_minor_pairs_both_ways.add((drug_a, drug_b))
            drugbank_minor_pairs_both_ways.add((drug_b, drug_a))
        elif row['severity'] == 1:
            drugbank_moderate_pairs_both_ways.add((drug_a, drug_b))
            drugbank_moderate_pairs_both_ways.add((drug_b, drug_a))
        elif row['severity'] == 2:
            drugbank_major_pairs_both_ways.add((drug_a, drug_b))
            drugbank_major_pairs_both_ways.add((drug_b, drug_a))
        else:
            print("severity is: " + str(row['severity']))
    
    intersection_major_pairs = set()
    intersection_moderate_pairs = set()
    intersection_minor_pairs = set()
    intersection_unknown_pairs = set()
    for index, row in drugcomb_df.iterrows():
        drug_pair = (row['drug_row'], row['drug_col'])
        reverse_drug_pair = (row['drug_col'], row['drug_row'])
        partial_row = [
                row['drug_row'],
                row['drug_col'],
                row['cell_line_name'],
                row['synergy_bliss'],
                row['synergy_hsa'],
                row['synergy_loewe'],
                row['synergy_zip'],
                row['S_max'],
                row['S_mean'],
                row['S_sum'],
        ]

        if drug_pair in drugbank_major_pairs_both_ways and reverse_drug_pair not in intersection_major_pairs:
            intersection_major_pairs.add(drug_pair)
            row_to_add = partial_row + ['Major']
            drugbank_syntox_df = pd.concat([drugbank_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        elif drug_pair in drugbank_moderate_pairs_both_ways and reverse_drug_pair not in intersection_moderate_pairs:
            intersection_moderate_pairs.add(drug_pair)
            row_to_add = partial_row + ['Moderate']
            drugbank_syntox_df = pd.concat([drugbank_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        elif drug_pair in drugbank_minor_pairs_both_ways and reverse_drug_pair not in intersection_minor_pairs:
            intersection_minor_pairs.add(drug_pair)
            row_to_add = partial_row + ['Minor']
            drugbank_syntox_df = pd.concat([drugbank_syntox_df, pd.DataFrame([row_to_add], columns=[
                'drug_row',
                'drug_col',
                'cell_line_name',
                'synergy_bliss',
                'synergy_hsa',
                'synergy_loewe',
                'synergy_zip',
                'S_max',
                'S_mean',
                'S_sum',
                'toxicity_category',
            ])])
        else:
            if drug_pair not in intersection_unknown_pairs and reverse_drug_pair not in intersection_unknown_pairs:
                intersection_unknown_pairs.add(drug_pair)
                row_to_add = partial_row + ['Unknown']
                drugbank_syntox_df = pd.concat([drugbank_syntox_df, pd.DataFrame([row_to_add], columns=[
                    'drug_row',
                    'drug_col',
                    'cell_line_name',
                    'synergy_bliss',
                    'synergy_hsa',
                    'synergy_loewe',
                    'synergy_zip',
                    'S_max',
                    'S_mean',
                    'S_sum',
                    'toxicity_category',
                ])])
            else:
                continue
    
    print("Major pairs in both DrugComb and in DrugBank: ", len(intersection_major_pairs))
    print("Moderate pairs in both DrugComb and in DrugBank: ", len(intersection_moderate_pairs))
    print("Minor pairs in both DrugComb and in DrugBank: ", len(intersection_minor_pairs))
    print("Unknown toxicity pairs in both DrugComb and in DrugBank: ", len(intersection_unknown_pairs))

    total_num_common_pairs = len(intersection_unknown_pairs) + len(intersection_major_pairs) + \
          len(intersection_moderate_pairs) + len(intersection_minor_pairs)

    print("Total common pairs: ", total_num_common_pairs)
    print("Total known pairs: ", total_num_common_pairs - len(intersection_unknown_pairs))

    # Save the known toxicity data to a csv file, all of these are known
    known_syntox_df = drugbank_syntox_df[drugbank_syntox_df['toxicity_category'] != 'Unknown']
    known_syntox_df.to_csv('data_processed/drugbank_syntox_known.csv', index=False)

    return known_syntox_df, intersection_major_pairs, intersection_moderate_pairs, intersection_minor_pairs, intersection_unknown_pairs


# Get the Jaccard similarity of two sets
# INPUT:
#   set1: (set) the first set
#   set2: (set) the second set
# OUTPUT:
#   (float) the Jaccard similarity of the two sets
def jaccard_similarity(set1, set2):
    # Get the intersection of the two sets
    intersection = len(set1.intersection(set2))
    # Get the union of the two sets
    union = len(set1.union(set2))
    # Return the Jaccard similarity
    return intersection / union


# Perform the Jonckheere-Terpstra test
# INPUT:
#   samples: (list) a list of lists, where each list is a sample
# OUTPUT:
#   z_stat: (float) the test statistic
#   p_value: (float) the p-value
def jonckheere_terpestra_test(samples):
    if not samples or len(samples) < 2:
        raise ValueError("At least two groups are required")
    
    jt_stat = 0 # initialize the test statistic
    n_i = [len(sample) for sample in samples] # get group sizes
    N = sum(n_i) # total number of samples

    for i in range(len(samples) - 1): # for each group
        for j in range(i + 1, len(samples)): # compare with all other groups
            for first_sample in samples[i]: # for each sample in the first group
                for second_sample in samples[j]: # for each sample in the second group
                    if first_sample < second_sample:
                        jt_stat += 1
                    elif first_sample == second_sample:
                        jt_stat += 0.5

    # Calculate mean under null hypothesis
    mean = ((N**2) - sum(size**2 for size in n_i)) / 4
    
    # Calculate variance under null hypothesis
    term1 = N**2 * (2*N + 3)
    term2 = sum(size**2 * (2*size + 3) for size in n_i)
    variance = (term1 - term2) / 72

    # Calculate standardized statistic
    z_stat = (jt_stat - mean) / math.sqrt(variance)

    # Calculate one-tail p-value
    p_value = 1.0 - norm.cdf(z_stat)
    
    return z_stat, p_value

# Calculates Cliff's Delta (non-parametric effect size for two groups).
# INPUT: 
#   x: (list or array-like) the first group of observations
#   y: (list or array-like) the second group of observations
# OUTPUT:
#   delta: (float) the Cliff's Delta value
def cliff_delta(x, y):
    """
    Calculates Cliff's Delta (non-parametric effect size for two groups).
    Delta quantifies the probability that a randomly selected observation
    from one group is larger/smaller than a randomly selected observation
    from the other group.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    n1 = len(x)
    n2 = len(y)

    # Initialize counter for concordant (c) and discordant (d) pairs
    # Cliff's Delta is (c - d) / (n1 * n2) where c is x > y and d is y > x.
    diff = x[:, None] - y
    c = np.sum(diff > 0)
    d = np.sum(diff < 0)
    delta = (c - d) / (n1 * n2)
    return delta