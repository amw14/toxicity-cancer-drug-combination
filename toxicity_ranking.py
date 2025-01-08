import math
import numpy as np
import pandas as pd
from scipy.stats import norm

# Get the drug comb data
# INPUT:
#   bliss: (bool) whether to filter out the antagonistic values for the bliss modality
#   loewe: (bool) whether to filter out the antagonistic values for the loewe modality
#   hsa: (bool) whether to filter out the antagonistic values for the hsa modality
#   zip: (bool) whether to filter out the antagonistic values for the zip modality
# OUTPUT:
#   drug_comb_data: (DataFrame) the drug comb data
def get_drug_comb_data(bliss=False, loewe=False, hsa=False, zip=False):
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
    if loewe:
        drugcomb_df['synergy_loewe'] = pd.to_numeric(drugcomb_df['synergy_loewe'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_loewe'])
    if hsa:
        drugcomb_df['synergy_hsa'] = pd.to_numeric(drugcomb_df['synergy_hsa'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_hsa'])
    if zip:
        drugcomb_df['synergy_zip'] = pd.to_numeric(drugcomb_df['synergy_zip'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_zip'])

    # Convert all values in drug_row and drug_col to lowercase
    drugcomb_df['drug_row'] = drugcomb_df['drug_row'].str.lower()
    drugcomb_df['drug_col'] = drugcomb_df['drug_col'].str.lower()
        
    print("Final shape of filtered drugcomb data: ", drugcomb_df.shape)

    return drugcomb_df


# Get the SIDER data
# INPUT:
#   None
# OUTPUT:
#   sider_cid_to_drugs_df: (DataFrame) the CID to drug name mapping
#   sider_all_side_effects_df: (DataFrame) the SIDER side effects data
def get_sider_data():
    # Read the SIDER data
    sider_cid_to_drugs_df = pd.read_csv('data/SIDER4.1/drug_names.tsv', sep='\t', index_col=False)
    sider_cid_to_drugs_df.columns = ["CID", "drug_name"]

    # Convert all drug_name to lowercase
    sider_cid_to_drugs_df['drug_name'] = sider_cid_to_drugs_df['drug_name'].str.lower()

    sider_all_side_effects_df = pd.read_csv('data/SIDER4.1/meddra_all_se.tsv', sep='\t', index_col=False)
    sider_all_side_effects_df.columns = ["CID_FLAT", "CID_STEREO", "UMLS_Label", "MedDRA_Concept", "UMLS_MedDRA", "Side_Effect"]
    return sider_cid_to_drugs_df, sider_all_side_effects_df


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

    return ddinter_df


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
    drug_syntox_df = pd.DataFrame(columns=['drug_row', 'drug_col', 'cell_line_name', 'synergy_zip', 'synergy_loewe', 'synergy_bliss', 'synergy_hsa', 'toxicity_category'])

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
        if drug_pair in ddinter_major_pairs_both_ways and reverse_drug_pair not in intersection_major_pairs:
            intersection_major_pairs.add(drug_pair)
            row_to_add = [row['drug_row'], row['drug_col'], row['cell_line_name'], row['synergy_zip'], row['synergy_loewe'], row['synergy_bliss'], row['synergy_hsa'], 'Major']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=['drug_row', 'drug_col', 'cell_line_name', 'synergy_zip', 'synergy_loewe', 'synergy_bliss', 'synergy_hsa', 'toxicity_category'])])
        elif drug_pair in ddinter_moderate_pairs_both_ways and reverse_drug_pair not in intersection_moderate_pairs:
            intersection_moderate_pairs.add(drug_pair)
            row_to_add = [row['drug_row'], row['drug_col'], row['cell_line_name'], row['synergy_zip'], row['synergy_loewe'], row['synergy_bliss'], row['synergy_hsa'], 'Moderate']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=['drug_row', 'drug_col', 'cell_line_name', 'synergy_zip', 'synergy_loewe', 'synergy_bliss', 'synergy_hsa', 'toxicity_category'])])
        elif drug_pair in ddinter_minor_pairs_both_ways and reverse_drug_pair not in intersection_minor_pairs:
            intersection_minor_pairs.add(drug_pair)
            row_to_add = [row['drug_row'], row['drug_col'], row['cell_line_name'], row['synergy_zip'], row['synergy_loewe'], row['synergy_bliss'], row['synergy_hsa'], 'Minor']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=['drug_row', 'drug_col', 'cell_line_name', 'synergy_zip', 'synergy_loewe', 'synergy_bliss', 'synergy_hsa', 'toxicity_category'])])
        elif drug_pair in ddinter_unknown_pairs_both_ways and reverse_drug_pair not in intersection_unknown_pairs:
            intersection_unknown_pairs.add(drug_pair)
            row_to_add = [row['drug_row'], row['drug_col'], row['cell_line_name'], row['synergy_zip'], row['synergy_loewe'], row['synergy_bliss'], row['synergy_hsa'], 'Unknown']
            drug_syntox_df = pd.concat([drug_syntox_df, pd.DataFrame([row_to_add], columns=['drug_row', 'drug_col', 'cell_line_name', 'synergy_zip', 'synergy_loewe', 'synergy_bliss', 'synergy_hsa', 'toxicity_category'])])
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

    return drug_syntox_df, intersection_major_pairs, intersection_moderate_pairs, intersection_minor_pairs, intersection_unknown_pairs



# Get the drug name to all unique side effects dictionary
# INPUT:
#   sider_cid_to_drugs_df: (DataFrame) the CID to drug name mapping
#   sider_all_side_effects_df: (DataFrame) the SIDER side effects data
# OUTPUT:
#   drug_name_to_side_effects: (dict) the drug name to all unique side effects dictionary
# Get the unique side effect set for each drug
def get_drug_to_side_effects(sider_cid_to_drugs_df, sider_all_side_effects_df):
    drug_name_to_side_effects = {}

    # Loop through each CID in the sider_cid_to_drugs_df
    for index, row in sider_cid_to_drugs_df.iterrows():
        # Get the CID
        CID = row['CID']
        # Get the drug name
        drug_name = row['drug_name']
        # Get the side effects for the CID
        side_effects = sider_all_side_effects_df[sider_all_side_effects_df['CID_FLAT'] == CID]['Side_Effect'].values
        # Store the side effects in the dictionary
        drug_name_to_side_effects[drug_name] = set(side_effects)

    return drug_name_to_side_effects


# Filter the drugcomb data -- remove the drugs that are not in the SIDER data
# INPUT:
#   drugcomb_df: (DataFrame) the drug comb data
#   sider_cid_to_drugs_df: (DataFrame) the CID to drug name mapping
# OUTPUT:
#   filtered_drug_comb_data: (DataFrame) the filtered drug comb data
#   common_drugs: (set) the set of common drugs in the drugcomb data and the SIDER data
#   unique_drug_pairs: (set) the set of unique drug pairs in the filtered drug comb data
def filter_drug_comb_data_by_sider(drugcomb_df, sider_cid_to_drugs_df):
    # Print the shape of the drugcomb data
    print("Original drugcomb data shape: ", drugcomb_df.shape)

    # Get the set of drugs in the drugcomb data and the SIDER data
    drugcomb_drugs = set(drugcomb_df['drug_row']).union(set(drugcomb_df['drug_col']))
    sider_drugs = set(sider_cid_to_drugs_df['drug_name'])
    common_drugs = sider_drugs.intersection(drugcomb_drugs)

    # Print how many drugs are common between drugcomb and sider
    print("Number of drugs in common between drugcomb and sider [lowercase enforced]: ", len(common_drugs))

    # Filter drugcomb data to only include drugs that are in the intersection
    filtered_drug_comb_data = drugcomb_df[(drugcomb_df['drug_row'].str.lower().isin(common_drugs)) & \
                                        (drugcomb_df['drug_col'].str.lower().isin(common_drugs))]
    
    # Print the shape of the filtered drugcomb data
    print("Filtered drugcomb data shape for both drugs being present in sider: ", filtered_drug_comb_data.shape)

    # How many unique drug pairs are there?
    unique_drug_pairs = set()

    for drug_row, drug_col in filtered_drug_comb_data[['drug_row', 'drug_col']].values:
        pair_one = (drug_row, drug_col)
        pair_two = (drug_col, drug_row)
        if pair_one not in unique_drug_pairs and pair_two not in unique_drug_pairs:
            unique_drug_pairs.add(pair_one)

    print("Number of unique drug pairs: ", len(unique_drug_pairs))

    return filtered_drug_comb_data, common_drugs, unique_drug_pairs


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


# Get the Jaccard similarity of all drug pairs
# INPUT:
#   unique_drug_pairs: (set) the set of unique drug pairs
#   sider_cid_to_drugs_df: (DataFrame) the CID to drug name mapping
#   sider_all_side_effects_df: (DataFrame) the SIDER side effects data
# OUTPUT:
#   drug_pair_to_jaccard: (dict) the drug pair to Jaccard similarity dictionary
#   drug_pair_to_side_effects: (dict) the drug pair to side effects dictionary
def drug_pair_to_jaccard_similarity(
    unique_drug_pairs,
    sider_cid_to_drugs_df,
    sider_all_side_effects_df,
):
    # Create a dictionary to store the Jaccard similarity of each drug pair
    drug_pair_to_jaccard = {}
    drug_pair_to_side_effects = {}

    # Loop through each drug pair in unique_drug_pairs
    for drug_pair in unique_drug_pairs:
        # Get the drugs in the pair
        drug1_name, drug2_name = drug_pair
        drug1_CID = sider_cid_to_drugs_df[sider_cid_to_drugs_df['drug_name'] == drug1_name]['CID'].values[0]
        drug2_CID = sider_cid_to_drugs_df[sider_cid_to_drugs_df['drug_name'] == drug2_name]['CID'].values[0]

        # Get the side effects of each drug
        side_effects1 = sider_all_side_effects_df[sider_all_side_effects_df['CID_FLAT'] == drug1_CID]['Side_Effect'].values
        side_effects2 = sider_all_side_effects_df[sider_all_side_effects_df['CID_FLAT'] == drug2_CID]['Side_Effect'].values

        # Calculate the Jaccard similarity of the side effects
        jaccard = jaccard_similarity(set(side_effects1), set(side_effects2))

        # Store the Jaccard similarity in the dictionary
        drug_pair_to_jaccard[drug_pair] = jaccard
        drug_pair_to_side_effects[drug_pair] = (side_effects1, side_effects2)

    return drug_pair_to_jaccard, drug_pair_to_side_effects


# Get rank based on similarity values
# INPUT:
#   drug_pair_to_similarity: (dict) the drug pair to some similarity value dictionary
# OUTPUT:
#   ranked_drug_pairs: (list) the ranked drug pairs based on Jaccard similarity
def rank_drug_pairs(drug_pair_to_similarity):
    # zip the drug pairs and their similarities
    drug_pairs = list(drug_pair_to_similarity.keys())
    similarity_values = list(drug_pair_to_similarity.values())
    zipped = list(zip(drug_pairs, similarity_values))

    # sort the zipped list by similarity value
    zipped.sort(key=lambda x: x[1], reverse=True)
    return zipped


def jonckheere_terpestra_test(samples):
    """
    Perform the Jonckheere-Terpstra test on the given samples.

    Parameters:
        samples: An array of arrays, where each inner array is a group containing the samples.

    Returns:
        A tuple containing the test statistic and the p-value.
    """
    
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
    p_value = 1 - norm.cdf(z_stat)
    
    return z_stat, p_value