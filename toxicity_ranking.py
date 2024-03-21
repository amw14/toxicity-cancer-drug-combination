import pandas as pd

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
    
    # Filter out antagonistic values
    if bliss:
        # drop any rows where synergy_bliss cannot be converted to a number
        drugcomb_df['synergy_bliss'] = pd.to_numeric(drugcomb_df['synergy_bliss'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_bliss'])
        drugcomb_df = drugcomb_df[drugcomb_df['synergy_bliss'] >= 0]
        #print("Shape after filtering out antagonistic bliss: ", drugcomb_df.shape)
    if loewe:
        drugcomb_df['synergy_loewe'] = pd.to_numeric(drugcomb_df['synergy_loewe'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_loewe'])
        drugcomb_df = drugcomb_df[drugcomb_df['synergy_loewe'] >= 0]
        #print("Shape after filtering out antagonistic loewe: ", drugcomb_df.shape)
    if hsa:
        drugcomb_df['synergy_hsa'] = pd.to_numeric(drugcomb_df['synergy_hsa'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_hsa'])
        drugcomb_df = drugcomb_df[drugcomb_df['synergy_hsa'] >= 0]
        #print("Shape after filtering out antagonistic hsa: ", drugcomb_df.shape)
    if zip:
        drugcomb_df['synergy_zip'] = pd.to_numeric(drugcomb_df['synergy_zip'], errors='coerce')
        drugcomb_df = drugcomb_df.dropna(subset=['synergy_zip'])
        drugcomb_df = drugcomb_df[drugcomb_df['synergy_zip'] >= 0]
        #print("Shape after filtering out antagonistic zip: ", drugcomb_df.shape)

    # Convert all values in drug_row and drug_col to lowercase
    drugcomb_df['drug_row'] = drugcomb_df['drug_row'].str.lower()
    drugcomb_df['drug_col'] = drugcomb_df['drug_col'].str.lower()
    
    # TODO: Drop columns that are not needed

    # TODO: Determine if want to do something with CSS or RI
        
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


# Filter the drugcomb data -- remove the drugs that are not in the SIDER data
# INPUT:
#   drugcomb_df: (DataFrame) the drug comb data
#   sider_cid_to_drugs_df: (DataFrame) the CID to drug name mapping
# OUTPUT:
#   filtered_drug_comb_data: (DataFrame) the filtered drug comb data
#   common_drugs: (set) the set of common drugs in the drugcomb data and the SIDER data
#   unique_drug_pairs: (set) the set of unique drug pairs in the filtered drug comb data
def filter_drug_comb_data(drugcomb_df, sider_cid_to_drugs_df):
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



if __name__ == "__main__":
    drugcomb_df = get_drug_comb_data(bliss=True, loewe=True, hsa=True, zip=True)
    sider_cid_to_drugs_df, sider_all_side_effects_df = get_sider_data()
    filtered_drug_comb_data, common_drugs, unique_drug_pairs = filter_drug_comb_data(drugcomb_df, sider_cid_to_drugs_df)

    