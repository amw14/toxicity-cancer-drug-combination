import pandas as pd

# Get the drug comb data
# INPUT:
#   None
# OUTPUT:
#   drug_comb_data: (DataFrame) the drug comb data
def get_drug_comb_data():
    # Read the drug comb data
    drugcomb_df = pd.read_csv('data/DrugComb/drugcomb_summary_v_1_5.csv', sep=',', index_col=False)

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
    
    # Drop the values if they are antagonistic
    # 'synergy_zip', 'synergy_loewe', 'synergy_hsa', 'synergy_bliss',

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

    sider_all_side_effects_df = pd.read_csv('data/SIDER4.1/meddra_all_se.tsv', sep='\t', index_col=False)
    sider_all_side_effects_df.columns = ["CID_FLAT", "CID_STEREO", "UMLS_Label", "MedDRA_Concept", "UMLS_MedDRA", "Side_Effect"]
    return sider_cid_to_drugs_df, sider_all_side_effects_df


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
