import pandas as pd

# Get the drug comb data
def get_drug_comb_data():
    # Read the drug comb data
    drug_comb_data = pd.read_csv('data/drug_combinations.csv')
    return drug_comb_data

# Get the SIDER data
def get_sider_data():
    # Read the SIDER data
    sider_data = pd.read_csv('data/meddra_all_se.csv')
    return sider_data

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
