# toxicity-cancer-drug-combination

## Preprocess the Databases ##

1. Run `test_toxicity_ranking.ipynb` and make sure all tests pass. This will ensure that the following have been pre-processed:
- DrugBank Drug-Drug Interactions Database
    - All NA values for drug names and toxicity severity have been removed
    - Stored at 'data_processed/drugbank_ddi.csv'
- DDInter Drug-Drug Interactions Database
    - All NA values for drug names and toxicity severity have been removed
    - Stored at 'data_processed/ddinter_data.csv'
- DrugComb Drug Combination Synergy Scores Database
    - All remaining combinations have all synergy scores (Bliss, HSA, Loewe, S_max, S_mean, S_sum, and ZIP)
    - All remaining combinations have both drugs in the pair named
    - All remaining combinations are limited to cancer drug combinations
    - Stored at 'data_processed/drugcomb_data.csv'
- Jaccard Similarity method is working
