# toxicity-cancer-drug-combination

## Preprocess the Databases ##

1. Run `test_preprocessing_functions.ipynb` and make sure all tests pass. This will ensure that the following have been pre-processed:
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
- DrugBank Drug Target Database
    - Parse XML file for drugs to their drug targets, compile into dataframe
    - Stored at 'data_processed/drugbank_drug_targets.csv'
- Reactome Pathway Database
    - Parse Reactome files for lowest pathways and all pathways
    - Lowest pathways stored at 'data_processed/reactome_lowest_pathways_homo_sapiens.csv'
    - All pathways stored at 'data_processed/reactome_all_pathways_homo_sapiens.csv'
- Jaccard Similarity method is working
- Jonckheere Terpestra Test is working
- Known Toxicity and Synergy Score INtersection between DrugComb and DrugBank
    - Stored at 'data_processed/drugbank_syntox_known.csv'
- Known Toxicity and Synergy Score Intersection between DrugComb and DDInter
    - Stored at 'data_processed/ddinter_syntox_known.csv'