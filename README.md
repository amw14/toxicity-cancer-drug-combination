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

2. Find the intersection between all databases
- Run ddinter_parse_toxsyntarpathstring_intersection.ipynb
    - Manually adds in some missing DrugBank Target information, stores all drug targets from the (known DDInter + DrugComb + DrugBank Target) intersection in 'data_processed/ddinter_syntoxtarg.csv'
    - Creates mapping of target proteins to UniProtIDs for mapping to Reactome Pathways/STRING, stored in 'ddinter_syntoxtarg_UniprotIDs.txt'
    - Finds intersection with Reactome lowest and all pathways, stores all drug targets from known DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways intersection in 'data_processed/ddinter_syntoxtarg_lowestpw.csv' and 'data_processed/ddinter_syntoxtarg_allpw.csv' respectively
    - Then finds intersection by mapping UniProtIDs to STRING IDs, storing all drug targets from DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING in 'data_processed/ddinter_syntoxtarg_lowestpw_string.csv' and 'data_processed/ddinter_syntoxtarg_allpw_string.csv' respectively
    - Then finds the drug combinations that contain all information between DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING with the average target distance, stored in 'data_processed/ddinter_processed_combos_syntoxtargallpw_string.csv'
    - Note, did not have to additionally filter by SMILES/Morgan Fingerprint representations because found all of these for DDInter
- Run drugbank_parse_toxsyntarpathstring_intersection.ipynb
    - Adds in some missing DrugBank Target information, from the academic licensed version of DrugBank, stores all drug targets from the (known DDInter + DrugComb + DrugBank Target) intersection in 'data_processed/drugbank_syntoxtarg.csv'
    - Creates mapping of target proteins to UniProtIDs for mapping to Reactome Pathways/STRING, stored in 'drugbank_syntoxtarg_UniprotIDs.txt'
    - Finds intersection with Reactome lowest and all pathways, stores all drug targets from known DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways intersection in 'data_processed/drugbank_syntoxtarg_lowestpw.csv' and 'data_processed/drugbank_syntoxtarg_allpw.csv' respectively
    - Then finds intersection by mapping UniProtIDs to STRING IDs, storing all drug targets from DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING in 'data_processed/drugbank_syntoxtarg_lowestpw_string.csv' and 'data_processed/drugbank_syntoxtarg_allpw_string.csv' respectively
    - Then finds the drug combinations that contain all information between DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING with the average target distance, stored in 'data_processed/drugbank_processed_combos_syntoxtargallpw_string.csv'
    - Note, did not have filter by SMILES/Morgan Fingerprint representations because found all of these for DDInter