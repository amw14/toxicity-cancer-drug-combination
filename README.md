# Characterizing clinical toxicity in cancer combination therapies 

[![DOI](https://zenodo.org/badge/727986478.svg)](https://doi.org/10.5281/zenodo.18159238)

This repository contains the code used for "Characterizing clinical toxicity in cancer combination therapies". The paper is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.13.648641v1).


## Preprocess the Databases ##

1. Run `test_preprocessing_functions.ipynb` and make sure all tests pass. This will ensure that the following have been pre-processed:
- DrugBank Drug-Drug Interactions Database
    - All NA values for drug names and toxicity severity have been removed
    - Stored at `data_processed/drugbank_ddi.csv`
- DDInter Drug-Drug Interactions Database
    - All NA values for drug names and toxicity severity have been removed
    - Stored at `data_processed/ddinter_data.csv`
- DrugComb Drug Combination Synergy Scores Database
    - All remaining combinations have all synergy scores (Bliss, HSA, Loewe, S_max, S_mean, S_sum, and ZIP)
    - All remaining combinations have both drugs in the pair named
    - All remaining combinations are limited to cancer drug combinations
    - Stored at `data_processed/drugcomb_data.csv`
- DrugBank Drug Target Database
    - Parse XML file for drugs to their drug targets, compile into dataframe
    - Stored at `data_processed/drugbank_drug_targets.csv`
- Reactome Pathway Database
    - Parse Reactome files for lowest pathways and all pathways
    - Lowest pathways stored at `data_processed/reactome_lowest_pathways_homo_sapiens.csv`
    - All pathways stored at `data_processed/reactome_all_pathways_homo_sapiens.csv`
- Jaccard Similarity method is working
- Jonckheere Terpestra Test is working
- Known Toxicity and Synergy Score INtersection between DrugComb and DrugBank
    - Stored at `data_processed/drugbank_syntox_known.csv`
- Known Toxicity and Synergy Score Intersection between DrugComb and DDInter
    - Stored at `data_processed/ddinter_syntox_known.csv`

2. Find the intersection between all databases
- Run ddinter_parse_toxsyntarpathstring_intersection.ipynb
    - Manually adds in some missing DrugBank Target information, stores all drug targets from the (known DDInter + DrugComb + DrugBank Target) intersection in `data_processed/ddinter_syntoxtarg.csv`
    - Creates mapping of target proteins to UniProtIDs for mapping to Reactome Pathways/STRING, stored in `ddinter_syntoxtarg_UniprotIDs.txt`
    - Finds intersection with Reactome lowest and all pathways, stores all drug targets from known DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways intersection in `data_processed/ddinter_syntoxtarg_lowestpw.csv` and `data_processed/ddinter_syntoxtarg_allpw.csv` respectively
    - Then finds intersection by mapping UniProtIDs to STRING IDs, storing all drug targets from DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING in `data_processed/ddinter_syntoxtarg_lowestpw_string.csv` and `data_processed/ddinter_syntoxtarg_allpw_string.csv` respectively
    - Then finds the drug combinations that contain all information between DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING with the average target distance, stored in `data_processed/ddinter_processed_combos_syntoxtargallpw_string.csv`
    - Note, did not have to additionally filter by SMILES/Morgan Fingerprint representations because found all of these for DDInter
- Run drugbank_parse_toxsyntarpathstring_intersection.ipynb
    - Adds in some missing DrugBank Target information, from the academic licensed version of DrugBank, stores all drug targets from the (known DDInter + DrugComb + DrugBank Target) intersection in `data_processed/drugbank_syntoxtarg.csv`
    - Creates mapping of target proteins to UniProtIDs for mapping to Reactome Pathways/STRING, stored in `drugbank_syntoxtarg_UniprotIDs.txt`
    - Finds intersection with Reactome lowest and all pathways, stores all drug targets from known DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways intersection in `data_processed/drugbank_syntoxtarg_lowestpw.csv` and `data_processed/drugbank_syntoxtarg_allpw.csv` respectively
    - Then finds intersection by mapping UniProtIDs to STRING IDs, storing all drug targets from DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING in `data_processed/drugbank_syntoxtarg_lowestpw_string.csv` and `data_processed/drugbank_syntoxtarg_allpw_string.csv` respectively
    - Then finds the drug combinations that contain all information between DDInter + DrugComb + DrugBank Target + Reactome Lowest/All Pathways + STRING with the average target distance, stored in `data_processed/drugbank_processed_combos_syntoxtargallpw_string.csv`
    - Note, did not have filter by SMILES/Morgan Fingerprint representations because found all of these for DDInter

## Synergy Score and Toxicity Analysis: Create Figure 2 ##

The first step of the analysis is to examine the relationship between synergy scores and toxicity. To do this and create the panels used in Figure 2 as well as Supplemental Figure S1, there are two files that need to be run:
1. DrugBank: `drugbank_syntox_analysis.ipynb`
2. DDInter: `ddinter_syntox_analysis.ipynb`

The figure files are stored in `results/figure2/` and the statistics are stored in:
- `results/synergy_score_tox_categ_distrib/kruskal_dunn_jonck_majmodmin_syntox_drugbank.tsv`
- `results/synergy_score_tox_categ_distrib/kruskal_dunn_jonck_majmodmin_syntox_ddinter.tsv`

There are supplementary histograms and statistics stored in `results/synergy_score_tox_categ_distrib` and `results/synergy_score_distrib`

## Drug Target and Pathway Analysis: Create Figure 3 and Figure 4 ##

The second step of the analysis is to examine which drug targets and pathways contribute to synergy scores and toxicity. To create Figure 3, Figure 4, and Supplemental Figure S2, the following files should be run:
1. DrugBank: `drugbank_targpw_analysis.ipynb`
2. DDInter: `ddinter_targpw_analysis.ipynb`

The figure files and statistical test results are stored in `results/jaccard_similarity_analysis/`

## Toxicity Score Analysis: Create Figure 5 ##

Figure 5 has three major components:
1. Assessing whether the average distance of targets on the STRING PPIN trends with toxicity levels -- run:
- `ddinter_avg_string_distance_analysis.ipynb`
- `drugbank_avg_string_distance_analysis.ipynb`
2. Assessing whether overlapping drug structure trends with toxicity levels -- run:
- `ddinter_deeptrasynergy_toxscore.ipynb`
- `drugbank_deeptrasynergy_toxscore.ipynb`
3. Assessing whether overlapping drug target protein neighborhoods on the STRING PPIN trends with toxicity levels -- run:
- `ddinter_graphsynergy_toxscore.ipynb`
- `drugbank_graphsynergy_toxscore.ipynb`

## Contact
If you have any questions, or find any issues with the code, please open an issue in this repository. We also welcome any contributions to the code - be sure to checkout the Contributing section below.

If you have questions or concerns with this project and do not want to create an issue, please contact [Alexandra Wong](mailto:alexandra_wong@brown.edu), [Cecile Meier-Scherling](cecile_meier-scherling@brown.edu), or [Lorin Crawford](mailto:lcrawford@microsoft.com).

## Relevant Citation (BibTeX)
If you use this code or method in your research, please consider citing the following:

```
@article {characterizing_toxicity,
	author = {Wong, Alexandra M and Meier-Scherling, Cecile and Crawford, Lorin},
	title = {Characterizing clinical toxicity in cancer combination therapies},
	elocation-id = {2025.04.13.648641},
	year = {2025},
	doi = {10.1101/2025.04.13.648641},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/04/19/2025.04.13.648641},
	eprint = {https://www.biorxiv.org/content/early/2025/04/19/2025.04.13.648641.full.pdf},
	journal = {bioRxiv}
}
```
## License

This project is available under the MIT License.
