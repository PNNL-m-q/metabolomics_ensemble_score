## Ensemble GC-MS Spectral Similarity Score

Publication: https://pubs.acs.org/doi/full/10.1021/jasms.5c00176

**Metadata**

1. compound_metadata.csv: Contains compound type information for a few compounds (e.g. whether they're amino acids, etc.)

1. score_metadata.csv: Contains all metadata information on each spectral similarity score

2. sample_metadata.csv: Contains all metadata information on each sample 

**Models**

1. model.RDS: The trained ensemble model with all scores

2. reduced_model.RDS: The trained ensemble model with the top 6 performing scores 

**Result_Data** 

1. BinSizes.csv: The number of candidate molecules per sample and retention index bin

2. FP_FN_Ranks.txt: Full model and reduced model predictions on the testing dataset

3. reduced_test_pred.RDS: An R object with the reduced model predictions on the testing dataset

4. test_pred.RDS: An R object with the full model predictions on the testing dataset

5. TP_Ranks.txt: Rankings of the true positive per sample and retention index bin for the top 6 scores, the full model, and the reduced model

*Note:* All other data used in this study is too large for a github repo and can be found here: https://data.pnnl.gov/group/nodes/dataset/33302

**Scripts**

1. build_dataset.R: Extracts all molecule information needed from this study after downloading https://data.pnnl.gov/group/nodes/dataset/33302

2. ensemble_model.R: Code to build the ensemble model after running build_dataset.R

3. false_positive_&_false_negative: Extracts all needed information about false positives and false negatives after running the ensemble model

4. top_N.R: Compares the true positive rankings of the built models and the top 6 scores 

**Visualization**

1. plots.R: Generates all visualizations of results for this study 

