
## Brief : Identify potential resistant variants from the solid media selection based on 
##           the frequency distribution of the variants with known phenotypes (variants in the DMS).            
## Inputs : 
##   - .csv with the frequency after selection of each variants (median of the replicates) (from variant_frequencies-solid-selection.py)
##   - .csv file with the results from the CaERG11 DMS (Bédard et al., 2024, supplementary table 3)
## Outputs : 
##   - Dataframe with the frequency after selection of each variant and its expected phenotype (sensitive, resistant) based on the DMS results
##   - Threshold values calculated from the frequency distribution of the known variants
##   - Dataframes of all the variants above the chosen frequency threshold and of the potentially resistant variants to validate individually
## Saved in an intermediate data folder : 
##   - Figure of the frequency distribution of the known variants with the different threshold values

#%% Import packages

import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from ast import literal_eval
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary


#%% Set-up 

# Version of the script (date; used to save files and plots)
version = "2026-04-29"

# Name of the experiment (used to save figures)
experiment = "solid-selection-analysis"

# Folder where to save the intermediate figures and dataframes
folder_figures = "../results/intermediate_data_selection-analysis/plots/"

# Folder where to save the final dataframes
folder_processed_data = "../results/"

# Color palette for the resistance phenotypes
colors_phenotypes = ["#e65354", "#7ea1fa"]
order_phenotypes = ["Resistant", "Sensitive"]

#%% Add the information from the DMS study (Bédard et al., 2024, Nat Micro) to the dataframe with the variant frequency after selection

# Variant frequencies after the selection
path = "../results/Variant_frequencies_after_selection_2026-04-29.csv"
df_freq = pd.read_csv(path, index_col = 0)
df_freq = df_freq[df_freq["variant_type"] != "wt"]

# Load the DMS results and filter for the voriconazole condition (which is misspelled in the published csv ahah)
path = "../data/dms_results.csv"
df_dms = pd.read_csv(path, index_col=0, sep = ";").reset_index().rename(columns={"Position" : "position"})
df_dms = df_dms[df_dms["Antifungal"] == "Voriconazoleco"]
df_dms["variant_id"] = df_dms["WT amino acid"] + df_dms["Variant"]     # Reformat the aa mutation code  

# Add the DMS info to the variant frequency dataframe
## Separate the df into categories to add the DMS results (DMS only has single mutants and synonymous variants)
df_single = df_freq[df_freq["hamming_dist_DNA"] == 1]
df_rest = df_freq[df_freq["hamming_dist_DNA"] != 1]

df_single["variant_id"] = [literal_eval(x)[0] for x in df_single["variant_id"]]
df_single = df_single.merge(df_dms[["variant_id", "Selection coefficient", "Category"]], on = "variant_id", how = "left")

## Concatenate the two dataframes back together
df_effect = pd.concat([df_single, df_rest])
df_effect["Category"] = df_effect["Category"].fillna(value = "Not in the DMS")

# Rename the "Category" column to better reflect what it actually is
df_effect = df_effect.rename(columns={"Category" : "expected_effect"})

# Convert the expected effect of the mutation to the expected resistance phenotype (beneficial -> resistant; neutral and deleterious -> sensitive)
phenotypes = {"Beneficial" : "Resistant", "Neutral" : "Sensitive", "Deleterious" : "Sensitive", "Not in the DMS" : "Not in the DMS"}
df_effect["expected_phenotype"] = [phenotypes[x] for x in df_effect["expected_effect"]]

# Get the number of previously characterized variants that were sequenced after the selection
## In this case, only the single mutants were considered as characterized in the DMS. That means that a variant 
## with a characterized mutation and a silent mutation (eg. variant [G464S, G465G] is considered as "Not in the DMS")
df_characterized = df_effect[df_effect["expected_effect"] != "Not in the DMS"]
print(f"There are {len(df_characterized)} variants that were previously characterized in the DMS")

# Save the dataframe
name = folder_processed_data + "Comparison_solid-selection_DMS_" + version + ".csv"
df_effect.to_csv(name)

#%% Identify threshold values based on the 95, 97.5 and 99th percentile for the frequency of the sensitive variants
# These thresholds will be used to identify potentially resistant variants to validate

# Load the dataframe
path = folder_processed_data + "Comparison_solid-selection_DMS_" + version + ".csv"
df_effect = pd.read_csv(path, index_col=0)

# Calculate different thresholds
thresholds = {}
for percentile in [0.95, 0.975, 0.99]:
    threshold_value = df_effect[df_effect["expected_phenotype"] == "Sensitive"]["median_log10_freq"].quantile(percentile)
    thresholds[percentile] = threshold_value

df_thresholds = pd.DataFrame.from_dict(thresholds, orient = "index", columns=["threshold"]).reset_index(names = "percentile")

# Save the dataframe
name = folder_processed_data + "Frequency_threshold_values_" + experiment + "_" + version + ".csv"
df_thresholds.to_csv(name)

#%% Check the frequency distribution for the variants with expected phenotypes (from the DMS)

# Load the dataframes
path = folder_processed_data + "Comparison_solid-selection_DMS_" + version + ".csv"
df_effect = pd.read_csv(path, index_col=0)

path = folder_processed_data + "Frequency_threshold_values_" + experiment + "_" + version + ".csv"
df_thresholds = pd.read_csv(path, index_col=1)

## Only keep the variants for which we know the expected phenotype
df_fig = df_effect[df_effect["expected_phenotype"] != "Not in the DMS"]

## Make a kde plot
fig, ax = plt.subplots()
#fig, ax = plt.subplots(figsize=(5, 3))
ax = sns.kdeplot(data = df_fig, x = "median_log10_freq", 
                hue = "expected_phenotype", palette = colors_phenotypes,
                hue_order = order_phenotypes, cut = 0, fill = True,
                alpha = 0.25, linewidth = 2)

## Add the threshold lines
plt.axvline(x = df_thresholds.loc[0.95]["threshold"], color = "#a4a4a4", linestyle = "--", label = "95")
plt.axvline(x = df_thresholds.loc[0.975]["threshold"], color = "#757575", linestyle = "--", label = "97.5")
plt.axvline(x = df_thresholds.loc[0.99]["threshold"], color = "black", linestyle = "--", label = "99")

# Save the figure (intermediate data)
fig_name = folder_figures + "Frequency_distribution_DMS_variants_" + experiment + "_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300,bbox_inches="tight")

## From the graph, a threshold at the 99th percentile seems the best to limit having sensitive variants in the variants to validate list

#%% Check if there is an enrichment of the resistant variants in the variants with a frequency above the threshold value

# Load the dataframes
path = folder_processed_data + "Comparison_solid-selection_DMS_" + version + ".csv"
df_effect = pd.read_csv(path, index_col=0)
df_known = df_effect[df_effect["expected_phenotype"] != "Not in the DMS"]

path = folder_processed_data + "Frequency_threshold_values_" + experiment + "_" + version + ".csv"
df_thresholds = pd.read_csv(path, index_col=1)

# Get the variants that have a relative frequency higher than the threshold value (99th percentile of the sensitive variants,
# based on the figure in the cell above)
threshold = df_thresholds.loc[0.99]["threshold"]
df_above = df_known[df_known["median_log10_freq"] >= threshold]
df_below = df_known[df_known["median_log10_freq"] < threshold]

# Calculate the different numbers needed to perform the Fisher's exact test based on this contingency table
##                            Phenotype
##                   | Resistant | Sensitive |
##           -------------------------------
##           | Above |           |           |
## Threshold  -------------------------------
##           | Below |           |           |
##           -------------------------------

nb_resistant_above = len(df_above[df_above["expected_phenotype"] == "Resistant"])
nb_sensitive_above = len(df_above[df_above["expected_phenotype"] == "Sensitive"])

nb_resistant_below = len(df_below[df_below["expected_phenotype"] == "Resistant"])
nb_sensitive_below = len(df_below[df_below["expected_phenotype"] == "Sensitive"])

# Perform Fisher's exact test
table = [[nb_resistant_above, nb_sensitive_above],
         [nb_resistant_below, nb_sensitive_below]]

statistic, pvalue = stats.fisher_exact(table, alternative = "greater")
if pvalue < 0.05 :
    print(f"There is a significant enrichment of variants with known resistance mutations in the variants with the highest frequency after selection (Fisher's exact test, p-value = {pvalue})")

else : 
    print(f"There is not a significant enrichment of variants with known resistance mutations in the variants with the highest frequency after selection (Fisher's exact test, p-value = {pvalue})")


#%% Apply the same threshold to the variants with mutations of unknown effect to identify interesting variants to validate

# Load the dataframes
path = folder_processed_data + "Comparison_solid-selection_DMS_" + version + ".csv"
df_effect = pd.read_csv(path, index_col=0)

path = folder_processed_data + "Frequency_threshold_values_" + experiment + "_" + version + ".csv"
df_thresholds = pd.read_csv(path, index_col=1)

# Get a dataframe of all the variants with a frequency above the threshold
threshold = df_thresholds.loc[0.99]["threshold"]
df_above = df_effect[df_effect["median_log10_freq"] >= threshold]
name = folder_processed_data + "All_variants_above_freq_threshold_after_selection_" + version + ".csv"
df_above.to_csv(name)
print(f"There are {len(df_above)} variants with a frequency above the threshold")

# Get a dataframe of potentially resistant variants to validate using growth curves (novel mutations)
df_to_validate = df_above[df_above["expected_phenotype"] == "Not in the DMS"]
print(f"Out of the variants with a frequency above the threshold, {len(df_to_validate)} were not characterized in the DMS")
df_to_validate = df_to_validate[df_to_validate["variant_type"] != "synonymous"]     # Remove variants with synonymous mutations
print(f"There are {len(df_to_validate)} potentially resistant variants to validate.")

name = folder_processed_data + "Potentially_resistant_variants_to_validate_" + version + ".xlsx"
df_to_validate.to_excel(name)

# %%
