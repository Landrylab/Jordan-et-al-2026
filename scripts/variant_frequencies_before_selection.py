
## Brief : Get the frequency before selection of the variants that were found in multiple replicates after selection          
## Inputs : 
##   - .csv with the read count of each variants before selection (output from the Aviti_Variant_calling_mutagenesis-library script)
##   - .csv with the total number of reads that passed the quality filters for each sample (output from the Aviti_Variant_calling_mutagenesis-library script)
##   - .csv with the variants that were found in multiple replicates after selection and their expected phenotype (output from the identify_potentially_resistant_variants script)
## Outputs : 
##   - Dataframe with the variant frequency before and after selection and their expected phenotype from the DMS
## Saved in an intermediate data folder : 
##   - Figures : correlation between replicates, correlation of the frequency before and after the screen

#%% Import packages

import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from ast import literal_eval
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary

#%% Set-up 

# Version of the script (date; used to save files and plots)
version = "2026-04-29"

# Name of the experiment (used to save figures)
experiment = "solid-selection-analysis"

# Folders where to save the intermediate figures
folder_figures = "../results/intermediate_data_selection-analysis/plots/"

# Folder where to save the final dataframe
folder_processed_data = "../results/"

# %% Load the dataframes

# Variants before selection
path = "../results/Variant_calling_all_samples_mutagenesis-library_2026-04-15.csv"
df_before = pd.read_csv(path, index_col=0)

# Total read counts for each sample
path = "../results/Total_reads_variant_calling_mutagenesis-library_2026-04-15.csv"
df_read_counts = pd.read_csv(path, index_col = 1).drop("Unnamed: 0", axis = 1)

# Variants in mulitple replicates after selection and their expected phenotype
path = "../results/Comparison_solid-selection_DMS_2026-04-29.csv"
df_after = pd.read_csv(path, index_col = 0)


#%% Calculate the relative frequency of each variant before the selection

for sample in ["F4-500", "F4-650", "F4-700", "F4-750"]:
    total_read_count = df_read_counts.loc[sample]["total_read_count"]

    df_sample_variants = df_before[df_before["sample_name"] == sample]
    df_sample_variants["freq"] = [x/total_read_count for x in df_sample_variants["read_count"]]

    # Put the samples back together in one dataframe
    if sample == "F4-500":
        df_before_freq = df_sample_variants

    else : 
        df_before_freq = pd.concat([df_before_freq, df_sample_variants])

# Get the log10 frequency
df_before_freq["log10_freq"] = np.log10(df_before_freq["freq"])

#%% Only keep the variants that are found in multiple replicates after the selection

# Reformat the variant_id column to fit with the df_after dataframe
df_before_single = df_before_freq[df_before_freq["hamming_dist_DNA"] == 1]
df_before_rest = df_before_freq[df_before_freq["hamming_dist_DNA"] != 1]

df_before_single["variant_id"] = [literal_eval(x)[0] for x in df_before_single["variant_id"]]
df_before_freq = pd.concat([df_before_single, df_before_rest])

# Keep only the variants found in at least two replicates after selection
list_variants_to_keep = df_after["variant_id"].to_list()
df_before_freq = df_before_freq[df_before_freq["variant_id"].isin(list_variants_to_keep)]

#%% Check the correlation of the variant frequency between replicates before selection for the variants found in multiple replicates after selection

replicates = [500, 650, 700, 750]

# Create a figure that will be a correlation matrix
fig, axes = plt.subplots(4, 4, figsize=(20, 20), sharex = True, sharey = True)

# Plot a scatterplot for each pair of replicates
for i in range(4):
    for j in range(4):
        ax = axes[i, j]
        rep1 = replicates[j]
        rep2 = replicates[i]

        df_rep1 = df_before_freq[(df_before_freq["condition"] == rep1) & (df_before_freq["variant_type"] != "wt")]
        df_rep2 = df_before_freq[(df_before_freq["condition"] == rep2) & (df_before_freq["variant_type"] != "wt")]

        ### Rename the column of the frequency for each replicate dataframe
        df_rep1 = df_rep1.rename(columns = {"log10_freq" : "log10_freq_rep1"})
        df_rep2 = df_rep2.rename(columns = {"log10_freq" : "log10_freq_rep2"})

        ### Get only the sequences present in both replicates
        # Merge the dataframes of the 2 replicates
        df_replicates = df_rep1.merge(df_rep2[["dna_mutation", "log10_freq_rep2"]], on = ["dna_mutation"], how = "inner")

        ### Spearman's correlation 
        rho, p_value = stats.spearmanr(df_replicates["log10_freq_rep1"], df_replicates["log10_freq_rep2"])
    
        ### Make the scatterplots
        sns.scatterplot(data = df_replicates, x = "log10_freq_rep1", y = "log10_freq_rep2",
                        edgecolor = "black", ax = ax)
    
        # Axis parameters
        ax.set_xlim(-6, -0)
        ax.set_ylim(-6, -0)

        # Label and legend parameters
        ax.set_xlabel(rep1, fontsize = 18)
        ax.set_ylabel(rep2, fontsize = 18)
        ax.tick_params(labelsize = 12)

        # Draw a diagonal line
        x_vals = np.array([-5.5, -0])
        y_vals = x_vals
        ax.plot(x_vals, y_vals, '--', color="#989898")

        # Add the rho value
        text = "ρ = " + str(round(rho,2))
        ax.text(-5.5, -0.5, text , fontsize = 18)

        # Add the number of sequences present in the 2 replicates
        text = "n = " + str(len(df_replicates))
        ax.text(-5.5, -1, text , fontsize = 18)
    
# Remove plots in the upper right part of the figure (same as the bottom left)
for i in range(4):
    fig.delaxes(axes[i, i])
    for j in range(i + 1, 4):
        fig.delaxes(axes[i, j])

fig.tight_layout()

# Save the figure
fig_name = folder_figures + "Log10_variant_freq_correlation_between_replicates_before_selection_for_variants_in_multiple_replicates_after_selection_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300, bbox_inches="tight")

#%% Get the median frequency before selection for each variant found in multiple replicates after selection

# Calculate the median frequency
df_median_freq_before = df_before_freq.groupby(["variant_id", "hamming_dist_DNA", "hamming_dist_AA", "variant_type", "fragment"])["log10_freq"].median().reset_index(name = "median_log10_freq")
print(f"The median frequency before selection was calculated for {len(df_median_freq_before)} variants")

# Merge with the df_after 
## Rename the frequency columns
df_median_freq_before = df_median_freq_before.rename(columns = {"median_log10_freq" : "median_log10_freq_before"})
df_after = df_after.rename(columns = {"median_log10_freq" : "median_log10_freq_after"})

## Merge
df_change = df_median_freq_before.merge(df_after[["variant_id", "median_log10_freq_after", "Selection coefficient", "expected_effect", "expected_phenotype"]], 
                                        on = "variant_id", how = "inner")

# Check for which variants I don't have the initial frequency (probably because the sequencing coverage was too low)
before_after_variants = df_change["variant_id"].to_list()
df_no_before = df_after[~df_after["variant_id"].isin(before_after_variants)]

# Save the dataframe
name = folder_processed_data + "Variant_frequencies_before_after_selection_" + version + ".csv"
df_change.to_csv(name)

#%% Check the correlation between the frequency before and after selection

# Load the dataframe
path = folder_processed_data + "Variant_frequencies_before_after_selection_" + version + ".csv"
df_change = pd.read_csv(path, index_col=0)

# Make a scatter plot
fig, ax = plt.subplots(figsize = (5, 5))

sns.scatterplot(data = df_change, x = "median_log10_freq_before", y = "median_log10_freq_after", ax = ax,
                    hue = "expected_phenotype", palette = ["#e65354", "#7ea1fa", "#cdcdcd"], hue_order = ["Resistant", "Sensitive", "Not in the DMS"],
                    edgecolor = "black")

# Save the figure
fig_name = folder_figures + "Correlation_frequency_before_after_selection_all_variants_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300, bbox_inches="tight")

# Same figure but with just the known variants
df_fig = df_change[df_change["expected_phenotype"] != "Not in the DMS"]

fig, ax = plt.subplots(figsize = (5, 5))

sns.scatterplot(data = df_fig, x = "median_log10_freq_before", y = "median_log10_freq_after", ax = ax,
                    hue = "expected_phenotype", palette = ["#e65354", "#7ea1fa"], hue_order = ["Resistant", "Sensitive"],
                    edgecolor = "black")

# Save the figure
fig_name = folder_figures + "Correlation_frequency_before_after_known_variants_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300, bbox_inches="tight")

#%% Interactive plot of the correlation between the frequency before and after selection

# Load the dataframe
path = folder_processed_data + "Variant_frequencies_before_after_selection_" + version + ".csv"
df_change = pd.read_csv(path, index_col=0)

# Make an interactive scatterplot
fig = px.scatter(df_change, x = "median_log10_freq_before", y = "median_log10_freq_after",
                 hover_name = "variant_id", color = "expected_phenotype", 
                 color_discrete_map = {"Resistant" : "#e65354",
                                       "Sensitive" : "#7ea1fa", 
                                       "Not in the DMS" : "#cdcdcd"},
                width = 700, height = 500)
fig.show()

# %%
