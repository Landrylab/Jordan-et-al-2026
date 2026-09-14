
## Brief : Get the frequency of each variant after the selection on solid media.            
## Inputs : 
##   - .csv with the read count of each variants (output from the Aviti_Variant_calling_solid-selection script)
##   - .csv with the total number of reads that passed the quality filters for each sample (output from the Aviti_Variant_calling_solid-selection script)
## Outputs : 
##   - Dataframe with the median variant frequencies after the solid selection
## Saved in an intermediate data folder : 
##   - Figures : correlation between replicates, distribution of frequency according to the mutation type (non-syn., syn., stop)


#%% Import packages

import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary

#%% Set-up 

# Version of the script (date; used to save files and plots)
version = "2026-04-29"

# Name of the experiment (used to save figures)
experiment = "solid-selection-analysis"

# Folder where to save the intermediate figures
if not os.path.exists("../results/intermediate_data_selection-analysis/plots/"):
    os.makedirs("../results/intermediate_data_selection-analysis/plots/")

folder_figures = "../results/intermediate_data_selection-analysis/plots/"

# Folder where to save the intermediate dataframes
if not os.path.exists("../results/intermediate_data_selection-analysis/dataframes/"):
    os.makedirs("../results/intermediate_data_selection-analysis/dataframes/")

folder_intermediate_data = "../results/intermediate_data_selection-analysis/dataframes/"

# Folder where to save the final dataframe
folder_processed_data = "../results/"

# %% Load the dataframes

# Variants
path = "../results/Variant_calling_all_samples_solid_selection_2026-04-15.csv"
df_variants = pd.read_csv(path, index_col=0)

# Read counts
path = "../results/Total_reads_variant_calling_solid_selection_2026-04-15.csv"
df_read_counts = pd.read_csv(path, index_col = 1).drop("Unnamed: 0", axis = 1)

#%% Calculate the relative frequency of each variant

for sample in ["F4-500", "F4-650", "F4-700", "F4-750"]:
    total_read_count = df_read_counts.loc[sample]["total_read_count"]

    df_sample_variants = df_variants[df_variants["sample_name"] == sample]
    df_sample_variants["freq"] = [x/total_read_count for x in df_sample_variants["read_count"]]

    # Sanity check : sum of the relative frequencies should equal 1
    print(df_sample_variants["freq"].sum())

    # Put the samples back together in one dataframe
    if sample == "F4-500":
        df_variants_freq = df_sample_variants

    else : 
        df_variants_freq = pd.concat([df_variants_freq, df_sample_variants])

# Get the log10 frequency
df_variants_freq["log10_freq"] = np.log10(df_variants_freq["freq"])

#%% Check the correlation in the variant frequency between replicates

replicates = [500, 650, 700, 750]

# Create a figure that will be a correlation matrix
fig, axes = plt.subplots(4, 4, figsize=(20, 20), sharex = True, sharey = True)

# Plot a scatterplot for each pair of replicates
for i in range(4):
    for j in range(4):
        ax = axes[i, j]
        rep1 = replicates[j]
        rep2 = replicates[i]

        df_rep1 = df_variants_freq[(df_variants_freq["condition"] == rep1) & (df_variants_freq["variant_type"] != "wt")]
        df_rep2 = df_variants_freq[(df_variants_freq["condition"] == rep2) & (df_variants_freq["variant_type"] != "wt")]

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
        ax.set_xlim(-5.5, -0)
        ax.set_ylim(-5.5, -0)

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
        ax.text(-5.25, -0.5, text , fontsize = 18)

        # Add the number of sequences present in the 2 replicates
        text = "n = " + str(len(df_replicates))
        ax.text(-5.25, -1, text , fontsize = 18)
    
# Remove plots in the upper right part of the figure (same as the bottom left)
for i in range(4):
    fig.delaxes(axes[i, i])
    for j in range(i + 1, 4):
        fig.delaxes(axes[i, j])

fig.tight_layout()

# Save the figure
fig_name = folder_figures + "Log10_variant_freq_correlation_between_replicates_after_selection_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300, bbox_inches="tight")

#%% Get the median frequency for variants sequenced in at least two replicates for the same fragment

# Get a list of the variants sequenced in multiple replicates
df_mult = df_variants_freq.groupby(by = ["nt_seq", "fragment"])["sample_name"].count().reset_index(name = "nb_replicates")
variants_in_mult_reps = df_mult[df_mult["nb_replicates"] >= 2]["nt_seq"].to_list()

# Filter the df to only keep sequences in two or more replicates
df_in_mult = df_variants_freq[df_variants_freq["nt_seq"].isin(variants_in_mult_reps)]

# Group by variant and calculate the median
## This step calculates the median of the replicates, but also for synonymous sequences (sequences with different DNA mutations that lead to the same aa mutations),
## but this does not combine variants with a silent mutations together
## Eg. variants with different DNA mutations leading to the F449L substitution (TTT449TTA and TTT449TTG) are grouped, 
## but variants with only F449L substitution are not grouped with variant with another silent mutation (eg. with variant [F449L, G465G])
df_freq = df_in_mult.groupby(["variant_id", "hamming_dist_DNA", "hamming_dist_AA", "variant_type", "fragment"])["log10_freq"].median().reset_index(name = "median_log10_freq")

# Get the number of variants in the df
print(f"There are {len(df_freq)} variants for which we calculated the median frequency after selection")

# Save the dataframe
name = folder_processed_data + "Variant_frequencies_after_selection_" + version + ".csv"
df_freq.to_csv(name)

#%% Check the distribution of the frequencies after selection by type of mutation

# Load the frequency dataframe
path = folder_processed_data + "Variant_frequencies_after_selection_" + version + ".csv"
df_freq = pd.read_csv(path, index_col=0)

# Make the boxplots
fig, ax = plt.subplots(figsize=(6, 4))
ax = sns.boxplot(data = df_freq[df_freq["variant_type"] != "wt"], x = "variant_type", y = "median_log10_freq",
                 color = "#cdcdcd", showfliers = False)

sns.stripplot(data = df_freq[df_freq["variant_type"] != "wt"], x = "variant_type", y = "median_log10_freq", ax = ax, dodge = True,
             color = "#cdcdcd", size = 3, edgecolor = "black", linewidth = 0.5, legend = False, alpha = 0.8, jitter=0.25)
    
ax.set_xlabel("Variant type")
ax.set_ylabel("Variant frequency (log10)")

# Save the figure
fig_name = folder_figures + "Boxplots_variant_frequency_by_type_after_selection_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300, bbox_inches="tight")

# %%
