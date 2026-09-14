
## Brief : Validation of potentially resistant mutants by comparing their growth in the presence of voriconazole to the wild-type strain         
## Inputs : 
##   - .xlsx file with the raw data from the Tecan plate reader (OD600 measurements every 15min for 48h)
##   - .xlsx file with the plate plan
## Outputs : 
##   - .csv with the AUC value for each variant and the results of a t-test to see if it is significantly different than that of the wild-type
## Saved in an intermediate data folder :
##   - Dataframes at each step of the processing

#%% Import packages

import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary


#%% Set-up 

# Name of the experiment (used to save figures)
experiment = "potential_resistant_validation_VCZ"

# Folder where to save the intermediate figures and dataframes
if not os.path.exists("../results/intermediate_data_growth-curves/plots/"):
    os.makedirs("../results/intermediate_data_growth-curves/plots/")

folder_figures = "../results/intermediate_data_growth-curves/plots/"

if not os.path.exists("../results/intermediate_data_growth-curves/dataframes/"):
    os.makedirs("../results/intermediate_data_growth-curves/dataframes/")

folder_intermediate_data = "../results/intermediate_data_growth-curves/dataframes/"

# Folder where to save the final dataframes
folder_processed_data = "../results/"

# %% Process the data from the Tecan

# Load the dataframes
## Plate plan
path = "../data/2026-04-14_validation_mutants_VCZ_plate2_plate-plan.xlsx"
plan = pd.read_excel(path, index_col=0, header=0)

## Tecan raw data
path = "../data/2026-04-14_Validation_mutants_VCZ_plate2.xlsx"
data = pd.read_excel(path, header=0, skiprows=2, skipfooter=34).reset_index(drop=True)

# Processing
## Drop the temperature column (it has no name in the Tecan output, but it is the second column)
data.drop(data.columns[1], axis=1, inplace=True)
## Rename the time column (it has no name in the Tecan output)
data.rename(columns={data.columns[0] : "time"}, inplace=True)
## Convert timepoints from seconds to hours
def format_timepoints(timepoint):
    no_s = timepoint.strip("s")
    timepoint_h = np.float64(no_s)/3600
    return timepoint_h
    
data["time_h"] = data["time"].apply(format_timepoints)
data.drop("time", axis=1, inplace=True)

# Link the plate plan info to the OD data
## Go from a wide to a long format
df_long = data.melt(id_vars = "time_h", var_name = "well", value_name = "OD")
df_long

## Merge on "well" to associate every OD measurement to the strain and concentration
## (inner join by default; all wells not in the plate plan will be removed)
df_merged = df_long.merge(right=plan, on=["well"])
df_merged

# Correct for the blank and variation in the initial cell inoculate density by substracting the median of the first 3 measurements
df_first_3_TP = df_merged[df_merged["time_h"] <= 0.5]
df_to_subtract = df_first_3_TP.groupby("well")["OD"].median().reset_index(name = "median_OD").set_index("well")

df_merged["OD_corrected"] = [x - df_to_subtract.loc[y]["median_OD"] for x, y in zip(df_merged["OD"], df_merged["well"])]

# Save the dataframe 
df_name = folder_intermediate_data + experiment + "_OD_data.csv"
df_merged.to_csv(df_name)

#%% Visualize the growth curves for each well (like what you can see on the plate reader software)

# Specify the well order so that it matches the position in the plate
well_order = []
for i in ["A", "B", "C", "D", "E", "F", "G", "H"]:
    for j in range(1, 13):
        well_order.append(i+str(j))

# Load the growth data
df_merged = pd.read_csv(folder_intermediate_data + experiment + "_OD_data.csv", index_col = 0)

# Make the figure
grid = sns.FacetGrid(data = df_merged, col = "well", col_wrap = 12, col_order = well_order)
grid.map(sns.lineplot, "time_h", "OD_corrected")
grid.set_titles(col_template="{col_name}", fontsize = 16)
grid.set_axis_labels("Time (h)", "OD", fontsize = 14)
fig_name = folder_figures + experiment + "_growth_curves_by_well.png"
plt.savefig(fig_name, format = "png", dpi = 300, bbox_inches = "tight")

#%% Visualize the growth curves with all the strains in the same plot

# Load the growth data
df_merged = pd.read_csv(folder_intermediate_data + experiment + "_OD_data.csv", index_col = 0)

# Plot the growth curves
fig, ax = plt.subplots()
ax = sns.lineplot(data = df_merged, x = "time_h", y = "OD_corrected", hue = "strain")
ax.axvline(x=24, color = "grey", linestyle = "--")        # Add lines at 24, 30 and 36h to see where to cut off for the growth parameters calculations
ax.axvline(x=30, color = "grey", linestyle = "--")
ax.axvline(x=36, color = "grey", linestyle = "--")
ax.set_xlabel("Time (h)", fontsize = 14)
ax.set_ylabel("OD", fontsize = 14)
sns.move_legend(ax, loc = "upper right", bbox_to_anchor = (1.25, 1))

# Save the figure
fig_name = folder_figures + experiment + "_growth_curves_all.png"
plt.savefig(fig_name, format = "png", dpi = 300, bbox_inches = "tight")

#%% Calculate the AUC

# Based on the figure from the cell above, 30h seems like a good cutoff point to calculate the AUC
#  (end of exponential phase for the wild-type strain)

# Load the growth data and plan
df_merged = pd.read_csv(folder_intermediate_data + experiment + "_OD_data.csv", index_col = 0)
plan = pd.read_excel("../data/2026-04-14_validation_mutants_VCZ_plate2_plate-plan.xlsx", index_col=0, header=0)

# Calculate the AUC and merge with the plate to add to the strains info
def get_auc(group):
    return np.trapezoid(group.OD_corrected)

auc = df_merged[df_merged["time_h"] <= 30].groupby("well")[["OD_corrected"]].apply(get_auc).reset_index(name="AUC_30")
df_auc = plan.merge(right=auc, on="well")

# Drop rows with no strain
df_auc = df_auc.dropna()

# Save the dataframe
name = folder_intermediate_data + experiment + "_AUC_30h.csv"
df_auc.to_csv(name)

#%% Check which variants are resistant by comparing growth to the wild-type

# Load the dataframe with the AUC
df_auc = pd.read_csv(folder_intermediate_data + experiment + "_AUC_30h.csv", index_col = 0)

# Filter out the blank wells
df_auc = df_auc[df_auc["strain"] != "Blank"]
    
# Normalize the AUC of the mutants to validate by the mean AUC for the wild-type strain
wt_mean = df_auc[df_auc["strain"] == "WT"]["AUC_30"].mean()
df_auc["norm_AUC"] = df_auc["AUC_30"]/wt_mean

# Test if the AUC of the mutant is different than that of the WT with a two sample t-test (two-sided)
## Initialize a dict to store the results of the t-test
t_test_results = {}
p_values = []

## Get the AUC values for the wild-type
wt_AUC = df_auc[df_auc["strain"] == "WT"]["norm_AUC"]

## Test every mutant against the wild-type
mutants = df_auc[df_auc["strain"] != "WT"]["strain"].unique()         # List of all the mutants to test
for mut in mutants:
    mutant_AUC = df_auc[df_auc["strain"] == mut]["norm_AUC"]
    t_statistic, p_value = stats.ttest_ind(mutant_AUC, wt_AUC, alternative = "two-sided")   # Perform the t-test

    t_test_results[mut] = {"mean_norm_AUC" : mutant_AUC.mean(),
                            "sd" : np.std(mutant_AUC),
                            "p_value" : p_value}
        
    p_values.append(p_value)
    
## Apply FDR correction since I am doing multiple comparisons
_, adjusted_pvalues, _, _ = multipletests(p_values, method='fdr_bh')

## Store the results in a dataframe
df_ttest = pd.DataFrame.from_dict(t_test_results, orient = "index", columns = ["mean_norm_AUC", "sd", "p_value"]).reset_index(names = "variant")
df_ttest["adjusted_p_value"] = adjusted_pvalues

# Classify the variants as resistant, deleterious or neutral based on the result of the t-test
def classify_variants(adjusted_p_value, auc):

    if adjusted_p_value < 0.05 and auc > 1 :             # Mean of the WT normalized AUC is 1
        category = "Beneficial" 
    elif adjusted_p_value < 0.05 and auc < 1 :
        category = "Deleterious"
    else :
        category = "Neutral"

    return category

## Apply the classification function to the variants
df_ttest["phenotype"] = [classify_variants(adjusted_p_value, auc) for adjusted_p_value, auc in zip(df_ttest["adjusted_p_value"], df_ttest["mean_norm_AUC"])]

# Save the t-test results dataframe
name = folder_processed_data + experiment + "_Variant_classification_t-test.csv"
df_ttest.to_csv(name)    

# Save the dataframe with the normalized AUC
name = folder_processed_data + experiment + "_normalized_AUC-30h.csv"
df_auc.to_csv(name)

#%% Visualize the distribution of the normalized AUC values for each strain

# Load the dataframe with the AUC
df_auc = pd.read_csv(folder_processed_data + experiment + "_normalized_AUC-30h.csv", index_col = 0)

# Load the dataframe with the variants classification
path = folder_processed_data + experiment + "_Variant_classification_t-test.csv"
variant_classification = pd.read_csv(path, index_col=0)

# Add the variant classification to the df with the AUC values
df_auc = df_auc.merge(variant_classification[["variant", "phenotype"]], how = "left", left_on = "strain", right_on = "variant")
df_auc["phenotype"] = df_auc["phenotype"].fillna("Wild-type")

strain_order = ["WT", "pRS31N", "G464S", "A409T", "T411M", "R414S", "Y415C", "A430T", "K451I", "P512L", "P515T", "A516T"]
        
color_palette = ["#8CB369", "#F4E285", "#BC4B51", "#cdcdcd"]
color_order = ["Beneficial", "Neutral", "Deleterious", "Wild-type"]


# Make the boxplot
sns.set_style("ticks")
fig, ax = plt.subplots(figsize = (16, 8))

ax = sns.boxplot(data = df_auc, x = "strain", y = "norm_AUC",
                hue = "phenotype", palette = color_palette, hue_order = color_order,
                order = strain_order, showfliers = False, linewidth = 2)

## Show every point
sns.stripplot(data = df_auc, x = "strain", y = "norm_AUC", order = strain_order, ax = ax,
            hue = "phenotype", palette = color_palette, hue_order = color_order, 
            size = 7, edgecolor = "black", linewidth = 1, legend = False)
    
## Axis labels, ticks and legend parameters
ax.set_xlabel("Variant", fontsize = 16)
ax.set_ylabel("Normalized AUC", fontsize = 16)
ax.set_ylim(0, 2.5)
ax.tick_params(labelsize = 14)
sns.move_legend(ax, "upper center", ncol = 4, title = None, frameon = False, bbox_to_anchor = (0.5, 1.07), fontsize = 14)
sns.despine()

## Save the figure as a png and svg file
fig_name = folder_figures + "Resistant_mutants_validation_AUC30_boxplots"
plt.savefig(f"{fig_name}.png", format='png', dpi=300)

# %%
