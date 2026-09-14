## Inputs :
### From the identify_potentially_resistant_variants script :
##   - .csv with the median frequency after the selection for the variants found in multiple replicates
##   - .csv with the threshold values (95, 97.5, 99 percentile of the frequency of the sensitive variants)
##   - .csv with the variants that have a frequency above the threshold (99th percentile)
### From the variant_frequencies_before_selection script :
##   - .csv with the median frequency before the selection for the variants found in multiple replicates after the selection
### From the potential_resistant_variants_growth_curves script :
##   - .csv with the normalized AUC data for the variants chosen for validation
##   - .csv with the results of the t-test comparing growth of the variants to the wild-type

## Outputs : 
##   - Figure 5 with panels C and D filled and panels A, B and E to be added in Inkscape
##   - Script to run in ChimeraX to generate pictures for panel E
##   - Figure S5

#%% Import packages

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np
from mpl_toolkits.axes_grid1 import inset_locator
import scipy.stats as stats

# Figure parameters
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary

# Version of the script (date; used to save files and plots)
version = "2026-08-03"

#%% Load the input dataframes

# Median frequency after selection with the expected phenotype
path = "../results/Comparison_solid-selection_DMS_2026-04-29.csv"
df_selection = pd.read_csv(path, index_col=0)

# Threshold values
path = "../results/Frequency_threshold_values_solid-selection-analysis_2026-04-29.csv"
df_thresholds = pd.read_csv(path, index_col=1)

# Variants above the threshold
path = "../results/All_variants_above_freq_threshold_after_selection_2026-04-29.csv"
df_variants_above = pd.read_csv(path, index_col=0)

# Frequency before selection for the variants found in multiple replicates after the selection
path = "../results/Variant_frequencies_before_after_selection_2026-04-29.csv"
df_before_after = pd.read_csv(path, index_col=1).drop(columns = ["Unnamed: 0"])

# AUC data for the variants chosen for validation
path = "../results/potential_resistant_validation_VCZ_normalized_AUC-30h.csv"
df_auc = pd.read_csv(path, index_col=0)

# Results of the t-test for the validation variants
path = "../results/potential_resistant_validation_VCZ_Variant_classification_t-test.csv"
df_ttest = pd.read_csv(path, index_col=1).drop(columns = ["Unnamed: 0"])

# Results of the DMS study (Bédard et al., 2024, supplementary table 3)
path = "../data/dms_results.csv"
df_dms = pd.read_csv(path, sep = ";")

#%% Prepare the dataframes used to make the figure 

# Frequency distributions
df_known = df_selection[df_selection["expected_phenotype"] != "Not in the DMS"]
df_unknown = df_selection[df_selection["expected_phenotype"] == "Not in the DMS"]

# Validation variants boxplots
## Add the variant classification to the auc dataframe to be able to color the boxplot based on the variant's phenotype 
df_auc = df_auc.merge(df_ttest["phenotype"], right_index = True, left_on = "strain", how = "left")
df_auc["phenotype"] = df_auc["phenotype"].fillna("WT")
categories = {"Beneficial" : "Resistant", "WT" : "Sensitive", "Neutral" : "Sensitive", "Deleterious" : "Sensitive"}
df_auc["Category"] = [categories[x] for x in df_auc["phenotype"]]

# Correlation of the frequency before and after
## Separate the known and unknown variants
df_known_before_after = df_before_after[(df_before_after["expected_phenotype"] != "Not in the DMS")]
df_unknown_before_after = df_before_after[(df_before_after["expected_phenotype"] == "Not in the DMS")]
## Separate the variants selection for validation
validation_variants = ["A409T", "T411M", "R414S", "Y415C", "A430T", "K451I", "P512L", "P515T", "A516T"]
df_validation = df_before_after.loc[validation_variants]


#%% Make figure 5

# Draw the figure and axes
## Create figure and gridspec
fig = plt.figure(constrained_layout = False, figsize = (6.8, 7.5))
gs = fig.add_gridspec(nrows = 75, ncols = 68, left=0, right=1, bottom=0, top=1)

## Draw axes for panels a and b (to be left blank)
ax1 = fig.add_subplot(gs[0:12, 0:32])
ax2 = fig.add_subplot(gs[0:12, 33:68])

## Draw axes for the frequency distribution histograms
ax3 = fig.add_subplot(gs[13:27, 0:23])
ax4 = fig.add_subplot(gs[30:42, 0:23])

## Draw the ax for the variant validation boxplots
ax5 = fig.add_subplot(gs[15:39, 29:67])

## Draw the ax for the protein structure visualization (to be left blank)
ax6 = fig.add_subplot(gs[50:70, 0:75])

## Make axes 1, 2 and 6 blank
ax1.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax1.spines[:].set_visible(False)

ax2.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax2.spines[:].set_visible(False)

ax6.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax6.spines[:].set_visible(False)

# Fill the other axes

###########################   Ax3 : frequency distribution of the known variants  ################################################

# Plot the distribution
sns.histplot(data = df_known, x = "median_log10_freq", ax = ax3, hue = "expected_phenotype", palette = ["#ff0000", "#0095ff"],
            hue_order = ["Resistant", "Sensitive"], fill = True, alpha = 0.25, linewidth = 0.7)

# Add a line for the threshold value
ax3.axvline(x = df_thresholds.loc[0.99]["threshold"], color = "black", linestyle = "--")

# Axis limits and tick parameters
ax3.set_xlim([-5.25, -1.5])
ax3.set_xticks([-5, -4, -3, -2])
ax3.tick_params("x", labelbottom = False)
ax3.set_ylim([0, 50])
ax3.set_yticks([0, 25, 50])
ax3.tick_params("y", labelsize = 8, pad = 2)

# Title and axis labels
ax3.set_title("Known variants", fontsize = 10, pad = 3)
ax3.set_xlabel("")
ax3.set_ylabel("Count", fontsize = 10, labelpad = 7)

# Legend
sns.move_legend(obj = ax3, loc = "upper left", title = None, frameon = False, bbox_to_anchor = (0, 1.05), 
                handletextpad = 0.6, handlelength = 1,  fontsize = 8, labelspacing = 0.2)

# Add the number of resistant and sensitive variants with a frequency over the threshold
nb_resistant = len(df_variants_above[df_variants_above["expected_phenotype"] == "Resistant"])
nb_sensitive = len(df_variants_above[df_variants_above["expected_phenotype"] == "Sensitive"])

ax3.text(-3, 20, f"{nb_resistant} resistant", fontsize = 10, color = "#e65354")
ax3.text(-2.9, 15, f"{nb_sensitive} sensitive", fontsize = 10, color = "#4792c6")

sns.despine(ax = ax3)

###########################   Ax4 : frequency distribution of the unknown variants  ################################################

df_unknown = df_selection[df_selection["expected_phenotype"] == "Not in the DMS"]

# Plot the frequency distribution
sns.histplot(data = df_unknown, x = "median_log10_freq", ax = ax4, hue = "expected_phenotype", palette = ["#8e6bb2"],
            hue_order = ["Not in the DMS"], fill = True, alpha = 0.5, linewidth = 0.7, legend = False)

# Add a line for the threshold value
ax4.axvline(x = df_thresholds.loc[0.99]["threshold"], color = "black", linestyle = "--")

# Axis limits and tick parameters
ax4.set_xlim([-5.25, -1.5])
ax4.set_xticks([-5, -4, -3, -2])
ax4.tick_params("x", labelsize = 8, pad = 2)
ax4.set_ylim([0, 205])
ax4.set_yticks([0, 50, 100, 150, 200])
ax4.tick_params("y", labelsize = 8, pad = 2)

# Title and axis labels
ax4.set_title("Unknown variants", fontsize = 10, pad = 3)
ax4.set_xlabel("Variant frequency (log10)", fontsize = 10, labelpad = 2)
ax4.set_ylabel("Count", fontsize = 10, labelpad = 1)

# Add the number of variants with a frequency over the threshold
df_not_in_the_DMS = df_variants_above[df_variants_above["expected_phenotype"] == "Not in the DMS"]
nb_nonsynonymous = len(df_not_in_the_DMS[df_not_in_the_DMS["variant_type"] == "non-synonymous"])
nb_synonymous = len(df_not_in_the_DMS[df_not_in_the_DMS["variant_type"] == "synonymous"])

ax4.text(-3, 60, f"{nb_nonsynonymous} missense", fontsize = 10, color = "#8e6bb2")
ax4.text(-2.9, 42, f"{nb_synonymous} silent", fontsize = 10, color = "#ba91e3")

sns.despine(ax = ax4)

###########################   Ax5 : AUC boxplots for the validation of potentially resistant variants  ################################################

# Order the variants to have the controls on the left and then in position order
variants_order = ["WT", "pRS31N", "G464S", "A409T", "T411M", "R414S", "Y415C", "A430T", "K451I", "P512L", "P515T", "A516T"]

# Make the boxplot
sns.boxplot(data = df_auc, x = "strain", y = "norm_AUC", ax = ax5, hue = "Category", palette = ["#fc5959", "#5aaae4"], 
            hue_order = ["Resistant", "Sensitive"], order = variants_order, linecolor = "black", showfliers = False, linewidth = 0.75)


# Show every point
sns.stripplot(data = df_auc, x = "strain", y = "norm_AUC", ax = ax5, size = 3, hue = "Category", palette = ["#fc5959", "#5aaae4"], 
              hue_order = ["Resistant", "Sensitive"], order = variants_order, edgecolor = "black", linewidth=0.75, legend = False)

# Axis labels, ticks and legend parameters
ax5.set_xlabel("Variant", fontsize = 10, labelpad = 3)
ax5.set_ylabel("Normalized AUC", fontsize = 10, labelpad = 2)
ax5.set_ylim(-0.1, 2.5)

ax5.set_xticks(ax5.get_xticks(), ax5.get_xticklabels(), rotation=45, rotation_mode = "anchor", ha = "right")
ax5.tick_params(labelsize = 8, pad = 2)

# Add a grid
ax5.grid(axis = "y")

# Add a grey background for the controls
ax5.axhspan(-0.1, 2.5, xmax=0.25, color = "#ebebeb")

# Add the results of the t-test
for variant in variants_order:
    if variant != "WT":
        x = variant
        max_y = df_auc[df_auc["strain"] == variant]["norm_AUC"].max()
        pvalue = df_ttest.loc[variant]["adjusted_p_value"]

        if pvalue < 0.001 :
            ax5.annotate("***", xy = (x, max_y), xycoords = "data", 
                        xytext = (0, 1.5), textcoords = "offset points", ha = "center")

        elif pvalue < 0.01 :
            ax5.annotate("**", xy = (x, max_y), xycoords = "data", xytext = (0, 1.5), textcoords = "offset points", ha = "center")

        elif pvalue < 0.05 :
            ax5.annotate("*", xy = (x, max_y), xycoords = "data", xytext = (0, 1.5), textcoords = "offset points", ha = "center")

# Place legend in the top right of the plot
sns.move_legend(obj = ax5, loc = "upper right", title = None, fontsize = 8, handletextpad = 0.5)

#######################################################################################################################################

# Add the panel letters
ax1.annotate('a', (-0.16, 1.1), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")
ax2.annotate('b', (-0.05, 1.1), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")
ax3.annotate('c', (-0.22, 1.1), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")
ax5.annotate('d', (-0.13, 1.07), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")
ax6.annotate('e', (0.01, 1.1), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")

# Save the figure
fig.savefig("figure5cd.svg", format = "svg", dpi=300, bbox_inches="tight")

# %% ChimeraX scripts for figure 5E and to calculate the distance between a residue and the itraconazole/heme molecule

# Get a list of all the positions with at least one substitution leading to resistance to voriconazole identified in the DMS study
## Keep only the voriconazole condition
df_dms = df_dms[df_dms["Antifungal"] == "Voriconazoleco"]
## Get the list
df_beneficial = df_dms[df_dms["Category"] == "Beneficial"]
list_positions = list(df_beneficial["Position"].unique())
## Format for ChimeraX
formatted_positions = ""
for i in range(len(list_positions)):
    pos = str(list_positions[i])
    if i == 0:
        formatted_positions += pos
    else : 
        formatted_positions += ", " + pos

# Script to make the panel in ChimeraX
with open(f'Figure5E_ChimeraX-script.txt', 'w') as script:
    script.write(f"open 5v5z fromDatabase pdb format mmcif\n") # Fetch the structure from the PDB
    script.write(f"set bgcolor white\nlighting depthCue false\nlighting soft\nlighting shadows false\ngraphics silhouettes true\n")  # Aesthetic parameters
    script.write(f"label delete\n")    # Remove the missing residues label
    script.write(f"select /A\nsurface (#!1 & sel)\ncartoon hide (#!1 & sel)\nhide (#!1 & sel) target a\ncolor sel gainsboro\n") # Show only the surface
    script.write(f"select ::name='1YN'\nshow sel target ab\ncolor sel dim gray\n") # Color the itraconazole molecule dark gray
    script.write(f"select ::name='HEM'\nshow (#!1-2 & sel) target ab\ncolor (#!1-2 & sel) black\n") # Color the heme molecule black
    script.write(f"select :373-520\ncolor sel #9f9f9f\n") # Color the random mutagenesis region in pale gray
    script.write(f"select :{formatted_positions}\ncolor sel #c24343\n") # Color the known variants in red
    script.write(f"select ~sel & ##selected\ntransparency (#!1-2 & sel) 80\n") # Make the rest of the surface more transparent
    script.write(f"select : 414, 415, 451, 515\ntransparency (#!1-2 & sel) 0\ncolor sel #8e6bb2\n") # Color the positions with the validated mutations in purple
    script.write(f"select : 373-520\ntransparency (#!1-2 & sel) 0\n") # Make the random mutagenesis region not transparent
    script.write(f"lighting flat\nlighting shadows true intensity 0.5\n")  # Adjust the lighting
    script.write(f"view matrix camera -0.98676,0.15527,0.046912,-33.774,-0.11494,-0.46528,-0.87767,-205.27,-0.11445,-0.87144,0.47697,128.31\n") # Place the camera in the chosen position for the first view
    script.write(f"2dlabels text 'R414' size 30 xpos 0.394 ypos 0.084\n2dlabels arrow start 0.426,0.121 end 0.450,0.162 weight 0.1 headStyle pointer\n") # Add label for R414
    script.write(f"2dlabels text 'Y415' size 30 xpos 0.482 ypos 0.073\n2dlabels arrow start 0.496,0.114 end 0.483,0.149 weight 0.1 headStyle pointer\n") # Add label for Y415
    script.write(f"2dlabels text 'P515' size 30 xpos 0.586 ypos 0.618\n2dlabels arrow start 0.601,0.615 end 0.581,0.564 weight 0.1 headStyle pointer\n") # Add label for P515
    script.write(f"select clear\n")     # Clear select before saving the image
    script.write(f"save '/Users/MathildeRC/Library/Mobile Documents/com~apple~CloudDocs/Complementation paper/data_analysis/figures/Figure5E_view1.png' width 2000 height 1265 supersample 3 transparentBackground true\n")   # Save the first view image
    script.write(f"turn y 180\n")   # Go to the second view
    script.write(f"2dlabels delete\n")  # Remove labels from view 1
    script.write(f"2dlabels text 'K451' size 30 xpos 0.364 ypos 0.234\n2dlabels arrow start 0.398,0.271 end 0.440,0.321 weight 0.1 headStyle pointer\n") # Add label for K451
    script.write(f"select clear\n")     # Clear select before saving the image
    script.write(f"save '/Users/MathildeRC/Library/Mobile Documents/com~apple~CloudDocs/Complementation paper/data_analysis/figures/Figure5E_view2.png' width 2000 height 1265 supersample 3 transparentBackground true\n")   # Save the image for view 2
    script.write(f"save '/Users/MathildeRC/Library/Mobile Documents/com~apple~CloudDocs/Complementation paper/data_analysis/figures/Figure5E.cxs'")   # Save the ChimeraX session

# Script to check all the atoms-atoms distance between the residue and itraconazole or the heme 
# All possible distances for each residue-molecule pairs are listed in the log from shortest to longest. 
# The log is saved in a html file
with open(f'Figure5E_Distance-to-active-site_ChimeraX-script.txt', 'w') as script:
    script.write(f"open 5v5z fromDatabase pdb format mmcif\n") # Fetch the structure from the PDB
    script.write(f"contacts :414 restrict :1YN dist 100 log true\n")    # Calculate all atoms-atoms distances between the R414 residue and itraconazole (1YN). Max distance of 100 to make sure I get all the possible distances
    script.write(f"contacts :414 restrict :HEM dist 100 log true\n")    # Distance to the heme (HEM)
    script.write(f"contacts :415 restrict :1YN dist 100 log true\n")
    script.write(f"contacts :415 restrict :HEM dist 100 log true\n")
    script.write(f"contacts :451 restrict :1YN dist 100 log true\n")
    script.write(f"contacts :451 restrict :HEM dist 100 log true\n")
    script.write(f"contacts :515 restrict :1YN dist 100 log true\n")
    script.write(f"contacts :515 restrict :HEM dist 100 log true\n")
    script.write(f"log save '/Users/MathildeRC/Library/Mobile Documents/com~apple~CloudDocs/Complementation paperdata_analysis/figures/Figure5E_distance-to-active-site_ChimeraX.html'\n")    # Save the log

# %% Figure S5 - Correlation of the frequency before and after selection

# Create the colormap for the effect of the known variants (center the colormap around 0 while keeping the asymmetry between min-max values)
color_palette = sns.color_palette("vlag", as_cmap=True)
max_s, min_s = df_known_before_after["Selection coefficient"].max(), df_known_before_after["Selection coefficient"].min()
norm = mcolors.TwoSlopeNorm(vcenter = 0, vmin = min_s, vmax = max_s)

### Plot
# Initialize the plot
fig, ax = plt.subplots(figsize = (5, 5))

# Make the scatterplot with the unknown variants
sns.scatterplot(data = df_unknown_before_after, x = "median_log10_freq_before", y = "median_log10_freq_after", ax = ax,
                    hue = "expected_phenotype", palette = ["#c0c0c0"], edgecolor = "black", linewidths = 0.2, zorder = 2, marker = "s")

# Make the variants that were individually reconstructed and validated black
sns.scatterplot(data = df_validation, x = "median_log10_freq_before", y = "median_log10_freq_after", ax = ax,
                    hue = "expected_phenotype", palette = ["#c0c0c0"], edgecolor = "black", linewidths = 1, zorder = 3, marker = "s", legend = False)

# Add the known variants on top
sns.scatterplot(data = df_known_before_after, x = "median_log10_freq_before", y = "median_log10_freq_after", ax = ax,
                    hue = "Selection coefficient", palette = color_palette, hue_norm = norm,
                    edgecolor = "black", linewidths = 0.2, legend = False, zorder = 4)

# Add a horizontal line at the threshold value used to select variants for validation
ax.axhline(y = df_thresholds.loc[0.99]["threshold"], color = "black", linestyle = "--", linewidth = 0.7, zorder = 1)

# Axis limits, tick parameters and labels
ax.set_xlim([-5.6, -1.5])
ax.set_xticks([-5, -4, -3, -2])
ax.set_ylim([-5.3, -1.5])
ax.set_yticks([-5, -4, -3, -2])
ax.tick_params("both", labelsize = 10)
ax.set_xlabel("Frequency before selection (log10)", fontsize = 12)
ax.set_ylabel("Frequency after selection (log10)", fontsize = 12)

# Colorbar 
## Get the colormap values in the right format
sm = plt.cm.ScalarMappable(cmap = color_palette, norm = norm)
sm.set_array([])
## Create an inset ax to place the color bar
ax_inset = inset_locator.inset_axes(ax, width = "7.5%", height = "65%", loc = "upper right", bbox_to_anchor = (0.1, -0.1, 1, 1),
                                    bbox_transform = ax.transAxes, borderpad = 0)
cbar = plt.colorbar(sm, cax = ax_inset)    # Add the colorbar to the plot
cbar.set_ticks([-0.6, -0.4, -0.2, 0, 0.2])
cbar.ax.set_yscale("linear")
cbar.set_label("DMS selection coefficient", fontsize = 10, labelpad = 7)

# Legend
sns.move_legend(obj = ax, loc = "lower right", title = None, frameon = False, bbox_to_anchor = (1.35, 0.1), 
                handletextpad = 0,  borderpad = 0.3, alignment = "left", fontsize = 10)

# Adding Spearman's rho
rho, pvalue = stats.spearmanr(df_before_after["median_log10_freq_before"], df_before_after["median_log10_freq_after"])
text_rho = "ρ = " + str(round(rho,2))
ax.text(-2.6, -4.8, text_rho , fontsize = 12)

# Add annotations for the variants validated to be resistant
labels = {}
for variant in validation_variants:
    x = df_before_after.loc[variant]["median_log10_freq_before"]
    y = df_before_after.loc[variant]["median_log10_freq_after"]
    labels[variant] = [x, y]

ax.annotate("R414S", xy = (labels["R414S"][0], labels["R414S"][1]), xycoords = "data",
            xytext = (10, -10), textcoords = "offset points", ha = "left", va = "center", 
            arrowprops = dict(arrowstyle = "-", shrinkA = 0, shrinkB = 4))
ax.annotate("Y415C", xy = (labels["Y415C"][0], labels["Y415C"][1]), xycoords = "data",
            xytext = (-5, 20), textcoords = "offset points", ha = "right", va = "center", 
            arrowprops = dict(arrowstyle = "-", shrinkA = 0, shrinkB = 4))
ax.annotate("K451I", xy = (labels["K451I"][0], labels["K451I"][1]), xycoords = "data",
            xytext = (0, 20), textcoords = "offset points", ha = "right", va = "center", 
            arrowprops = dict(arrowstyle = "-", shrinkA = 0, shrinkB = 4))
ax.annotate("P515T", xy = (labels["P515T"][0], labels["P515T"][1]), xycoords = "data",
            xytext = (-10, 0), textcoords = "offset points", ha = "right", va = "center", 
            arrowprops = dict(arrowstyle = "-", shrinkA = 0, shrinkB = 4))

# Save the figure
fig_name = "FigureS5"
plt.savefig(f"{fig_name}.png", format='png', dpi=300, bbox_inches="tight")
plt.savefig(f"{fig_name}.svg", format='svg', dpi=300, bbox_inches="tight")


# %%
