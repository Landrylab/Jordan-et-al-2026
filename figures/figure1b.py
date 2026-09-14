
# This script generates heatmaps of pairwise comparison of protein sequence identity for orthologs of Erg3, Erg6 and Erg11 (figure 1B)
# Sequence identity matrices are from the pid_calculation.Rmd script


#%% Import packages

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Figure parameters
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary

#%% Load the sequence identity matrices

# Load the sequence identity matrices (output from the pid_calculation.Rmd script)
df_erg3 = pd.read_csv("../results/erg3_seqid.csv", index_col = 0)
df_erg6 = pd.read_csv("../results/erg6_seqid.csv", index_col = 0)
df_erg11 = pd.read_csv("../results/erg11_seqid.csv", index_col = 0)

# Reorder the dataframes to match the cladogram (figure 1a)
order_erg3 = ["CtroERG3", "CalbERG3", "CparERG3", "CaurERG3", "ScerERG3", "NglaERG3"]
df_erg3 = df_erg3.reindex(index = order_erg3, columns = order_erg3)

order_erg6 = ["CtroERG6", "CalbERG6", "CparERG6", "CaurERG6", "ScerERG6", "NglaERG6"]
df_erg6 = df_erg6.reindex(index = order_erg6, columns = order_erg6)

order_erg11 = ["CtroERG11", "CalbERG11", "CparERG11", "CaurERG11", "ScerERG11", "NglaERG11", "AfumCYP51A", "AfumCYP51B"]
df_erg11 = df_erg11.reindex(index = order_erg11, columns = order_erg11)

# %% Make figure 1

labels = ["C. tropicalis", "C. albicans", "C. parapsilosis", "C. auris", "S. cerevisiae", "N. glabratus"]
labels_erg11 = ["C. tropicalis", "C. albicans", "C. parapsilosis", "C. auris", "S. cerevisiae", "N. glabratus", "A. fumigatus $\mathrm{Cyp51A}$", "A. fumigatus $\mathrm{Cyp51B}$"]

# Draw the figure and axes
## Create figure and gridspec
fig = plt.figure(constrained_layout = False, figsize = (6.5, 5.5))
gs = fig.add_gridspec(nrows = 55, ncols = 65, left=0, right=1, bottom=0, top=1)

## Draw ax for the cladogram 
ax1 = fig.add_subplot(gs[0:22, 0:20])

## Draw ax for the plasmid and strain schematics
ax2 = fig.add_subplot(gs[24:56, 0:66])

## Draw ax for the sequence identity heatmaps
ax3 = fig.add_subplot(gs[0:16, 26:37])  # Erg3
ax4 = fig.add_subplot(gs[0:16, 38:49])  # Erg6
ax5 = fig.add_subplot(gs[0:20, 50:65])  # Erg11
ax6 = fig.add_subplot(gs[18:20, 32:44]) # Colorbar

# Make the ax1 and ax3 blank spaces to fill in Inkscape
ax1.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax1.spines[:].set_visible(False)
ax1.annotate('a',(0, 0.95), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")

ax2.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax2.spines[:].set_visible(False)
ax2.annotate('c', (0, 0.95), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")

# Fill axes with the sequence identity heatmaps
## Erg3
sns.heatmap(df_erg3, ax = ax3, cmap = "YlGnBu_r", vmin = 0, vmax = 100, square = True, linewidth = 0.5, linecolor = "black",
            clip_on = False, xticklabels = labels, yticklabels = False,
            cbar_ax = ax6, cbar_kws = {"orientation" : "horizontal", "ticks" : [0, 25, 50, 75, 100]})
ax3.tick_params(axis = "both", length = 0, pad = 1)
ax3.set_xticks(ax3.get_xticks(), ax3.get_xticklabels(), rotation=45, ha = "right", rotation_mode = "anchor", fontstyle = "italic", fontsize = 6)
ax3.set_title("Erg3", fontsize = 8, fontweight='bold')

## Erg6
sns.heatmap(df_erg6, ax = ax4, cmap = "YlGnBu_r", vmin = 0, vmax = 100, square = True, linewidth = 0.5, linecolor = "black",
            clip_on = False, xticklabels = labels, yticklabels = False, cbar = False)
ax4.tick_params(axis = "both", length = 0, pad = 1)
ax4.set_xticks(ax4.get_xticks(), ax4.get_xticklabels(), rotation=45, ha = "right", rotation_mode = "anchor", fontstyle = "italic", fontsize = 6)
ax4.set_title("Erg6", fontsize = 8, fontweight='bold')

## Erg11
sns.heatmap(df_erg11, ax = ax5, cmap = "YlGnBu_r", vmin = 0, vmax = 100, square = True, linewidth = 0.5, linecolor = "black",
            clip_on = False, xticklabels = labels_erg11, yticklabels = labels_erg11, cbar = False)
ax5.tick_params(axis = "x", length = 0, pad = 1)
ax5.tick_params(axis = "y", length = 0, pad = 2)
ax5.set_xticks(ax5.get_xticks(), ax5.get_xticklabels(), rotation=45, ha = "right", rotation_mode = "anchor", fontstyle = "italic", fontsize = 6)
ax5.yaxis.tick_right()
ax5.set_yticks(ax5.get_yticks(), ax5.get_yticklabels(), rotation=0, fontstyle = "italic", fontsize = 6)
ax5.set_title("Erg11", fontsize = 8, fontweight='bold')

## Colorbar parameters
ax6.set_xticks(ax6.get_xticks(), ax6.get_xticklabels(), fontsize = 6)
ax6.tick_params(axis = "x", length = 1, pad = 1)
ax6.set_xlabel("Identity (%)", fontsize = 7, labelpad = 1)
ax6.spines["outline"].set(visible = True, edgecolor = "black", linewidth = 0.5)

## Add panel label
ax3.annotate('b', (-0.2, 1.1), xycoords='axes fraction', annotation_clip=False, fontsize=12, fontweight='bold', fontname = "Arial")

# Save the figure
fig.savefig("figure1b.svg", format = "svg", dpi=300, bbox_inches="tight")

# %%
