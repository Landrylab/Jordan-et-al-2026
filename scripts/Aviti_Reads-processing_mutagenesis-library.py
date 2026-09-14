#%%

# This script takes demultiplexed data containing RC primers and prepares the reads for variant calling in part 2. Quality control along the way with graphs and dataframes to show the read counts

### Pipeline

#1. Quality control with fastQC on .fastq.gz files
#       - Figure : number of reads obtained
#2. From demultiplexed .fastq.gz files, pandaseq will merge the forward and reverse reads and remove the remaining primers
#       - Figure : number of reads retrieved after merging
#3. From .fasta files, vsearch will aggregate the reads
#       - Figure : number of sequences with one, two and more than two reads
#       - Figure : number of reads associated with sequences with one, two and more than two reads
#5. From .fasta files, needle will align the reads to the wild-type (wt) sequence
#       - Figure : number of reads lost at each filtering step
#       - Figure : final number of reads that will be used for downstream analyses

# Outputs the aligned reads (Needle files) and a dataframe with the total reads retrieved after processing
# Files generated after each processing step prior to alignment with Needle are saved in an intermediate data folder.

#%% Set up

# Import the necessary packages
import subprocess
import openpyxl       # Necessary to read excel files
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Version of the script (date)
version = "2026-04-09"

# Name of the experiment (to save figures)
experiment = "mutagenesis-library"

# Paths to the different config files
sample_sheet_path = "../data/mutagenesis_library/sample-sheet.xlsx"
wt_seq = "../data/mutagenesis_library/F4_wt_seq.fasta"        # WT sequence after trimming

# Path to the folder containing the sequencing reads
folder_reads = "../data/mutagenesis_library/"

# Paths to the folders where to save the intermediate dataframes and quality control figures
if not os.path.exists("../results/intermediate_data_mutagenesis-library/plots/"):
    os.makedirs("../results/intermediate_data_mutagenesis-library/plots/")

if not os.path.exists("../results/intermediate_data_mutagenesis-library/dataframes/fastqc/"):
    os.makedirs("../results/intermediate_data_mutagenesis-library/dataframes/fastqc/")

if not os.path.exists("../results/intermediate_data_mutagenesis-library/dataframes/merged_reads/"):
    os.makedirs("../results/intermediate_data_mutagenesis-library/dataframes/merged_reads/")
    
if not os.path.exists("../results/intermediate_data_mutagenesis-library/dataframes/aggregated_reads/"):
     os.makedirs("../results/intermediate_data_mutagenesis-library/dataframes/aggregated_reads/")

if not os.path.exists("../results/intermediate_data_mutagenesis-library/dataframes/read_counts/"):
     os.makedirs("../results/intermediate_data_mutagenesis-library/dataframes/read_counts/")

folder_qc_figures = "../results/intermediate_data_mutagenesis-library/plots/"
folder_fastqc = "../results/intermediate_data_mutagenesis-library/dataframes/fastqc/"
folder_merged = "../results/intermediate_data_mutagenesis-library/dataframes/merged_reads/"
folder_agg = "../results/intermediate_data_mutagenesis-library/dataframes/aggregated_reads/"
folder_read_counts = "../results/intermediate_data_mutagenesis-library/dataframes/read_counts/"

# Path to the folder where to save the aligned reads and dataframe with the total reads processed
if not os.path.exists("../results/"):
     os.makedirs("../results/")

folder_aligned = "../results/"
folder_processed = "../results/"

# Pandaseq parameters
l = "430"         # Maximum length of the assembled sequence (length before trimming)
L = "530"         # Minimum length of the assembled sequence (length before trimming)
o = "95"          # Maximum length of the overlap between the R1 and R2 reads (length before trimming)
O = "145"         # Minimum length of the overlap between the R1 and R2 reads

# %%

# Load the sample sheet
info = pd.read_excel(sample_sheet_path, header=0, index_col=0, keep_default_na=False)    # keep_default_na = False -> so that the "None" condition is not convert to NaN
print(info.head())

# %% FastQC

def quality_fastqc(input_fastq_file):
    subprocess.check_output('fastqc '+input_fastq_file +" --outdir "+folder_fastqc, shell=True)
    return

# Run FastQC on all sequencing files (files with a fastqc.gz extension)
ext = ('fastq.gz')
 
for files in os.listdir(folder_reads):
    if files.endswith(ext):
        quality_fastqc(folder_reads + files)

##################  Make sure that the sequencing quality is good by checking the FastQC reports before continuing ############################################
        
# %% Count the total number of sequencing reads received

before_merge_depth_dict = {}
for Sample in list(info.index):
    fullpath = folder_reads +info.loc[Sample]['reads_file_For']
    call ="gunzip -c " + fullpath + " |wc -l"     # Decompress the .gz files
    result = subprocess.check_output(call, shell=True)  
    numseqs = int(result)/4.0  # Divide by 4 because there are 4 lines in the file per read
    before_merge_depth_dict[str(info.loc[Sample]["name"])] = numseqs

# Put number of reads for each sample in a dataframe and save it
reads_before = pd.DataFrame.from_dict(before_merge_depth_dict, orient="index", columns = ["reads_before"])
name = folder_read_counts + "Reads_received_" + experiment + "_" + version  
reads_before.to_csv(f"{name}.csv")

# Check the number of reads received
reads_before

# Sanity check : check that the number of reads received makes sense by ctrl F the run ID in the raw reads
#                file for some samples (should be equal to the number of reads). Can also check if it matches
#                what is written in the FastQC report

# %% Visualize the number of reads received

fig, ax = plt.subplots(figsize=(4, 8))

ax.bar(info["name"], reads_before["reads_before"], color = "#9bd179", edgecolor = "black", linewidth = 0.5)
#plt.axhline(y= expected_read_count, color = "black", linestyle = "--", linewidth = 1)
ax.set_ylim(0, 600000)  
ax.set_title("Reads received", fontsize=18)   
ax.set_xlabel('Sample', fontsize=14)
ax.set_ylabel('Read count (thousand reads)', fontsize=14)
ax.tick_params(axis = "both", labelsize=12)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=30, ha='right')

# Format the y axis ticks
def y_tick_formatter(x, pos):
    return f"{x/1000}"
ax.yaxis.set_major_formatter(FuncFormatter(y_tick_formatter))

# Add annotations : number of reads 
for p in ax.patches:
    ax.annotate(f"{(p.get_height()/1000):.0f}k", (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='center', xytext =(0, 7), textcoords='offset points')

#Save the figure
name = folder_qc_figures + "Reads_before_merging_" + experiment + "_" + version 
plt.savefig(f"{name}.png", format='png', dpi=300,bbox_inches="tight")

#%% Merge and trim reads with Pandaseq

def merge_trim_reads(Sample):
    # This function generates and calls a Pandaseq call to merge R1 and R2
    # and trim the primers for each demultiplexed samples

    # Path to the sequencing reads
    filepath_for = folder_reads + info.loc[Sample]['reads_file_For']
    filepath_rev = folder_reads + info.loc[Sample]['reads_file_Rev']

    # Path to the folder where to save the merged reads and pandaseq stats
    filepath_out = folder_merged + str(info.loc[Sample]['name'])+'.fasta'
    stats_out = folder_merged + str(info.loc[Sample]['name']) + "_stats.txt"

    # Sequence of the primers to trim
    for_seq = info.loc[Sample]["RC_for_seq"]
    rev_seq = info.loc[Sample]["RC_rev_seq"]

    # Other parameters (same as Bédard et al., 2024): 
    # -B : allow input to lack a barcode
    # -k : number of sequence locations per kmer (increased from the default value of 2)
    # -N : eliminate all sequences with Ns
    # -t 0.5 : score threshold
    # -T 6 : use 6 threads
    # -w filepath_out : write output to filepath_out
    # -g stats_out : write stats to stats_out
    panda_seq_call = 'pandaseq -f '+filepath_for+' -r '+filepath_rev+' -p '+for_seq+' -q '+rev_seq+' -L '+L+' -l '+l+' -O '+O+' -o '+o+' -k 4 -B -N -t 0.5 -T 6 -w '+ filepath_out + ' -g ' + stats_out
    print(panda_seq_call)
    subprocess.check_output(panda_seq_call, shell=True)

    return filepath_out

# Running Pandaseq for each sample
for Sample in list(info.index):
    merge_trim_reads(Sample)


########################  Make sure that the merging looks ok before continuing  ######################################
# Look at the _stats.txt file : 
# - Check that the best overlap ("BESTOLP") is what it should be (for the length before trimming)
# - Check that most reads are OK, if not check in which category the reads fall and adjust the parameters in consequence
# Align a couple of the merged sequences from the .fasta file to the wt sequence on Benchling to make sure that the length is ok.

#%% Count the number of reads after merging

def count_merged(Sample):
    
    # Define filepath based on sample name
    filepath = folder_merged + str(info.loc[Sample]['name']) + '.fasta'
    print(filepath)

    # Counter for the number of reads
    after_merge_depth = 0
    with open(filepath, 'r') as source:
        for line in source:
            if line.startswith('>'):    # Increment by one for each header
                after_merge_depth += 1

    after_merge_depth_dict[str(info.loc[Sample]["name"])] = after_merge_depth
    return after_merge_depth

# Empty containers to hold the post-merge number of reads for each sample
after_merge_depth_dict = {}

# Count merged reads for each sample
for Sample in list(info.index):
    count_merged(Sample)

# Put post-merge number of reads for each sample in a dataframe 
reads_after_merge = pd.DataFrame.from_dict(after_merge_depth_dict, orient="index", columns = ["reads_after_merge"])
reads_after_merge

# Merge with the reads_before dataframe
reads_before = pd.read_csv(folder_read_counts + "Reads_received_" + experiment + "_" + version + ".csv", index_col = 0)
df_merge = pd.merge(reads_before, reads_after_merge, left_index=True, right_index=True)
df_merge["lost_merge"] = df_merge["reads_before"] - df_merge["reads_after_merge"]
df_merge["percent_lost_merge"] = (df_merge["lost_merge"]/df_merge["reads_before"])*100

# Save the dataframe
name = folder_read_counts + "Reads_after_merging_" + experiment + "_" + version  
df_merge.to_csv(f"{name}.csv")

#%% Visualize the number of reads lost during merging

fig, ax = plt.subplots(figsize=(4, 8))

ax.bar(x=info["name"], height=reads_before["reads_before"], color = "#d45366", edgecolor = "black", label = "Lost during merging")
ax.bar(x=info["name"], height=reads_after_merge["reads_after_merge"], color = "#9bd179", edgecolor = "black", label = "Good")
#plt.axhline(y=expected_read_count, color = "black", linestyle = "--", linewidth = 1, label = "Expected read count")
ax.set_ylim(0, 600000)
ax.set_title("Quality control - Merging", fontsize=18)
ax.set_xlabel('Sample', fontsize=14)
ax.set_ylabel('Read count (thousands reads)', fontsize=14)
ax.tick_params(axis = "both", labelsize=12)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=30, ha='right')

# Format the y-axis ticks
def y_tick_formatter(x, pos):
    return f"{x/1000}"
ax.yaxis.set_major_formatter(FuncFormatter(y_tick_formatter))

#Save the figure
name = folder_qc_figures + "Reads_after_merging_" + experiment + "_" + version
plt.savefig(f"{name}.png", format='png', dpi=300, bbox_inches="tight")

#%% Aggregate reads with vsearch

def aggregate(Sample):
    # This function generates and runs a vsearch call to aggregate the merged reads.
       
    filepath_merged = folder_merged +str(info.loc[Sample]['name'])+'.fasta'  
    aggregate_path = folder_agg +str(info.loc[Sample]['name'])+ '.fasta'          
    
    vsearch_call = 'vsearch --derep_fulllength '+ filepath_merged +' --relabel seq --output '+aggregate_path+' --sizeout'
    subprocess.check_output(vsearch_call, shell=True)
    # Generates and runs a vsearch call with the following parameters:
    #    --derep_fulllength trim_path    vsearch program to be used and input. derep_fulllength aggregates perfectly identical 
    #                                    sequences in the input fasta, dramatically reducing the number of sequences that then 
    #                                    need to be aligned
    #    --relabel seq                   changes header so that all sequence names begin by seq
    #    --output                        path where output will be written
    #    --sizeout                       append the size of the sequence cluster to the fasta header
    
for Sample in list(info.index):
    aggregate(Sample)

#%%  Count the unique sequences and singletons in aggregated reads

# Containers to hold the number of unique sequences, singletons and sequences with 2 reads
after_agg_unique_dict = {}
after_agg_single_dict = {}
after_agg_two_dict = {}

for Sample in list(info.index):
    file_to_count=folder_agg + str(info.loc[Sample]['name'])+'.fasta'
    
    # Bash line to count all unique reads per library
    cmdline_unique=f'grep -c ">" {file_to_count}'
    unique_seq_number=subprocess.getoutput(cmdline_unique)
    # Write number of unique sequences to dict
    after_agg_unique_dict[str(info.loc[Sample]["name"])] = int(unique_seq_number)
    
    # Bash line to count all sequences with one read (singletons)
    cmdline_singleton = f'grep -c "size=1;" {file_to_count}'
    singleton_number=subprocess.check_output(cmdline_singleton, shell=True)
    #Write number of singleton sequences to dict
    after_agg_single_dict[str(info.loc[Sample]["name"])] = int(singleton_number)
    
    # Bash line to count all sequences with 2 reads (almost singletons)
    cmdline_two = f'grep -c "size=2;" {file_to_count}'
    two_number=subprocess.check_output(cmdline_two, shell=True)
    #Write number of singleton sequences to dict
    after_agg_two_dict[str(info.loc[Sample]["name"])] = int(two_number)

# Delete temporary variables
del(file_to_count, cmdline_unique, unique_seq_number, cmdline_singleton, singleton_number, Sample )

# Put number of unique sequences/singletons/sequences with two reads for each sample in a dataframe 
seq_agg = pd.DataFrame.from_dict(after_agg_unique_dict, orient="index", columns = ["nbr_seq_unique"])
seq_agg_single = pd.DataFrame.from_dict(after_agg_single_dict, orient="index", columns = ["nbr_seq_single"])
seq_agg_two = pd.DataFrame.from_dict(after_agg_two_dict, orient="index", columns = ["nbr_seq_two"])

# Add to the df_merge dataframe
df_merge = pd.read_csv(folder_read_counts + "Reads_after_merging_" + experiment + "_" + version + ".csv", index_col = 0)
df_agg = pd.merge(df_merge, seq_agg, left_index= True, right_index= True)
df_agg = pd.merge(df_agg, seq_agg_single, left_index= True, right_index= True)
df_agg = pd.merge(df_agg, seq_agg_two, left_index= True, right_index= True)

# Number of reads for the sequences with 2 reads
df_agg["reads_seq_two"] = df_agg["nbr_seq_two"]*2

# Calculate the percentage of merged reads that the singletons and sequences with only two reads represent
df_agg["percent_merged_single"] = (df_agg["nbr_seq_single"]/df_agg["reads_after_merge"])*100
df_agg["percent_merged_two"] = (df_agg["reads_seq_two"]/df_agg["reads_after_merge"])*100

# Save the dataframe
name = folder_read_counts + "Reads_after_aggregating_" + experiment+"_" + version  
df_agg.to_csv(f"{name}.csv")

# %% Visualize the number of unique sequences that are singletons and sequences with two reads

fig, ax = plt.subplots(figsize=(4, 8))

ax.bar(info["name"], seq_agg["nbr_seq_unique"], color = "#cdcdcd", edgecolor = "black", label = "Singletons")
ax.bar(info["name"], (seq_agg["nbr_seq_unique"] -  seq_agg_single["nbr_seq_single"]), color = "#a6cee3", edgecolor = "black", label = "Sequences with 2 reads")
ax.bar(info["name"], (seq_agg["nbr_seq_unique"] -  seq_agg_single["nbr_seq_single"] - seq_agg_two["nbr_seq_two"]), color = "#1f78b4", edgecolor = "black", label = "Sequences with > 2 reads")

ax.set_title("Unique sequences", fontsize=18)
ax.set_xlabel('Sample', fontsize=14)
ax.set_ylabel('Count', fontsize=14) 
ax.tick_params(axis = "both", labelsize=12)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=30, ha='right')
ax.legend(loc = "upper left", frameon = False, fontsize = 14)

name = folder_qc_figures + "Unique_sequences_" + experiment + "_" + version 
plt.savefig(f"{name}.png", format='png', dpi=300,bbox_inches="tight")

#%% Visualize the proportion of the reads representing the sequences with only 1 or 2 reads

fig, ax = plt.subplots(figsize=(4, 8))

ax.bar(info["name"], df_agg["reads_after_merge"], color = "#cdcdcd", edgecolor = "black", label = "Singletons")
ax.bar(info["name"], (df_agg["reads_after_merge"] - df_agg["nbr_seq_single"]), color = "#a6cee3", edgecolor = "black", label = "Sequences with 2 reads")
ax.bar(info["name"], (df_agg["reads_after_merge"] - df_agg["nbr_seq_single"] - df_agg["reads_seq_two"]), color = "#1f78b4", edgecolor = "black", label = "Sequences with > 2 reads")
#plt.axhline(y=expected_read_count, color = "black", linestyle = "--", linewidth = 1, label = "Expected read count")
ax.set_ylim(0, 600000)     
ax.set_title("Number of singleton reads", fontsize=18)
ax.set_xlabel('Sample', fontsize=14)
ax.set_ylabel('Read count (thousands reads)', fontsize=14)     
ax.tick_params(axis = "both", labelsize=12)
ax.legend(loc = "upper left", frameon = False, fontsize = 14)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=30, ha='right')

# Format the y-axis ticks
def y_tick_formatter(x, pos):
    return f"{x/1000}"
ax.yaxis.set_major_formatter(FuncFormatter(y_tick_formatter))

name = folder_qc_figures + "Singleton_prop_" + experiment + "_" + version  
plt.savefig(f"{name}.png", format='png', dpi=300,bbox_inches="tight")

#%% Visualize the number of reads that passed all the filters and are used for the alignment with Needle

fig, ax = plt.subplots(figsize=(4, 8))

ax.bar(info["name"], df_agg["reads_before"], color = "#d45366", edgecolor = "black", label = "Merging")
ax.bar(info["name"], df_agg["reads_after_merge"], color = "#ece16f", edgecolor = "black", label = "Singletons")
ax.bar(info["name"], df_agg["reads_after_merge"] - df_agg["nbr_seq_single"], color = "#9bd179", edgecolor = "black", label = "Aligned reads")
#plt.axhline(y=expected_read_count, color = "black", linestyle = "--", linewidth = 1, label = "Expected read count")
ax.set_ylim(0, 600000) 
ax.set_title("Quality control - All filtering steps", fontsize=18)
ax.set_xlabel('Sample', fontsize=14)
ax.set_ylabel('Read count (thousands reads)', fontsize=14)   
ax.tick_params(axis = "both", labelsize=12)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=30, ha='right')
ax.legend(loc = "upper right", frameon = False, fontsize = 14)

# Format the y-axis ticks
def y_tick_formatter(x, pos):
    return f"{x/1000}"
ax.yaxis.set_major_formatter(FuncFormatter(y_tick_formatter))

name = folder_qc_figures + "All_filtering_" + experiment + "_" + version
plt.savefig(f"{name}.png", format='png', dpi=300,bbox_inches="tight")

#%% Align to WT sequence with Needle

# Remove singletons before alignment
def filter_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        keep_entry = False
        entry = []
        for line in infile:
            if line.startswith('>'):
                # Extract size value using regex
                match = re.search(r'size=(\d+)', line)
                size = int(match.group(1)) if match else 0               
                if entry and keep_entry:
                    outfile.writelines(entry)
                keep_entry = size > 1  # Keep sequences where size > 1
                entry = [line] if keep_entry else []
            else:
                if keep_entry:
                    entry.append(line)   
        if entry and keep_entry:
            outfile.writelines(entry)
            
for Sample in list(info.index):
    # Define input and output file paths
    input_fasta = folder_agg + str(info.loc[Sample]['name']) + '.fasta'
    output_fasta = folder_agg + str(info.loc[Sample]['name']) + '_nosingle.fasta'
    # Run the filtering function
    filter_fasta(input_fasta, output_fasta)
    
# Define and run the needle call
for Sample in list(info.index):
    
    # Define file paths
    file_in= folder_agg +str(info.loc[Sample]['name'])+ '_nosingle.fasta'
    file_out= folder_aligned +str(info.loc[Sample]['name'])+ '_mutagenesis-library.needle'
    
    # Write command line to run the needle alignment tool
    needle_call = 'needle -auto -gapopen 100 -asequence '+ wt_seq + ' -bsequence '+ file_in +' -aformat3 markx10 -outfile '+file_out   
    subprocess.check_output(needle_call, shell = True)

    # generate and run a needle call with the following parameters:
    #    -auto                        do not run in interactive mode
    #    -gapopen 100                 increase penalty for gap in the alignment. This helps
    #                                 with instances where there are no indels in the read but
    #                                 needle still opens a gap 
    #    -asequence ref_fasta_path    reference sequence fasta file
    #    -bsequence filepath          query sequences in fasta format
    #    -aformat3 markx10            choose the ouput format
    #    -ourfile needle_out          define output file path


#%% Count the number of reads that were aligned with Needle

def count_aligned(Sample):
    
    # Define filepath based on sample name
    filepath = folder_aligned + str(info.loc[Sample]['name']) + '_mutagenesis-library.needle'
    print(filepath)

    # Counter for the number of reads
    align_depth = 0
    with open(filepath, 'r') as source:
        for line in source:
            if line.startswith('>>>'):    # Detect header
                align_name = line         # Get alignment name (contains the read count of the query sequence)
                match = re.search(r'size=(\d+)', align_name)
                size = int(match.group(1)) if match else 0
                align_depth += size
                
    align_depth_dict[str(info.loc[Sample]["name"])] = align_depth
    return align_depth_dict

# Empty containers to hold the number of reads that were aligned for each sample
align_depth_dict = {}

# Count merged reads for each sample
for Sample in list(info.index):
    count_aligned(Sample)

# Put post-merge number of reads for each sample in a dataframe 
reads_aligned = pd.DataFrame.from_dict(align_depth_dict, orient="index", columns = ["reads_aligned"])
reads_aligned

# Merge with the dataframe with the number of reads after each processing step
df_agg = pd.read_csv(folder_read_counts + "Reads_after_aggregating_" + experiment+"_" + version + ".csv", index_col = 0)
df_final = df_agg.merge(reads_aligned, left_index=True, right_index=True)

# Calculate the proportion of the total reads received that passed all the processing steps
df_final["percent_reads_aligned"] = (df_final["reads_aligned"]/df_final["reads_before"])*100

# Save the dataframe
name =  folder_processed + "Total_reads_processed_" + experiment + "_" + version  
df_final.to_csv(f"{name}.csv")

#### Sanity check : the number of reads aligned should be equal to the number of reads after merging - number of singleton reads

#%% Visualize the final read count (total of reads aligned)

fig, ax = plt.subplots(figsize=(4, 8))

ax.bar(info["name"], df_final["reads_aligned"], color = "#9bd179", edgecolor = "black", linewidth = 0.5)
#plt.axhline(y=expected_read_count, color = "black", linestyle = "--", linewidth = 1, label = "Expected read count")
ax.set_ylim(0, 600000)        
ax.set_title("Reads processed", fontsize=18)
ax.set_xlabel('Sample', fontsize=14)
ax.set_ylabel('Read count (thousands reads)', fontsize=14) 
ax.tick_params(axis = "both", labelsize=12)
ax.legend(fontsize=14)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=30, ha='right')

# Format the y-axis ticks
def y_tick_formatter(x, pos):
    return f"{x/1000}"
ax.yaxis.set_major_formatter(FuncFormatter(y_tick_formatter))

## Add annotations : total of reads aligned 
for p in ax.patches:
    ax.annotate(f"{(p.get_height()/1000):.1f}k", (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='center', xytext =(0, 7), textcoords='offset points')

name = folder_qc_figures + "Total_reads_processed_" + experiment + "_" + version 
plt.savefig(f"{name}.png", format='png', dpi=300,bbox_inches="tight")

# %%
