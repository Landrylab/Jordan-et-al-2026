
## From needle alignments, calls mutations in a sequence.
## Inputs :
##   - Needle files (output from the Aviti_Reads_processing script)
##   - Sample sheet with information for each sample (same as the one used in the Aviti_Reads_processing script)
##   - Total_reads_processed csv file (output from the Aviti_Reads_processing script)
## Outputs : 
##   - Dataframe with all the variants for all samples
##   - Other dataframes and figures generated along the way are saved in an intermediate data folder

#%% Import packages

import os,sys,re
import openpyxl             # Necessary to open .xlsx files
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'       # To ensure proper conversion of text labels for downstream editing in Inkscape if necessary

#%% Set up - Ideally this is the only cell that has to be modified

# Version of the script (date)
version = "2026-04-15"

# Name of the experiment (to save figures)
experiment = "mutagenesis-library"

# Paths to the config files
sample_sheet_path = "../data/mutagenesis_library/sample-sheet.xlsx"
wt_seq_path = "../data/mutagenesis_library/F4_wt_seq.fasta"        # WT sequence after trimming

# Path to the needle alignments files
folder_aligned = "../results/"

# Path to the total number of reads processed for each sample (output of the Aviti_reads_processing_solid-selection.py script)
read_counts_path = "../results/Total_reads_processed_mutagenesis-library_2026-04-09.csv"

# Path to the folder where to save the final dataframes
folder_processed_data = "../results/"

# Path to the folder where to save intermediate figures and dataframes
folder_figures = "../results/intermediate_data_mutagenesis-library/plots/"
folder_dataframes = "../results/intermediate_data_mutagenesis-library/dataframes/"

# Codons_mutated : all codons in this case since it is random mutagenesis (can be changed if the library was constructed with NNK oligos for example)
codons_mutated = ["GCT","GCC","TGC","TGT","GAT","GAC","GAG","GAA","TTC","TTT",
                  "GGT","GGA","CAC","CAT","ATC","ATT","AAG","AAA","TTG","TTA",
                  "ATG","AAC","AAT","CCT","CCA","CAG","CAA","AGA","CGT","TCT",
                  "TCC","ACC","ACT","GTT","GTC","TGG","TAC","TAT","TAA","GTA",
                  "GCA","GTG","ATA","GCG","CTC","CTA","CGG","TAG","TGA","CCC",
                  "GGG","TCG","AGG","CGA","CGC","AGT","CTG","ACG","TCA","AGC",
                  "GGC","CTT","ACA","CCG"]

# Genetic code (standard)
genetic_code = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }

#%% Load the input dataframes (sample sheet + read counts)

# Sample sheet
info = pd.read_excel(sample_sheet_path, header=0, index_col=0)

# Read counts 
df_read_counts = pd.read_csv(read_counts_path, header=0, index_col=0)

############################ Make sure that both are ok before continuing   ##############################

#%%
############################################################################################################
##  Next we define a bunch of functions that will all be run at the end of the script                     ##  
############################################################################################################

##  get_wt_seq() : retrieves the wild-type sequence for the fragment corresponding to the input sample,
##                       and prints the sequence to be able to check that it is ok

##  find_mutations() : finds the differences between the aligned query and reference sequences.    
##      - Calls parse_needle_output(path) to extract the reference and query sequences from the needle alignments.
##      - **FILTERING STEP** : removes sequences with mutations in nucleotides not covered by the random mutagenesis target
##                             region and sequences that cover less than 80% of the reference sequence. Filtering 
##                             parameters can be changed if needed.

##  get_n_align() : get the total number of reads that passed the alignment filter (used to calculate 
##                              variant frequency in downstream analysis)

##  process_mutations() : separates sequences with indels from sequences with only SNPs (and the wild-type). Puts 
##                          sequences with an invalid mutation format in a separate list. This list should be
##                          normally be empty and is there to check for exceptions that were not taken into account
##                          in the previous functions.

##  get_mutated_codons_df() : get the codon which has been mutated from the nucleotide position of the mutation. Formats
##                          everything into a dataframe. Also puts sequences with an invalid mutation format in a list 
##                          (should be empty)

##  get_aa_substitution() : translates the wt and mutated codons to amino acids. Classifies each mutation (one mutation per row)
##                           as synonymous, non-synonymous or stop.

##  group_by_sequence() : groups each mutation from the same sequence together. Returns a df with one row per sequence instead of
##                      one row per mutation

##  get_hamming_dist() : calculates the Hamming distance (difference) between two strings. Used to count the number of mutations
##                          in the DNA and amino acid sequences of each variant compared to the wild-type

##  get_indels_df() : processes the list of sequences with indels from the process_mutations() function. Sequences with indels are not
##                      a focus for this project, so they are not integrated in the same dataframe as the other variants, but they can
##                      still be interesting more for validation purposes

#%% Defining functions - get_wt_seq()
def get_wt_seq(wt_seq_path):

    wt_seq = ""

    with open(wt_seq_path, 'r') as source:

        for line in source:

            if line.startswith('>') == False:  # skips the header
                wt_seq += line.strip('\n')
    
    print(f"The trimmed wild-type sequence is : {wt_seq}")
    
    return wt_seq

#%% Defining functions - parse_needle_output()

def parse_needle_output(path):

    # Initialize counter for the number of alignments
    n_aligns = 0

    # Initialize dictionary that will hold the aligned sequences
    align_seqs_dict = {}
          
    with open(path, 'r') as source:
        
        # Initialize strings for data processing
        current_align = ''
        current_sseq = ''           # Subject sequence (reference)
        current_qseq = ''           # Query sequence (mutated sequence)

        # Initialize a counter to tell if the subject sequence (first sequence of the alignment in the needle file) is complete
        sseq_done = 0

        # Loop through each line of the needle file
        for line in source:

            # Detect header of the alignment -> means that we moved from one alignment to the next (previous alignment is done processing)
            if line.startswith('>>>') == True:

                # Increment alignment counter
                n_aligns +=1

                # Get alignment name (seq_wt; xxx nt vs seq1; size=x; xxx nt)
                align_name = line.strip('>>>')
                
                # If this is not the first alignment 
                if n_aligns != 1:
                    
                    # Add the information on the previous alignment to the dict (done processing the previous alignment)
                    align_seqs_dict[current_align] = [current_sseq, current_qseq]
                    
                    # Update the name for the new entry (start processing the new alignment) and reset temporary variables
                    current_align = align_name
                    current_sseq = ''
                    current_qseq = ''
                    sseq_done = 0

                # If it is the first alignment, no need to store the previous alignment, just start processing the first alignment
                else:
                    current_align = align_name

            # Detect lines that correspond to the sequences (not the header and not useless lines)
            elif line.startswith(';') == False and line.startswith('>') == False and line.startswith('\n') == False and line.startswith('#') == False:

                # If the subject sequence is complete, add the content of the line to the query sequence
                if sseq_done == 1:
                    current_qseq += line.strip('\n')

                # If the subject sequence is not complete, continue to update it
                else:
                    current_sseq += line.strip('\n')

            # Detect the end of the last alignement in the file, add it to the dict
            elif line.startswith('#--') == True:
                align_seqs_dict[align_name] = [current_sseq, current_qseq]

            # Detect the end of the first sequence of the alignment (subject sequence), mark subject sequence as done
            else:
                if sseq_done == 0 and current_sseq != '':
                    sseq_done = 1

    return align_seqs_dict, n_aligns

#%% Defining functions - find_mutations()

def find_mutations(path, nt_before, nt_after):

    # Initialize a dictionary that will hold the nt sequence : mutation information
    allele_dict = {}

    # Initialize a dictionary that will hold the sequences that do not pass the filtering steps 
    unexpected_seqs = {}

    # Parse the needle file and get a dictionary containing the alignments
    align_dict, align_count = parse_needle_output(path)

    # Loop through each alignment
    for entry in list(align_dict.keys()):

        # Initialize the list of mutations found in the query sequence
        read_var_list = []

        # Retrieve the reference and query sequences
        ref_seq = align_dict[entry][0]
        query_seq = align_dict[entry][1]

        # Initialize counters and variables needed to account for the presence of indels
        ## Value used to adjust the index for the presence of insertions in the query sequence
        gap_adjust = 0
        ## Temporary variable to hold the sequence of an indel as a string. Annotation added to read_var_list when the gap ends
        temp_var = None
        ## Start position (in the reference sequence) of the indel annotation, with adjustment for gap presence
        indel_start = 0
        ## Reference sequence with gaps removed (equivalent to the wild-type sequence of the fragment)
        ref_seq_no_gaps = ref_seq.replace('-','')
        ## Query sequences when gaps are removed (what was actually sequenced)
        query_seq_no_gaps = query_seq.replace('-','')


        # Loop through each nucleotide of the alignment
        for nt in range(0, len(ref_seq)):

            # Detect a deletion mutant
            if query_seq[nt] == '-':

                # Logic for detection and annotation of deletion mutants
                # 
                # Suppose we have a deletion mutant with this alignment : 
                #
                # 1 2 3 4 5 6 7 8
                # A T C G A T C G       Reference sequence
                # A T C G - - C G       Query sequence
                # 1 2 3 4 5 6 7 8
                #
                # In that case, the deletion gaps are indexed because the index is based on the reference. 
                # The value of "nt" is therefore equal to the index of the reference

                # If it is the start of a gap in the alignment
                if indel_start == 0:

                    temp_var = '|del|'+ ref_seq[nt]

                    # Start position of the deletion (+1 to convert from 0-based index), compensated for previous insertions if any
                    indel_start = ((nt+1) - gap_adjust)

                # If it is a continuation of a gap in the alignment
                else:
                    temp_var += ref_seq[nt]

            # Detect an insertion mutant
            elif ref_seq[nt] == '-':

                # Logic for detection and annotation of insertion mutants
                # 
                # Suppose we have a deletion mutant with this alignment : 
                #
                # 1 2 3 4     5 6 7 8
                # A T C G - - A T C G       Reference sequence
                # A T C G A A A T C G       Query sequence
                # 1 2 3 4 5 6 7 8 9 10
                #
                # In that case, the insertion gaps are not indexed because the index is based on the reference.
                # The value of nt is therefore not equal to the index of the reference, so we must adjust it when we calculate the position.

                # If it is the start of a gap in the alignment
                if indel_start == 0:

                    temp_var = '|ins|'+ query_seq[nt]
                    indel_start = ((nt+1) - gap_adjust)
                    gap_adjust += 1
                
                # If it is a continuation of a gap in the alignment
                else:
                    temp_var += query_seq[nt]
                    gap_adjust += 1

            # Detect a mutation (mismatch between the reference and query sequence)
            elif query_seq[nt] != ref_seq[nt]:
                # Calculate the position : +1 to convert from 0-based index and adjust for the presence of insertions in the sequence (if there is any)
                position = str(((nt+1) - gap_adjust))
                # Creates an annotation for the mismatch and appends to the list of mutations for that mutant
                variant = ref_seq[nt]+'|' + position + '|' + query_seq[nt]
                read_var_list.append(variant)

            # Detect the end of a gap (no mismatch but gap is open)
            else:
                 if indel_start != 0:
                    # Now that the gap is ended, add the indel to the list of mutations for that mutant
                    read_var_list.append(str((indel_start))+temp_var)

                    # Reset temporary variables for the next indel entry
                    temp_var = None
                    indel_start = 0
        
        # If the sequence ends with an indel (gap is not closed after looping through every nt), add the indel to the list of mutations for that mutant now
        if indel_start != 0:
            read_var_list.append(str((indel_start))+temp_var)

        # Filtering steps 

        ## 1. If there are mutations in the nt positions in the first and last codons that are not actually covered by the fragment,
        ##  remove them
        if query_seq[0:nt_before] != ref_seq[0:nt_before]:
            unexpected_seqs[entry] = query_seq_no_gaps, read_var_list

        elif query_seq[-nt_after:] != ref_seq[-nt_after:]:
            unexpected_seqs[entry] = query_seq_no_gaps, read_var_list

        ## 2. The alignment must cover at least 80% of the reference sequence 
        ## and there must be less than 25 differences between the query and the reference sequence (insertions, deletions or mutations).
        elif len(query_seq_no_gaps) >= len(ref_seq_no_gaps)*0.8 and len(read_var_list) < 25:             
            allele_dict[entry] = query_seq_no_gaps, read_var_list
        
        else :
            unexpected_seqs[entry] = query_seq_no_gaps, read_var_list

    # Show the number of alignment processed
    print(f"{align_count} alignments have been processed for this sample")

    return allele_dict, unexpected_seqs, align_count

#%% Defining functions - get_n_align()

def get_n_align(allele_dict):

    # Initialize the counter for the number of reads
    n_align = 0

    # Loop through every sequence
    for seq in list(allele_dict.keys()):

        # Get the number of reads for that sequence ("size") and add it to the count
        var_info = seq.split(',')
        var_count =int(var_info[1].split(';')[1].strip('size='))
        n_align += var_count

    # Add the number of reads for that sample to a dictionary
    dict_reads_aligned[str(info.loc[Sample]["name"])] = n_align
    
    # Print the number of reads that passed the alignment filters
    return n_align

#%% Defining functions - process_mutations()

def process_mutations(wt_seq, allele_dict):

    # Get the length of the reference sequence
    length_ref = len(wt_seq)

    # Lists to store the sequences with SNPs and sequences with indels
    sequences = []
    sequences_indels = []

    # List to store the sequences with invalid formats
    invalid_sequences = []

    # Process mutations
    for key, (nt_seq, list_mutations) in allele_dict.items():
        
        # Initialize a list of the valid mutated positions for that sequence
        positions = []

        # Initialize a variable to categorize the sequence as : wt or mut
        category = ""

        # Separate the sequences with indels from the sequences with only SNPs
        if len(nt_seq) != length_ref:              
            sequences_indels.append((key, nt_seq, list_mutations))

        # Categorize sequences as "wt" if the mutation list is empty
        elif list_mutations == []:
            category = "wt"
            sequences.append((key, nt_seq, list_mutations, category))
        
        # For other sequences, check if the format of all the mutations in the mutation list is ok
        else : 
            for mut in list_mutations:
                try:
                # Split the mutation by '|', [mutation_type, position, new_base]
                    split_mut = mut.split('|')

                    # Ensure that the split resulted in exactly 3 parts
                    if len(split_mut) == 3:
                        position = int(split_mut[1])  # Extract position
                        positions.append(position)
                    else:
                        print(f"Invalid mutation format: {mut}")
                        invalid_sequences.append((key, nt_seq, list_mutations))

                # If the position is not an integer       
                except ValueError:
                    print(f"Skipping mutation due to invalid position value: {mut}")
                    invalid_sequences.append((key, nt_seq, list_mutations))
                    continue # Skip this mutation

                # If there are other unexpected errors
                except IndexError:
                    print(f"Skipping mutation due to unexpected format: {mut}")
                    invalid_sequences.append((key, nt_seq, list_mutations))
                    continue  # Skip this mutation

        # Proceed with sequences with valid mutations
        if positions and (key, nt_seq, list_mutations) not in invalid_sequences :
            category = "mut"
            sequences.append((key, nt_seq, list_mutations, category))

    # Return the lists of the processed and invalid sequences
    return sequences, sequences_indels, invalid_sequences


#%% Defining functions - get_mutated_codons_df()

def get_mutated_codons_df(sequences_list, wt_seq, sample_name, aa_pos_start):
    
    # Define wild-type codons based on the reference sequence
    sequence_length = len(wt_seq)
    wt_codons = [wt_seq[i:i+3] for i in range(0, sequence_length, 3)]

    # Initialize a list to store all the processed sequences (to transform into a df at the end)
    mut_sequences_data = []

    # Initialize a list to store invalid sequences
    invalid_sequences = []

    # Loop through all the sequences
    for key, nt_seq, list_mutations, category in sequences_list:
        
        # Extract sequence name (seqXXXX) and size (read count)
        seq_match = re.search(r'(seq\d+)', key)
        size_match = re.search(r'size=(\d+)', key)

        seq_name = seq_match.group(1) if seq_match else key  # Default to full key if not found
        read_count = int(size_match.group(1)) if size_match else None  # Extracted sequence read count

        # Add the wild-type sequence to the sequence list
        if category == "wt":
            mut_sequences_data.append({
                "nt_seq" : nt_seq,
                "seq_name" : seq_name + "_" + sample_name,
                "read_count" : read_count,
                "position" : None,
                "wt_codon" : None,
                "mut_codon" : None,
                "seq_type" : "wt"
            })

        # Initialize a dictionary that will contain the mutated position and codon for each sequence
        codon_mutations = {}

        # Loop through all the mutations for one sequence
        for mut in list_mutations:
            try:
                wt_base, position, new_base = mut.split('|')
                position = int(position) # Convert position to integer

                # Find the codon index (0-based)
                codon_pos = (position - 1) // 3  

                # For sequences with a valid codon position (not outside of the range of the fragment or with an invalid format)
                if 0 <= codon_pos < len(wt_codons):
                    # Convert position to 0-based within codon
                    pos_in_codon = (position - 1) % 3

                    # If it is the first mutation in that codon, add the wild-type codon in the dict of mutated codon for that sequence
                    if codon_pos not in codon_mutations:                    
                        codon_mutations[codon_pos] = list(wt_codons[codon_pos])  
                    
                    # Apply the mutation to the wild-type codon stored in the dict
                    # (If it is not the first mutation in that codon, applies the mutation to the already mutated codon)
                    codon_mutations[codon_pos][pos_in_codon] = new_base
                
                # If invalid codon position
                else : 
                    invalid_sequences.append(seq_name)
                    
            except (ValueError, IndexError):
                invalid_sequences.append(seq_name)
                continue

        # Store the mutated codons for all the valid sequences 
        for codon_pos, mutated_codon in codon_mutations.items():

            # From the position in the fragment, get the codon position in the whole sequence
            pos = codon_pos + aa_pos_start                                                           

            mut_sequences_data.append({
                "nt_seq" : nt_seq,
                "seq_name" : seq_name + "_" + sample_name,
                "read_count" : read_count,
                "position" : int(pos),
                "wt_codon" : wt_codons[codon_pos],
                "mut_codon" : "".join(mutated_codon),
                "seq_type" : "mut"})

    # Convert to DataFrame and remove the sequences with a mutation in the last codon (if incomplete)
    df = pd.DataFrame(mut_sequences_data)

    return df, invalid_sequences

#%% Defining functions - get_aa_substitutions()
def get_aa_substitutions(input_df, genetic_code):

    # Initialize a list to put the data for each row of the input df
    df_data = []
    for index, row in input_df.iterrows():

        # Get the full amino acid sequence for each mutated sequence
        nt_seq = row["nt_seq"]
        aa_seq = ""
        for i in range(0, len(nt_seq), 3) :
            codon = nt_seq[i:i+3]
            aa = genetic_code[codon]
            aa_seq += aa
        
        # Get the info on the mutation
        wt_codon = row["wt_codon"]
        mut_codon = row["mut_codon"]
        position = str(row["position"]).replace(".0", "")      # Removes the decimal which is not supposed to be there
        seq_type = row["seq_type"]

        # Process the wild-type sequence
        if seq_type == "wt":
            wt_aa = None
            mut_aa = None
            key_codon = None
            key_aa = None
            variant_id = None
            mut_type = "wt"

        # Process the sequences with mutations
        else : 
            # Get the mutation key (ex. AAA123AAB) for each mutated codon
            key_codon = wt_codon + str(position) + mut_codon

            # Get the mutation key (ex. Y132F) and classify the mutation (non-synonymous, synonymous, stop)
            wt_aa = genetic_code[wt_codon]
            mut_aa = genetic_code[mut_codon]
            if wt_aa != mut_aa : 
                key_aa = wt_aa + str(position) + mut_aa
                variant_id = key_aa
                if mut_aa == "*" :
                    mut_type = "stop"
                else :
                    mut_type = "non-synonymous"
            else : 
                key_aa = None       # No aa mutation key for synonymous mutations 
                variant_id = wt_aa + str(position) + mut_aa
                mut_type = "synonymous"

        # Store the data for that mutation
        df_data.append({
            "nt_seq" : nt_seq,
            "seq_name": row["seq_name"],
            "read_count": row["read_count"],
            "position": position,
            "wt_codon" : wt_codon,
            "mut_codon": mut_codon,
            "dna_mutation" : key_codon,
            "aa_seq" : aa_seq,
            "wt_aa" : wt_aa,
            "mut_aa" : mut_aa,
            "aa_mutation" : key_aa,
            "variant_id" : variant_id,
            "mut_type" : mut_type})

    # Convert list to DataFrame
    return pd.DataFrame(df_data)

#%% Defining functions - group_by_sequence()

def group_by_sequence(input_df):

    # Initialize an empty list to store the data
    df_data = []

    # List of all the sequences
    sequences = list(input_df["nt_seq"].unique())

    for seq in sequences:

        df_seq = input_df[input_df["nt_seq"] == seq]

        # Get the info for the sequence
        seq_name = df_seq.iloc[0]["seq_name"]   
        read_count = df_seq.iloc[0]["read_count"]    
        aa_seq = df_seq.iloc[0]["aa_seq"]  

        ## Initiliaze list of mutations for the sequence
        mut_positions = []
        dna_mutations = []
        aa_mutations = []
        variant_id = []

        ## Initialize the variant_type variable
        variant_type = ""

        ## Process the sequence
        for _, row in df_seq.iterrows():

            # Add the mutated position and mutation code to the list for that sequence
            mut_positions.append(row["position"])
            dna_mutations.append(row["dna_mutation"])
            aa_mut = row["aa_mutation"]
            if aa_mut != None:
                aa_mutations.append(aa_mut)
            variant_id.append(row["variant_id"])
            
            # Get the mutation type
            mut_type = row["mut_type"]

            # If one of the mutation in the sequence leads to a stop codon, categorize the sequence as "stop"
            if mut_type == "stop":
                variant_type = "stop"

            # If this is the wild-type sequence
            elif mut_type == "wt":
                variant_type = "wt"

        # If the sequence contains only synonymous mutations and is not the wild-type
        if aa_mutations == [] and variant_type != "wt":
            variant_type = "synonymous" 

        # If the sequence contains at least one non-synonymous mutation and no stop codons:
        elif variant_type == "":
            variant_type = "non-synonymous"
  
        # Add the sequence information to make the dataframe
        df_data.append({
            "nt_seq" : seq,
            "seq_name" : seq_name,
            "read_count" : read_count,
            "position" : mut_positions,
            "dna_mutation" : dna_mutations, 
            "aa_seq" : aa_seq,
            "aa_mutation" : aa_mutations,
            "variant_id" : variant_id,
            "variant_type" : variant_type})
    
    df = pd.DataFrame(df_data)

    return df

#%% Defining functions - get_hamming_dist()

def get_hamming_dist(stringA, stringB):
    i = 0
    count = 0
    while i < len(stringA) : 
        if stringA[i] != stringB[i] : 
            count +=1
        i += 1
    return count

#%% Defining functions - get_indels_df()

def get_indels_df(seq_list, aa_pos_start, sample_name):
    
    # Initialize a list to store the info for each sequence 
    indels_data = []

    # Loop through each sequence
    for key, nt_seq, list_mutations in seq_list:
        
        # Extract sequence name (seqXXXX) and size (read count)
        seq_match = re.search(r'(seq\d+)', key)
        size_match = re.search(r'size=(\d+)', key)

        seq_name = seq_match.group(1) if seq_match else key  # Default to full key if not found
        size = int(size_match.group(1)) if size_match else None  # Extracted sequence size

        # Loop through the mutations for that sequence
        for muta in list_mutations : 
            a, b, c = muta.split("|")

            # For indels
            if b in ["ins", "del"]:

                # Get the codon position in the fragment and convert to the position in the full sequence
                position = int(a)
                codon_pos = (position - 1) // 3
                pos_in_sequence = codon_pos + aa_pos_start

                # Get the length of the indel and if it causes a frameshift
                indel_length = len(c)     
                frameshift = "frameshift"
                if indel_length % 3 == 0 : 
                    frameshift = "no_frameshift"

        # Add the sequence info
        indels_data.append({
            "nt_seq" : nt_seq,
            "seq_name": seq_name + "_" + sample_name,
            "read_count": size,
            "indel_position": pos_in_sequence,  
            "list_indels_muts" : list_mutations, 
            "frameshift" : frameshift})

    # Convert to DataFrame and remove the sequences with a mutation in the last codon (if incomplete)
    df = pd.DataFrame(indels_data)
    return df

#%% Defining functions - get_heatmap_format_df()

def get_heatmap_format_df(input_df, list_codons, start_pos, end_pos):

    # Initialize the dataframe
    df_heatmap = pd.DataFrame(columns = list(range(start_pos, end_pos+1)), index = codons_mutated, dtype = float)
    
    # Fill the dataframe
    for pos in range(start_pos, end_pos + 1):
        pos_data = input_df[input_df["position"] == pos]        # Make a subset of the df for each position
        pos_data = pos_data.groupby(["mutated", "position"], as_index = False)["log2_read_count"].median()      # This is another formatting step, the median() function doesn't do anything because there is only one value per group
        pos_data.set_index("mutated", inplace = True)
        for mut in pos_data.index:        # Add the info to the heatmap dataframe
            df_heatmap.at[mut, pos] = pos_data.at[mut, "log2_read_count"]

    df_heatmap = df_heatmap.fillna(0)

#%% Putting everything together!

# This step outputs a dataframe containing all the sequences with their mutations for each sample

# Initialize a dictionnary to store the number of reads that passed the alignment filter (for the get_n_align() function)
dict_reads_aligned = {}

# Process each sample
for Sample in list(info.index):
    
    sample_name = str(info.loc[Sample]['name'])

    # Print things to follow along the analysis
    print(f"Now processing alignments from sample {sample_name}")

    # Retrieve the wild-type sequence for the fragment from the fasta file
    wt_seq = get_wt_seq(wt_seq_path)

    # Path to the file with the aligned sequences
    needle_file = folder_aligned + sample_name + '_mutagenesis-library.needle'

    # Get the number of nucleotides in the first and last codon that are not covered by the fragment
    nt_before = info.loc[Sample]["nt_before"]
    nt_after = info.loc[Sample]["nt_after"]

    # Create dictionnary of all the sequences with a list of the mutations
    mutated_seqs, unexpected_seqs, nb_align = find_mutations(needle_file, nt_before, nt_after) 
    print(f"There were {len(unexpected_seqs)} unexpected sequences for sample {sample_name}")

    # Get the number of reads of the sequences that passed the quality filters
    # and compare to the total number of reads that were aligned
    nb_reads_ok = get_n_align(mutated_seqs)  
    total_reads_aligned = df_read_counts.loc[sample_name]["reads_aligned"]
    percent_reads_ok = 100*(nb_reads_ok/total_reads_aligned)
    print(f"{nb_reads_ok} reads passed the quality filters ({percent_reads_ok:.3f}% of the total number of aligned reads)")

    # Separate the sequences with indels
    mut_sequences, indel_sequences, invalid_sequences1 = process_mutations(wt_seq, mutated_seqs)
    print(f"There are {len(mut_sequences)} mutated sequences and {len(indel_sequences)} sequences with indels in sample {sample_name}")
    print(f"There are {len(invalid_sequences1)} sequences with an invalid mutation code format for the 'process_mutations' step")        # This number should be 0

    # For the mutated sequences, get which codon is mutated for all mutations. At this point, there is one line
    # per mutation so sequences with multiple mutations are spread across different rows
    aa_start = info.loc[Sample]['aa_pos_start']
    df_mut_sequences, invalid_sequences2 = get_mutated_codons_df(mut_sequences, wt_seq, sample_name, aa_start)
    print(f"There are {len(invalid_sequences2)} sequences with an invalid mutation code format for the 'get_mutated_codons_df' step")        # This number should be 0

    # Translate the codons to amino acids and get the mutation key (ex. Y132F) and mutation type (synonymous, non-synonymous, stop)
    df_mut_sequences_aa = get_aa_substitutions(df_mut_sequences, genetic_code)

    # Group all the mutations from the same sequence (get a df with one row per sequence instead of one row per mutation)
    df_final = group_by_sequence(df_mut_sequences_aa)

    # Calculate the hamming distance between the mutated sequences and the wild-type (DNA and aa sequences)
    df_final["hamming_dist_DNA"] = [get_hamming_dist(wt_seq, x) for x in df_final["nt_seq"]]

    wt_seq_aa = df_final[df_final["variant_type"] == "wt"].loc[0]["aa_seq"]
    df_final["hamming_dist_AA"] = [get_hamming_dist(wt_seq_aa, x) for x in df_final["aa_seq"]]

    # Add the sample information
    df_final["fragment"] = info.loc[Sample]["fragment"]
    df_final["sample_name"] = info.loc[Sample]["name"]
    df_final["condition"] = info.loc[Sample]["condition"]
    
    # Process the indels and add the sample information 
    df_indels = get_indels_df(indel_sequences, aa_start, sample_name)
    df_indels["fragment"] = info.loc[Sample]["fragment"]
    df_indels["sample_name"] = info.loc[Sample]["name"]
    df_indels["condition"] = info.loc[Sample]["condition"]

    # Concatenate the dataframes (variants + indels) of each sample into one
    if Sample == 1:
        df_all = df_final
        df_all_indels = df_indels
    
    else :
        df_all = pd.concat([df_all, df_final])
        df_all_indels = pd.concat([df_all_indels, df_indels])


    # For the single mutants, make a read count heatmap
    df_single = df_final[df_final["hamming_dist_DNA"] == 1]
    df_single["mut_codon"] = [x[0][-3:] for x in df_single["dna_mutation"]]
    df_single["pos"] = [x[0] for x in df_single["position"]]
    df_single["log2_read_count"] = np.log2(df_single["read_count"])
    
    ## Pivot the dataframe to have a df in heatmap format
    df_heatmap = df_single.pivot_table(index = "mut_codon", columns = "pos", values = "log2_read_count")
    df_heatmap.fillna(value = 0, inplace = True)

    ## Make the heatmap
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize = (20, 20))
    ax = sns.heatmap(df_heatmap, square = True, cmap = "crest", vmin = 1,            # Mutations present in the library
                mask = df_heatmap == 0, linewidths = 0.002,
                cbar_kws= {"pad": 0.01,'label': 'Log2 read count', "shrink" : 0.4})
    color = matplotlib.colors.ListedColormap("#dedcdc")
    ax = sns.heatmap(df_heatmap, mask = df_heatmap != 0, cmap = color,
                     linewidths = 0.002, cbar = False)
    
    ## Title and labels 
    ax.set_xlabel("Position", fontsize = 16)
    ax.set_ylabel("Codons", fontsize = 16)
    ax.set_title(f"Read counts for sample {sample_name}", fontsize = 18)

    ## Save the heatmap 
    fig_name = folder_figures + "Heatmap_read_count_DNA_single_mutants_" + sample_name + "_" + experiment + "_" + version
    plt.savefig(f"{fig_name}.png", format = "png", dpi = 300, bbox_inches = "tight")

    print(f"Done processing sample {sample_name}!")
    print(f"#########################################")


# Save the final dataframe (with all the sequences for all samples)
name = folder_processed_data + "Variant_calling_all_samples_" + experiment + "_" + version + ".csv"
df_all.to_csv(name)

# Make a dataframe with the total number of reads that passed the quality filters for each sample and save it 
df_reads_ok = pd.DataFrame.from_dict(data = dict_reads_aligned, orient = "index", columns = ["total_read_count"]).reset_index(names = "sample")
name = folder_processed_data + "Total_reads_variant_calling_" + experiment + "_" + version + ".csv"
df_reads_ok.to_csv(name)

# Save the indels dataframe
name = folder_dataframes + "Variant_calling_indels_all_samples_" + experiment + "_" + version + ".csv"
df_all_indels.to_csv(name)


#%% Check the coverage of each sequence

# Load the dataframes
## Variants
path = "../results/Variant_calling_all_samples_mutagenesis-library_2026-04-15.csv"
df_variants = pd.read_csv(path, index_col=0)
## Read counts
path = "../results/Total_reads_variant_calling_mutagenesis-library_2026-04-15.csv"
df_read_counts = pd.read_csv(path, index_col = 1).drop("Unnamed: 0", axis = 1)

# Calculate the proportion of wild-type reads in each library
df_wt = df_variants[df_variants["variant_type"] == "wt"]
df_wt["total_reads"] = [df_read_counts.loc[x]["total_read_count"] for x in df_wt["sample_name"]]
df_wt["%_total_reads"] = df_wt["read_count"] / df_wt["total_reads"]
## Save the dataframe
name = folder_dataframes + "Percent_wt_reads_" + experiment + "_" + version + ".csv"
df_wt.to_csv(name)

# Count the mean (including and excluding the wild-type sequence) and median read count for each library
df_mean_rc = df_variants.groupby("sample_name")["read_count"].mean().reset_index(name = "mean_read_count")
df_mean_rc_no_wt = df_variants[df_variants["variant_type"] != "wt"].groupby("sample_name")["read_count"].mean().reset_index(name = "mean_read_count_no_wt")
df_median_rc = df_variants.groupby("sample_name")["read_count"].median().reset_index(name = "median_read_count")
## Merge into one dataframe and save it
df_read_count_stats = df_mean_rc.merge(df_mean_rc_no_wt, on = "sample_name")
df_read_count_stats = df_read_count_stats.merge(df_median_rc, on = "sample_name")

name = folder_dataframes + "Read_count_stats_" + experiment + "_" + version + ".csv"
df_read_count_stats.to_csv(name)

#%% Look at the read count distribution for each library 

# Load the dataframe
path = "../results/Variant_calling_all_samples_mutagenesis-library_2026-04-15.csv"
df_variants = pd.read_csv(path, index_col=0)

# Calculate the log2 read count
df_variants["log10_rc"] = np.log10(df_variants["read_count"])

# Make the violin plot
fig, ax = plt.subplots(figsize = (5, 4))
sns.violinplot(data = df_variants, x = "sample_name", y = "log10_rc", inner = "box")

ax.set_xlabel("Library")
ax.set_ylabel("Read count (log10)")

# Save the figure
fig_name = folder_figures + "Read_count_distribution_" + experiment + "_" + version
plt.savefig(f"{fig_name}.png", format='png', dpi=300,bbox_inches="tight")

# %% Check the library diversity

# Load the dataframe
path = "../results/Variant_calling_all_samples_mutagenesis-library_2026-04-15.csv"
df_variants = pd.read_csv(path, index_col=0)

# Get the number of sequences with one, two, etc. mutations (DNA) in each library
df_diversity_DNA = df_variants.groupby(["sample_name", "hamming_dist_DNA"])["nt_seq"].count().reset_index(name = "nb_of_sequences")

# Do the same for a.a mutations
df_diversity_AA = df_variants.groupby(["sample_name", "hamming_dist_AA"])["aa_seq"].count().reset_index(name = "nb_of_sequences")

#%% Look at the overlap between libraries

# Load the dataframe
path = "../results/Variant_calling_all_samples_mutagenesis-library_2026-04-15.csv"
df_variants = pd.read_csv(path, index_col=0)

# Count in how many libraries each DNA sequence is in 
df_overlap_DNA = df_variants.groupby("nt_seq")["sample_name"].count().reset_index(name = "nb_of_libraries")
print(f"There are {len(df_overlap_DNA[df_overlap_DNA["nb_of_libraries"] >= 2])} sequences present in at least two libraries")
 
# Same thing for amino acid sequences
df_overlap_AA = df_variants.groupby("aa_seq")["sample_name"].count().reset_index(name = "nb_of_libraries")
print(f"There are {len(df_overlap_AA[df_overlap_AA["nb_of_libraries"] >= 2])} variants present in at least two libraries")
