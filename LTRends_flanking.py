# From each insertion 8 regions:
# R1 = Genomic DNA flanking LTR-1
# R2 = LTR-1-left 
# R3 = LTR-1-right
# R4 = TE coding region flanking LTR-1
# R5 = TE coding region flanking LTR-2
# R6 = LTR-2-left 
# R7 = LTR-2-right
# R8 = Genomic DNA flanking LTR-2

import os
import pandas as pd
from Bio import SeqIO

def LTRe_fDNA(output, path, output_blast, minlenTE, maxlenTE):

    # Part 1: create a .bed file with the LTR positions of each TE
    df = pd.read_csv(output_blast, sep='\t', header=None)

    # Sort the DataFrame by Chromosome and LTR starting position
    df_sorted = df.sort_values(by=[1, 2], ascending=[True, True])
    df_LTR = pd.DataFrame(columns=["Chr", "Start", "End", "Score", "Strand"])

    # Iterate through df_sorted and check for the right distance between two LTRs
    for i in range(len(df_sorted) - 1):
        if df_sorted.iloc[i + 1, 1] == df_sorted.iloc[i, 1]:        # Check if they are on the same Chromosome
            if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] < maxlenTE:
                if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] > minlenTE:        # Check if they are within the specified range
                    new_row1 = df_sorted.iloc[i, [1, 2, 3, 4, 5]].tolist()
                    new_row2 = df_sorted.iloc[i + 1, [1, 2, 3, 4, 5]].tolist()
                    df_LTR.loc[len(df_LTR)] = new_row1  # Add LTR left
                    df_LTR.loc[len(df_LTR)] = new_row2  # Add LTR right
    
    df_LTR["Strand"] = df_LTR["Strand"].replace("plus", "+")
    df_LTR["Strand"] = df_LTR["Strand"].replace("minus", "-")

    #Part 2: use the start and end position of each LTR pair to extract the sequence of the 8 regions of interest
    x = 10 # Genomic region
    y = 10 # LTR region
    c = 10 # Coding region

    df_regions = pd.DataFrame(columns=["Chr", "Start", "End", "Region", "Score", "Strand"])

    for index, row in df_LTR.iterrows():

        if index % 2 == 0:
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - x - 1, row['Start'] - 1, 'R1', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
        
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - 1, row['Start'] + y - 1, 'R2', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - y, row['End'], 'R3', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'], row['End'] + c, 'R4', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])

        else:
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - x - 1, row['Start'] - 1, 'R5', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
        
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - 1, row['Start'] + y - 1, 'R6', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - y, row['End'], 'R7', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])
            
            df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'], row['End'] + c, 'R8', 0, row['Strand']]],
                                                              columns=["Chr", "Start", "End", "Region", "Score", "Strand"])])


    output_df_regions = os.path.join(path, output + "_regions" + ".bed")
    df_regions.to_csv(output_df_regions, sep='\t', header=False, index=False)
    print(f"Boundary regions of the LTR saved to {output_df_regions}")

    return output_df_regions


def LTRe_fDNA_out(regions4_fasta, output, path):

    sequences = list(SeqIO.parse(regions4_fasta, "fasta"))
    tr = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'a':'t', 't': 'a', 'c': 'g', 'g': 'c', 'n':'n'}


    if len(sequences) % 8 == 0:
        str_regions = []
        string_seq = ""
        counter_seq = 1
        neg_counter = 1
        for sequence in sequences:
            if sequence.id.endswith("(+)"):
                if counter_seq == 1:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "||"
                elif counter_seq == 2:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "..."
                elif counter_seq == 3:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "|"
                elif counter_seq == 4:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "......"
                elif counter_seq == 5:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "|"
                elif counter_seq == 6:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "..."
                elif counter_seq == 7:
                    string_seq = string_seq + (str(sequence.seq).upper()) + "||"
                elif counter_seq == 8:
                    string_seq = string_seq + (str(sequence.seq).upper())
                    str_regions.append(string_seq)
                    counter_seq = 0
                    string_seq = ""
                counter_seq = counter_seq + 1
            else:
                if neg_counter == 1:
                    result = [tr[base] for base in reversed(sequence.seq)]  #bedtools provided the reverse complement but on the genomic regions
                    R1 = ''.join(result)                                    # it is not needed, so by performing it again I get the original sequence
                    R1 = R1.upper()
                elif neg_counter == 2:
                    R7 = str(sequence.seq).upper()
                elif neg_counter == 3:
                    R6 = str(sequence.seq).upper()
                elif neg_counter == 4:
                    R5 = str(sequence.seq).upper()
                elif neg_counter == 5:
                    R4 = str(sequence.seq).upper()
                elif neg_counter == 6:
                    R3 = str(sequence.seq).upper()
                elif neg_counter == 7:
                    R2 = str(sequence.seq).upper()
                elif neg_counter == 8:
                    result2 = [tr[base] for base in reversed(sequence.seq)]
                    R8 = ''.join(result2)
                    R8 = R8.upper()                    
                    R = R1 + "||" + R2 + "..." + R3 + "|" + R4 + "......" + R5 + "|" + R6 + "..." + R7 + "||" + R8
                    str_regions.append(R)
                    neg_counter = 0
                    string_seq = ""
                neg_counter = neg_counter + 1

        output_str_regions = os.path.join(path, output + "_regions" + ".txt")

        with open(output_str_regions, "w") as file:
            for item in str_regions:
                file.write(f"{item}\n")

        print(f"The genomic flanking, the LTR and the coding regions were saved to {output_str_regions}")

    else:
        print("The number of sequences in the regions.bed file should be a multiple of 4")

    return output_str_regions

# This function works, but if there are deletion inside the LTR the nucleotides are shifted
# A better approach is to take the Start and End of each LTR instead of the start and end of the TE
# def LTRe_fDNAv1(output_bed, LTR, output):
#     x = 10 # Genomic region
#     y = 10 # LTR region
#     c = 10 # Coding region

#     sequence = list(SeqIO.parse(LTR, "fasta"))
#     first_sequence = sequence[0]
#     LTR_length = len(first_sequence)

#     column_names = ["Chr", "Start", "End", "Label", "Score", "Strand"]
#     df = pd.read_csv(output_bed, sep='\t', header=None, names=column_names)
#     df_regions = pd.DataFrame(columns=["Chr", "Start", "End", "Region", "Score", "Strand"])
    
#     for _, row in df.iterrows():

#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] - x - 1, row['Start'] - 1, 'R1']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'], row['Start'] + y, 'R2']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] + LTR_length - y, row['Start'] + LTR_length, 'R3']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['Start'] + LTR_length + 1, row['Start'] + LTR_length + c + 1, 'R4']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - LTR_length - c - 1, row['End'] - LTR_length - 1, 'R5']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - LTR_length , row['End'] - LTR_length + y, 'R6']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] - y, row['End'], 'R7']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
        
#         df_regions = pd.concat([df_regions, pd.DataFrame([[row['Chr'], row['End'] + 1, row['End'] + x + 1, 'R8']],
#                                                           columns=["Chr", "Start", "End", "Region"])])
    
#     path = os.path.dirname(LTR)
#     output_df_regions = os.path.join(path, output + "_regions" + ".bed")
#     df_regions.to_csv(output_df_regions, sep='\t', header=False, index=False)
#     print(f"Boundary regions of the LTR saved to {output_df_regions}")

#     return output_df_regions