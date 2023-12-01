import pandas as pd
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import math

def TE_position(output, output_blast, minlenTE, maxlenTE, path):
    df = pd.read_csv(output_blast, sep='\t', header=None)
    output_bed = os.path.join(path, output + ".bed")
    output_summary = os.path.join(path, output + "_summary" + ".txt")
    output_TE_len = os.path.join(path, output + "_TE_len" + ".png")

    # Invert values for the minus strand
    df2 = df.copy(deep=True)
    for i in range(len(df2)-1):
        if df2.iloc[i, 5] == "minus":
            x = df2.iloc[i, 2]
            y = df2.iloc[i, 3]
            df2.iloc[i, 2] = y
            df2.iloc[i, 3] = x
    blast2 = os.path.join(path, output +'_blast_rev.txt')
    df2.to_csv(blast2, sep='\t', header=False, index=False)

    df3 = pd.read_csv(blast2, sep='\t', header=None)
    # Sort the DataFrame by Chromosome and LTR starting position
    df_sorted = df3.sort_values(by=[1, 2], ascending=[True, True])
    summary_df = pd.DataFrame(columns=[0, 1, 2, 3, 4])
    
    # Iterate through df_sorted and check for the right distance between two LTRs
    for i in range(len(df_sorted) - 1):
        if df_sorted.iloc[i + 1, 1] == df_sorted.iloc[i, 1]:        # Check if they are on the same Chromosome
            if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] < maxlenTE:
                if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] > minlenTE:        # Check if they are within the specified range
                    if df_sorted.iloc[i + 1, 5] == df_sorted.iloc[i, 5]:
                        if df_sorted.iloc[i + 1, 5] == df_sorted.iloc[i, 5]:                # Check if they are on the same strand
                            new_row = df_sorted.iloc[i, [0, 1, 2]].tolist() + df_sorted.iloc[i + 1, [3, 5]].tolist()
                            summary_df.loc[len(summary_df)] = new_row

    summary_df.columns = ['TE', 'Chromosome', 'Start', 'End', 'Strand']
    summary_df["Strand"] = summary_df["Strand"].replace("plus", "+")
    summary_df["Strand"] = summary_df["Strand"].replace("minus", "-")

    summary_df['Length'] = summary_df['End'] - summary_df['Start']
    #summary_df = summary_df.replace('LTR','', regex=True)
    #summary_df = summary_df.replace('_','', regex=True)
    summary_df = summary_df.sort_values(by=['Length'], ascending=[False])
    summary_df.to_csv(output_summary, sep='\t', header=False, index=False)
    print(f"Number of TEs: ", len(summary_df))

    # Genrate an image of the TE length distribution
    min_bin = math.floor(minlenTE / 100) * 100
    max_bin = math.ceil(maxlenTE / 100) * 100
    plt.figure(figsize=(8, 6))

    # Create a histogram with specified bins
    plt.hist(summary_df['Length'], bins=range(min_bin, max_bin + 101, 100), edgecolor='k')

    # Set labels and title
    plt.xlabel('Bins (Size 100)')
    plt.ylabel('Count')
    plt.title('TE length distribution')

    plt.savefig(output_TE_len)
    print(f"TE length distribution saved in {output_TE_len}")

    bed = pd.DataFrame({'Chromosome': summary_df['Chromosome'], 'Start': summary_df['Start'], 'End': summary_df['End'], 'Label': "L", 'Score': 0, 'Strand': summary_df['Strand']})
    bed_sorted = bed.sort_values(["Chromosome", "Start"], ascending=True)
    bed_sorted.to_csv(output_bed, sep='\t', header=False, index=False)

    print(f".bed results saved to {output_bed}")
    
    return output_bed, bed_sorted, blast2


def getfasta(genome, output_bed, output, path, out_extra):
    output_fasta = os.path.join(path, output + out_extra +".fasta")
    command = ["bedtools", "getfasta", "-fi", genome, "-bed", output_bed, "-s", "-fo", output_fasta]

    try:
        subprocess.run(command, check=True)
        print(".bed extraction executed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error:", e)

    print(f"Sequences saved to {output_fasta}")
    return output_fasta
