import pandas as pd
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import math

def TE_position(output, output_blast, minlenTE, maxlenTE):
    df = pd.read_csv(output_blast, sep='\t', header=None)
    path = os.path.dirname(output_blast)
    output_bed = os.path.join(path, output + ".bed")
    output_summary = os.path.join(path, output + "_summary" + ".txt")
    output_TE_len = os.path.join(path, output + "_TE_len" + ".png")


    # Sort the DataFrame by Chromosome and LTR starting position
    df_sorted = df.sort_values(by=[1, 2], ascending=[True, True])
    summary_df = pd.DataFrame(columns=[0, 1, 2, 3])

    # Iterate through df_sorted and check for the right distance between two LTRs
    for i in range(len(df_sorted) - 1):
        if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] < maxlenTE:
            if df_sorted.iloc[i + 1, 3] - df_sorted.iloc[i, 2] > minlenTE:
                new_row = df_sorted.iloc[i, [0, 1, 2]].tolist() + [df_sorted.iloc[i + 1, 3]]
                summary_df.loc[len(summary_df)] = new_row

    summary_df.columns = ['TE', 'Chromosome', 'Start', 'End']
    summary_df['Length'] = summary_df['End'] - summary_df['Start']
    summary_df = summary_df.replace('LTR','', regex=True)
    summary_df = summary_df.replace('_','', regex=True)
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


    bed = summary_df[['Chromosome', 'Start', 'End']]
    bed.to_csv(output_bed, sep='\t', header=False, index=False)

    print(f".bed results saved to {output_bed}")
    return output_bed, bed

def TEfasta(genome, output_bed, output):
    path = os.path.dirname(output_bed)
    output_fasta = os.path.join(path, output + ".fasta")
    command = ["bedtools", "getfasta", "-fi", genome, "-bed", output_bed, "-fo", output_fasta]

    try:
        subprocess.run(command, check=True)
        print(".bed extraction executed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error:", e)

    print(f"TE sequences saved to {output_fasta}")
    return output_fasta
