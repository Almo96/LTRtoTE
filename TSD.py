import subprocess
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

def TSDbed(bed, tsd, genome, output):
    # Positions to the left and the right of the TE
    df_left = bed.copy()
    df_left['Start-TSD'] = df_left['Start'] - tsd
    df_left['Start'] = df_left['Start'] - 1
    df_left = df_left.drop('End', axis=1)
    df_left = df_left.reindex(columns=['Chromosome', 'Start-TSD', 'Start'])

    df_right = bed.copy()
    df_right['End+TSD'] = df_right['End'] + tsd
    df_right['End'] = df_right['End'] + 1
    df_right = df_right.drop('Start', axis=1)

    # Create the temporary .bed files
    path = os.path.dirname(genome)
    output_left = os.path.join(path, output + "TSD_left" + ".bed")
    output_right = os.path.join(path, output + "TSD_right" + ".bed")

    df_left.to_csv(output_left, sep='\t', header=False, index=False)
    df_right.to_csv(output_right, sep='\t', header=False, index=False)

    # Use bedtools to search to retrieve the fasta sequences
    left_fasta = os.path.join(path, output + "TSD_left" + ".fasta")
    right_fasta = os.path.join(path, output + "TSD_right" + ".fasta")

    command_l = ["bedtools", "getfasta", "-fi", genome, "-bed", output_left, "-fo", left_fasta]
    command_r = ["bedtools", "getfasta", "-fi", genome, "-bed", output_right, "-fo", right_fasta]

    try:
        subprocess.run(command_l, check=True)
        print("fasta sequence left generated")
    except subprocess.CalledProcessError as e:
        print("Error:", e)

    try:
        subprocess.run(command_r, check=True)
        print("fasta sequence right generated")
    except subprocess.CalledProcessError as e:
        print("Error:", e)
    
    return output_left, output_right, left_fasta, right_fasta

def remove_temp(output_left, output_right, left_fasta, right_fasta):
    command_remove = ["rm", output_left, output_right, left_fasta, right_fasta]

    try:
        subprocess.run(command_remove, check=True)
        print("Temporary files removed")
    except subprocess.CalledProcessError as e:
        print("Error:", e)

def checkTSD(left_fasta, right_fasta, output, genome):
    # Load sequences from the two FASTA files
    sequences1 = list(SeqIO.parse(left_fasta, "fasta"))
    sequences2 = list(SeqIO.parse(right_fasta, "fasta"))
    
    if len(sequences1) != len(sequences2):
        print("Error: The number of sequences in the two files is not the same.")
        return
    
    data = {
        "Sequence ID": [],
        "Common TSD": [],
        "TSD Length": []
    }
    
    for seq1, seq2 in zip(sequences1, sequences2):
        seq1_str = str(seq1.seq).upper()
        seq2_str = str(seq2.seq).upper()

        common_tsd = ""
        
        while len(seq1_str) != 0:
            if seq1_str == seq2_str:
                common_tsd = seq1_str
                data["Sequence ID"].append(seq1.id)
                data["Common TSD"].append(common_tsd)
                data["TSD Length"].append(len(common_tsd))
                break
            else:
                seq1_str = seq1_str[:-1]  # Remove one letter from the left of seq1
                seq2_str = seq2_str[1:]   # Remove one letter from the right of seq2
        
        if not common_tsd:
            data["Sequence ID"].append(seq1.id)
            data["Common TSD"].append("No tsd")
            data["TSD Length"].append(0)
    
    path = os.path.dirname(genome)
    output_TSD = os.path.join(path, output + "TSD" + ".txt")

    df_TSD = pd.DataFrame(data)
    df_TSD.to_csv(output_TSD, sep='\t', header=False, index=False)

    print(f"TSD analysis saved to {output_TSD}")

    # Generate a figure with the TSD length distribution
    bins = list(range(0, max(df_TSD['TSD Length']) + 3))
    plt.figure(figsize=(8, 6))

    # Use plt.hist to calculate the histogram data (counts)
    hist, _ = np.histogram(df_TSD['TSD Length'], bins=bins)
    plt.bar(bins[:-1], hist, width=1.0, align='center', alpha=0.7)

    plt.xlabel('TSD Length')
    plt.ylabel('Count')
    plt.title('TSD Length Distribution')
    plt.grid(axis='y', alpha=0.75)

    plt.xticks(bins)
    plt.yticks(range(int(max(hist)) + 1))

    plt.savefig(os.path.join(path, output + "_TSD_distribution" + ".png"))
    print(f"TSD length distribution saved in {output_TSD}")
