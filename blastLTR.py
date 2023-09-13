import subprocess
import os

def blast(LTR, genome, output, path):
    output_blast = os.path.join(path, output + "_blast" +".txt")


    # Define the BLAST command
    blast_command = [
        "blastn",
        "-query", LTR,
        "-subject", genome,
        "-outfmt", "6 qseqid sseqid sstart send evalue sstrand",
        "-out", output_blast,
    ]

    # Run the BLAST command
    subprocess.run(blast_command)

    print(f"BLAST results saved to {output_blast}")
    if os.path.getsize(output_blast) == 0:
        TE_presence = False
    else:
        TE_presence = True

    return output_blast, TE_presence