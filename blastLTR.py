import subprocess
import os

def blast(LTR, genome, output):
    path = os.path.dirname(LTR)
    output_blast = os.path.join(path, output + "_blast" +".txt")


    # Define the BLAST command
    blast_command = [
        "blastn",
        "-query", LTR,
        "-subject", genome,
        "-outfmt", "6 qseqid sseqid sstart send evalue",
        "-out", output_blast,
    ]

    # Run the BLAST command
    subprocess.run(blast_command)

    print(f"BLAST results saved to {output_blast}")
    return output_blast