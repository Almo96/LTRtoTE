# Author Almorò Scarpa
# Program to find LTR-retrotransposons in a genome from their LTR

import argparse
import os
import InputOutput
import blastLTR
import TEfromLTR
import TSD

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script with mandatory and optional fields.")
    
    # Add the mandatory positional arguments
    parser.add_argument("LTR", help="LTR file (.fasta or .fa)")
    parser.add_argument("genome", help="Genome file (.fasta or .fa)")
    parser.add_argument("output", help="Output file (.fasta or .fa)")
    parser.add_argument("minlenTE", help="Minimum lenght of the TE (integer)")
    parser.add_argument("maxlenTE", help="Maximum lenght of the TE (integer)")

    # Add the optional argument with a default value of 20
    parser.add_argument("--tsd", help="Optional TSD field (integer)", default=20, type=int)
    args = parser.parse_args()

    LTR = args.LTR
    genome = args.genome
    output = args.output
    minlenTE = int(args.minlenTE)
    maxlenTE = int(args.maxlenTE)
    tsd = int(args.tsd)

    InputOutput.controlIO(LTR, genome, output, minlenTE, maxlenTE, tsd)

    blast_path = blastLTR.blast(LTR, genome, output)
    
    output_bed, bed = TEfromLTR.TE_position(output, blast_path, minlenTE, maxlenTE)

    output_fasta = TEfromLTR.TEfasta(genome, output_bed, output)

    output_left, output_right, left_fasta, right_fasta = TSD.TSDbed(bed, tsd, genome, output)

    TSD.checkTSD(left_fasta, right_fasta, output, genome)

    TSD.remove_temp(output_left, output_right, left_fasta, right_fasta)
