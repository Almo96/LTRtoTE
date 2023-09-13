# Author Almorò Scarpa
# Program to find LTR-retrotransposons in a genome from their LTR

import argparse
import os
import InputOutput
import blastLTR
import TEfromLTR
import TSD
import LTRends_flanking

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script with mandatory and optional fields.")
    
    # Add the mandatory positional arguments
    parser.add_argument("LTR", help="LTR file (.fasta or .fa)")
    parser.add_argument("genome", help="Genome file (.fasta or .fa)")
    parser.add_argument("output_path", help="Output path (.fasta or .fa)")
    parser.add_argument("output", help="Output file (.fasta or .fa)")
    parser.add_argument("minlenTE", help="Minimum lenght of the TE (integer)")
    parser.add_argument("maxlenTE", help="Maximum lenght of the TE (integer)")

    # Add the optional argument with a default value of 20
    parser.add_argument("--tsd", help="Optional TSD field (integer)", default=20, type=int)
    args = parser.parse_args()

    LTR = args.LTR
    genome = args.genome
    output = args.output
    output_path = args.output_path
    minlenTE = int(args.minlenTE)
    maxlenTE = int(args.maxlenTE)
    tsd = int(args.tsd)

    InputOutput.controlIO(LTR, genome, output, output_path, minlenTE, maxlenTE, tsd)

    output_blast, TE_presence = blastLTR.blast(LTR, genome, output, output_path)
    
    if TE_presence == True:

        output_bed, bed, blast2 = TEfromLTR.TE_position(output, output_blast, minlenTE, maxlenTE, output_path)

        output_fasta = TEfromLTR.getfasta(genome, output_bed, output, output_path, "")

        #output_left, output_right, left_fasta, right_fasta = TSD.TSDbed(bed, tsd, genome, output)

        #TSD.checkTSD(left_fasta, right_fasta, output, genome)

        regions4_bed = LTRends_flanking.LTRe_fDNA(output, output_path, blast2, minlenTE, maxlenTE)

        regions4_fasta = TEfromLTR.getfasta(genome, regions4_bed, output, output_path, "_regions")

        output_regions = LTRends_flanking.LTRe_fDNA_out(regions4_fasta, output, output_path)

        TSD.checkTSD(output_regions, output, output_path)
        
    else:
        
        print("No TEs have been found")
