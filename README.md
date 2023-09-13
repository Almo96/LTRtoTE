LTRtoTE
================

LTR retrotransposons are selfish genetic elements that are abundant in
most eukaryotic species. LTRtoTE is a software designed for utilizing an
LTR sequence to retrieve full-length transposable element (TE)
insertions. It is user-friendly, requiring only the LTR sequence and a
genome to operate. LTRtoTE will furnish information regarding the number
of transposable elements in the genome, their sequences, and the target
site duplications (TSD) generated by the insertion of these elements.

**Input files:**

- genome.fasta
- LTR.fasta

**Output files:**

- output_blast.txt blast results for the LTR in the genome
- output_blast_rev.txt blast result with the positions of the minus strand swapped, necessary for bedtools
- output_summary.txt TE start and stop positions and length
- output.bed all the TE positions: chr start stop
- output.fasta the sequence of every TE
- output.TElen.png histogram of the length distribution with bins of
  100bp
- output_regions.bed: chr start stop Region Emptycol Strand
- output_regions.fasta: sequence of the regions
- output_regions.txt: R1||R2...R3|R4......R5|R6...R7|R8
- output_tsd.txt: sequence of the target site duplications for each insertion
- output_tsd_len.txt: length of each TSD if present


Informations on the terminal are provided during the execution showing
the location of all the files generated and the number of TEs found with
the parameters provided.

How to run the code:

``` bash
python LTRtoTE.py path/to/LTR_T.fasta path/to/genome path/to/output output_name minlenTE maxlenTE
```

- **output_path** the path where all the output files will be saved
- **output_name** is the name the name for the outputs, NOT A PATH, on
  the terminal you will see where all the files will be saved.
- **minlenTE** & **maxlenTE** minimum and maximum length of the TE (it
  will discard larger or smaller TEs/fragments)
- **tsd** optional parameter, default is 20, it specifies the flanking
  region to check for the TSD, most TSD are 2 or 4 bp, but if you want
  you can increase the range

Example:

``` bash
python LTRtoTE.py ../genomes/LTR_T.fasta ../genomes/genome ../genomes/output Dmel732 4500 6000 --tsd 50
```

**Requirements:**

*Python modules:*

- pandas
- numpy
- matplotlib
- biopython

*Bioinformatic tools:*

- blast
- bedtools

**Notes:**

- I suggest you to put the two input files in a folder, all the output
  will be generated in the same folder
- I suggest you to call the LTR file LTR_TEname, where TEname is the
  name of your TE, this is not required, but the code will automatically
  remove LTR\_ and keep the TE name for some files.
- There are some other .bed and .fasta files generated during the
  process that contain information about the flanking regions, if you
  want them just delete the TSD.remove_temp() function from TSDtoTE.py.


**Regions:**
In the code and the output_regions.txt there are 8 regions. for a region
we refer to the following:

- R1: Genomic sequence flanking the LTR
- R2: Start sequence of the first LTR
- R3: End sequence of the first LTR
- R4: Coding region of the next to the first LTR
- R5: Coding region of the next to the second LTR
- R6: Start sequence of the second LTR
- R7: End sequence of the second LTR
- R8: Genomic sequence of the second LTR



Wow you read up to this point, ok if you are curious to know **how it
works:**

- Blast the LTR
- Create a db with the results and check for LTR in a distance range
  specified by the user
- Use the start and end of every LTR pair to create a bed file for the
  TE
- Use bedtools to get the sequence for every insertion and save it in a
  .fasta file
- Create a figure with the TE length distribution
- New bed files with -x bases on the start and +x bases on the left
  (don’t worry I also remove one base from the other side that would be
  coming from the LTR) in this way we have flanking reions
- Check for TSD in the flanking regions and save them in a db
- Generate a figure with the distribution of the TSD length
