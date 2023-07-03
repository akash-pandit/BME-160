#### BME 160 Final Project
## ORF-Based Genome Comparison Tool

This is a command line tool designed to run on a linux kernel that compares a genome file to 2 other genome files, determining which of the 2 genome files the first file is closest to. This 'closeness' is determined by comparing each ORF (Open Reading Frame) in each file, and whichever of the latter 2 files has more matching ORFs with the first file is considered the file the first file is closer to.

Codon: any collection of 3 nucleotides. Examples: ATG, TAA, CTC, TAG, CAT, TGC, AAA

Open Reading Frames: a span of codons between a start codon and end codon (ATG...TAA). Example: ATGGCCTGTTAA.

Start codons used: ATG; Stop codons used: TAA, TGA, TAG

#### compare.py
The main file, defines classes for comparing the files themselves & for reading from the command line with argparse

#### sequenceAnalysis.py
A collection of modules built in prior labs for this course, Lab 5's ORFfinder class & the professor provided FastAReader class are imported into this project

#### testTarget/C1/C2.fa
randomly generated sequences of nucleotides (2 million characters per) used for testing due to relatively small size compared to other genome files (E. Coli genomes are 4.6 million bases, Mouse genomes are 2.7 billion bases, Human genome is 3 billion bases)


Target File - the file to be compared

C1 File - the first of the 2 comparison files

C2 File - the second of the 2 comparison files

## Usage
This program is ideally used to compare files that are already very similar genetically. The -f, -c1, and -c2 flags are all required, while the -g flag is optional and will display your results in a matplotlib figure. This is how one would execute the program:

```py compare.py -f [target file] -c1 [first comparison file] -c2 [second comparison file] -g```

Make sure your compare.py and sequenceAnalysis.py files are in the same directory. If in the same directory, genome files can be accessed with just the file name, otherwise the full file path should be provided.

If you make any errors when executing the script, user friendly error handling will tell you clearly what you may have done wrong.

