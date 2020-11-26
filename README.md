# Working with Whole Genome Data

# Program
MATLAB

# Data Import

I used FASTA format data of **Chromosome 1** among human DNA data.
Prior to analysis, the gz-compressed file was extracted as a 'fa' file through the '7zip' program and saved it in the Matlab folder.

- Data Source: ensemblGenomes (http://ensemblgenomes.org/)


# Processing

## Memory Mapping (Create MEMMAPFILE Object)

Since the data of the Chromosome 1 file is very large, **a memory mapping process** was performed so that it can work without directly loading it into memory.
To simplify the indexing operation, the newline (newline) characters that existed in the original FASTA file were removed, and the data was written into a new file by processing the data as an 8-bit integer (unit8 class).
Then, a memmapfile object 'chr1' that maps a file to memory was constructed so that we can easily access to the desired data. 
Finally, the mapped file data was accessed by indexing operation and the sequence was printed by using int2nt or seqdisp.

# Data Analysis

## 1. Calculation of GC Content per block and plot result

As the first analysis, I calculated the relative ratios of G and C for each block and expressed the results as a line graph.
As a result, the GC ratio was relatively high in the region between 0 Mbp and 10 Mbp. Also, the graph was cut around 140 Mbp so we can suggest that other bases other than A, G, C, and T existed in this region.


## 2. Finding Regions of High GC Content

In the second analysis, I determined the areas with high relative proportions of G and C.
(*   I defined the region with high GC-content in Chromosome 1 as the section with a ratio of more than 0.5. )


## 3. Finding Regions of Low GC Content

In the third analysis, in contrast to the second analysis, I found the regions with low relative proportions of G and C.
(*   I defined the region with low GC-content in Chromosome 1 as the section with a ratio of less than 0.4. )

##  4. Search for sequences other than A, T, C, G

As a final analysis, I found the types of nucleotide sequences other than A, T, C, and G present in chromosome 1 and the ratio of the sequence content relative to the entire sequence. Finally, I expressed the results as a line graph.
It was found that the type of the sequences other than A, T, C, and G was N, and the relative ratio of 'N' was largely around 140 Mbp as analyzed in 'Analysis 1' above.
