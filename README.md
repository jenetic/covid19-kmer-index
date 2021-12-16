# COVID19-kmer-index
By: Jenny Lam, Riley Kong, and Andrew Hong
Mentor: Hojoon Lee PhD

Using Python, we created a comprehensive k-mer (DNA substrings with k nucleotides) index of the SARS-CoV-2 genome that contains all possible substition mutations along with their corresponding unique mutant 10-mers and 20-mers. Our mutation table uses the SARS-CoV-2 reference genome (accession number: NC_045512.2) from the National Center for Biotechnology Information. For each position in the reference genome, we applied all possible subsititutions to generate an altered genomic sequence, which was used to generate mutant animo acids, 10-mers, and 20-mers. We made sure to account for the ORF7a/ORF7b overlap, as well as the ribosomal frameshift at C13468 in the SARS-CoV-2 genome. Our k-mer index allows for easy lookup of all possible SARS-CoV-2 mutations, as well as identification of mutations in mutated genomes by matching up the unique 20-mers, as every possible mutation has a unique 20-mer associated with it.

The Python script generates a .csv file which contains the following columns:
- **Start**: start position of chromosome
  - the table uses a 0-based coordinate system, which counts between nucleotides
- **End**: end position of chromosome
- **REF**: original reference nucleotide at that position
- **ALT: altered nucleotide
  - For each reference nucleotide, there are 3 altered nucleotides (3 rows per reference nucleotide) to cover every possible nucleotide substitution. 
- **WT AA (wildtype animo acid)**: amino acid that results from the nucleotides of the reference genome
- **Mut AA (mutant animo acid)**: amino acid that results from nucleotide substitition described in ALT column
- **Gene**: gene corresponding to chromosome position
- **AA pos (amino acid position)**: position number of animo acid
- **10-mers**: list of all possible 10-mers containing mutation
- **20-mers**: list of all possible 20-mers containing mutation
- **Unique 10-mers**: list of unique 10-mers containing mutation (10-mers that don't appear in any other mutation)
- **Unique 20-mers**: list of unique 20-mers containing mutation (20-mers that don't appear in any other mutation)
