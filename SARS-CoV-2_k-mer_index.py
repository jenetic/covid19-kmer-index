"""
Process Reference Genome
"""

# Get sequence of reference genome
!git clone https://github.com/nageshsinghc4/COVID-19-coronavirus.git

# Reference genome
with open('COVID-19-coronavirus/cov2.fasta') as f:
  genome = ''.join([line[:-1] for line in f.readlines()[1:]])

"""Genes, amino acids, nucleotides"""

genes = [
  ['ORF1ab', 265, 21555],
  # ['ORF1ab', 265, 13468],
  # ['ORF1a', 265, 13483],
  ['S', 21562, 25384],
  ['ORF3a', 25392, 26220],
  ['E', 26244, 26472],
  ['M', 26522, 27191],
  ['ORF6', 27201, 27387],
  ['ORF7a', 27393, 27759], 
  ['ORF7b', 27755, 27887], # Overlap
  ['ORF8', 27893, 28259],
  ['N', 28273, 29533],
  ['ORF10', 29557, 29674]
]
ORF7b_start, ORF7b_end = 27755, 27887

aminoacids = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

nucleotides = ['A', 'C', 'G', 'T']

"""Wildtype k-mer count for reference genome (for unique 10/20-mers)"""

tenmer_count = {}
# Count wildtype 10-mers
for i in range(len(genome) - 10 + 1):
  seq = genome[i:i + 10]
  rev_comp = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
  kmer = min(seq, rev_comp)
  if kmer in tenmer_count:
    tenmer_count[kmer] += 1
  else:
    tenmer_count[kmer] = 1

twentymer_count = {}
# Count wildtype 20-mers
for i in range(len(genome) - 20 + 1):
  seq = genome[i:i + 20]
  rev_comp = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
  kmer = min(seq, rev_comp)
  if kmer in twentymer_count:
    twentymer_count[kmer] += 1
  else:
    twentymer_count[kmer] = 1
  
# Get wildtype 20-mers for ref genome (for identifying mutant k-mers in diff genomes)
wildtype_kmers = set()
k = 20

for i in range(len(genome) - k + 1):
  seq = genome[i:i + k]
  rev_comp = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
  kmer = min(seq, rev_comp)
  wildtype_kmers.add(kmer)

"""Mutation Table"""

# Returns list of possible mutant k-mers & adds to k-mer counter given pos in genome (i), k, and k-mer count dict
def find_kmers_count(i, k, kmer_count_dict):
  kmers = []
  for j in range(max(0, i-(k-1)), min(len(genome)-k,i)+1):
    seq = genome[j:i]+alt+genome[i+1:j+k]
    rev_comp = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
    kmer = min(seq, rev_comp)
    kmers.append(kmer)
    if kmer in kmer_count_dict:
      kmer_count_dict[kmer] += 1
    else:
      kmer_count_dict[kmer] = 1
  return kmers

# Same thing as above but without the counter (so k-mers at 13467-13468 won't duplicate in counter)
def find_kmers(i, k):
  kmers = []
  for j in range(max(0, i-(k-1)), min(len(genome)-k,i)+1):
    seq = genome[j:i]+alt+genome[i+1:j+k]
    rev_comp = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
    kmer = min(seq, rev_comp)
    kmers.append(kmer)
  return kmers

# Mutation Table
out = []

synon = 0
nonsynon = 0

overlap = False

for i in range(len(genome) + 4):
  
  # To make copy of 27755-27758 due to ORF7a/ORF7b overlap
  if i == 27759:
    overlap = True
  if overlap:
    i += -4
  
  coding = False

  for j, gene in enumerate(genes):
    gene_start, gene_end = gene[1], gene[2]  
    # Check if i is in a gene
    if gene_start <= i and i <= gene_end-1:
      coding = True
  
      # Get codon start and end pos
      if i >= 13468 and i <= 21555: # Frameshift 1 nucleotide back of reading frame btwn 13467 and 21555
        startpos = (gene_start+2) + ((i-(gene_start+2))//3)*3 # +2 to gene_start for frameshift
        endpos = startpos+3
        geneindex = j
      elif overlap and i <= ORF7b_end:
        startpos = (ORF7b_start) + ((i-(ORF7b_start))//3)*3
        endpos = startpos+3
      else:
        startpos = gene_start + ((i-gene_start)//3)*3
        endpos = startpos+3
        geneindex = j
      break     
  
  for alt in nucleotides:
    if genome[i] == alt:
      continue

    row = []
    row.append('chrV')
    row.append(str(i)) # Position in genome
    row.append(str(i+1)) # Position in genome
    row.append(genome[i]) # Nucleotide
    row.append(alt) # Mutant nucleotide

    # If i is in gene
    if coding:

      # Append reference amino acids & amino acids as result of mutation
      ref_codon = genome[startpos:endpos]
      mut_codon = genome[startpos:i] + alt + genome[i+1:endpos]
      row.append(aminoacids[ref_codon])
      row.append(aminoacids[mut_codon])

      # Gene name
      if overlap and i <= 27887:
        row.append('ORF7b')
      else:
        row.append(genes[geneindex][0]) 
      
      # Amino acid position
      if i >= 13468 and i <= 21555:
        row.append(((i+1)-genes[geneindex][1])//3 + 1) # Shifted back 1
      elif overlap and i <= ORF7b_end:
        row.append((i-ORF7b_start)//3 + 1)
      else:
        row.append((i-genes[geneindex][1])//3 + 1) # genes[geneindex][1] = start of gene
      
      # Counter of synonymous and nonsyonymous mutations 
      if aminoacids[ref_codon] == aminoacids[mut_codon]:
        synon += 1
      else:
        nonsynon += 1
    
    # if i isn't in gene
    else:
      row.append('NA')
      row.append('NA')
      row.append('NA')
      row.append('NA')

    # mutant 10-mers and 20-mers
    if overlap and i <= 27758: # So overlapping area btwn 27755-27758 not duplicated in counter
      row.append(','.join(find_kmers(i, 10))) # row[9]
      row.append(','.join(find_kmers(i, 20,))) # row[10]
    else:  
      row.append(','.join(find_kmers_count(i, 10, tenmer_count))) # row[9]
      row.append(','.join(find_kmers_count(i, 20, twentymer_count))) # row[10]

    # Add row to output
    out.append(row)

  # Makes copy of nucleotide 13467-13468 due to ribosomal frameshift
  if i == 13467:
    for alt in nucleotides:
      if genome[i] != alt:
        row = []
        row.append('chrV')
        row.append(str(i))
        row.append(str(i+1))
        row.append(genome[i])
        row.append(alt)
        
        startpos = (gene_start+2) + ((i-(gene_start+2))//3)*3 # +2 to gene_start for frameshift
        endpos = startpos+3 
        ref_codon = genome[startpos:endpos]
        mut_codon = genome[startpos:i] + alt + genome[i+1:endpos]
        row.append(aminoacids[ref_codon])
        row.append(aminoacids[mut_codon])

        row.append(genes[geneindex][0]) # Gene name
        row.append(((i+1)-genes[geneindex][1])//3+1) # Shifted amino acid position
        row.append(','.join(find_kmers(i, 10)))
        row.append(','.join(find_kmers(i, 20)))
        out.append(row)
  
for row in out:
  # Unique 10-mers (among wildtype and mutant k-mers combined)
  unique_tenmers = []
  kmers10 = row[9].split(',')
  for kmer in kmers10:
    if tenmer_count[kmer] == 1:
      unique_tenmers.append(kmer)
  row.append(','.join(unique_tenmers)) # row[11]

  # Unique 20-mers (among wildtype and mutant k-mers combined)
  unique_twentymers = []
  kmers20 = row[10].split(',')
  for kmer in kmers20:
    if twentymer_count[kmer] == 1:
      unique_twentymers.append(kmer)
  row.append(','.join(unique_twentymers)) # row[12]

"""Export mutation table to CSV"""

import csv

with open('output.csv', 'w') as f:
  writer = csv.writer(f)
  writer.writerows(out)