# Compare the sequenced library to the expected gRNAs
# Conclusion: gRNAs were synthesized correctly, any gRNAs found to be "not-present" may be due to Illumina sequencing errors

import _config, _lib

def get_library_seq_fn(gene):
  return '/cluster/nishar/mm10_fasta/ES_TF_data/MERA_new_round/gRNA_' + _lib.get_name(gene) + '_library.fastq'

genes = ['msh2', 'sdhd', 'brca2']

for gene in genes:
  print gene
  with open(get_library_seq_fn(gene)) as f:
    lbs = [s.strip() for s in f.readlines()]

  with open(_lib.get_grna_fn(gene)) as f:
    grnas = [s.strip() for s in f.readlines()]
  grnas = set(grnas)
  counts = dict.fromkeys(grnas, 0)

  starts = [26, 25, 27]
  lens = [21, 20, 19]
  for l in lbs:
    for st in starts:
      for le in lens:
        s = l[st : st + le]
        if s in grnas:
          counts[s] += 1

  for i in range(100):
    print i, sum([1 for key in counts if counts[key] == i])