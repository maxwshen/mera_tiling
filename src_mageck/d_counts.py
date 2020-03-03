# Summarize data into counts csv files

import _config, _lib
import sys, os, fnmatch, datetime, subprocess, math
import numpy as np
import pandas as pd

from mylib import util
from collections import defaultdict

# Default params
inp_dir = _config.OUT_PLACE + 'c2_combine_dfs/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'

PSEUDOCOUNT = 1


# Functions

def get_gene_locs(gene):
  starts = []
  ends = []
  headers = []
  if '12k' in gene:
    inp_fn = '/cluster/mshen/prj/mera_tiling/subprjs/2018-01-11/data/expdesign2/Apobec3b.gRNA.filt_2bp.library.locs.txt'
  elif '7sk' in gene:
    inp_fn = '/data/gl/g5/bh0085/mera_tiling/data/expdesign2/Hg38tf.gRNA.filt_2bp.library.locs.txt'
  else:
    inp_fn = '/cluster/mshen/prj/mera_tiling/subprjs/2018-01-11/data/expdesign2/DNArepair.gRNA.filt_2bp.library.locs.txt'
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if 'Hg38tf' not in inp_fn and 'DNArepair' not in inp_fn:
        words = line.split()
        starts.append(words[0])
        ends.append(words[1])
      else:
        headers.append( line.strip() )
  return starts, ends, headers

def process_genes(gene_nm):
  if 'GFP' in gene_nm:
    return 'GFP_Control'
  elif 'Control' in gene_nm:
    return 'Non_Targeting_Control'
  else:
    return gene_nm.split('/')[0]

def counts(inp_fn, gene):
  starts, ends, headers = get_gene_locs(gene)
  counts = pd.read_csv(inp_fn, index_col = 0)

  old_cols = counts.columns
  for col in counts.columns:
    print '%s: %s' % (col, sum(counts[col]))

  if '12k' in gene:
    counts['sgRNA'] = [int(s) if s.isdigit() else s for s in starts]
    counts['end'] = [int(s) if s.isdigit() else s for s in ends]
    counts = counts[['sgRNA', 'end'] + list(old_cols)]
  else:
    counts['sgRNA'] = [s.replace(' ', '_') for s in headers]
    counts['gene'] = [process_genes(s) for s in headers]
    counts = counts[['sgRNA', 'gene'] + list(old_cols)]

  if 'h7sk' not in gene:
    # handle technical replicates
    counts['%stot_Bulk_after_dSa_APA' % (gene)] = counts['%s-1_Bulk_after_dSa_APA' % (gene)] + counts['%s-2_Bulk_after_dSa_APA' % (gene)]
    counts['%stot_neg' % (gene)] = counts['%s-1_neg' % (gene)] + counts['%s-2_neg' % (gene)]
  else:
    counts['%stot_Bulk_after_h7sk-NFKB' % (gene)] = counts['%s_1Bulk_after_h7sk-NFKB' % (gene)] + counts['%s_2Bulk_after_h7sk-NFKB' % (gene)]
    counts['%stot_GFP_loss' % (gene)] = counts['%s_1GFP_loss' % (gene)] + counts['%s_2GFP_loss' % (gene)]

  counts.to_csv(out_dir + gene + '_counts.csv')
  counts.to_csv(out_dir + gene + '_counts.tsv', sep = '\t', index = False)

  return


@util.time_dec
def main():
  print NAME  
  util.ensure_dir_exists(out_dir)

  # Function calls
  for gene in _config.d.GENES:
  # for gene in ['DLD-1_APO-P2A-GFP_12k_rep3', 'DLD-1_APO-P2A-GFP_12k_rep4']:
    print '\t', gene
    inp_fn = inp_dir + '%s.csv' % (gene)
    counts(inp_fn, gene)

  return 


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main()
