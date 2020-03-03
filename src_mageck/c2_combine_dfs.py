##

import _config, _lib
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/cluster/mshen/')
import numpy as np
import pandas as pd
from Bio import SeqIO

from mylib import util


# Default params
NAME = util.get_fn(__file__)
inp_place = _config.OUT_PLACE + 'c_countgrna/'
out_dir = _config.OUT_PLACE + NAME + '/'

# Functions
def combine_dfs(gene):
  inp_dir = inp_place + gene + '/'

  mdf = None
  num_dfs = 0

  for fn in os.listdir(inp_dir):
    if 'unmatched' not in fn:
      df = pd.read_csv(inp_dir + fn, index_col = 0)

      if mdf is None:
        mdf = df
      else:
        # cast to dict, add, then cast to df?
        mdf = mdf + df

      num_dfs += 1

  print 'Combined %s dfs' % (num_dfs)
  mdf.to_csv(out_dir + '%s.csv' % (gene))
  return

@util.time_dec
def main():
  print NAME  
  util.ensure_dir_exists(out_dir)

  # Function calls
  for gene in _config.d.GENES:
  # for gene in ['DLD-1_APO-P2A-GFP_12k_rep3', 'DLD-1_APO-P2A-GFP_12k_rep4',]:
    print '\t', gene
    combine_dfs(gene)

  return


if __name__ == '__main__':
  main()
