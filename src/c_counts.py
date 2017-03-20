# Summarize data into counts csv files

import _config, _lib
import sys, os, fnmatch, datetime, subprocess, math
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from mylib import util
from collections import defaultdict

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'b_countgrna/'
NAME = util.get_fn(__file__)

PSEUDOCOUNT = 1


# Functions

def get_gene_locs(gene):
  starts = []
  ends = []
  with open(_lib.get_grna_locs_fn(gene)) as f:
    for i, line in enumerate(f):
      words = line.split()
      starts.append(words[0])
      ends.append(words[1])
  return starts, ends


def counts(inp_fn, out_place, gene):
  starts, ends = get_gene_locs(gene)
  counts = pd.DataFrame.from_csv(inp_fn)
  counts['start'] = [int(s) if s.isdigit() else s for s in starts]
  counts['end'] = [int(s) if s.isdigit() else s for s in ends]

  with open(out_place + gene + '_counts.csv', 'w') as f:
    counts.sort_values(by = 'start', ascending = False).to_csv(f)

  return


@util.time_dec
def main(inp_dir, out_place, run = True):
  print NAME  
  util.ensure_dir_exists(out_place)
  if not run:
    print '\tskipped'
    return out_place

  # Function calls
  for gene in _config.d.GENES:
    print '\t', gene
    inp_fn = inp_dir + gene +'/' + _config.d.FN
    counts(inp_fn, out_place, gene)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
