# Finds differential enrichment across 2050 TFs for mouseORF

import _config, _lib
import sys, os, fnmatch, datetime, subprocess, math
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from mylib import util
from collections import defaultdict

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'c_counts/'
NAME = util.get_fn(__file__)

# Functions

def mORFenrich(counts):
  with open(_config.DATA_DIR + 'expdesign2/mouseORFs_unique.txt') as f:
    keys = [s.strip() for s in f.readlines()]

  for exp in _config.d.g_morf_exps:
    print exp
    datas = []
    for cond in ['pos', 'neg']:
      nm = exp + cond
      data = defaultdict(list)
      for i in range(len(counts['start'])):
        tf = counts['start'][i].split('_')[0]
        data[tf].append(counts[nm][i] + _config.PSEUDOCOUNT)
      datas.append(data)

    pos_data, neg_data = datas[0], datas[1]
    pos_means, neg_means = {}, {}
    enrichment = {}

    for key in pos_data:
      enrichment[key] = abs(np.mean( np.array(pos_data[key]) / np.array(neg_data[key]) ))

    diff_TFs = sorted(enrichment, key = enrichment.get, reverse = True)

    import code; code.interact(local=dict(globals(), **locals()))

  return

@util.time_dec
def main(inp_dir, out_place, run = True):
  print NAME  
  util.ensure_dir_exists(out_place)
  if not run:
    print '\tskipped'
    return out_place

  # Function calls
  if 'mouseORF' not in _config.d.GENES:
    print 'Error: MouseORF not an experiment.'
    sys.exit(0)

  counts = pd.read_csv(inp_dir + 'mouseORF_counts.csv')
  mORFenrich(counts)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
