# Summarize data into counts csv files

import _config, _lib
import sys, os, fnmatch, datetime, subprocess, math
import numpy as np
import pandas as pd
import scipy.stats

import seaborn as sns
import matplotlib.pyplot as plt

from mylib import util
from collections import defaultdict

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'c_counts/'
NAME = util.get_fn(__file__)

SIG_VALUE = 0.05

# Functions
def make_negative(counts, grp, out_fn):
  ratios = []
  for i in range(len(counts)):
    if _lib.is_neg_control( counts['start'][i] ):
      ratio = math.log( float(counts[grp[2]][i] + _config.PSEUDOCOUNT) / float(counts[grp[0]][i] + _config.PSEUDOCOUNT), 2)
      ratios.append( ratio )

  mean, std = np.mean(ratios), np.std(ratios)
  print '\t\tmean:', mean, '\tstd:', std

  filtered_data = [['', grp[0], grp[1], grp[2], 'start', 'end']]
  for i in range(len(counts)):
    r = math.log( float(counts[grp[2]][i] + _config.PSEUDOCOUNT) / float(counts[grp[0]][i] + _config.PSEUDOCOUNT), 2)
    cdf = scipy.stats.norm(mean, std).cdf(r)

    if 1 - cdf < SIG_VALUE:
      d = [ counts['Unnamed: 0'][i], counts[grp[0]][i], counts[grp[1]][i], counts[grp[2]][i], counts['start'][i], counts['end'][i] ]
      filtered_data.append(d)

  print '\t\tSignificant finds:', len(filtered_data), 'out of', len(counts)

  df = pd.DataFrame(filtered_data)
  df.to_csv(out_fn, header = False, index = False)

  return mean, std

def significance(inp_fn, out_place, gene):
  counts = pd.read_csv(inp_fn)

  our_groups = [s for s in _config.d.GROUPS if gene in s[0]]
  for group in our_groups:
    print '\t\t', group[2]
    trimmed_group = [s.replace(gene, '') for s in group]
    out_fn = out_place + group[2] + '_sig_counts.csv'
    make_negative(counts, trimmed_group, out_fn)
  
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
    inp_fn = inp_dir + gene + '_counts.csv'
    significance(inp_fn, out_place, gene)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
