# Condense into fasta with unique lines. Counts are in header

import _config, _lib
import sys, os, fnmatch, datetime, subprocess
import numpy as np, pandas as pd
from Bio import SeqIO
import seaborn as sns
from mylib import util
import matplotlib.pyplot as plt

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a_demultiplex/'
NAME = util.get_fn(__file__)

# Functions

def filter_quality(quals):
  le = 22
  st = 20
  if np.mean(quals[st : st + le]) <= _config.QUALITY_CUTOFF:
    return False
  return True

def condense(inp_fn, out_dir):
  lc = util.line_count(inp_fn) / 4
  r = SeqIO.parse(inp_fn, 'fastq')
  timer = util.Timer(total = lc)

  qs = [dict() for x in range(100)]
  reads = dict()
  rejected = 0
  while True:
    try:
      rx = r.next()
      s = str(rx.seq)
      quals = rx.letter_annotations['phred_quality']

      for i in range(len(quals)):
        if quals[i] not in qs[i]:
          qs[i][quals[i]] = 0
        qs[i][quals[i]] += 1

      if filter_quality(quals):
        if s not in reads:
          reads[s] = 0
        reads[s] += 1
      else:
        rejected += 1
      timer.update(print_progress = True)
    except StopIteration:
      break

  if lc > 0:
    print '\tQuality: Rejected', rejected, 'out of', lc, float(rejected) / lc, 'percent'

  out_fn = out_dir + _config.d.FN
  with open(out_fn, 'w') as f:
    for key in reads:
      f.write('>' + str(reads[key]) + '\n' + key + '\n')

  print '\tUnique: Kept', util.line_count(out_fn), 'out of', lc * 2, float(util.line_count(out_fn) * 100) / lc * 2, 'percent'

  out_fig_fn = out_dir + 'quals.pdf'
  df = pd.DataFrame(qs)
  df = df.T
  df = df / df.apply(np.nansum)   # normalize
  df.to_csv(out_dir + 'quals.csv')

  sns.set_style('whitegrid')
  ax = sns.heatmap(data = df, vmin = 0, vmax = 1)
  ax.invert_yaxis()
  plt.ylabel('PHRED Quality Score')
  plt.xlabel('Positions')
  plt.title('Sequencing Quality')
  plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=8)
  plt.savefig(out_fig_fn)
  plt.clf()
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
    exps = [nm.replace(gene, '') for nm in _config.SPLITS if gene in nm]
    for exp in exps:
      print '\t', gene, exp
      inp_fn = inp_dir + gene + exp + '/' + _config.d.FN
      out_dir = out_place + gene + exp + '/'
      util.ensure_dir_exists(out_dir)
      condense(inp_fn, out_dir)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
