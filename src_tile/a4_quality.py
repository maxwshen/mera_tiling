# Condense into fasta with unique lines. Counts are in header

import _config, _lib
import sys, os, fnmatch, datetime, subprocess
import numpy as np, pandas as pd
from Bio import SeqIO
import seaborn as sns
from mylib import util
import matplotlib.pyplot as plt

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a2_condense/'
NAME = util.get_fn(__file__)
NUM_SPLITS = 16

# Functions

def quality(inp_dir, out_dir, gene, exp):
  q = pd.read_csv(inp_dir + 'split0/' + gene + exp + '/quals.csv', index_col = 0)
  for i in range(1, NUM_SPLITS):
    q += pd.read_csv(inp_dir + 'split' + str(i) + '/' + gene + exp + '/quals.csv', index_col = 0)

  q /= NUM_SPLITS
  df = q
  df.to_csv(out_dir + 'quals.csv')
  sns.set_style('whitegrid')
  ax = sns.heatmap(data = df, vmin = 0, vmax = 1)
  ax.invert_yaxis()
  plt.ylabel('PHRED Quality Score')
  plt.xlabel('Positions')
  plt.title(gene + exp + ' | Sequencing Quality')
  plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=8)
  plt.savefig(out_dir + '/quals.pdf')
  plt.clf()  

  return

def pdfunite(out_dir):
  print 'Joining pdfs...'
  s = ['pdfunite']
  for exp in _config.SPLITS:
    s.append( out_dir + exp + '/quals.pdf' )
  s.append(out_dir + 'quals.pdf')
  subprocess.call(' '.join(s), shell = True)
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
      out_dir = out_place + gene + exp + '/'
      util.ensure_dir_exists(out_dir)
      quality(inp_dir, out_dir, gene, exp)
  pdfunite(out_place)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
