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
DEFAULT_INP_DIR = _config.OUT_PLACE + 'c_counts/'
NAME = util.get_fn(__file__)

# Functions

def plot_controls(counts, typ, tg, nm, out_fn):
  if typ == 'pos':
    ctrls = _config.POS_CONTROLS
    ttl = nm + ' positive (GFPneg) control'
    colr = 'red'
  elif typ == 'neg':
    ctrls = _config.NEG_CONTROLS
    ttl = nm + ' negative (wildtype) control'
    colr = 'green'
  else:
    print 'bad type of control:', typ
  xlbl = 'Log enrichment of bulk/GFP- (base 2)'

  neg = [s for s in tg if 'neg' in s.lower()][0]
  bulk = [s for s in tg if 'bulk' in s.lower()][0]

  data = []
  for i in range(len(counts)):
    # import pdb; pdb.set_trace()
    if _lib.is_pos_control(counts['start'][i]) and typ == 'pos':
      val = float(counts[bulk][i] + _config.PSEUDOCOUNT) / ( counts[neg][i] + _config.PSEUDOCOUNT )
      data.append( math.log(val, 2) )
    if _lib.is_neg_control(counts['start'][i]) and typ == 'neg':
      val = float(counts[bulk][i] + _config.PSEUDOCOUNT) / (counts[neg][i] + _config.PSEUDOCOUNT)
      data.append( math.log(val, 2) ) 
  
  if len(data) == 0:
    print 'No controls found - check _config.CONTROLS list to see if text is included'
    print _config.POS_CONTROLS, _config.NEG_CONTROLS
    return

  sns.distplot(data, color = colr, kde = False, bins = np.arange(min(data), max(data) + 0.25, 0.25))
  plt.title(ttl)
  plt.xlim(min(-0.5, min(data)), max(0.5, max(data)))
  plt.xlabel(xlbl)
  plt.ylabel('Quantity')
  plt.savefig(out_fn)
  plt.clf()
  return


def plotter(inp_fn, out_place, gene):
  counts = pd.read_csv(inp_fn)

  group = _config.d.GROUPS[_config.d.GENES.index(gene)]
  nm = gene

  pos_fn = out_place + nm + '_poscontrol.pdf'
  neg_fn = out_place + nm + '_negcontrol.pdf'
  plot_controls(counts, 'pos', group, nm, pos_fn)
  plot_controls(counts, 'neg', group, nm, neg_fn)

  # for key in _config.d.GROUPS:
  #   if gene in key:
  #     nm = key
  #     group = _config.d.GROUPS[key]

  #     pos_fn = out_place + nm + '_poscontrol.pdf'
  #     neg_fn = out_place + nm + '_negcontrol.pdf'
  #     plot_controls(counts, 'pos', group, nm, pos_fn)
  #     plot_controls(counts, 'neg', group, nm, neg_fn)

  # subprocess.call('pdfunite ' + pos_fn + ' ' + neg_fn + ' ' + out_place + 'controls.pdf')


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
    plotter(inp_fn, out_place, gene)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
