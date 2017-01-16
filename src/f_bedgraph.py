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

def write_bulk(counts, out_place, gene):
  print '\t\tbulk'
  regions = []

  timer = util.Timer(total = len(counts))
  for i in range(len(counts)):
    for exps in [s for s in _config.SPLITS if gene in s]:
      exp = exps.replace(gene, '')
      # print exps, gene, exp
      if counts[exp][i] >= _config.GRNA_COUNT_THRESHOLD:
        regions.append( [counts['start'][i], counts['end'][i], 1] )
    timer.update()

  regions = combine_regions(regions)
  out_fn = out_place + gene + '_bulk.bg'
  with open(out_fn, 'w') as f:
    f.write(_lib.bg_header('bulk', gene) + '\n')
    for r in regions:
      f.write(_config.CHROMS[gene] + '\t' + '\t'.join([str(s) for s in r]) + '\n')
  return

def write_enrichment(counts, out_place, gene, exp, bulk_nm):
  print '\t\t', exp
  regions = []
  
  xp = exp.replace(gene, '')
  bnm = bulk_nm.replace(gene, '')

  # import pdb; pdb.set_trace()

  timer = util.Timer(total = len(counts))
  for i in range(len(counts)):
    if counts[xp][i] >= _config.GRNA_COUNT_THRESHOLD:
      enrich = float(counts[xp][i] + _config.PSEUDOCOUNT) / float(counts[bnm][i] + _config.PSEUDOCOUNT)
      regions.append( [counts['start'][i], counts['end'][i], enrich] )
    timer.update()

  regions = combine_avg_regions(regions)
  out_fn = out_place + gene + '_' + xp + '.bg'
  with open(out_fn, 'w') as f:
    f.write(_lib.bg_header(xp, gene) + '\n')
    for r in regions:
      if not _lib.is_control(r[0]):
        # if it is, it's a control, so don't write
        f.write(_config.CHROMS[gene] + '\t' + '\t'.join([str(s) for s in r]) + '\n')
  return

def combine_regions(regions):
  # 0-th element: start
  # 1-st element: end
  # 2-nd element: value (only combine regions w/ same value)
  # USE FOR BULK ONLY

  regions = sorted(regions, key = lambda x: x[0])
  # import pdb; pdb.set_trace()
  new_regions = []
  curr_region = None
  for i in range(len(regions)):
    if curr_region == None:
      curr_region = regions[i]
    else:
      if regions[i][2] == curr_region[2] and regions[i][0] <= curr_region[1] and regions[i][1] > curr_region[1]:
        # print curr_region, '\t', regions[i]
        curr_region[1] = regions[i][1]
      else:
        new_regions.append(curr_region)
        curr_region = regions[i]

  print '\t\tKept', len(new_regions), 'out of', len(regions), 'regions'
  return new_regions

def remove_control(regions):
  new_regions = []
  for r in regions:
    if not _lib.is_control(r[0]):
      new_regions.append(r)
  return new_regions

def combine_avg_regions(regions):
  # USE FOR MED/NEG only
  #
  # regions defined as a list of:
  #   [counts['start'][i], counts['end'][i], enrich]
  #
  # If 2 regions overlap, and they aren't the same enrichment value, then 3 regions will be produced.
  # In this way, it's not guaranteed that the number of output regions is less than the number of input regions.

  regions = remove_control(regions)
  regions = sorted(regions, key = lambda x: x[0])
  if len(regions) == 0:
    print '\tNo regions found'
    return []

  keys = range(int(regions[0][0]), int(regions[-1][1]) + 2)
  
  # b[i][0] stores enrichment total (at position i)
  # b[i][1] stores the number of regions that position i is involved in
  b = defaultdict(list)
  for k in keys:
    b[k] = [0, 0]
  for r in regions:
    for i in range(int(r[0]), int(r[1]) + 1):
      b[i][0] += float(r[2])
      b[i][1] += 1

  # Scan through all positions in b, combining regions that have the same enrichment, and separating the overlapping regions with different enrichment
  new_regions = []
  curr_start = 0
  curr_val = 0
  for i in range(int(regions[0][0]), int(regions[-1][1]) + 2):
    if curr_start == 0:
      if b[i][1] != 0:
        # Only track regions that have reads in them
        curr_start = i
        curr_val = b[i][0]
    else:
      if b[i][0] == curr_val:
        pass
      else:
        # Once we reach a position with a different enrichment total, split it out
        new_regions.append([curr_start, i - 1, float(curr_val) / float(b[i - 1][1])])
        if b[i][1] == 0:
          curr_start = 0
          curr_val = 0
          curr_num = 0
        else:
          curr_start = i
          curr_val = b[i][0]
  # if curr_start != 0:
    # new_regions.append([curr_start, i - 1, float(curr_val) / float(curr_num)])


  print '\t\tSeparated into', len(new_regions), 'regions from an original', len(regions), 'regions'
  return new_regions

def bedgraph(inp_fn, out_place, gene):
  counts = pd.read_csv(inp_fn)

  for grp in _config.d.GROUPS:
    if gene in grp[0]:
      bulk_nm, med_nm, neg_nm = grp[0], grp[1], grp[2]
      write_enrichment(counts, out_place, gene, med_nm, bulk_nm)
      write_enrichment(counts, out_place, gene, neg_nm, bulk_nm)
      write_bulk(counts, out_place, gene)
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
    bedgraph(inp_fn, out_place, gene)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
