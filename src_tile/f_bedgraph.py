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
      if gene == 'Apobec3b' and '7k' in exps:
        continue
      # exp = exps.replace(gene, '')
      exp = exps
      # print exps, gene, exp
      if counts[exp][i] >= _config.GRNA_COUNT_THRESHOLD:
        s, e = reorder_start_end(counts, i)
        if 'control' not in s:
          regions.append( [s, e, 1] )
    timer.update()

  regions = combine_regions(regions)
  out_fn = out_place + gene + '_bulk.bg'
  with open(out_fn, 'w') as f:
    f.write(_lib.bg_header('bulk', gene) + '\n')
    for r in regions:
      f.write(_config.CHROMS[gene.lower()] + '\t' + '\t'.join([str(s) for s in r]) + '\n')
  return

def reorder_start_end(counts, i):
  s, e = counts['start'][i], counts['end'][i]
  if s > e:
    s, e = e, s
  return s, e

def write_combined_sum(counts, out_place, gene, grp, nm):
  print '\t\tcombinedsum', nm
  
  combined_counts = counts[grp[0]]
  for i in range(1, len(grp)):
    combined_counts += counts[grp[i]]

  regions = []
  timer = util.Timer(total = len(combined_counts))
  for i in range(len(combined_counts)):
    if combined_counts[i] >= _config.GRNA_COUNT_THRESHOLD:
      ct = float(combined_counts[i])
      s, e = reorder_start_end(counts, i)
      regions.append( [s, e, ct] )
    timer.update()

  # regions = combine_avg_regions(regions)
  out_fn = out_place + gene + '_' + nm + '_sum.bg'
  with open(out_fn, 'w') as f:
    f.write(_lib.bg_header(nm, gene) + '\n')
    for r in regions:
      if not _lib.is_control(r[0]):
        # if it is, it's a control, so don't write
        f.write(_config.CHROMS[gene.lower()] + '\t' + '\t'.join([str(s) for s in r]) + '\n')
  return

def write_enrichment(counts, out_place, gene, exp, bulk_nm):
  print '\t\t', exp
  regions = []
  
  # xp = exp.replace(gene, '')
  # bnm = bulk_nm.replace(gene, '')
  xp = exp
  bnm = bulk_nm

  timer = util.Timer(total = len(counts))
  for i in range(len(counts)):
    if counts[xp][i] >= _config.GRNA_COUNT_THRESHOLD:
      enrich = float(counts[xp][i] + _config.PSEUDOCOUNT) / float(counts[bnm][i] + _config.PSEUDOCOUNT)
      s, e = reorder_start_end(counts, i)
      regions.append( [s, e, enrich] )
    timer.update()

  # regions = combine_avg_regions(regions)
  out_fn = out_place + gene + '_' + xp + '.bg'
  with open(out_fn, 'w') as f:
    f.write(_lib.bg_header(xp, gene) + '\n')
    for r in regions:
      if not _lib.is_control(r[0]):
        # if it is, it's a control, so don't write
        f.write(_config.CHROMS[gene.lower()] + '\t' + '\t'.join([str(s) for s in r]) + '\n')
  return

def combine_regions(regions):
  # 0-th element: start
  # 1-st element: end
  # 2-nd element: value (only combine regions w/ same value)
  # USE FOR BULK ONLY

  regions = sorted(regions, key = lambda x: x[0])
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

def assign_names(group):
  bulk_nm, med_nm, pos_nm, neg_nm = '', '', '', ''
  for g in group:
    if 'bulk' in g or 'Bulk' in g:
      bulk_nm = g
    if 'med' in g:
      med_nm = g
    if 'pos' in g:
      pos_nm = g
    if 'neg' in g:
      neg_nm = g
  return bulk_nm, pos_nm, med_nm, neg_nm

def bedgraph(inp_fn, out_place, gene):
  counts = pd.read_csv(inp_fn)

  for nm in _config.d.F_COMBINED.keys():
    if gene in nm:
      write_combined_sum(counts, out_place, gene, _config.d.F_COMBINED[nm], nm)

  for grp in _config.d.F_BEDGRAPH_GROUPS:
    if gene in grp[0]:
      bulk_nm, pos_nm, med_nm, neg_nm = assign_names(grp)

      if bulk_nm == '':
        print 'ERROR: No bulk name in ', group

      if pos_nm != '':
        write_enrichment(counts, out_place, gene, pos_nm, bulk_nm)
      if neg_nm != '':
        write_enrichment(counts, out_place, gene, neg_nm, bulk_nm)
      if med_nm != '':
        write_enrichment(counts, out_place, gene, med_nm, bulk_nm)
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
  # for gene in _config.d.GENES:
  for gene in ['Apobec3b']:
    print '\t', gene
    inp_fn = inp_dir + gene + '_counts.csv'
    bedgraph(inp_fn, out_place, gene)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
