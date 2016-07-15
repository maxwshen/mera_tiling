# counts the number of copies of each gRNA present in reads
#   Primer:
#     Allows 1bp insertion or deletion before gRNA
#     allows any mismatches before gRNA  
#   gRNA sequence
#     requires exact match to gRNA
#   

import _config, _lib
import sys, os, fnmatch, datetime, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO

from mylib import util


# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a2_condense/'
NAME = util.get_fn(__file__)

# Functions

def get_grnas(gene):
  grna_fn = _lib.get_grna_fn(gene)
  with open(grna_fn) as f:
    return [s.strip() for s in f.readlines()]
  print 'ERROR: No gRNA file for', gene
  return []

def find_grna(name, grnas, seq, df, ct):
  # Finds which gRNA perfectly matches the sequence
  # Assumes that a sequence can only have 1 gRNA (no overlap)
  #
  # Ensure gRNAs is a set, so 'in' is O(1)
  # 'in' is O(n) for lists
  num_found = 0
  
  if _config.REPRODUCE:
    # Nisha only considers one starting point
    starts = [21]
  else:
    starts = [21, 20, 22]
  lens = [21, 20, 19]
  seqs = []
  added_once = False
  for st in starts:
    for l in lens:
      s = seq[st : st + l]
      if s in grnas:
        if not added_once:
          df[name][s] += ct
          added_once = True
        num_found += 1
  if num_found > 1:
    print seq[5:], 'matched', num_found, 'gRNAs'

  if num_found > 0:
    return True
  return False


def count_grna(inp_dir, out_dir, gene):
  grnas = get_grnas(gene)
  if len(set(grnas)) != len(grnas):
    print '\tERROR: Duplicate gRNAs in file'

  exps = [nm.replace(gene, '') for nm in _config.SPLITS if gene in nm]
  df = pd.DataFrame(index = grnas, columns = exps, dtype = int)
  df[:][:] = 0

  grnas = set(grnas)
  for name in exps:
    print '\t\t', name
    reads_fn = inp_dir + gene + name + '/' + _config.d.FN

    num_reads_matched = 0
    r = SeqIO.parse(reads_fn, 'fasta')
    nl = util.line_count(reads_fn) / 2
    timer = util.Timer(total = nl, print_interval = 5)
    while True:
      try:
        rx = r.next()
        ct = int(rx.description.strip())
        if find_grna(name, grnas, str(rx.seq), df, ct):
          num_reads_matched += 1
        timer.update()
      except StopIteration:
        break

    pct_reads_matched = float(num_reads_matched) / nl
    print num_reads_matched, '/', nl
    print pct_reads_matched * 100, '%', ' unique reads matched'

  out_fn = out_dir + _config.d.FN
  with open(out_fn, 'w') as f:
    df.to_csv(f)
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
    out_dir = out_place + gene + '/'
    util.ensure_dir_exists(out_dir)
    count_grna(inp_dir, out_dir, gene)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
