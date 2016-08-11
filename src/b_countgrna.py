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
    print '\tWarning: Settings are for Nisha reproduction - only exact matches to primer are allowed, undercounting by 1.5-2%'
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


###############################################
################ MISCELLANEOUS ################
###############################################



def find_grna_hamming(name, grnas, seq, df, ct):
  # Exploring the question: Why only 5% of reads match gRNAs?
  # ANS: 5% ~ HDR rate
  s = seq[21 : 21 + 21]
  hammings = []
  grnas_list = list(grnas)
  for g in grnas_list:
    if len(g) == 21:
      hammings.append(sum([1 for i in range(len(s)) if s[i] != g[i]]))
    else:
      hammings.append(100)

  # print min(hammings), ' ',
  if min(hammings) > 6:
    print seq[21: 21 + 21], '\n', grnas_list[hammings.index(min(hammings))], min(hammings)
  if min(hammings) <= 1:
    return True
  return False

def align_grna(seq, grnas):
  # Exploring the question: Why only 5% of reads match gRNAs?
  # ANS: 5% ~ HDR rate
  def score_alignment(ans):
    w = ans.split()
    score = sum([1 for i in range(len(w[0])) if w[0][i] != w[1][i]])
    return score

  sa_tool = '/cluster/mshen/tools/seq-align/bin/needleman_wunsch'
  sa_options = '--match 10 --mismatch -8 --gapopen -20 --gapextend -1'

  aligns = []
  scores = []
  timer = util.Timer(total = len(grnas), print_interval = 20)
  for g in grnas:
    ans = subprocess.check_output(' '.join([sa_tool, 
          sa_options, 
          seq[21:21+21], 
          g]),
          shell = True)
    aligns.append(ans)
    scores.append(score_alignment(ans))
    timer.update()

  print min(scores), '\n', aligns[scores.index(min(scores))]
  return
