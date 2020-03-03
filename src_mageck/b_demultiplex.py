# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/cluster/mshen/')
import numpy as np
from collections import defaultdict
from mylib import util, compbio
import pandas as pd
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO


# Default params
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

##
# Functions
##
def match(seq, description):
  for (bc, index) in _config.BARCODES:
    bc_score = match_barcode(seq, bc)
    idx_score = match_index(description, index)

    if idx_score is True:
      if bc_score is None:
        return _config.SPLITS_FROM_BARCODE[(bc, index)]
      elif bc_score <= _config.BARCODE_MM:
        return _config.SPLITS_FROM_BARCODE[(bc, index)]
  return 'other'

def match_index(description, index):
  read_idx = description.split(':')[-1].strip()
  return bool(index == read_idx)

def match_barcode(seq, bc):
  if bc != 'NONE':
    s = seq[:len(bc)]
    score = sum([1 for i in range(len(bc)) if s[i] != bc[i]])
  else:
    score = None
  return score

##
# primary
##
def demultiplex(split):
  inp_fn = inp_dir + '%s.fq' % (split)
  for name in _config.SPLITS + ['other']:
    util.ensure_dir_exists(out_dir + name)
    util.exists_empty_fn(out_dir + name + '/%s.fa' % (split))

  lc = util.line_count(inp_fn) / 4
  r = SeqIO.parse(inp_fn, 'fastq')
  num_bad_q, num_tot = 0, 0
  timer = util.Timer(total = lc)
  while True:
    try:
      rx = r.next()
      tag = match(str(rx.seq), str(rx.description))
      if tag != 'other' and '4-C' not in tag and 'Dam' not in tag:
        bc_len = len(_config.BARCODE_FROM_SPLITS[tag])
      else:
        bc_len = 0
      
      qs = compbio.SeqIO_fastq_qual_string(rx)[bc_len:]
      quals = [ord(s)-33 for s in qs]
      if np.mean(quals) < 23:
        num_bad_q += 1
        continue

      num_tot += 1

      with open(out_dir + tag + '/%s.fa' % (split), 'a') as f:
        f.write('>' + rx.description + '\n' + str(rx.seq[bc_len:]) + '\n')
      timer.update()
    except StopIteration:
      break

  print 'Rejected %s fraction of reads' % (num_bad_q / num_tot)

  return

##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print 'Generating qsub scripts...'
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0
  for idx in range(0, 57):
    command = 'python %s.py %s' % (NAME, idx)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, idx)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -m e -wd %s %s' % (_config.SRC_DIR, sh_fn))

  # Save commands
  with open(qsubs_dir + '_commands.txt', 'w') as f:
    f.write('\n'.join(qsub_commands))

  print 'Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir)
  return

##
# Main
##
@util.time_dec
def main(split = ''):
  print NAME  

  if split == '':
    gen_qsubs()
    return

  demultiplex(split)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1])
  else:
    main()