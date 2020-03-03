# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/cluster/mshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd

# Default params
inp_dir = _config.DATA_DIR
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

##
# Functions
##
def divide():
  for fn in os.listdir(_config.OUT_PLACE + 'b_demultiplex/'):
    # if fn[:4] != 'DLD1':
      # continue

    if fn not in ['U2OS_APO-P2A-GFP_7k_rep3-1_neg.fa', 'U2OS_APO-P2A-GFP_7k_rep4_Bulk_before_Cas9.fa', 'DLD-1_APO-P2A-GFP_7k_rep3-1_Bulk_after_Cas9.fa', 'U2OS7k+h7sk_NFKBpool_rep3_1GFP_loss.fa', 'U2OS7k+h7sk_NFKBpool_rep3_2GFP_loss.fa']:
      continue

    if '.fa' not in fn:
      continue

    inp_fn = _config.OUT_PLACE + 'b_demultiplex/%s' % (fn)
    inp_fn_numlines = util.line_count(inp_fn)

    num_splits = 60
    split_size = int(inp_fn_numlines / num_splits)
    if num_splits * split_size < inp_fn_numlines:
      split_size += 1
    while split_size % 2 != 0:
      split_size += 1
    # print 'Using split size %s' % (split_size)

    split_num = 0
    for idx in range(1, inp_fn_numlines, split_size):
      start = idx
      end = start + split_size  
      out_fn = _config.OUT_PLACE + 'b_demultiplex/%s/%s.fa' % (fn.replace('.fa', ''), split_num)
      command = 'tail -n +%s %s | head -n %s > %s' % (start, inp_fn, end - start, out_fn)
      split_num += 1
      print command

  return

##
# Main
##
@util.time_dec
def main():
  print NAME  
  
  # Function calls
  divide()

  return


if __name__ == '__main__':
  main()