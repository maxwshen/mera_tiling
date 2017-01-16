# If parallelize is used for a_demultiplex, use this combine function to gather split output together

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

from mylib import util
from mylib import compbio


# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a2_condense/'
NAME = util.get_fn(__file__)

def combine(inp_dir):
  num_splits = 16

  for nm in _config.SPLITS:
    print '\tCombining', nm, '...'
    for rfn in [_config.d.FN]:
      locs = [inp_dir + '/' + 'split' + str(s) + '/' + nm + '/' + rfn for s in range(num_splits)]
      util.ensure_dir_exists(inp_dir + '/' + nm)
      out_fn = inp_dir + '/' + nm + '/' + rfn
      util.exists_empty_fn(out_fn)
      subprocess.call('cat ' + ' '.join(locs) + ' >> ' + out_fn, 
        shell = True)
      print '\t\toutlines:', util.line_count(out_fn)

  return


def main(inp_dir, run = True):
  print NAME  
  if not run:
    print '\tskipped'
    return

  combine(inp_dir)
  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR)