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
DEFAULT_INP_DIR = _config.OUT_PLACE + 'b_countgrna/'
NAME = util.get_fn(__file__)

def combine(inp_dir):
  num_splits = 16

  for gene in _config.d.GENES:
    print gene
    cs = pd.read_csv(inp_dir + 'split0' + '/' + gene + '/' + _config.d.FN)
    for s in range(1, num_splits):
      cs += pd.read_csv(inp_dir + 'split' + str(s) + '/' + gene + '/' + _config.d.FN)

    util.ensure_dir_exists(inp_dir + gene + '/')
    cs.to_csv(inp_dir + gene + '/' + _config.d.FN, index = False)

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