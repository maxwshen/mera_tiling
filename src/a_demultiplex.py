# Separate into groups based on a barcode file

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np
from Bio import SeqIO

from mylib import util
from mylib import compbio


# Default params
DEFAULT_INP_DIR = _config.DATA_DIR + _config.DATA_FOLD
NAME = util.get_fn(__file__)


def match(seq):
  keys = _config.BARCODES
  for k in keys:
    s = seq[:len(k)]
    dists.append(sum([1 for i in range(len(k)) if s[i] != k[i]]))
  
  if min(dists) <= _config.BARCODE_MM:
    if dists.count(min(dists)) == 1:
      # import pdb; pdb.set_trace()
      return _config.SPLITS_FROM_BARCODE[keys[dists.index(min(dists))]]
    if dists.count(min(dists)) > 1:
      return 'other'
  return 'other'

def demultiplex(inp_fn, out_dir):
  for name in _config.SPLITS + ['other']:
    util.ensure_dir_exists(out_dir + name)
    util.exists_empty_fn(out_dir + name + '/' + _config.d.FN)

  lc = util.line_count(inp_fn) / 4
  r = SeqIO.parse(inp_fn, 'fastq')
  timer = util.Timer(total = lc)
  while True:
    try:
      rx = r.next()
      tag = match(str(rx.seq))
      if tag != 'other':
        bc_len = _config.BARCODE_FROM_SPLITS[tag]
      else:
        bc_len = 5
      
      with open(out_dir + tag + '/' +  _config.d.FN, 'a') as f:
        f.write('@' + rx.description + '\n' + str(rx.seq[bc_len:]) + '\n+\n' + compbio.SeqIO_fastq_qual_string(rx)[bc_len:] + '\n')
      timer.update(print_progress = True)
    except StopIteration:
      break
  
  return


@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  inp_fn = inp_dir + _config.d.FN 
  demultiplex(inp_fn, out_dir)
  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')