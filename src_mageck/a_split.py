# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/cluster/mshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import seaborn as sns

# Default params
inp_dir = _config.DATA_DIR
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

##
# Functions
##
def divide():
  inp_fn = inp_dir + 'SHE2518.fq'

  # num_lines = util.line_count(inp_fn)
  batch_size = int(20e6)
  num_lines = 1138721480
  curr_num = 0
  commands = []
  for start in range(1, num_lines + 1, batch_size):
    end = start + batch_size

    out_fn = out_dir + '%s.fq' % (curr_num)
    command = 'tail -n +%s %s | head -n %s > %s' % (start, inp_fn, end - start, out_fn)

    curr_num += 1

    commands.append(command)

  print '\n'.join(commands)

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