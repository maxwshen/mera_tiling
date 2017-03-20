# Check if gRNAs are present in other barcodes
# I wrote zb_otherbarcodes.py to search for the occurence of positive control (GFP) gRNAs in all barcodes, including others. This script uses grep to find occurences of each of the 77 control gRNAs in the "other" file that doesn't match the 12 expected brca2/sdhd barcodes. I then use the remaining expected barcodes to categorize all the hits. I find that nearly all positive control gRNAs are in barcode 25 / Y / GATGCGTA. 

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

from mylib import util
from mylib import compbio

# Default params
NAME = util.get_fn(__file__)

otherbc = ['ATCGATCGTTC', 'ATCGATGAGAG', 'ATCGATGGTAA', 'ATCGATTTCCT', 'TCGATACGTA', 'TCGATCAGGT', 'TCGATCTCTG', 'TCGATGTCTG', 'CGATAACCG', 'CGATCCTGT', 'CGATTGAGC', 'GATATCCG', 'GATGCGTA', 'GATTCTGC', 'ATAGTTA', 'ATCCAAT', 'ATGGACT', 'TCTGGA', 'TTCGTC', 'TTAGTG', 'GCTGT', 'TGGAA', 'AACTC', 'AAGCT']


def get_gfp_ctrls():
  num_gfp_control = 77
  with open('/cluster/mshen/prj/mera_tiling/data/expdesign2/Brca2_gRNA.txt') as f:
    txt = f.readlines()
  grnas = [s.strip() for s in txt[-num_gfp_control:]]
  print 'Found', len(grnas), 'gRNAs'
  return grnas

def grep_match(grna):
  # fn = '\' /cluster/mshen/prj/mera_tiling/out/2016-08-11/a_demultiplex/other/SHE1245_AHF7FYAFXX_S1_L001_R1_001.fastq'
  fn = '\' /cluster/mshen/prj/mera_tiling/out/' + _config.DATA_FOLD + 'a_demultiplex/other/' + _config.d.FN
  try:
    matches = subprocess.check_output('grep \'' + grna + fn, shell = True)
  except:
    matches = ''

  counts = np.zeros((1, len(otherbc)))
  m = [s.strip() for s in matches.split()]
  for seq in m:
    for i in range(len(otherbc)):
      bc = otherbc[i]
      if seq[:len(bc)] == bc:
        counts[0][i] += 1
        break

  return counts

def grep_match_expected(grna):
  sdhdm = subprocess.check_output('grep \'' + grna + '\' /cluster/mshen/prj/mera_tiling/out/' + _config.DATA_FOLD + 'a_demultiplex/sdhd*/' + _config.d.FN + ' | wc -l', shell = True)
  brcam = subprocess.check_output('grep \'' + grna + '\' /cluster/mshen/prj/mera_tiling/out/' + _config.DATA_FOLD + 'a_demultiplex/brca*/' + _config.d.FN + ' | wc -l', shell = True)
  return int(brcam), int(sdhdm)

def zb_otherbarcodes():
  g = get_gfp_ctrls()
  
  total_counts = np.zeros((1, len(otherbc)))
  brca_total, sdhd_total = 0, 0
  timer = util.Timer(total = len(g))
  for grna in g:
    # counts = grep_match(grna)
    # print grna, counts
    # total_counts += counts
    # print total_counts

    brca, sdhd = grep_match_expected(grna)
    brca_total += brca
    sdhd_total += sdhd
    print brca_total, sdhd_total

    timer.update()

  d = {}
  for i in range(len(otherbc)):
    d[otherbc[i]] = total_counts[0][i]

  order = ['TCGATACGTA', 'ATCGATCGTTC', 'AACTC', 'TCTGGA', 'ATGGACT', 'GATTCTGC', 'CGATAACCG', 'TCGATCTCTG', 'ATCGATGGTAA', 'AAGCT', 'TTCGTC', 'ATCCAAT', 'GATGCGTA', 'CGATTGAGC', 'TCGATGTCTG', 'ATCGATTTCCT', 'GCTGT', 'TTAGTG', 'ATAGTTA', 'GATATCCG', 'CGATCCTGT', 'TCGATCAGGT', 'ATCGATGAGAG', 'TGGAA']
  print brca_total, '\n', sdhd_total
  for o in order:
    print o, d[o]

  return

if __name__ == '__main__':
  zb_otherbarcodes()