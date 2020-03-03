# Summarize data into counts csv files
from __future__ import division
import _config, _lib
import sys, os, fnmatch, datetime, subprocess, math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from mylib import util
from collections import defaultdict
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'c_counts/'
NAME = util.get_fn(__file__)

# Functions
def plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn):
  log_fold = []
  tot_pos_col = sum(counts[pos_col])
  tot_neg_col = sum(counts[neg_col])
  for idx, row in counts.iterrows():
    normalized_pos = row[pos_col] / tot_pos_col
    normalized_neg = row[neg_col] / tot_neg_col
    pos_log = np.log(normalized_pos + 1e-6)
    neg_log = np.log(normalized_neg + 1e-6)
    log_fold.append( pos_log - neg_log )

  sorted_names = [x for _,x in sorted(zip(log_fold, counts[nm_col]))]
  log_fold = sorted(log_fold)
  xdata = range(len(log_fold))
  
  log_fold_dict = dict()
  for idx, log_score in enumerate(log_fold):
    log_fold_dict[sorted_names[idx]] = log_score

  data_pos_control = defaultdict(list)
  data_neg_control = defaultdict(list)
  data_rest = defaultdict(list)

  for idx, nm in enumerate(sorted_names):
    if 'pos_control' in nm or 'GFP_' in nm:
      data_pos_control['x'].append( xdata[idx] )
      data_pos_control['y'].append( log_fold[idx] )
    elif 'neg_control' in nm or 'Non-Targeting Control' in nm:
      data_neg_control['x'].append( xdata[idx] )
      data_neg_control['y'].append( log_fold[idx] )
    else:
      data_rest['x'].append( xdata[idx] )
      data_rest['y'].append( log_fold[idx] )

  print 'Plotting to %s' % (out_fn)
  with PdfPages(out_fn, 'w') as pdf:
    # Plot ranked log enrichment scatter plot
    plt.plot(data_rest['x'], data_rest['y'], '.k', ms = 5, alpha = 0.25)
    plt.plot(data_neg_control['x'], data_neg_control['y'], '.g')
    plt.plot(data_pos_control['x'], data_pos_control['y'], '.r')
    plt.xlabel('gRNAs')
    plt.title(out_fn.split('/')[-1].replace('.pdf', ''))
    plt.ylabel('Log-fold change (bulk over neg)')
    plt.xlim([-500, len(log_fold) + 500])
    red_patch = mpatches.Patch(color='red', label='GFP- Control')
    green_patch = mpatches.Patch(color='green', label='GFP+ Control')
    plt.legend(handles=[red_patch, green_patch], loc = 4)
    pdf.savefig()
    plt.clf()

    # Plot histogram of controls
    plt.subplot(2, 1, 1)
    sns.violinplot(data_pos_control['y'], color = 'r', cut = 0, bw = 0.2, inner = 'stick')
    plt.title('Log-fold enrichment of controls')
    plt.ylabel('GFP- Control')
    plt.xlim([min(log_fold), max(log_fold)])

    plt.subplot(2, 1, 2)
    sns.violinplot(data_neg_control['y'], color = 'g', cut = 0, bw = 0.2, inner = 'stick')
    plt.xlabel('Log-fold change (bulk over neg)')
    plt.ylabel('GFP+ Control')
    plt.xlim([min(log_fold), max(log_fold)])

    pdf.savefig()
    plt.clf()

  return log_fold_dict


def plot_reproducibility(r1d, r2d, out_fn):
  minmin = min(min(r1d.values()), min(r2d.values()))
  maxmax = max(max(r1d.values()), max(r2d.values()))
  
  data_pos_control = defaultdict(list)
  data_neg_control = defaultdict(list)
  data_rest = defaultdict(list)

  for nm in r1d:
    if 'pos_control' in nm or 'GFP_' in nm:
      data_pos_control['x'].append( r1d[nm] )
      data_pos_control['y'].append( r2d[nm] )
    elif 'neg_control' in nm or 'Non-Targeting Control' in nm:
      data_neg_control['x'].append( r1d[nm] )
      data_neg_control['y'].append( r2d[nm] )
    else:
      data_rest['x'].append( r1d[nm] )
      data_rest['y'].append( r2d[nm] )

  r1all = data_rest['x'] + data_pos_control['x'] + data_neg_control['x']
  r2all = data_rest['y'] + data_pos_control['y'] + data_neg_control['y']
  rsq = spearmanr(r1all, r2all)[0]**2

  print 'Plotting to %s' % (out_fn)
  plt.plot(data_rest['x'], data_rest['y'], '.k', alpha = 0.25)
  plt.plot(data_pos_control['x'], data_pos_control['y'], '.r')
  plt.plot(data_neg_control['x'], data_neg_control['y'], '.g')
  plt.plot( [minmin, minmin], [maxmax, maxmax], 'k-')
  plt.title(out_fn.split('/')[-1].replace('.pdf', '') + ' | Rsq: ' + str(rsq))
  plt.xlabel('Replicate 1: Log fold change (bulk/neg)')
  plt.ylabel('Replicate 2: Log fold change (bulk/neg)')
  red_patch = mpatches.Patch(color='red', label='GFP- Control')
  green_patch = mpatches.Patch(color='green', label='GFP+ Control')
  plt.legend(handles=[red_patch, green_patch], loc = 4)
  plt.savefig(out_fn)
  plt.clf()
  return


def plot_label_noise(counts, out_fn):

  pos_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep1_Bulk'
  neg_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep1_GFP_neg'
  log_fold_r1_u2os = []
  tot_pos_col = sum(counts[pos_col])
  tot_neg_col = sum(counts[neg_col])
  for idx, row in counts.iterrows():
    if row[pos_col] >= 10 and row[neg_col] >= 10:
      normalized_pos = row[pos_col] / tot_pos_col
      normalized_neg = row[neg_col] / tot_neg_col
      pos_log = np.log(normalized_pos + 1e-6)
      neg_log = np.log(normalized_neg + 1e-6)
      log_fold_r1_u2os.append( pos_log - neg_log )

  pos_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep2_Bulk'
  neg_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep2_GFP_neg'
  log_fold_r2_u2os = []
  tot_pos_col = sum(counts[pos_col])
  tot_neg_col = sum(counts[neg_col])
  for idx, row in counts.iterrows():
    if row[pos_col] >= 10 and row[neg_col] >= 10:
      normalized_pos = row[pos_col] / tot_pos_col
      normalized_neg = row[neg_col] / tot_neg_col
      pos_log = np.log(normalized_pos + 1e-6)
      neg_log = np.log(normalized_neg + 1e-6)
      log_fold_r2_u2os.append( pos_log - neg_log )

  pos_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep1_Bulk'
  neg_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep1_GFP_neg'
  log_fold_r1_dld1 = []
  tot_pos_col = sum(counts[pos_col])
  tot_neg_col = sum(counts[neg_col])
  for idx, row in counts.iterrows():
    if row[pos_col] >= 10 and row[neg_col] >= 10:
      normalized_pos = row[pos_col] / tot_pos_col
      normalized_neg = row[neg_col] / tot_neg_col
      pos_log = np.log(normalized_pos + 1e-6)
      neg_log = np.log(normalized_neg + 1e-6)
      log_fold_r1_dld1.append( pos_log - neg_log )

  pos_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep2_Bulk'
  neg_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep2_GFP_neg'
  log_fold_r2_dld1 = []
  tot_pos_col = sum(counts[pos_col])
  tot_neg_col = sum(counts[neg_col])
  for idx, row in counts.iterrows():
    if row[pos_col] >= 10 and row[neg_col] >= 10:
      normalized_pos = row[pos_col] / tot_pos_col
      normalized_neg = row[neg_col] / tot_neg_col
      pos_log = np.log(normalized_pos + 1e-6)
      neg_log = np.log(normalized_neg + 1e-6)
      log_fold_r2_dld1.append( pos_log - neg_log )


  tot_u2os = len(log_fold_r1_u2os) + len(log_fold_r2_u2os)
  tot_dld1 = len(log_fold_r1_dld1) + len(log_fold_r2_dld1)

  for threshold in [0.5, 1, 2]:
    print sum([-threshold <= s <= threshold for s in log_fold_r1_u2os + log_fold_r2_u2os]) / float(tot_u2os)
    print sum([-threshold <= s <= threshold for s in log_fold_r1_dld1 + log_fold_r2_dld1]) / float(tot_dld1)
    print threshold, '\n'

  import code; code.interact(local=dict(globals(), **locals()))
  return


def plotter(inp_dir, out_place):
  gene = 'Apobec3b'
  inp_fn = inp_dir + gene + '_counts.csv'
  counts = pd.read_csv(inp_fn)
  group = _config.d.GROUPS[_config.d.GENES.index(gene)]
  nm = gene
  nm_col = 'start'

  pos_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep1_Bulk'
  neg_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep1_GFP_neg'
  out_fn = out_place + 'Apobec3b_U2OS_rep1.pdf'
  r1d_u2os = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  pos_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep2_Bulk'
  neg_col = 'Apobec3b_U2OS-P2A-GFP_12k_rep2_GFP_neg'
  out_fn = out_place + 'Apobec3b_U2OS_rep2.pdf'
  r2d_u2os = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  out_fn = out_place + 'reproducibility_Apobec3b_U2OS.pdf'
  plot_reproducibility(r1d_u2os, r2d_u2os, out_fn)

  pos_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep1_Bulk'
  neg_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep1_GFP_neg'
  out_fn = out_place + 'Apobec3b_DLD1_rep1.pdf'
  r1d_dld1 = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  pos_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep2_Bulk'
  neg_col = 'Apobec3b_DLD1-P2A-GFP_12k_rep2_GFP_neg'
  out_fn = out_place + 'Apobec3b_DLD1_rep2.pdf'
  r2d_dld1 = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  out_fn = out_place + 'reproducibility_Apobec3b_DLD1.pdf'
  plot_reproducibility(r1d_dld1, r2d_dld1, out_fn)

  out_fn = out_place + 'label_noise_U2OS_vs_DLD1.pdf'
  plot_label_noise(counts, out_fn)

  ###################################################
  gene = 'Hg38tf'
  inp_fn = inp_dir + gene + '_counts.csv'
  counts = pd.read_csv(inp_fn)
  group = _config.d.GROUPS[_config.d.GENES.index(gene)]
  nm = gene
  nm_col = 'sequence info'

  pos_col = 'Apobec3b_U2OS-P2A-GFP_7k_rep1_Bulk'
  neg_col = 'Apobec3b_U2OS-P2A-GFP_7k_rep1_GFP_neg'
  out_fn = out_place + 'Apobec3b_TFs_U2OS_rep1.pdf'
  r1d = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  pos_col = 'Apobec3b_U2OS-P2A-GFP_7k_rep2_Bulk'
  neg_col = 'Apobec3b_U2OS-P2A-GFP_7k_rep2_GFP_neg'
  out_fn = out_place + 'Apobec3b_TFs_U2OS_rep2.pdf'
  r2d = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  out_fn = out_place + 'reproducibility_Apobec3b_TFs_U2OS.pdf'
  plot_reproducibility(r1d, r2d, out_fn)

  pos_col = 'Apobec3b_DLD1-P2A-GFP_7k_rep1_Bulk'
  neg_col = 'Apobec3b_DLD1-P2A-GFP_7k_rep1_GFP_neg'
  out_fn = out_place + 'Apobec3b_TFs_DLD1_rep1.pdf'
  r1d = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  pos_col = 'Apobec3b_DLD1-P2A-GFP_7k_rep2_Bulk'
  neg_col = 'Apobec3b_DLD1-P2A-GFP_7k_rep2_GFP_neg'
  out_fn = out_place + 'Apobec3b_TFs_DLD1_rep2.pdf'
  r2d = plot_enrichment(counts, pos_col, neg_col, nm_col, out_fn)

  out_fn = out_place + 'reproducibility_Apobec3b_TFs_DLD1.pdf'
  plot_reproducibility(r1d, r2d, out_fn)

  import code; code.interact(local=dict(globals(), **locals()))
  return


@util.time_dec
def main(inp_dir, out_place, run = True):
  print NAME  
  util.ensure_dir_exists(out_place)
  if not run:
    print '\tskipped'
    return out_place

  # Function calls
  plotter(inp_dir, out_place)

  return out_place


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
