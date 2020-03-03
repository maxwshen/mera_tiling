# Dumb parallelization
# Designed for scripts that take a directory as input.
#
# We split up a directory into many subdirectories, then launch 
# a script on each subdirectory

import _config, _parallel_config
import sys, os, fnmatch, importlib, subprocess


from mylib import util

NAME = util.get_fn(__file__)

P_SCRIPT = 'a2_condense'
SPLITS = 16
SPLIT_TYPE = 'line'
REGEX_FILTER = '*fastq*'
LINES_DIVISOR = 4

def split_folds(fold):
  for i in range(SPLITS):
    util.ensure_dir_exists(fold + 'split' + str(i))
  return

def split_by_lines(inp_dir):
  # Splits a folder into groups by lines within each file
  # Used for scripts that operate line-by-line on all files

  for nm in _config.SPLITS:
    for fn in os.listdir(inp_dir + nm + '/'):
      if fnmatch.fnmatch(fn, REGEX_FILTER):
        nl = util.line_count(inp_dir + nm + '/' + fn)
        jump = nl / SPLITS
        jump = (jump / LINES_DIVISOR) * LINES_DIVISOR

        for i in range(SPLITS):
          if i < SPLITS - 1:
            arg = str(jump *  i + 1) + ',' + \
            str(jump * (i + 1)) + 'p;' + \
            str(jump * (i + 1) + 1) + 'q'
          else:
            arg = str(jump *  i + 1) + ',' + \
            str(nl) + 'p'
            
          # sed grabs a range of lines in a file
          subprocess.call('sed -n \'' + arg + '\' ' + inp_dir + nm + '/' + fn + ' > ' + inp_dir + 'split' + str(i) + '/' + nm + '/' + fn,
            shell = True)
  return


def split_by_files(inp_dir):
  # Splits a folder into groups of files
  # Used for scripts that process each file independently
  pass


def start_threads(inp_dir, out_dir, script):
  # Main meat of this 
  inps = [inp_dir + 'split' + str(s) + '/' for s in range(SPLITS)]
  outs = [out_dir + 'split' + str(s) + '/' for s in range(SPLITS)]
  logs = [out_dir + 'split' + str(s) + '/stdout.out' for s in range(SPLITS)]
  errs = [out_dir + 'split' + str(s) + '/stderr.out' for s in range(SPLITS)]

  outfs, errfs = [open(s, 'w') for s in logs], [open(s, 'w') for s in errs]

  ps = []
  for i in range(len(inps)):
    p = subprocess.Popen('python ' + 
      ' '.join([script + '.py', inps[i], outs[i]]),
      stdout = outfs[i], 
      stderr = errfs[i],
      shell = True)
    ps.append(p)
    print '\tstarted split', i, '\tpid:', p.pid

  exit_codes = [p.wait() for p in ps]
  print exit_codes
  print 'All processes done'

  return

def combine_outputs(out_dir):
  # Concatenates all split outputs together 
  # into the main output directory 
  out_splits = [out_dir + 'split' + str(s) + '/' for s in range(SPLITS)]

  fns = set()
  for s in out_splits:
    for fn in os.listdir(s):
      if fnmatch.fnmatch(s + fn, REGEX_FILTER):
        fns.add(fn)

  for fn in fns:
    util.exists_empty_fn(out_dir + fn)
    print '\tCombining', fn, '...'
    locs = [s + fn for s in out_splits]
    subprocess.call('cat ' + ' '.join(locs) + ' > ' + out_dir + fn, 
      shell = True)

  return


def main(inp_dir, out_dir, script):
  print '\tParallelizing', script, 'with', SPLITS, 'splits'

  # Make split folders
  split_folds(inp_dir)
  split_folds(out_dir)

  skip = True
  if not skip:
    # Split the input
    print '\tSplitting by', SPLIT_TYPE
    if SPLIT_TYPE == 'line':
      split_by_lines(inp_dir)
      pass
    elif SPLIT_TYPE == 'file':
      split_by_files(inp_dir)
    else:
      print 'ERROR: Invalid split type'
      return

  # Parallelize on split input
  start_threads(inp_dir, out_dir, script)

  # Combine split output into overall results
  combine_outputs(out_dir)

  print 'Done'
  return out_dir


if __name__ == '__main__':
  if len(sys.argv) != 4:
    mod = importlib.import_module(P_SCRIPT)
    inp_dir = getattr(mod, 'DEFAULT_INP_DIR')
    out_dir = _config.OUT_PLACE + getattr(mod, 'NAME') + '/'
    script = P_SCRIPT
  else:
    inp_dir = sys.argv[1]
    out_dir = sys.argv[2]
    script = sys.argv[3]

  main(inp_dir, out_dir, script)
