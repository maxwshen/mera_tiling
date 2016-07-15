import sys, os, datetime

import _config, _clean

from mylib import util

# Import all your steps
import a_demultiplex, a_parallelize, a2_condense, b_countgrna, c_bedgraph

##############################################################
##############################################################
util.shell_cp(_config.SRC_DIR + '_config.py', _config.RESULTS_PLACE)

start = datetime.datetime.now()
print start
##############################################################
##############################################################

runbool = False
runbool = True

a_res_dir = a_parallelize.main(_config.DATA_DIR + _config.DATA_FOLD, 
  _config.OUT_PLACE + 'a_demultiplex/', 'a_demultiplex', 
  run = runbool)


a2_res_dir = a2_condense.main(a_res_dir, 
  _config.OUT_PLACE + 'a2_condense/',
  run = runbool)


b_res_dir = b_countgrna.main(a2_res_dir, 
  _config.OUT_PLACE + 'b_countgrna/',
  run = runbool)


c_res_dir = c_bedgraph.main(b_res_dir, 
  _config.OUT_PLACE + 'c_bedgraph/',
  run = runbool)


##############################################################
##############################################################


print '\n\nDone.\n'
end = datetime.datetime.now()
print end, '\nTotal Time:', end - start

_clean.main()