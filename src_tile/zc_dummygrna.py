import _config
import subprocess

dummygrna = 'GAGGCGTCTGGGTGGCTCTTG'

fold = '/cluster/mshen/prj/mera_tiling/out/2017-03-06/a_demultiplex'

for exp in _config.SPLITS:
  # print exp, subprocess.check_output('grep "' + dummygrna + '" ' + fold + '/' + exp + '/' + _config.d.FN + ' | wc -l', shell = True)
  pass


Gs = 'G' * 30
for exp in _config.SPLITS:
  print exp, subprocess.check_output('grep "' + Gs + '" ' + fold + '/' + exp + '/' + _config.d.FN + ' | wc -l', shell = True)