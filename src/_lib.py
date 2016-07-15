# Stores project-specific functions

# Stores project-specific functions

import _config

test = 1

def get_grna_fn(gene):
  return _config.DATA_DIR + 'expdesign2/' + \
    gene.capitalize() + '_gRNA.txt'

def get_grna_locs_fn(gene):
  return _config.DATA_DIR + 'expdesign2/' + \
    gene.capitalize() + '.gRNA.filt_2bp.library.locs.txt'

def bg_header(exp, gene):
  if exp == 'bulk':
    ans = 'track type=bedGraph name="' + gene +  '_' +  exp + '" description="' + gene + '_' + exp + '" visibility=' + _config.BEDGRAPH_VISIBILITY[exp] + ' color=' + _config.BEDGRAPH_COLORS[exp]
  if exp in ['med', 'neg']:
    ans = 'track type=bedGraph name="' + gene +  '_gfp' +  exp + '" description="' + gene + '_gfp' + exp + '_ratio" visibility=' + _config.BEDGRAPH_VISIBILITY[exp] + ' color=' + _config.BEDGRAPH_COLORS[exp]
  return ans