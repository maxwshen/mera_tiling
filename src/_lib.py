# Stores project-specific functions

# Stores project-specific functions

import _config

def get_name(gene):
  # Standardize the capitalization of gene names
  if gene.lower() == 'mouseorf':
    name = 'mouseORF'
  else:
    name = gene.lower().capitalize()
  return name

def get_grna_fn(gene):
  # Finds the grna file given a gene
  return _config.DATA_DIR + 'expdesign2/' + \
    get_name(gene) + '_gRNA.txt'

def get_grna_locs_fn(gene):
  # Finds the gRNA location file given a gene
  return _config.DATA_DIR + 'expdesign2/' + \
    get_name(gene) + '.gRNA.filt_2bp.library.locs.txt'

def bg_header(exp, gene):
  # Creates a bedgraph file header
  if exp == 'bulk':
    ans = 'track type=bedGraph name="' + gene +  '_' +  exp + '" description="' + gene + '_' + exp + '" visibility=' + _config.BEDGRAPH_VISIBILITY[exp] + ' color=' + _config.BEDGRAPH_COLORS[exp]
  if exp in ['med', 'neg']:
    ans = 'track type=bedGraph name="' + gene +  '_gfp' +  exp + '" description="' + gene + '_gfp' + exp + '_ratio" visibility=' + _config.BEDGRAPH_VISIBILITY[exp] + ' color=' + _config.BEDGRAPH_COLORS[exp]
  return ans

def find_control(st):
  # Given a start location
  # returns which kind of control it is, if any
  for s in _config.POS_CONTROLS:
    if s in st:
      return 'pos'
  for s in _config.NEG_CONTROLS:
    if s in st:
      return 'neg'
  return 'none'