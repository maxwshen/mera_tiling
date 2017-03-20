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
  ok = False
  if 'bulk' in exp:
    ans = 'track type=bedGraph name="' + gene +  '_' +  exp + '" description="' + gene + '_' + exp + '" visibility=' + _config.BEDGRAPH_VISIBILITY['bulk'] + ' color=' + _config.BEDGRAPH_COLORS['bulk']
    ok = True
  if 'med' in exp:
    ans = 'track type=bedGraph name="' + gene +  '_gfp' +  exp + '" description="' + gene + '_gfp' + exp + '_ratio" visibility=' + _config.BEDGRAPH_VISIBILITY['med'] + ' color=' + _config.BEDGRAPH_COLORS['med']
    ok = True
  else:
    # if 'neg' in exp:
    ans = 'track type=bedGraph name="' + gene +  '_gfp' +  exp + '" description="' + gene + '_gfp' + exp + '_ratio" visibility=' + _config.BEDGRAPH_VISIBILITY['neg'] + ' color=' + _config.BEDGRAPH_COLORS['neg']
    ok = True
  if not ok:
    print 'Error: exp name does not contain bulk/med/neg'
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

def is_control(obj):
  if isinstance(obj, str):
    return not obj.isdigit()
  return False


def is_pos_control(obj):
  if not is_control(obj):
    return False
  else:
    if sum([1 for s in _config.POS_CONTROLS if s in obj]) > 0:
      return True
    else:
      return False

def is_neg_control(obj):
  if not is_control(obj):
    return False
  else:
    if sum([1 for s in _config.NEG_CONTROLS if s in obj]) > 0:
      return True
    else:
      return False