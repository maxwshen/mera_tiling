import sys

PRJ_DIR = '/data/gl/g5/bh0085/mera_tiling/'  
SRC_DIR = PRJ_DIR + 'src/'

# toy = True
toy = False
if toy:
  PRJ_DIR += 'toy/'
#######################################################
# Note: Directories should end in / always
#######################################################
DATA_DIR = PRJ_DIR + 'data/'
OUT_PLACE = PRJ_DIR + 'out/'
RESULTS_PLACE = PRJ_DIR + 'results/'
#######################################################
#######################################################

CLEAN = False       # Values = 'ask', True, False

# which data are we using? import that data's parameters
# DATA_FOLD = '2016-06-30/'
# DATA_FOLD = '2016-06-30-rep1/'
# DATA_FOLD = '2016-06-30-rep1-repr/'
REPRODUCE = False
# DATA_FOLD = '2016-07-06/'
# DATA_FOLD = '2016-08-11/'
# DATA_FOLD = '2016-08-26/'
# DATA_FOLD = '2016-10-17/'
# DATA_FOLD = '2017-01-16/'
# DATA_FOLD = '2017-03-06/'
# DATA_FOLD = '2017-04-21/'
# DATA_FOLD = '2017-06-05/'
DATA_FOLD = '2017-10-13/'


sys.path.insert(0, DATA_DIR + DATA_FOLD)
import _dataconfig as d
print 'Using data folder:\n', DATA_DIR + DATA_FOLD
OUT_PLACE += DATA_FOLD
RESULTS_PLACE += DATA_FOLD

#######################################################
# Project-specific parameters
#######################################################


BARCODE_MM = 1
QUALITY_CUTOFF = 20

PSEUDOCOUNT = 1


# Read in barcodes. 1-to-1 mapping between 2 columns
BARCODES = []
SPLITS = []
with open(DATA_DIR + DATA_FOLD + 'barcodes.txt') as f:
  for i, line in enumerate(f):
    words = line.split()
    if len(words) != 2:
      print 'ERROR: Bad barcodes file'
    SPLITS.append(words[0])
    BARCODES.append(words[1].upper())
SPLITS_FROM_BARCODE, BARCODE_FROM_SPLITS = dict(), dict()
for i in range(len(BARCODES)):
  SPLITS_FROM_BARCODE[BARCODES[i]] = SPLITS[i]
  BARCODE_FROM_SPLITS[SPLITS[i]] = BARCODES[i]
BARCODES, SPLITS = sorted(BARCODES), sorted(SPLITS)

GRNA_COUNT_THRESHOLD = 5

CHROMS = {'myh9': 'chr17', 'msh2': 'chr17', 'brca2': 'chr5', 'hoxa1': 'chr6', 'mouseORF': 'chr5', 'sdhd': 'chr9', 'pou5f1': 'chr17', 'tdgf1': 'chr9', 'p53': 'chr11', 'sox2': 'chr3', 'cdx2': 'chr5', 'apobec3': 'chr15', 'apobec3b': 'chr22'}

BEDGRAPH_COLORS = {'bulk': '0,0,0', 'med': '0,255,0', 'neg': '255,0,0'}
BEDGRAPH_VISIBILITY = {'bulk': 'dense', 'med': 'full', 'neg': 'full'}

POS_CONTROLS = ['eGFP_control', 'eGFP_offtarget', 'control_', 'pos_control_', 'GFP_']
NEG_CONTROLS = ['Neg_control', 'Neg_Control', 'neg_control', 'Non-Targeting Control']