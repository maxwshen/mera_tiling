# Config parameters
# imported by src/_config

FN = 'LIB025602_GEN00084827.fastq'
GENES = ['CDX2', 'HOXa1', 'mouseORF']

START = len('TGTGGAAAGGACGAAACACC')
BP_RELAX = 1

GROUPS = [['CDX2_bulk_gDNA', 'CDX2_Endoderm_R1_GFPpos', 'CDX2_Endoderm_R1_GFPneg', 'CDX2_Endoderm_R2_GFPpos', 'CDX2_Endoderm_R2_GFPneg'],
          ['HOXa1_bulk_gDNA', 'HOXa1_Endoderm_R1_GFPpos', 'HOXa1_Endoderm_R1_GFPneg', 'HOXa1_Endoderm_R2_GFPpos', 'HOXa1_Endoderm_R2_GFPneg', 'HOXa1_Undirected_R1_GFPpos', 'HOXa1_Undirected_R1_GFPneg', 'HOXa1_Undirected_R2_GFPpos', 'HOXa1_Undirected_R2_GFPneg'],
          ['CDX2_ORF_bulk_gDNA', 'CDX2_Endoderm_R1_ORF_GFPpos', 'CDX2_Endoderm_R1_ORF_GFPneg', 'CDX2_Endoderm_R2_ORF_GFPpos', 'CDX2_Endoderm_R2_ORF_GFPneg', 'CDX2_Undirected_R1_ORF_GFPpos', 'CDX2_Undirected_R1_ORF_GFPneg', 'CDX2_Undirected_R2_ORF_GFPpos', 'CDX2_Undirected_R2_ORF_GFPneg', 'HOXa1_ORF_bulk_gDNA', 'HOXa1_Endoderm_R1_ORF_GFPpos', 'HOXa1_Endoderm_R1_ORF_GFPneg', 'HOXa1_Undirected_R1_ORF_GFPpos', 'HOXa1_Undirected_R1_ORF_GFPneg', 'HOXa1_Undirected_R2_ORF_GFPpos', 'HOXa1_Undirected_R2_ORF_GFPneg']]

F_BEDGRAPH_GROUPS = [['CDX2_bulk_gDNA', 'CDX2_Endoderm_R1_GFPpos', 'CDX2_Endoderm_R1_GFPneg'], ['CDX2_bulk_gDNA', 'CDX2_Endoderm_R2_GFPpos', 'CDX2_Endoderm_R2_GFPneg'], ['HOXa1_bulk_gDNA', 'HOXa1_Endoderm_R1_GFPpos', 'HOXa1_Endoderm_R1_GFPneg'], ['HOXa1_bulk_gDNA', 'HOXa1_Endoderm_R2_GFPpos', 'HOXa1_Endoderm_R2_GFPneg'], ['HOXa1_bulk_gDNA', 'HOXa1_Undirected_R1_GFPpos', 'HOXa1_Undirected_R1_GFPneg'], ['HOXa1_bulk_gDNA', 'HOXa1_Undirected_R2_GFPpos', 'HOXa1_Undirected_R2_GFPneg']]

# Chose groups that have good representation of reads
'''
- currently: bedgraph for 
CDX2: Endoderm_R1_GFP pos/neg
CDX2: Endoderm_R2_GFP pos/neg
HOXa1: Endoderm_R1_GFP pos/neg
HOXa1: Endoderm_R2_GFP pos/neg
HOXa1: Undirected_R1_GFP pos/neg
HOXa1: Undirected_R2_GFP pos/neg
'''

g_morf_exps = ['CDX2_Endoderm_R1_ORF_GFP', 'CDX2_Endoderm_R2_ORF_GFP', 'CDX2_Undirected_R1_ORF_GFP', 'CDX2_Undirected_R2_ORF_GFP', 'HOXa1_Endoderm_R1_ORF_GFP', 'HOXa1_Undirected_R1_ORF_GFP', 'HOXa1_Undirected_R2_ORF_GFP']