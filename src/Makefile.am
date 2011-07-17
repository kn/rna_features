bin_PROGRAMS = RNAFeatures

RNAFeatures_SOURCES = RNAFeatures.c cmdline.c cmdline.h

RNAFeatures_LDADD = ../lib_rna_features/libRNAFeatures.a ../lib_phylip/libPHYLIP.a -L../librna2 -lRNA

INCLUDES = -I../lib_rna_features -I../librna2 -I../lib_phylip