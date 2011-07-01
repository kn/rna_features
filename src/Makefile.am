bin_PROGRAMS = RNAFeatures

RNAFeatures_SOURCES = RNAFeatures.c cmdline.c cmdline.h

RNAFeatures_LDADD = ../lib_rna_features/libRNAFeatures.a -L../librna -lRNA -lm ../lib_phylip/libPHYLIP.a

INCLUDES = -I../lib_rna_features -I../librna -I../lib_phylip