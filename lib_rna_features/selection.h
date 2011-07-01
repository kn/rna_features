#include "lrfutils.h"

extern void classify(int numBps, bp_info* bps, int len, int* classes);
extern void separate_sequences_by_class(int numSeq, int len, int* classes, char** seqs,
										int *len1, char** seqC1, int *len2, char** seqC2,
										int *len3, char** seqC3);