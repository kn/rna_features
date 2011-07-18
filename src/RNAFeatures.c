/*********************************************************************                
 *                                                                   *
 *                          RNAFeatures.c                            *
 *                                                                   *
 *	                   Author: Katsuya Noguchi                       *
 *                                                                   *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cmdline.h"
#include "config.h"
#include "cotransfold.h"
#include "dnapars.h"
#include "fold.h"
#include "fold_vars.h"
#include "selection.h"
#include "boltzmann.h"

#define PRIVATE static
#define TRUE 1
#define FALSE 0

#define BUF_SIZE 1024

PRIVATE void formatSequence(char *seq);
PRIVATE void usage(void);
PRIVATE void help(void);
PRIVATE void version(void);

/********************************************************************
 *                                                                  *
 * main -- main program                                             *
 *                                                                  *
 ********************************************************************/

int main (int argc, char *argv[]) {
	int is_plain = TRUE, i = 0, num_char = 0, len = 0;
	double mfe;
	float mfe_prob;
	char *seq = NULL;
	char *structure = NULL;
	FILE *seqFile=stdin;
	struct gengetopt_args_info args;
	
	double *cis = (double*) allocate(sizeof(double), "cis");
	double *trans =(double*) allocate(sizeof(double), "trans");
	
	double *Q, *X, *Y;
	Q = (double*) allocate(sizeof(double), "Q");
	X = (double*) allocate(sizeof(double), "X");
	Y = (double*) allocate(sizeof(double), "Y");
	
	bp_info *bps;
	int *numBps = (int*)allocate(sizeof(int), "numBps");
	int* classes;

	if (cmdline_parser(argc, argv, &args) != 0) {
		usage();
		exit(EXIT_FAILURE);
	}
	if (args.help_given){
		help();
		exit(EXIT_SUCCESS);
	}
	if (args.version_given){
		version();
		exit(EXIT_SUCCESS);
	}
	if (args.inputs_num < 1){
		perror("This program requires following one arguments:\n-filename containing a sequence\ne.g. coTransFold tRNA_seq\n");
		exit(1);
	}

	if (args.weight_energy_given) {
		is_plain = FALSE;
	}
	
	// Global RNA package variables.
	  do_backtrack = 1; 
	  dangles = 2;

	// Open input files.
	seqFile = fopen(args.inputs[0], "r");
	if (seqFile == NULL){
		fprintf(stderr, "ERROR: Can't open input file %s\n", args.inputs[0]);
		exit(1);
	}

	// Read Sequence and Structure from input files.
	while (fgetc(seqFile) != EOF)
		num_char++;
	seq = (char*)allocate(sizeof(char) * (num_char + 1), "seq");
	
	rewind(seqFile);
	while (i < num_char)
		seq[i++] = fgetc(seqFile);
	seq[i] = '\0';

	// Close input files.
	fclose(seqFile);
	
	formatSequence(seq);
	
	len = strlen(seq);
	structure = (char*)allocate(sizeof(char) * len, "structure");
	classes = (int*)allocate(sizeof(int) * len, "classes");
	bps = (bp_info*)allocate(sizeof(bp_info) * len, "bps");
	mfe = fold(seq, structure);
	
	find_base_pairs(structure, numBps, bps);
	
	cotransfold(seq, *numBps, bps, cis, trans, is_plain);
	
	classify(*numBps, bps, len, classes);
	
	mfe_prob = boltzmann(seq, NULL, mfe, Q, X, Y);

	// Prints out result to console
	printf("Length: %d\nBase pairs and closing internal loop:\n", len);
	for (i = 0; i < *numBps; i++)
		printf("(%d, %d) : (%d, %d)\n", bps[i].bp.i, bps[i].bp.j, bps[i].loop.i, bps[i].loop.j);
	printf("Sequence:\n%s\nStructure:\n%s\nMFE: %f\n", seq, structure, mfe);
	printf("CIS: %f\nTRANS: %f\n", *cis, *trans);		
	printf("Classes (1 = class I, 2 = class II, 3 = class III, 0 = no class, -1 = discarded):\n");
	for (i = 0; i < len; i++)
		printf("%d ", classes[i]);
	printf("\n");
	
	free(cis);
	free(trans);
	free(bps);
	free(numBps);
	free(classes);	
	free(seq);
	free(structure);
	
	/* Test Parsimony */
	char** sequences = (char**)malloc(sizeof(char*) * 5);
	int lclasses[13] = {1, 1, -1, 2, 2, -1, 3, 3, -1, 1, 2, 3, 0};
	double *avg1 = (double*)malloc(sizeof(double));
	double *avg2 = (double*)malloc(sizeof(double));
	double *avg3 = (double*)malloc(sizeof(double));
	
	for (i = 0; i < 5; i++)
		sequences[i] = (char*)malloc(sizeof(char) * 13);
	strcpy(sequences[0], "AACGUGGCCAAAU");
	strcpy(sequences[1], "AAGGUCGCCAAAC");
	strcpy(sequences[2], "CAUUUCGUCACAA");
	strcpy(sequences[3], "GGUAUUUCGGCCU");
	strcpy(sequences[4], "GGGAUCUCGGCCC");
	
	get_avg_sub_rate_for_classes(5, 13, sequences, lclasses,
								avg1, avg2, avg3);
	printf("\nParsimony score:\nC1: %f\nC2: %f\nC3: %f\n", *avg1, *avg2, *avg3);
	
	printf("\nEb: %f\nVb: %f\nMFE probability: %f\n", *X / *Q, *Y / *Q - pow(*X / *Q, 2), mfe_prob);
}

PRIVATE void formatSequence(char *seq) {
	int i = 0;
	int j = 0;
	int len = strlen(seq);
	
	for (i = 0; i < len; i++) {
		seq[i] = toupper(seq[i]);
		if (seq[i] == 'T')
			seq[i] = 'U';
		if (seq[i] == 'A' || seq[i] == 'C' || seq[i] == 'G' || seq[i] == 'U')
			seq[j++] = seq[i];
	}
	seq[j] = '\0';
}

/********************************************************************
 *                                                                  *
 * usage, help, version - shows information and exits               *
 *                                                                  *
 ********************************************************************/

PRIVATE void usage(void) {
	help();
}

PRIVATE void help(void) {
	cmdline_parser_print_version();

	printf("\nUsage: %s [OPTIONS]... [FILES]\n\n", CMDLINE_PARSER_PACKAGE);
	printf("%s\n","  -h, --help              Print help and exit");
	printf("%s\n","  -V, --version           Print version and exit\n");
}

PRIVATE void version(void) {
	printf("RNAFeatures version " PACKAGE_VERSION "\n");
	exit(EXIT_SUCCESS);
}