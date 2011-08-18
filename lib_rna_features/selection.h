/**
 *  \file selection.h
 *  \brief Provides functions to compute statistic measures of structural influences on RNA gene selection.
 */

#ifndef SELECTION_H
#define SELECTION_H

#include "lrfutils.h"

/**
 *  \brief Classify base pairs into class 1, 2 and 3
 * 
 *  Sufficient space must be allocated for 'classes' before calling this function.
 * 
 *  \param numBps	A number of base pairs in a given consensus structure
 *  \param bps		An array of base pairs of consensus structure
 *	\param len		The length of a sequence
 *	\param classes	An array of integers representing a class of each of the bases.
 */
extern void classify_base_pairs(int numBps, bp_info* bps, int len, int* classes);

/**
 *	\brief Calculates the average substitution rate for bases in each classes.
 *
 *	Sufficient space must be allocated for 'avg1', 'avg2' and 'avg3' before calling
 *	this function.
 *
 *	\param numSeq	A number of sequences in an alignment
 *	\param len		The length of an alignment
 *	\param seqs		The sequences in an alignment
 *	\param classes	Classes of bases of consensus structure
 *	\param avg1		The average substitution rate of class 1 bases
 *	\param avg2		The average substitution rate of class 2 bases
 *	\param avg3		The average substitution rate of class 3 bases
 */
extern void get_avg_sub_rate_for_classes(int numSeq, int len, char** seqs, int *classes,
										double *avg1, double *avg2, double *avg3);
										
#endif