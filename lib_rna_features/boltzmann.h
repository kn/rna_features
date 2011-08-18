/**
 *  \file boltzmann.h
 * 
 *  \brief Calculates statistic measures of higher moments of Boltzmann distribution
 *	of single RNA sequences.
 */

#ifndef BOLTZMANN_H
#define BOLTZMANN_H

#include "data_structures.h"
#define FLT_OR_DBL double

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  a flag indicating that auxilary arrays are needed throughout the computations which are necessary for stochastic backtracking
 */
extern  int st_back;

/**
 *  \brief Computes the statistic measures of higher moments of Boltzmann distribution
 *	for a given RNA sequence.
 * 
 *  \param sequence		The RNA sequence to be computed
 *	\param mfe			The MFE of the RNA sequence
 *	\param Q
 *	\param X
 *	\param Y
 *
 *  \returns          The log probablity of a given MFE
 */
float   boltzmann(const char *sequence, double mfe, FLT_OR_DBL *Q, FLT_OR_DBL *X, FLT_OR_DBL *Y);

/**
 *  \brief Free arrays from pf_fold()
 */
void    free_pf_arrays(void);

#endif
