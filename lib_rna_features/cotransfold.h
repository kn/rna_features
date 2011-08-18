/**
 *  \file cotransfold.h
 *  \brief Provides a function to calculates Cis and Trans of a given RNA sequence and its structure.
 */

#ifndef COTRANSFOLD_H
#define COTRANSFOLD_H

#include "lrfutils.h"

/**
 *  \brief Calculates plain or energy weighted Cis and Trans.
 * 
 *  Sufficient space must be allocated for 'cis' and 'trans' before calling this function.
 * 
 *	\param seq		A sequence
 *  \param numBps	A number of base pairs in a given secondary structure
 *  \param bps		An array of base pairs of secondary structure
 *	\param cis		A poniter to which a Cis value will be stored
 *	\param trans	A poniter to which a Trans value will be stored
 *	\param is_plain	An integer indicating plain (1) or energy_weighted (0)
 */
extern void cotransfold(char *seq, int numBps, bp_info *bps, double *cis, double *trans, int is_plain);

#endif