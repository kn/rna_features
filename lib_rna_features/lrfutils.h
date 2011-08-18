/**
 *  \file lrfutils.h
 *  \brief Provides utility functions to help other programs in RNAFeatures.
 */

#ifndef LRFUTILS_H
#define LRFUTILS_H

typedef struct {
	int i;
	int j;
} tuple;

typedef struct {
	tuple bp;
	tuple loop;
} bp_info;

/**
 *  \brief Allocates memory of a given size. Outputs error message if it fails.
 * 
 *  \param size		A size of memory to be allocated
 *  \param name		A name of the pointer for which a memory is allocated
 */
extern void *allocate(unsigned int size, char* name);

/**
 *	\brief Determine base pairs and types of loops they close, given its parenthesis structure.
 *
 *	\param structure	A parenthesis structure of a sequence
 *	\param numBps		A number of base pairs in the structure
 */
extern void find_base_pairs(const char* structure, int* numBps, bp_info* bps);
#endif