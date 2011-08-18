/*
 *  Unit Test for RNAFeatures
 */

#include <string.h>
#include <stdlib.h>
#include "CUnit/Basic.h"

#include "cotransfold.h"
#include "selection.h"
#include "boltzmann.h"

#define EPSILON 0.000001

char *seq, *nopair_seq, *onlytrans_seq, *onlycis_seq, *noadjpair_seq;
char *structure, *nopair_structure, *onlytrans_structure, *onlycis_structure, *noadjpair_structure;
bp_info *bps, *nopair_bps, *onlytrans_bps, *onlycis_bps, *noadjpair_bps; 
int *numBps, *nopair_numBps, *onlytrans_numBps, *onlycis_numBps, *noadjpair_numBps;
int *classes, *nopair_classes;

/* The suite initialization function.
 * Opens the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */
int init_suite1(void)
{
	seq = (char*)malloc(sizeof(char) * 74);
	structure = (char*)malloc(sizeof(char) * 74);
	bps = (bp_info*)malloc(sizeof(bp_info) * 24);
	numBps = (int*)malloc(sizeof(int));
	classes = (int*)malloc(sizeof(int) * 73);
	strcpy(seq,       "GCCUUGUUGGCGCAAUCGGUAGCGCGUAUGACUCUUAAUCAUAAGGUUAGGGGUUCGAGCCCCCUACAGGGCU");
	strcpy(structure, "(((((((..((((........))))...((((.((((....))))))))(((((....)))))..))))))).");
	if (seq == NULL || structure == NULL || bps == NULL || numBps == NULL)
		return 1;
		
	nopair_seq = (char*)malloc(sizeof(char) * 74);
	nopair_structure = (char*)malloc(sizeof(char) * 74);
	nopair_bps = (bp_info*)malloc(sizeof(bp_info) * 1);
	nopair_numBps = (int*)malloc(sizeof(int));
	nopair_classes = (int*)malloc(sizeof(int) * 73);
	strcpy(nopair_seq,       "GCCUUGUUGGCGCAAUCGGUAGCGCGUAUGACUCUUAAUCAUAAGGUUAGGGGUUCGAGCCCCCUACAGGGCU");
	strcpy(nopair_structure, ".........................................................................");
	if (nopair_seq == NULL || nopair_structure == NULL || nopair_bps == NULL || nopair_numBps == NULL)
		return 1;
		
	onlytrans_seq = (char*)malloc(sizeof(char) * 21);
	onlytrans_structure = (char*)malloc(sizeof(char) * 21);
	onlytrans_bps = (bp_info*)malloc(sizeof(bp_info) * 9);
	onlytrans_numBps = (int*)malloc(sizeof(int));
	strcpy(onlytrans_seq,       "GGGGGGGGGUUUUUUUUUUU");
	strcpy(onlytrans_structure, "(((((((((.))))))))).");
	if (onlytrans_seq == NULL || onlytrans_structure == NULL || onlytrans_bps == NULL || onlytrans_numBps == NULL)
		return 1;
		
	onlycis_seq = (char*)malloc(sizeof(char) * 47);
	onlycis_structure = (char*)malloc(sizeof(char) * 47);
	onlycis_bps = (bp_info*)malloc(sizeof(bp_info) * 6);
	onlycis_numBps = (int*)malloc(sizeof(int));
	strcpy(onlycis_seq,       "UUUUUUUUUAAAAAAAAACCGGUGUUAACCCCCCCCCGGGGGGGGG");
	strcpy(onlycis_structure, ".................((()))((())).................");
	if (onlycis_seq == NULL || onlycis_structure == NULL || onlycis_bps == NULL || onlycis_numBps == NULL)
		return 1;
		
	noadjpair_seq = (char*)malloc(sizeof(char) * 47);
	noadjpair_structure = (char*)malloc(sizeof(char) * 47);
	noadjpair_bps = (bp_info*)malloc(sizeof(bp_info) * 6);
	noadjpair_numBps = (int*)malloc(sizeof(int));
	strcpy(noadjpair_seq,       "UUUUUUUUUAAAAAAAAAUGCCCCCCCCCGGGGGGGGG");
	strcpy(noadjpair_structure, ".................()().................");
	if (noadjpair_seq == NULL || noadjpair_structure == NULL || noadjpair_bps == NULL || noadjpair_numBps == NULL)
		return 1;
		
	return 0;
}

/* The suite cleanup function.
 * Closes the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */
int clean_suite1(void)
{
	free(seq); free(structure); free(bps); free(numBps);
	free(nopair_seq); free(nopair_structure); free(nopair_bps); free(nopair_numBps);
	free(onlytrans_seq); free(onlytrans_structure); free(onlytrans_bps); free(onlytrans_numBps);
	free(onlycis_seq); free(onlycis_structure); free(onlycis_bps); free(onlycis_numBps);
	free(noadjpair_seq); free(noadjpair_structure); free(noadjpair_bps); free(noadjpair_numBps);
	return 0;
}

/* 
 * Unit test of find_base_pairs().
 */
void testFIND_BASE_PAIRS(void)
{
	int i;
	int exp_bps_i[24] = {0, 1, 2, 3, 4, 5, 6, 9, 10, 11,
						12, 28, 29, 30, 31, 33, 34, 35, 36, 49,
						50, 51, 52, 53};
	int exp_bps_j[24] = {71, 70, 69, 68, 67, 66, 65, 24, 23, 22,
						21, 48, 47, 46, 45, 44, 43, 42, 41, 62,
						61, 60, 59, 58};
	find_base_pairs(structure, numBps, bps);
	CU_ASSERT(24 == *numBps);
	for (i = 0; i < 24; i++) {
		CU_ASSERT(exp_bps_i[i] == bps[i].bp.i);
		CU_ASSERT(exp_bps_j[i] == bps[i].bp.j);
	}
	
	find_base_pairs(nopair_structure, nopair_numBps, nopair_bps);
	CU_ASSERT(0 == *nopair_numBps);
	
	find_base_pairs(onlytrans_structure, onlytrans_numBps, onlytrans_bps);
	CU_ASSERT(9 == *onlytrans_numBps);
	
	find_base_pairs(onlycis_structure, onlycis_numBps, onlycis_bps);
	CU_ASSERT(6 == *onlycis_numBps);
	
	find_base_pairs(noadjpair_structure, noadjpair_numBps, noadjpair_bps);
	CU_ASSERT(2 == *noadjpair_numBps);
}

/* 
 * Unit test of cotransfold().
 */
void testCOTRANSFOLD(void)
{
	double *cis = (double*)malloc(sizeof(double));
	double *trans = (double*)malloc(sizeof(double));

   	cotransfold(seq, *numBps, bps, cis, trans, 1);
   	CU_ASSERT(fabs(*cis - 0.068216) < EPSILON);
	CU_ASSERT(fabs(*trans - 0.325186) < EPSILON);
	
	cotransfold(seq, *numBps, bps, cis, trans, 0);
   	CU_ASSERT(fabs(*cis - 0.955024) < EPSILON);
	CU_ASSERT(fabs(*trans - 4.552597) < EPSILON);
	
	cotransfold(nopair_seq, *nopair_numBps, nopair_bps, cis, trans, 1);
   	CU_ASSERT(fabs(*cis - 0.0) < EPSILON);
	CU_ASSERT(fabs(*trans - 0.0) < EPSILON);
	
	cotransfold(onlytrans_seq, *onlytrans_numBps, onlytrans_bps, cis, trans, 1);
   	CU_ASSERT(fabs(*cis - 0.0) < EPSILON);
	CU_ASSERT(fabs(*trans - 0.434294) < EPSILON);
	
	cotransfold(onlycis_seq, *onlycis_numBps, onlycis_bps, cis, trans, 1);
   	CU_ASSERT(fabs(*cis - 0.0) < EPSILON);
	CU_ASSERT(fabs(*trans - 0.0) < EPSILON);
	
	cotransfold(noadjpair_seq, *noadjpair_numBps, noadjpair_bps, cis, trans, 1);
   	CU_ASSERT(fabs(*cis - 0.0) < EPSILON);
	CU_ASSERT(fabs(*trans - 0.0) < EPSILON);
	
	free(cis); free(trans);
}

/*
 * Unit test of classify_base_pairs().
 */
void testCLASSIFY_BASE_PAIRS(void) {
	int i;
	int exp_classes[73] = {1, 2, 3, 3, 3, 2, 1, 0, 0, 1,
						   2, 2, 1, 0, 0, 0, 0, 0, 0, 0,
						   0, 1, 2, 2, 1, 0, 0, 0, 1, 2,
						   -1, -1, 0, -1, -1, 2, 1, 0, 0, 0,
						   0, 1, 2, -1, -1, -1, -1, 2, 1, 1,
						   2, 3, 2, 1, 0, 0, 0, 0, 1, 2,
						   3, 2, 1, 0, 0, 1, 2, 3, 3, 3,
						   2, 1, 0};
	classify_base_pairs(*numBps, bps, 73, classes);
	for (i = 0; i < 73; i++)
		CU_ASSERT(exp_classes[i] == classes[i]);
	classify_base_pairs(*nopair_numBps, nopair_bps, 73, nopair_classes);
	for (i = 0; i < 73; i++)
		CU_ASSERT(0 == nopair_classes[i]);
}

/*
 * Unit test of boltzmann().
 */
void testBOLTZMANN(void) {
	/**
	 * To be implemented
	 */
}

/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
int main()
{
   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite("Suite_1", init_suite1, clean_suite1);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   /* NOTE - ORDER IS IMPORTANT - MUST TEST fread() AFTER fprintf() */
   if ((NULL == CU_add_test(pSuite, "test of find_base_pairs()", testFIND_BASE_PAIRS)) ||
       (NULL == CU_add_test(pSuite, "test of cotransfold()", testCOTRANSFOLD)) ||
	   (NULL == CU_add_test(pSuite, "test of classify_base_pairs()", testCLASSIFY_BASE_PAIRS)) ||
	   (NULL == CU_add_test(pSuite, "test of boltzmann()", testBOLTZMANN)))
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}