/*! \file  BlaFormat.c
 *
 *  \brief Subroutines for matrix format conversion
 *
 *  \note  This file contains Level-1 (Bla) functions. It requires:
 *         AuxArray.c, AuxMemory.c, AuxThreads.c, and BlaSparseCSR.c,
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include "faspxx_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT faspxx_format_dcoo_dcsr (const dCOOmat *A, dCSRmat *B)
 *
 * \brief Transform a DBL matrix from its IJ format to its CSR format.
 *
 * \param A   Pointer to dCOOmat matrix
 * \param B   Pointer to dCSRmat matrix
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
SHORT faspxx_format_dcoo_dcsr(const dCOOmat* A, dCSRmat* B)
{
    const INT m = A->row, n = A->col, nnz = A->nnz;
    INT       iind, jind, i;

    faspxx_dcsr_alloc(m, n, nnz, B);
    INT* ia = B->IA;

    INT* ind = (INT*)faspxx_mem_calloc(m + 1, sizeof(INT));
    memset(ind, 0, sizeof(INT) * (m + 1));             // initialize ind
    for (i = 0; i < nnz; ++i) ind[A->rowind[i] + 1]++; // count nnz in each row

    ia[0] = 0; // first index starting from zero
    for (i = 1; i <= m; ++i) {
        ia[i]  = ia[i - 1] + ind[i]; // set row_idx
        ind[i] = ia[i];
    }

    // loop over nnz and set col_idx and val
    for (i = 0; i < nnz; ++i) {
        iind         = A->rowind[i];
        jind         = ind[iind];
        B->JA[jind]  = A->colind[i];
        B->val[jind] = A->val[i];
        ind[iind]    = ++jind;
    }

    faspxx_mem_free(ind);
    ind = NULL;

    return SUCCESS;
}

/**
 * \fn SHORT faspxx_format_dcsr_dcoo (const dCSRmat *A, dCOOmat *B)
 *
 * \brief Transform a DBL matrix from its CSR format to its IJ format.
 *
 * \param A   Pointer to dCSRmat matrix
 * \param B   Pointer to dCOOmat matrix
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
SHORT faspxx_format_dcsr_dcoo(const dCSRmat* A, dCOOmat* B)
{
    const INT m = A->row, nnz = A->nnz;
    INT       i, j;

    B->rowind = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
    B->colind = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
    B->val    = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));

#ifdef _OPENMP
#pragma omp parallel for if (m > OPENMP_HOLDS) private(i, j)
#endif
    for (i = 0; i < m; ++i) {
        for (j = A->IA[i]; j < A->IA[i + 1]; ++j) B->rowind[j] = i;
    }

    memcpy(B->colind, A->JA, nnz * sizeof(INT));
    memcpy(B->val, A->val, nnz * sizeof(DBL));

    return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
