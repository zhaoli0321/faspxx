/*! \file  BlaSparseCOO.c
 *
 *  \brief Sparse matrix operations for dCOOmat matrices
 *
 *  \note  This file contains Level-1 (Bla) functions. It requires:
 *         AuxMemory.c and AuxThreads.c
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
 * \fn dCOOmat faspxx_dcoo_create (const INT m, const INT n, const INT nnz)
 *
 * \brief Create IJ sparse matrix data memory space
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   The new dCOOmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
dCOOmat faspxx_dcoo_create(const INT m, const INT n, const INT nnz)
{
    dCOOmat A;

    A.rowind = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
    A.colind = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
    A.val    = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));

    A.row = m;
    A.col = n;
    A.nnz = nnz;

    return A;
}

/**
 * \fn void faspxx_dcoo_alloc (const INT m, const INT n, const INT nnz, dCOOmat *A)
 *
 * \brief Allocate COO sparse matrix memory space
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCSRmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcoo_alloc(const INT m, const INT n, const INT nnz, dCOOmat* A)
{

    if (nnz > 0) {
        A->rowind = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
        A->colind = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
        A->val    = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));
    } else {
        A->rowind = NULL;
        A->colind = NULL;
        A->val    = NULL;
    }

    A->row = m;
    A->col = n;
    A->nnz = nnz;

    return;
}

/**
 * \fn void faspxx_dcoo_free (dCOOmat *A)
 *
 * \brief Free IJ sparse matrix data memory space
 *
 * \param A   Pointer to the dCOOmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcoo_free(dCOOmat* A)
{
    if (A == NULL) return;

    faspxx_mem_free(A->rowind);
    A->rowind = NULL;
    faspxx_mem_free(A->colind);
    A->colind = NULL;
    faspxx_mem_free(A->val);
    A->val = NULL;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
