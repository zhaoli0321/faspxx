/*! \file  BlaSparseCSR.c
 *
 *  \brief Sparse matrix operations for dCSRmat matrices
 *
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
 * \fn dCSRmat faspxx_dcsr_create (const INT m, const INT n, const INT nnz)
 *
 * \brief Create CSR sparse matrix data memory space
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   the new dCSRmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
dCSRmat faspxx_dcsr_create(const INT m, const INT n, const INT nnz)
{
    dCSRmat A;

    if (m > 0) {
        A.IA = (INT*)faspxx_mem_calloc(m + 1, sizeof(INT));
    } else {
        A.IA = NULL;
    }

    if (n > 0) {
        A.JA = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
    } else {
        A.JA = NULL;
    }

    if (nnz > 0) {
        A.val = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));
    } else {
        A.val = NULL;
    }

    A.row = m;
    A.col = n;
    A.nnz = nnz;

    return A;
}

/**
 * \fn void faspxx_dcsr_alloc (const INT m, const INT n, const INT nnz, dCSRmat *A)
 *
 * \brief Allocate CSR sparse matrix memory space
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCSRmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_alloc(const INT m, const INT n, const INT nnz, dCSRmat* A)
{
    if (m <= 0 || n <= 0) {
        printf("### ERROR: Matrix dim %d, %d must be positive! [%s]\n", m, n,
               __FUNCTION__);
        return;
    }

    if (m > 0) {
        A->IA = (INT*)faspxx_mem_calloc(m + 1, sizeof(INT));
    } else {
        A->IA = NULL;
    }

    if (nnz > 0) {
        A->JA  = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
        A->val = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));
    } else {
        A->JA  = NULL;
        A->val = NULL;
    }

    A->row = m;
    A->col = n;
    A->nnz = nnz;

    return;
}

/**
 * \fn void faspxx_dcsr_free (dCSRmat *A)
 *
 * \brief Free CSR sparse matrix data memory space
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_free(dCSRmat* A)
{
    if (A == NULL) return;

    faspxx_mem_free(A->IA);
    A->IA = NULL;
    faspxx_mem_free(A->JA);
    A->JA = NULL;
    faspxx_mem_free(A->val);
    A->val = NULL;
    A->col = 0;
    A->row = 0;
    A->nnz = 0;
    A      = NULL;
}

/**
 * \fn void faspxx_dcsr_sort (dCSRmat *A)
 *
 * \brief Sort each row of A in ascending order w.r.t. column indices.
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_sort(dCSRmat* A)
{
    const INT n = A->col;
    INT       i, j, start, row_length;

    // temp memory for sorting rows of A
    INT *index, *ja;
    DBL* a;

    index = (INT*)faspxx_mem_calloc(n, sizeof(INT));
    ja    = (INT*)faspxx_mem_calloc(n, sizeof(INT));
    a     = (DBL*)faspxx_mem_calloc(n, sizeof(DBL));

    for (i = 0; i < n; ++i) {
        start      = A->IA[i];
        row_length = A->IA[i + 1] - start;

        for (j = 0; j < row_length; ++j) index[j] = j;

        faspxx_aux_iQuickSortIndex(&(A->JA[start]), 0, row_length - 1, index);

        for (j = 0; j < row_length; ++j) {
            ja[j] = A->JA[start + index[j]];
            a[j]  = A->val[start + index[j]];
        }

        for (j = 0; j < row_length; ++j) {
            A->JA[start + j]  = ja[j];
            A->val[start + j] = a[j];
        }
    }

    // clean up memory
    faspxx_mem_free(index);
    index = NULL;
    faspxx_mem_free(ja);
    ja = NULL;
    faspxx_mem_free(a);
    a = NULL;
}

/**
 * \fn void faspxx_dcsr_getdiag (INT n, const dCSRmat *A, dvector *diag)
 *
 * \brief Get first n diagonal entries of a CSR matrix A
 *
 * \param n     Number of diagonal entries to get (if n=0, then get all diagonal
 * entries) \param A     Pointer to dCSRmat CSR matrix \param diag  Pointer to the
 * diagonal as a dvector
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_getdiag(INT n, const dCSRmat* A, dvector* diag)
{
    INT i, k, j, ibegin, iend;

    SHORT nthreads = 1, use_openmp = FALSE;

    if (n == 0 || n > A->row || n > A->col) n = MIN(A->row, A->col);

#ifdef _OPENMP
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    faspxx_dvec_alloc(n, diag);

    if (use_openmp) {
        INT mybegin, myend, myid;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, k, j)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                ibegin = A->IA[i];
                iend   = A->IA[i + 1];
                for (k = ibegin; k < iend; ++k) {
                    j = A->JA[k];
                    if ((j - i) == 0) {
                        diag->val[i] = A->val[k];
                        break;
                    } // end if
                }     // end for k
            }         // end for i
        }
    } else {
        for (i = 0; i < n; ++i) {
            ibegin = A->IA[i];
            iend   = A->IA[i + 1];
            for (k = ibegin; k < iend; ++k) {
                j = A->JA[k];
                if ((j - i) == 0) {
                    diag->val[i] = A->val[k];
                    break;
                } // end if
            }     // end for k
        }         // end for i
    }
}

/*!
 * \fn void faspxx_dcsr_diagpref (dCSRmat *A)
 *
 * \brief Re-order the column and data arrays of a CSR matrix,
 *        so that the first entry in each row is the diagonal
 *
 * \param A   Pointer to the matrix to be re-ordered
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note Reordering is done in place.
 *
 * Modified by Chensong Zhang on Dec/21/2012
 */
void faspxx_dcsr_diagpref(dCSRmat* A)
{
    const INT num_rowsA = A->row;
    DBL*      A_data    = A->val;
    INT*      A_i       = A->IA;
    INT*      A_j       = A->JA;

    // Local variable
    INT i, j;
    INT tempi, row_size;
    DBL tempd;

#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend, ibegin, iend;
    INT nthreads = faspxx_get_num_threads();
#endif

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

#ifdef _OPENMP
    if (num_rowsA > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, i, j, ibegin, iend, tempi, tempd, mybegin, myend)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, num_rowsA, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                ibegin = A_i[i];
                iend   = A_i[i + 1];
                // check whether the first entry is already diagonal
                if (A_j[ibegin] != i) {
                    for (j = ibegin + 1; j < iend; j++) {
                        if (A_j[j] == i) {
#if DEBUG_MODE > 2
                            printf("### DEBUG: Switch entry_%d with entry_0\n", j);
#endif
                            tempi       = A_j[ibegin];
                            A_j[ibegin] = A_j[j];
                            A_j[j]      = tempi;

                            tempd          = A_data[ibegin];
                            A_data[ibegin] = A_data[j];
                            A_data[j]      = tempd;
                            break;
                        }
                    }
                    if (j == iend) {
                        printf("### ERROR: Diagonal entry %d is zero!\n", i);
                        faspxx_chkerr(ERROR_MISC, __FUNCTION__);
                    }
                }
            }
        }
    } else {
#endif
        for (i = 0; i < num_rowsA; i++) {
            row_size = A_i[i + 1] - A_i[i];
            // check whether the first entry is already diagonal
            if (A_j[0] != i) {
                for (j = 1; j < row_size; j++) {
                    if (A_j[j] == i) {
#if DEBUG_MODE > 2
                        printf("### DEBUG: Switch entry_%d with entry_0\n", j);
#endif
                        tempi  = A_j[0];
                        A_j[0] = A_j[j];
                        A_j[j] = tempi;

                        tempd     = A_data[0];
                        A_data[0] = A_data[j];
                        A_data[j] = tempd;

                        break;
                    }
                }
                if (j == row_size) {
                    printf("### ERROR: Diagonal entry %d is zero!\n", i);
                    faspxx_chkerr(ERROR_MISC, __FUNCTION__);
                }
            }
            A_j += row_size;
            A_data += row_size;
        }
#ifdef _OPENMP
    }
#endif

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif
}

/**
 * \fn void faspxx_dcsr_cp (const dCSRmat *A, dCSRmat *B)
 *
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param B   Pointer to the dCSRmat matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_cp(const dCSRmat* A, dCSRmat* B)
{
    B->row = A->row;
    B->col = A->col;
    B->nnz = A->nnz;

    faspxx_iarray_cp(A->row + 1, A->IA, B->IA);
    faspxx_iarray_cp(A->nnz, A->JA, B->JA);
    faspxx_darray_cp(A->nnz, A->val, B->val);
}

/**
 * \fn void faspxx_dcsr_trans (const dCSRmat *A, dCSRmat *AT)
 *
 * \brief Find transpose of dCSRmat matrix A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param AT  Pointer to the transpose of dCSRmat matrix A (output)
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
INT faspxx_dcsr_trans(const dCSRmat* A, dCSRmat* AT)
{
    const INT n = A->row, m = A->col, nnz = A->nnz;

    // Local variables
    INT i, j, k, p;

    AT->row = m;
    AT->col = n;
    AT->nnz = nnz;

    AT->IA = (INT*)faspxx_mem_calloc(m + 1, sizeof(INT));

    AT->JA = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));

    if (A->val) {
        AT->val = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));

    } else {
        AT->val = NULL;
    }

    // first pass: find the Number of nonzeros in the first m-1 columns of A
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1

    // faspxx_iarray_set(m+1, AT->IA, 0);
    memset(AT->IA, 0, sizeof(INT) * (m + 1));

    for (j = 0; j < nnz; ++j) {
        i = A->JA[j]; // column Number of A = row Number of A'
        if (i < m - 1) AT->IA[i + 2]++;
    }

    for (i = 2; i <= m; ++i) AT->IA[i] += AT->IA[i - 1];

    // second pass: form A'
    if (A->val) {
        for (i = 0; i < n; ++i) {
            INT ibegin = A->IA[i], iend = A->IA[i + 1];
            for (p = ibegin; p < iend; p++) {
                j          = A->JA[p] + 1;
                k          = AT->IA[j];
                AT->JA[k]  = i;
                AT->val[k] = A->val[p];
                AT->IA[j]  = k + 1;
            } // end for p
        }     // end for i
    } else {
        for (i = 0; i < n; ++i) {
            INT ibegin = A->IA[i], iend1 = A->IA[i + 1];
            for (p = ibegin; p < iend1; p++) {
                j         = A->JA[p] + 1;
                k         = AT->IA[j];
                AT->JA[k] = i;
                AT->IA[j] = k + 1;
            } // end for p
        }     // end of i
    }         // end if

    return SUCCESS;
}

/**
 * \fn void faspxx_dcsr_compress (const dCSRmat *A, dCSRmat *B, const DBL dtol)
 *
 * \brief Compress a CSR matrix A and store in CSR matrix B by
 *        dropping small entries abs(aij)<=dtol
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param B     Pointer to dCSRmat CSR matrix
 * \param dtol  Drop tolerance
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_compress(const dCSRmat* A, dCSRmat* B, const DBL dtol)
{
    INT i, j, k;
    INT ibegin, iend1;

    SHORT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP
    if (B->nnz > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    INT* index = (INT*)faspxx_mem_calloc(A->nnz, sizeof(INT));

    B->row = A->row;
    B->col = A->col;

    B->IA = (INT*)faspxx_mem_calloc(A->row + 1, sizeof(INT));

    B->IA[0] = A->IA[0];

    // first pass: determine the size of B
    k = 0;
    for (i = 0; i < A->row; ++i) {
        ibegin = A->IA[i];
        iend1  = A->IA[i + 1];
        for (j = ibegin; j < iend1; ++j)
            if (ABS(A->val[j]) > dtol) {
                index[k] = j;
                ++k;
            } /* end of j */
        B->IA[i + 1] = k;
    } /* end of i */
    B->nnz = k;
    B->JA  = (INT*)faspxx_mem_calloc(B->nnz, sizeof(INT));
    B->val = (DBL*)faspxx_mem_calloc(B->nnz, sizeof(DBL));

    // second pass: generate the index and element to B
    if (use_openmp) {
        INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, B->nnz, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) {
                B->JA[i]  = A->JA[index[i]];
                B->val[i] = A->val[index[i]];
            }
        }
    } else {
        for (i = 0; i < B->nnz; ++i) {
            B->JA[i]  = A->JA[index[i]];
            B->val[i] = A->val[index[i]];
        }
    }

    faspxx_mem_free(index);
    index = NULL;
}

/**
 * \fn SHORT faspxx_dcsr_compress_inplace (dCSRmat *A, const DBL dtol)
 *
 * \brief Compress a CSR matrix A IN PLACE by
 *        dropping small entries abs(aij)<=dtol
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param dtol  Drop tolerance
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note This routine can be modified for filtering.
 */
SHORT faspxx_dcsr_compress_inplace(dCSRmat* A, const DBL dtol)
{
    const INT row = A->row;
    const INT nnz = A->nnz;

    INT   i, j, k;
    INT   ibegin, iend = A->IA[0];
    SHORT status = SUCCESS;
    k            = 0;
    for (i = 0; i < row; ++i) {
        ibegin = iend;
        iend   = A->IA[i + 1];
        for (j = ibegin; j < iend; ++j)
            if (ABS(A->val[j]) > dtol || i == A->JA[j]) {
                A->JA[k]  = A->JA[j];
                A->val[k] = A->val[j];
                ++k;
            } /* end of j */
        A->IA[i + 1] = k;
    } /* end of i */

    if (k <= nnz) {
        A->nnz = k;
        A->JA  = (INT*)faspxx_mem_realloc(A->JA, k * sizeof(INT));
        A->val = (DBL*)faspxx_mem_realloc(A->val, k * sizeof(DBL));
    } else {
        printf("### WARNING: Size of compressed matrix is bigger than original!\n");
        status = ERROR_UNKNOWN;
    }

    return (status);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
