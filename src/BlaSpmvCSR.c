/*! \file  BlaSpmvCSR.c
 *
 *  \brief Linear algebraic operations for dCSRmat matrices
 *
 *  \note  This file contains Level-2 (Bla) functions. It requires:
 *         AuxArray.c, AuxMemory.c, AuxThreads.c, BlaSparseCSR.c, BlaSparseUtil.c,
 *         and BlaArray.c
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include "faspxx_functs.h"

extern unsigned long total_alloc_mem;
extern unsigned long total_alloc_count;

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT faspxx_blas_dcsr_add (const dCSRmat *A, const DBL alpha,
 *                               const dCSRmat *B, const DBL beta, dCSRmat *C)
 *
 * \brief compute C = alpha*A + beta*B in CSR format
 *
 * \param A      Pointer to dCSRmat matrix
 * \param alpha  DBL factor alpha
 * \param B      Pointer to dCSRmat matrix
 * \param beta   DBL factor beta
 * \param C      Pointer to dCSRmat matrix
 *
 * \return       SUCCESS if succeed, ERROR if not
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
SHORT faspxx_blas_dcsr_add(const dCSRmat* A, const DBL alpha, const dCSRmat* B,
                           const DBL beta, dCSRmat* C)
{
    INT i, j, k, l;
    INT count = 0, added, countrow;

    SHORT status = SUCCESS, use_openmp = FALSE;

#ifdef _OPENMP
    INT mybegin, myend, myid, nthreads;
    if (A->nnz > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (A->row != B->row || A->col != B->col) {
        printf("### ERROR: Matrix sizes do not match!\n");
        status = ERROR_MAT_SIZE;
        goto FINISHED;
    }

    if (A == NULL && B == NULL) {
        C->row = 0;
        C->col = 0;
        C->nnz = 0;
        status = SUCCESS;
        goto FINISHED;
    }

    if (A->nnz == 0 && B->nnz == 0) {
        C->row = A->row;
        C->col = A->col;
        C->nnz = A->nnz;
        status = SUCCESS;
        goto FINISHED;
    }

    // empty matrix A
    if (A->nnz == 0 || A == NULL) {
        faspxx_dcsr_alloc(B->row, B->col, B->nnz, C);
        memcpy(C->IA, B->IA, (B->row + 1) * sizeof(INT));
        memcpy(C->JA, B->JA, (B->nnz) * sizeof(INT));

        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i)
            {
                myid = omp_get_thread_num();
                faspxx_get_start_end(myid, nthreads, A->nnz, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) C->val[i] = B->val[i] * beta;
            }
#endif
        } else {
            for (i = 0; i < A->nnz; ++i) C->val[i] = B->val[i] * beta;
        }

        status = SUCCESS;
        goto FINISHED;
    }

    // empty matrix B
    if (B->nnz == 0 || B == NULL) {
        faspxx_dcsr_alloc(A->row, A->col, A->nnz, C);
        memcpy(C->IA, A->IA, (A->row + 1) * sizeof(INT));
        memcpy(C->JA, A->JA, (A->nnz) * sizeof(INT));

        if (use_openmp) {
#ifdef _OPENMP
            INT mybegin, myend, myid;
#pragma omp parallel private(myid, mybegin, myend, i)
            {
                myid = omp_get_thread_num();
                faspxx_get_start_end(myid, nthreads, A->nnz, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) C->val[i] = A->val[i] * alpha;
            }
#endif
        } else {
            for (i = 0; i < A->nnz; ++i) C->val[i] = A->val[i] * alpha;
        }

        status = SUCCESS;
        goto FINISHED;
    }

    C->row = A->row;
    C->col = A->col;

    C->IA = (INT*)faspxx_mem_calloc(C->row + 1, sizeof(INT));

    // allocate work space for C->JA and C->val
    C->JA = (INT*)faspxx_mem_calloc(A->nnz + B->nnz, sizeof(INT));

    C->val = (DBL*)faspxx_mem_calloc(A->nnz + B->nnz, sizeof(DBL));

    // initial C->IA
    memset(C->IA, 0, sizeof(INT) * (C->row + 1));
    memset(C->JA, -1, sizeof(INT) * (A->nnz + B->nnz));

    for (i = 0; i < A->row; ++i) {
        countrow = 0;
        for (j = A->IA[i]; j < A->IA[i + 1]; ++j) {
            C->val[count] = alpha * A->val[j];
            C->JA[count]  = A->JA[j];
            C->IA[i + 1]++;
            count++;
            countrow++;
        } // end for js

        for (k = B->IA[i]; k < B->IA[i + 1]; ++k) {
            added = 0;

            for (l = C->IA[i]; l < C->IA[i] + countrow + 1; l++) {
                if (B->JA[k] == C->JA[l]) {
                    C->val[l] = C->val[l] + beta * B->val[k];
                    added     = 1;
                    break;
                }
            } // end for l

            if (added == 0) {
                C->val[count] = beta * B->val[k];
                C->JA[count]  = B->JA[k];
                C->IA[i + 1]++;
                count++;
            }

        } // end for k

        C->IA[i + 1] += C->IA[i];
    }

    C->nnz = count;
    C->JA  = (INT*)faspxx_mem_realloc(C->JA, (count) * sizeof(INT));
    C->val = (DBL*)faspxx_mem_realloc(C->val, (count) * sizeof(DBL));

FINISHED:
    return status;
}

/**
 * \fn void faspxx_blas_dcsr_axm (dCSRmat *A, const DBL alpha)
 *
 * \brief Multiply a sparse matrix A in CSR format by a scalar alpha.
 *
 * \param A      Pointer to dCSRmat matrix A
 * \param alpha  DBL factor alpha
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_axm(dCSRmat* A, const DBL alpha)
{
    const INT nnz = A->nnz;

    // A direct calculation can be written as:
    faspxx_blas_darray_ax(nnz, alpha, A->val);
}

/**
 * \fn void faspxx_blas_dcsr_mxv (const dCSRmat *A, const DBL *x, DBL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_mxv(const dCSRmat* A, const DBL* x, DBL* y)
{
    const INT  m  = A->row;
    const INT *ia = A->IA, *ja = A->JA;
    const DBL* aj = A->val;

    INT          i, k, begin_row, end_row, nnz_row;
    register DBL temp;

    SHORT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP
    if (m > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
        INT myid, mybegin, myend;

#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row,    \
                                     nnz_row, k)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                nnz_row   = end_row - begin_row;
                switch (nnz_row) {
                    case 3:
                        k = begin_row;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        break;
                    case 4:
                        k = begin_row;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        break;
                    case 5:
                        k = begin_row;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        break;
                    case 6:
                        k = begin_row;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        break;
                    case 7:
                        k = begin_row;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        k++;
                        temp += aj[k] * x[ja[k]];
                        break;
                    default:
                        for (k = begin_row; k < end_row; ++k) {
                            temp += aj[k] * x[ja[k]];
                        }
                        break;
                }
                y[i] = temp;
            }
        }
    }

    else {
        for (i = 0; i < m; ++i) {
            temp      = 0.0;
            begin_row = ia[i];
            end_row   = ia[i + 1];
            nnz_row   = end_row - begin_row;
            switch (nnz_row) {
                case 3:
                    k = begin_row;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    break;
                case 4:
                    k = begin_row;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    break;
                case 5:
                    k = begin_row;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    break;
                case 6:
                    k = begin_row;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    break;
                case 7:
                    k = begin_row;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    k++;
                    temp += aj[k] * x[ja[k]];
                    break;
                default:
                    for (k = begin_row; k < end_row; ++k) {
                        temp += aj[k] * x[ja[k]];
                    }
                    break;
            }
            y[i] = temp;
        }
    }
}

/**
 * \fn void faspxx_blas_dcsr_mxv_agg (const dCSRmat *A, const DBL *x, DBL *y)
 *
 * \brief Matrix-vector multiplication y = A*x (nonzeros of A = 1)
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_mxv_agg(const dCSRmat* A, const DBL* x, DBL* y)
{
    const INT    m  = A->row;
    const INT *  ia = A->IA, *ja = A->JA;
    INT          i, k, begin_row, end_row;
    register DBL temp;

#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend;
    INT nthreads = faspxx_get_num_threads();
#endif

#ifdef _OPENMP
    if (m > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, i, mybegin, myend, temp, begin_row, end_row, k)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                y[i] = temp;
            }
        }
    } else {
#endif
        for (i = 0; i < m; ++i) {
            temp      = 0.0;
            begin_row = ia[i];
            end_row   = ia[i + 1];
            for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
            y[i] = temp;
        }
#ifdef _OPENMP
    }
#endif
}

/**
 * \fn void faspxx_blas_dcsr_aAxpy (const DBL alpha, const dCSRmat *A,
 *                                const DBL *x, DBL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  DBL factor alpha
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_aAxpy(const DBL alpha, const dCSRmat* A, const DBL* x, DBL* y)
{
    const INT    m  = A->row;
    const INT *  ia = A->IA, *ja = A->JA;
    const DBL*   aj = A->val;
    INT          i, k, begin_row, end_row;
    register DBL temp;

    SHORT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP
    if (m > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (alpha == 1.0) {
        if (use_openmp) {
            INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k)
#endif
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    temp      = 0.0;
                    begin_row = ia[i];
                    end_row   = ia[i + 1];
                    for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
                    y[i] += temp;
                }
            }
        } else {
            for (i = 0; i < m; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
                y[i] += temp;
            }
        }
    }

    else if (alpha == -1.0) {
        if (use_openmp) {
            INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, temp, i, begin_row, end_row, k)
#endif
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    temp      = 0.0;
                    begin_row = ia[i];
                    end_row   = ia[i + 1];
                    for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
                    y[i] -= temp;
                }
            }
        } else {
            for (i = 0; i < m; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
                y[i] -= temp;
            }
        }
    }

    else {
        if (use_openmp) {
            INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k)
#endif
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    temp      = 0.0;
                    begin_row = ia[i];
                    end_row   = ia[i + 1];
                    for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
                    y[i] += temp * alpha;
                }
            }
        } else {
            for (i = 0; i < m; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
                y[i] += temp * alpha;
            }
        }
    }
}

/**
 * \fn void faspxx_blas_dcsr_aAxpy_agg (const DBL alpha, const dCSRmat *A,
 *                                    const DBL *x, DBL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y (nonzeros of A = 1)
 *
 * \param alpha  DBL factor alpha
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_aAxpy_agg(const DBL alpha, const dCSRmat* A, const DBL* x, DBL* y)
{
    const INT  m  = A->row;
    const INT *ia = A->IA, *ja = A->JA;

    INT          i, k, begin_row, end_row;
    register DBL temp;

    if (alpha == 1.0) {
#ifdef _OPENMP
        if (m > OPENMP_HOLDS) {
            INT myid, mybegin, myend;
            INT nthreads = faspxx_get_num_threads();
#pragma omp parallel for private(myid, i, mybegin, myend, begin_row, end_row, temp, k)
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    temp      = 0.0;
                    begin_row = ia[i];
                    end_row   = ia[i + 1];
                    for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                    y[i] += temp;
                }
            }
        } else {
#endif
            for (i = 0; i < m; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                y[i] += temp;
            }
#ifdef _OPENMP
        }
#endif
    } else if (alpha == -1.0) {
#ifdef _OPENMP
        if (m > OPENMP_HOLDS) {
            INT myid, mybegin, myend;
            INT nthreads = faspxx_get_num_threads();
#pragma omp parallel for private(myid, i, mybegin, myend, begin_row, end_row, temp, k)
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    temp      = 0.0;
                    begin_row = ia[i];
                    end_row   = ia[i + 1];
                    for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                    y[i] -= temp;
                }
            }
        } else {
#endif
            for (i = 0; i < m; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                y[i] -= temp;
            }
#ifdef _OPENMP
        }
#endif
    }

    else {
#ifdef _OPENMP
        if (m > OPENMP_HOLDS) {
            INT myid, mybegin, myend;
            INT nthreads = faspxx_get_num_threads();
#pragma omp parallel for private(myid, i, mybegin, myend, begin_row, end_row, temp, k)
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, m, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    temp      = 0.0;
                    begin_row = ia[i];
                    end_row   = ia[i + 1];
                    for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                    y[i] += temp * alpha;
                }
            }
        } else {
#endif
            for (i = 0; i < m; ++i) {
                temp      = 0.0;
                begin_row = ia[i];
                end_row   = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) temp += x[ja[k]];
                y[i] += temp * alpha;
            }
#ifdef _OPENMP
        }
#endif
    }
}

/**
 * \fn DBL faspxx_blas_dcsr_vmv (const dCSRmat *A, const DBL *x, const DBL *y)
 *
 * \brief vector-Matrix-vector multiplication alpha = y'*A*x
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
DBL faspxx_blas_dcsr_vmv(const dCSRmat* A, const DBL* x, const DBL* y)
{
    register DBL value = 0.0;
    const INT    m     = A->row;
    const INT *  ia = A->IA, *ja = A->JA;
    const DBL*   aj = A->val;
    INT          i, k, begin_row, end_row;
    register DBL temp;

    SHORT use_openmp = FALSE;

#ifdef _OPENMP
    if (m > OPENMP_HOLDS) {
        use_openmp = TRUE;
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : value) private(i, temp, begin_row, end_row, k)
#endif
        for (i = 0; i < m; ++i) {
            temp      = 0.0;
            begin_row = ia[i];
            end_row   = ia[i + 1];
            for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
            value += y[i] * temp;
        }
    } else {
        for (i = 0; i < m; ++i) {
            temp      = 0.0;
            begin_row = ia[i];
            end_row   = ia[i + 1];
            for (k = begin_row; k < end_row; ++k) temp += aj[k] * x[ja[k]];
            value += y[i] * temp;
        }
    }
    return value;
}

/**
 * \fn void faspxx_blas_dcsr_mxm (const dCSRmat *A, const dCSRmat *B, dCSRmat *C)
 *
 * \brief Sparse matrix multiplication C=A*B
 *
 * \param A   Pointer to the dCSRmat matrix A
 * \param B   Pointer to the dCSRmat matrix B
 * \param C   Pointer to dCSRmat matrix equal to A*B
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_mxm(const dCSRmat* A, const dCSRmat* B, dCSRmat* C)
{
    INT i, j, k, l, count;

    INT* JD = (INT*)faspxx_mem_calloc(B->col, sizeof(INT));

    C->row = A->row;
    C->col = B->col;
    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)faspxx_mem_calloc(C->row + 1, sizeof(INT));

    for (i = 0; i < B->col; ++i) JD[i] = -1;

    // step 1: Find first the structure IA of C
    for (i = 0; i < C->row; ++i) {
        count = 0;

        for (k = A->IA[i]; k < A->IA[i + 1]; ++k) {
            for (j = B->IA[A->JA[k]]; j < B->IA[A->JA[k] + 1]; ++j) {
                for (l = 0; l < count; l++) {
                    if (JD[l] == B->JA[j]) break;
                }

                if (l == count) {
                    JD[count] = B->JA[j];
                    count++;
                }
            }
        }
        C->IA[i + 1] = count;
        for (j = 0; j < count; ++j) {
            JD[j] = -1;
        }
    }

    for (i = 0; i < C->row; ++i) C->IA[i + 1] += C->IA[i];

    // step 2: Find the structure JA of C
    INT countJD;

    C->JA = (INT*)faspxx_mem_calloc(C->IA[C->row], sizeof(INT));

    for (i = 0; i < C->row; ++i) {
        countJD = 0;
        count   = C->IA[i];
        for (k = A->IA[i]; k < A->IA[i + 1]; ++k) {
            for (j = B->IA[A->JA[k]]; j < B->IA[A->JA[k] + 1]; ++j) {
                for (l = 0; l < countJD; l++) {
                    if (JD[l] == B->JA[j]) break;
                }

                if (l == countJD) {
                    C->JA[count] = B->JA[j];
                    JD[countJD]  = B->JA[j];
                    count++;
                    countJD++;
                }
            }
        }

        // for (j=0;j<countJD;++j) JD[j]=-1;
        faspxx_iarray_set(countJD, JD, -1);
    }

    faspxx_mem_free(JD);
    JD = NULL;

    // step 3: Find the structure A of C
    C->val = (DBL*)faspxx_mem_calloc(C->IA[C->row], sizeof(DBL));

    for (i = 0; i < C->row; ++i) {
        for (j = C->IA[i]; j < C->IA[i + 1]; ++j) {
            C->val[j] = 0;
            for (k = A->IA[i]; k < A->IA[i + 1]; ++k) {
                for (l = B->IA[A->JA[k]]; l < B->IA[A->JA[k] + 1]; l++) {
                    if (B->JA[l] == C->JA[j]) {
                        C->val[j] += A->val[k] * B->val[l];
                    } // end if
                }     // end for l
            }         // end for k
        }             // end for j
    }                 // end for i

    C->nnz = C->IA[C->row] - C->IA[0];
}

/**
 * \fn void faspxx_blas_dcsr_rap (const dCSRmat *R, const dCSRmat *A,
 *                              const dCSRmat *P, dCSRmat *RAP)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param R   Pointer to the dCSRmat matrix R
 * \param A   Pointer to the dCSRmat matrix A
 * \param P   Pointer to the dCSRmat matrix P
 * \param RAP Pointer to dCSRmat matrix equal to R*A*P
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package.
 *       Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void faspxx_blas_dcsr_rap(const dCSRmat* R, const dCSRmat* A, const dCSRmat* P,
                          dCSRmat* RAP)
{
    const INT  n_coarse = R->row;
    const INT* R_i      = R->IA;
    const INT* R_j      = R->JA;
    const DBL* R_data   = R->val;

    const INT  n_fine = A->row;
    const INT* A_i    = A->IA;
    const INT* A_j    = A->JA;
    const DBL* A_data = A->val;

    const INT* P_i    = P->IA;
    const INT* P_j    = P->JA;
    const DBL* P_data = P->val;

    INT  RAP_size;
    INT* RAP_i    = NULL;
    INT* RAP_j    = NULL;
    DBL* RAP_data = NULL;

#ifdef _OPENMP
    INT* P_marker = NULL;
    INT* A_marker = NULL;
#endif

    INT* Ps_marker = NULL;
    INT* As_marker = NULL;

    INT ic, i1, i2, i3, jj1, jj2, jj3;
    INT jj_counter, jj_row_begining;
    DBL r_entry, r_a_product, r_a_p_product;

    INT nthreads = 1;

#ifdef _OPENMP
    INT myid, mybegin, myend, Ctemp;
    nthreads = faspxx_get_num_threads();
#endif

    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads   = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length    = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc        = minus_one_length + coarse_add_nthreads + nthreads;

    Ps_marker = (INT*)faspxx_mem_calloc(total_calloc, sizeof(INT));
    As_marker = Ps_marker + coarse_mul_nthreads;

    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
    RAP_i = (INT*)faspxx_mem_calloc(n_coarse + 1, sizeof(INT));

    faspxx_iarray_set(minus_one_length, Ps_marker, -1);

#ifdef _OPENMP
    INT* RAP_temp = As_marker + fine_mul_nthreads;
    INT* part_end = RAP_temp + coarse_add_nthreads;

    if (n_coarse > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, Ctemp, P_marker, A_marker,      \
                                     jj_counter, ic, jj_row_begining, jj1, i1, jj2,    \
                                     i2, jj3, i3)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker   = Ps_marker + myid * n_coarse;
            A_marker   = As_marker + myid * n_fine;
            jj_counter = 0;
            for (ic = mybegin; ic < myend; ic++) {
                P_marker[ic]    = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter++;

                for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic) {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3] = jj_counter;
                                    jj_counter++;
                                }
                            }
                        }
                    }
                }

                RAP_temp[ic + myid] = jj_row_begining;
            }
            RAP_temp[myend + myid] = jj_counter;

            part_end[myid] = myend + myid + 1;
        }
        faspxx_iarray_cp(part_end[0], RAP_temp, RAP_i);
        jj_counter = part_end[0];
        Ctemp      = 0;
        for (i1 = 1; i1 < nthreads; i1++) {
            Ctemp += RAP_temp[part_end[i1 - 1] - 1];
            for (jj1 = part_end[i1 - 1] + 1; jj1 < part_end[i1]; jj1++) {
                RAP_i[jj_counter] = RAP_temp[jj1] + Ctemp;
                jj_counter++;
            }
        }
        RAP_size = RAP_i[n_coarse];
    }

    else {
#endif
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic++) {
            Ps_marker[ic]   = jj_counter;
            jj_row_begining = jj_counter;
            jj_counter++;

            for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                i1 = R_j[jj1];

                for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic) {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3] = jj_counter;
                                jj_counter++;
                            }
                        }
                    }
                }
            }

            RAP_i[ic] = jj_row_begining;
        }

        RAP_i[n_coarse] = jj_counter;
        RAP_size        = jj_counter;
#ifdef _OPENMP
    }
#endif

    RAP_j    = (INT*)faspxx_mem_calloc(RAP_size, sizeof(INT));
    RAP_data = (DBL*)faspxx_mem_calloc(RAP_size, sizeof(DBL));

    faspxx_iarray_set(minus_one_length, Ps_marker, -1);

#ifdef _OPENMP
    if (n_coarse > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, P_marker, A_marker, jj_counter, \
                                     ic, jj_row_begining, jj1, r_entry, i1, jj2,       \
                                     r_a_product, i2, jj3, r_a_p_product, i3)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker   = Ps_marker + myid * n_coarse;
            A_marker   = As_marker + myid * n_fine;
            jj_counter = RAP_i[mybegin];
            for (ic = mybegin; ic < myend; ic++) {
                P_marker[ic]         = jj_counter;
                jj_row_begining      = jj_counter;
                RAP_j[jj_counter]    = ic;
                RAP_data[jj_counter] = 0.0;
                jj_counter++;
                for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                    r_entry = R_data[jj1];

                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                        r_a_product = r_entry * A_data[jj2];

                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic) {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                                r_a_p_product = r_a_product * P_data[jj3];

                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3]         = jj_counter;
                                    RAP_data[jj_counter] = r_a_p_product;
                                    RAP_j[jj_counter]    = i3;
                                    jj_counter++;
                                } else {
                                    RAP_data[P_marker[i3]] += r_a_p_product;
                                }
                            }
                        } else {
                            for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                                i3            = P_j[jj3];
                                r_a_p_product = r_a_product * P_data[jj3];
                                RAP_data[P_marker[i3]] += r_a_p_product;
                            }
                        }
                    }
                }
            }
        }
    } else {
#endif
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic++) {
            Ps_marker[ic]        = jj_counter;
            jj_row_begining      = jj_counter;
            RAP_j[jj_counter]    = ic;
            RAP_data[jj_counter] = 0.0;
            jj_counter++;

            for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                r_entry = R_data[jj1];

                i1 = R_j[jj1];
                for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                    r_a_product = r_entry * A_data[jj2];

                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic) {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                            r_a_p_product = r_a_product * P_data[jj3];

                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3]        = jj_counter;
                                RAP_data[jj_counter] = r_a_p_product;
                                RAP_j[jj_counter]    = i3;
                                jj_counter++;
                            } else {
                                RAP_data[Ps_marker[i3]] += r_a_p_product;
                            }
                        }
                    } else {
                        for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                            i3            = P_j[jj3];
                            r_a_p_product = r_a_product * P_data[jj3];
                            RAP_data[Ps_marker[i3]] += r_a_p_product;
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
    }
#endif

    RAP->row = n_coarse;
    RAP->col = n_coarse;
    RAP->nnz = RAP_size;
    RAP->IA  = RAP_i;
    RAP->JA  = RAP_j;
    RAP->val = RAP_data;

    faspxx_mem_free(Ps_marker);
    Ps_marker = NULL;
}

/**
 * \fn void faspxx_blas_dcsr_rap_agg (const dCSRmat *R, const dCSRmat *A,
 *                                  const dCSRmat *P, dCSRmat *RAP)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P  (nonzeros of R, P = 1)
 *
 * \param R   Pointer to the dCSRmat matrix R
 * \param A   Pointer to the dCSRmat matrix A
 * \param P   Pointer to the dCSRmat matrix P
 * \param RAP Pointer to dCSRmat matrix equal to R*A*P
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_dcsr_rap_agg(const dCSRmat* R, const dCSRmat* A, const dCSRmat* P,
                              dCSRmat* RAP)
{
    const INT  n_coarse = R->row;
    const INT* R_i      = R->IA;
    const INT* R_j      = R->JA;

    const INT  n_fine = A->row;
    const INT* A_i    = A->IA;
    const INT* A_j    = A->JA;
    const DBL* A_data = A->val;

    const INT* P_i = P->IA;
    const INT* P_j = P->JA;

    INT  RAP_size;
    INT* RAP_i    = NULL;
    INT* RAP_j    = NULL;
    DBL* RAP_data = NULL;

#ifdef _OPENMP
    INT* P_marker = NULL;
    INT* A_marker = NULL;
#endif

    INT* Ps_marker = NULL;
    INT* As_marker = NULL;

    INT ic, i1, i2, i3, jj1, jj2, jj3;
    INT jj_counter, jj_row_begining;

    INT nthreads = 1;

#ifdef _OPENMP
    INT myid, mybegin, myend, Ctemp;
    nthreads = faspxx_get_num_threads();
#endif

    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads   = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length    = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc        = minus_one_length + coarse_add_nthreads + nthreads;

    Ps_marker = (INT*)faspxx_mem_calloc(total_calloc, sizeof(INT));
    As_marker = Ps_marker + coarse_mul_nthreads;

    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
    RAP_i = (INT*)faspxx_mem_calloc(n_coarse + 1, sizeof(INT));

    faspxx_iarray_set(minus_one_length, Ps_marker, -1);

#ifdef _OPENMP
    INT* RAP_temp = As_marker + fine_mul_nthreads;
    INT* part_end = RAP_temp + coarse_add_nthreads;

    if (n_coarse > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, Ctemp, P_marker, A_marker,      \
                                     jj_counter, ic, jj_row_begining, jj1, i1, jj2,    \
                                     i2, jj3, i3)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker   = Ps_marker + myid * n_coarse;
            A_marker   = As_marker + myid * n_fine;
            jj_counter = 0;
            for (ic = mybegin; ic < myend; ic++) {
                P_marker[ic]    = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter++;

                for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic) {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3] = jj_counter;
                                    jj_counter++;
                                }
                            }
                        }
                    }
                }

                RAP_temp[ic + myid] = jj_row_begining;
            }
            RAP_temp[myend + myid] = jj_counter;

            part_end[myid] = myend + myid + 1;
        }
        faspxx_iarray_cp(part_end[0], RAP_temp, RAP_i);
        jj_counter = part_end[0];
        Ctemp      = 0;
        for (i1 = 1; i1 < nthreads; i1++) {
            Ctemp += RAP_temp[part_end[i1 - 1] - 1];
            for (jj1 = part_end[i1 - 1] + 1; jj1 < part_end[i1]; jj1++) {
                RAP_i[jj_counter] = RAP_temp[jj1] + Ctemp;
                jj_counter++;
            }
        }
        RAP_size = RAP_i[n_coarse];
    }

    else {
#endif
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic++) {
            Ps_marker[ic]   = jj_counter;
            jj_row_begining = jj_counter;
            jj_counter++;

            for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                i1 = R_j[jj1];

                for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic) {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3] = jj_counter;
                                jj_counter++;
                            }
                        }
                    }
                }
            }

            RAP_i[ic] = jj_row_begining;
        }

        RAP_i[n_coarse] = jj_counter;
        RAP_size        = jj_counter;
#ifdef _OPENMP
    }
#endif

    RAP_j    = (INT*)faspxx_mem_calloc(RAP_size, sizeof(INT));
    RAP_data = (DBL*)faspxx_mem_calloc(RAP_size, sizeof(DBL));

    faspxx_iarray_set(minus_one_length, Ps_marker, -1);

#ifdef _OPENMP
    if (n_coarse > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, P_marker, A_marker, jj_counter, \
                                     ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker   = Ps_marker + myid * n_coarse;
            A_marker   = As_marker + myid * n_fine;
            jj_counter = RAP_i[mybegin];
            for (ic = mybegin; ic < myend; ic++) {
                P_marker[ic]         = jj_counter;
                jj_row_begining      = jj_counter;
                RAP_j[jj_counter]    = ic;
                RAP_data[jj_counter] = 0.0;
                jj_counter++;
                for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {

                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {

                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic) {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {

                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3]         = jj_counter;
                                    RAP_data[jj_counter] = A_data[jj2];
                                    RAP_j[jj_counter]    = i3;
                                    jj_counter++;
                                } else {
                                    RAP_data[P_marker[i3]] += A_data[jj2];
                                }
                            }
                        } else {
                            for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                                i3 = P_j[jj3];
                                RAP_data[P_marker[i3]] += A_data[jj2];
                            }
                        }
                    }
                }
            }
        }
    } else {
#endif
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic++) {
            Ps_marker[ic]        = jj_counter;
            jj_row_begining      = jj_counter;
            RAP_j[jj_counter]    = ic;
            RAP_data[jj_counter] = 0.0;
            jj_counter++;

            for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
                i1 = R_j[jj1];
                for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic) {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3]        = jj_counter;
                                RAP_data[jj_counter] = A_data[jj2];
                                RAP_j[jj_counter]    = i3;
                                jj_counter++;
                            } else {
                                RAP_data[Ps_marker[i3]] += A_data[jj2];
                            }
                        }
                    } else {
                        for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                            i3 = P_j[jj3];
                            RAP_data[Ps_marker[i3]] += A_data[jj2];
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
    }
#endif

    RAP->row = n_coarse;
    RAP->col = n_coarse;
    RAP->nnz = RAP_size;
    RAP->IA  = RAP_i;
    RAP->JA  = RAP_j;
    RAP->val = RAP_data;

    faspxx_mem_free(Ps_marker);
    Ps_marker = NULL;
}

/**
 * \fn void faspxx_blas_dcsr_rap_agg1 (const dCSRmat *R, const dCSRmat *A,
 *                                   const dCSRmat *P, dCSRmat *B)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P (nonzeros of R, P = 1)
 *
 * \param R   Pointer to the dCSRmat matrix R
 * \param A   Pointer to the dCSRmat matrix A
 * \param P   Pointer to the dCSRmat matrix P
 * \param B   Pointer to dCSRmat matrix equal to R*A*P
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package.
 *       Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void faspxx_blas_dcsr_rap_agg1(const dCSRmat* R, const dCSRmat* A, const dCSRmat* P,
                               dCSRmat* B)
{
    const INT  row = R->row, col = P->col;
    const INT *ir = R->IA, *ia = A->IA, *ip = P->IA;
    const INT *jr = R->JA, *ja = A->JA, *jp = P->JA;
    const DBL* aj = A->val;

    INT *iac, *jac;
    DBL* acj;

    INT* index  = (INT*)faspxx_mem_calloc(A->col, sizeof(INT));
    INT* iindex = (INT*)faspxx_mem_calloc(col, sizeof(INT));

    INT nB = A->nnz;
    INT i, i1, j, jj, k, length;
    INT begin_row, end_row, begin_rowA, end_rowA, begin_rowR, end_rowR;
    INT istart, iistart, count;

    // for (i=0; i<A->col; ++i) index[i] = -2;
    faspxx_iarray_set(A->col, index, -2);

    // memcpy(iindex,index,col*sizeof(INT));
    faspxx_iarray_cp(col, index, iindex);

    jac = (INT*)faspxx_mem_calloc(nB, sizeof(INT));

    iac = (INT*)faspxx_mem_calloc(row + 1, sizeof(INT));

    DBL* temp = (DBL*)faspxx_mem_calloc(A->col, sizeof(DBL));

    iac[0] = 0;

    // First loop: form sparsity partern of R*A*P
    for (i = 0; i < row; ++i) {
        // reset istart and length at the begining of each loop
        istart = -1;
        length = 0;
        i1     = i + 1;

        // go across the rows in R
        begin_rowR = ir[i];
        end_rowR   = ir[i1];
        for (jj = begin_rowR; jj < end_rowR; ++jj) {
            j = jr[jj];
            // for each column in A
            begin_rowA = ia[j];
            end_rowA   = ia[j + 1];
            for (k = begin_rowA; k < end_rowA; ++k) {
                if (index[ja[k]] == -2) {
                    index[ja[k]] = istart;
                    istart       = ja[k];
                    ++length;
                }
            }
        }

        // book-keeping [reseting length and setting iistart]
        count   = length;
        iistart = -1;
        length  = 0;

        // use each column that would have resulted from R*A
        for (j = 0; j < count; ++j) {
            jj        = istart;
            istart    = index[istart];
            index[jj] = -2;

            // go across the row of P
            begin_row = ip[jj];
            end_row   = ip[jj + 1];
            for (k = begin_row; k < end_row; ++k) {
                // pull out the appropriate columns of P
                if (iindex[jp[k]] == -2) {
                    iindex[jp[k]] = iistart;
                    iistart       = jp[k];
                    ++length;
                }
            } // end for k
        }     // end for j

        // set B->IA
        iac[i1] = iac[i] + length;

        if (iac[i1] > nB) { // Memory not enough!!!
            nB  = nB * 2;
            jac = (INT*)faspxx_mem_realloc(jac, nB * sizeof(INT));
        }

        // put the correct columns of p into the column list of the products
        begin_row = iac[i];
        end_row   = iac[i1];
        for (j = begin_row; j < end_row; ++j) {
            // put the value in B->JA
            jac[j] = iistart;
            // set istart to the next value
            iistart = iindex[iistart];
            // set the iindex spot to 0
            iindex[jac[j]] = -2;
        } // end j

    } // end i: First loop

    jac = (INT*)faspxx_mem_realloc(jac, (iac[row]) * sizeof(INT));

    acj = (DBL*)faspxx_mem_calloc(iac[row], sizeof(DBL));

    INT* BTindex = (INT*)faspxx_mem_calloc(col, sizeof(INT));

    // Second loop: compute entries of R*A*P
    for (i = 0; i < row; ++i) {
        i1 = i + 1;

        // each col of B
        begin_row = iac[i];
        end_row   = iac[i1];
        for (j = begin_row; j < end_row; ++j) {
            BTindex[jac[j]] = j;
        }

        // reset istart and length at the beginning of each loop
        istart = -1;
        length = 0;

        // go across the rows in R
        begin_rowR = ir[i];
        end_rowR   = ir[i1];
        for (jj = begin_rowR; jj < end_rowR; ++jj) {
            j = jr[jj];

            // for each column in A
            begin_rowA = ia[j];
            end_rowA   = ia[j + 1];
            for (k = begin_rowA; k < end_rowA; ++k) {
                if (index[ja[k]] == -2) {
                    index[ja[k]] = istart;
                    istart       = ja[k];
                    ++length;
                }
                temp[ja[k]] += aj[k];
            }
        }

        // book-keeping [resetting length and setting iistart]
        // use each column that would have resulted from R*A
        for (j = 0; j < length; ++j) {
            jj        = istart;
            istart    = index[istart];
            index[jj] = -2;

            // go across the row of P
            begin_row = ip[jj];
            end_row   = ip[jj + 1];
            for (k = begin_row; k < end_row; ++k) {
                // pull out the appropriate columns of P
                acj[BTindex[jp[k]]] += temp[jj];
            }
            temp[jj] = 0.0;
        }

    } // end for i: Second loop

    // setup coarse matrix B
    B->row = row;
    B->col = col;
    B->IA  = iac;
    B->JA  = jac;
    B->val = acj;
    B->nnz = B->IA[B->row] - B->IA[0];

    faspxx_mem_free(temp);
    temp = NULL;
    faspxx_mem_free(index);
    index = NULL;
    faspxx_mem_free(iindex);
    iindex = NULL;
    faspxx_mem_free(BTindex);
    BTindex = NULL;
}

/**
 * \fn void faspxx_blas_dcsr_ptap (const dCSRmat *Pt, const dCSRmat *A,
 *                               const dCSRmat *P, dCSRmat *Ac)
 *
 * \brief Triple sparse matrix multiplication B=P'*A*P
 *
 * \param Pt  Pointer to the restriction matrix
 * \param A   Pointer to the fine coefficient matrix
 * \param P   Pointer to the prolongation matrix
 * \param Ac  Pointer to the coarse coefficient matrix (output)
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note Driver to compute triple matrix product P'*A*P using ltz CSR format.
 *       In ltx format: ia[0]=1, ja[0] and a[0] are used as usual. When called
 *       from Fortran, ia[0], ja[0] and a[0] will be just ia(1),ja(1),a(1).
 *       For the indices,
 *           ia_ltz[k] = ia_usual[k]+1,
 *           ja_ltz[k] = ja_usual[k]+1,
 *            a_ltz[k] =  a_usual[k].
 */
void faspxx_blas_dcsr_ptap(const dCSRmat* Pt, const dCSRmat* A, const dCSRmat* P,
                           dCSRmat* Ac)
{
    const INT nc = Pt->row, n = Pt->col, nnzP = P->nnz, nnzA = A->nnz;
    INT       i, maxrpout;

    // shift A from usual to ltz format
#ifdef _OPENMP
#pragma omp parallel for if (n > OPENMP_HOLDS)
#endif
    for (i = 0; i <= n; ++i) {
        A->IA[i]++;
        P->IA[i]++;
    }

#ifdef _OPENMP
#pragma omp parallel for if (nnzA > OPENMP_HOLDS)
#endif
    for (i = 0; i < nnzA; ++i) {
        A->JA[i]++;
    }

#ifdef _OPENMP
#pragma omp parallel for if (nc > OPENMP_HOLDS)
#endif
    for (i = 0; i <= nc; ++i) {
        Pt->IA[i]++;
    }

#ifdef _OPENMP
#pragma omp parallel for if (nnzP > OPENMP_HOLDS)
#endif
    for (i = 0; i < nnzP; ++i) {
        P->JA[i]++;
        Pt->JA[i]++;
    }

    // compute P' A P
    dCSRmat PtAP =
        faspxx_blas_dcsr_rap2(Pt->IA, Pt->JA, Pt->val, A->IA, A->JA, A->val, Pt->IA,
                              Pt->JA, Pt->val, n, nc, &maxrpout, P->IA, P->JA);

    Ac->row = PtAP.row;
    Ac->col = PtAP.col;
    Ac->nnz = PtAP.nnz;
    Ac->IA  = PtAP.IA;
    Ac->JA  = PtAP.JA;
    Ac->val = PtAP.val;

    // shift A back from ltz format
#ifdef _OPENMP
#pragma omp parallel for if (Ac->row > OPENMP_HOLDS)
#endif
    for (i = 0; i <= Ac->row; ++i) Ac->IA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (Ac->nnz > OPENMP_HOLDS)
#endif
    for (i = 0; i < Ac->nnz; ++i) Ac->JA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (n > OPENMP_HOLDS)
#endif
    for (i = 0; i <= n; ++i) A->IA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (nnzA > OPENMP_HOLDS)
#endif
    for (i = 0; i < nnzA; ++i) A->JA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (n > OPENMP_HOLDS)
#endif
    for (i = 0; i <= n; ++i) P->IA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (nnzP > OPENMP_HOLDS)
#endif
    for (i = 0; i < nnzP; ++i) P->JA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (nc > OPENMP_HOLDS)
#endif
    for (i = 0; i <= nc; ++i) Pt->IA[i]--;

#ifdef _OPENMP
#pragma omp parallel for if (nnzP > OPENMP_HOLDS)
#endif
    for (i = 0; i < nnzP; ++i) Pt->JA[i]--;

    return;
}

/*!
 * \fn dCSRmat faspxx_blas_dcsr_rap2 (INT *ir, INT *jr, DBL *r,
 *                                  INT *ia, INT *ja, DBL *a,
 *                                  INT *ipt, INT *jpt, DBL *pt,
 *                                  INT n, INT nc,
 *                                  INT *maxrpout, INT *ipin, INT *jpin)
 *
 * \brief Compute R*A*P
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note It uses dCSRmat only. The functions called from here are in sparse_util.c.
 *       Not used for the moment!
 */
dCSRmat faspxx_blas_dcsr_rap2(INT* ir, INT* jr, DBL* r, INT* ia, INT* ja, DBL* a,
                              INT* ipt, INT* jpt, DBL* pt, INT n, INT nc, INT* maxrpout,
                              INT* ipin, INT* jpin)
{
    dCSRmat ac;
    INT     n1, nc1, nnzp, maxrp;
    INT *   ip = NULL, *jp = NULL;

    /*
     if ipin is null, this
     means that we need to do the transpose of p here; otherwise,
     these are considered to be input
     */
    maxrp = 0;
    nnzp  = ipt[nc] - 1;
    n1    = n + 1;

    if (!ipin) {
        ip = (INT*)calloc(n1, sizeof(INT));
        jp = (INT*)calloc(nnzp, sizeof(INT));
        /* these must be null anyway, so no need to assign null
         ipin=NULL;
         jpin=NULL;
         */
    } else {
        ip = ipin;
        jp = jpin;
    }

    faspxx_sparse_iit_(ipt, jpt, &nc, &n, ip, jp);

    /* triple matrix product: R * A * transpose(P^T)=R*A*P.*/
    /* A is square n by n*/
    /* Note: to compute R*A* P the input are R, A and P^T */
    /* we need to transpose now the structure of P, because the input is P^T */
    /* end of transpose of the boolean corresponding to P */
    /* ic are the addresses of the rows of the output */
    nc1   = nc + 1;
    ac.IA = (INT*)calloc(nc1, sizeof(INT));

    /*
     First call is with jc=null so that we find the number of
     nonzeros in the result
     */
    ac.JA = NULL;
    faspxx_sparse_rapms_(ir, jr, ia, ja, ip, jp, &n, &nc, ac.IA, ac.JA, &maxrp);
    ac.nnz = ac.IA[nc] - 1;
    ac.JA  = (INT*)calloc(ac.nnz, sizeof(INT));

    /*
     second call is to fill the column indexes array jc.
     */
    faspxx_sparse_rapms_(ir, jr, ia, ja, ip, jp, &n, &nc, ac.IA, ac.JA, &maxrp);
    if (!ipin) {
        if (ip) free(ip);
        if (jp) free(jp);
    }
    ac.val = (DBL*)calloc(ac.nnz, sizeof(DBL));
    /* this is the compute with the entries */
    faspxx_sparse_rapcmp_(ir, jr, r, ia, ja, a, ipt, jpt, pt, &n, &nc, ac.IA, ac.JA,
                          ac.val, &maxrp);
    ac.row = nc;
    ac.col = nc;

    /*=========================================================*/
    *maxrpout = maxrp;

    return ac;
}

/**
 * \fn void faspxx_blas_dcsr_rap4 (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B,
 *                               INT *icor_ysk)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param R   pointer to the dCSRmat matrix
 * \param A   pointer to the dCSRmat matrix
 * \param P   pointer to the dCSRmat matrix
 * \param B   pointer to dCSRmat matrix equal to R*A*P
 * \param icor_ysk pointer to the array
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package.
 *            Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void faspxx_blas_dcsr_rap4(dCSRmat* R, dCSRmat* A, dCSRmat* P, dCSRmat* B,
                           INT* icor_ysk)
{
    SHORT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP
    if (R->row > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
        const INT row = R->row, col = P->col;
        INT *     ir = R->IA, *ia = A->IA, *ip = P->IA;
        INT *     jr = R->JA, *ja = A->JA, *jp = P->JA;
        DBL *     rj = R->val, *aj = A->val, *pj = P->val;
        INT       istart, iistart;
        INT       end_row, end_rowA, end_rowR;
        INT       i, j, jj, k, length, myid, mybegin, myend;
        INT       jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
        INT*      index   = NULL;
        INT*      iindex  = NULL;
        INT*      BTindex = NULL;
        DBL*      temp    = NULL;
        INT       FiveMyid, min_A, min_P, A_pos, P_pos, FiveIc;
        INT       minus_one_length_A = icor_ysk[5 * nthreads];
        INT       minus_one_length_P = icor_ysk[5 * nthreads + 1];
        INT       minus_one_length   = minus_one_length_A + minus_one_length_P;

        INT* iindexs =
            (INT*)faspxx_mem_calloc(minus_one_length + minus_one_length_P, sizeof(INT));

#if DEBUG_MODE > 1
        total_alloc_mem += minus_one_length * sizeof(INT);
#endif
        INT* indexs   = iindexs + minus_one_length_P;
        INT* BTindexs = indexs + minus_one_length_A;

        INT* iac = (INT*)faspxx_mem_calloc(row + 1, sizeof(INT));

#if DEBUG_MODE > 1
        total_alloc_mem += (row + 1) * sizeof(INT);
#endif

        INT* part_end = (INT*)faspxx_mem_calloc(2 * nthreads + row, sizeof(INT));

#if DEBUG_MODE > 1
        total_alloc_mem += (2 * nthreads + row) * sizeof(INT);
#endif

        INT*  iac_temp     = part_end + nthreads;
        INT** iindex_array = (INT**)faspxx_mem_calloc(nthreads, sizeof(INT*));
        INT** index_array  = (INT**)faspxx_mem_calloc(nthreads, sizeof(INT*));

        faspxx_iarray_set(minus_one_length, iindexs, -2);

#ifdef _OPENMP
#pragma omp parallel for private(myid, FiveMyid, mybegin, myend, min_A, min_P, index,  \
                                     iindex, A_pos, P_pos, ic, FiveIc, jj_counter,     \
                                     jj_row_begining, end_rowR, jj1, i1, end_rowA,     \
                                     jj2, i2, end_row, jj3, i3)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            FiveMyid = myid * 5;
            mybegin  = icor_ysk[FiveMyid];
            if (myid == nthreads - 1) {
                myend = row;
            } else {
                myend = icor_ysk[FiveMyid + 5];
            }
            min_A = icor_ysk[FiveMyid + 2];
            min_P = icor_ysk[FiveMyid + 4];
            A_pos = 0;
            P_pos = 0;
            for (ic = myid - 1; ic >= 0; ic--) {
                FiveIc = ic * 5;
                A_pos += icor_ysk[FiveIc + 1];
                P_pos += icor_ysk[FiveIc + 3];
            }
            iindex_array[myid] = iindex = iindexs + P_pos - min_P;
            index_array[myid] = index = indexs + A_pos - min_A;
            jj_counter                = 0;
            for (ic = mybegin; ic < myend; ic++) {
                iindex[ic]      = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter++;
                end_rowR = ir[ic + 1];
                for (jj1 = ir[ic]; jj1 < end_rowR; jj1++) {
                    i1       = jr[jj1];
                    end_rowA = ia[i1 + 1];
                    for (jj2 = ia[i1]; jj2 < end_rowA; jj2++) {
                        i2 = ja[jj2];
                        if (index[i2] != ic) {
                            index[i2] = ic;
                            end_row   = ip[i2 + 1];
                            for (jj3 = ip[i2]; jj3 < end_row; jj3++) {
                                i3 = jp[jj3];
                                if (iindex[i3] < jj_row_begining) {
                                    iindex[i3] = jj_counter;
                                    jj_counter++;
                                }
                            }
                        }
                    }
                }
                iac_temp[ic + myid] = jj_row_begining;
            }
            iac_temp[myend + myid] = jj_counter;
            part_end[myid]         = myend + myid + 1;
        }
        faspxx_iarray_cp(part_end[0], iac_temp, iac);
        jj_counter = part_end[0];
        INT Ctemp  = 0;
        for (i1 = 1; i1 < nthreads; i1++) {
            Ctemp += iac_temp[part_end[i1 - 1] - 1];
            for (jj1 = part_end[i1 - 1] + 1; jj1 < part_end[i1]; jj1++) {
                iac[jj_counter] = iac_temp[jj1] + Ctemp;
                jj_counter++;
            }
        }
        INT* jac = (INT*)faspxx_mem_calloc(iac[row], sizeof(INT));
#if DEBUG_MODE > 1
        total_alloc_mem += iac[row] * sizeof(INT);
#endif
        faspxx_iarray_set(minus_one_length, iindexs, -2);
#ifdef _OPENMP
#pragma omp parallel for private(myid, index, iindex, FiveMyid, mybegin, myend, i,     \
                                     istart, length, i1, end_rowR, jj, j, end_rowA, k, \
                                     iistart, end_row)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            iindex   = iindex_array[myid];
            index    = index_array[myid];
            FiveMyid = myid * 5;
            mybegin  = icor_ysk[FiveMyid];
            if (myid == nthreads - 1) {
                myend = row;
            } else {
                myend = icor_ysk[FiveMyid + 5];
            }
            for (i = mybegin; i < myend; ++i) {
                istart = -1;
                length = 0;
                i1     = i + 1;
                // go across the rows in R
                end_rowR = ir[i1];
                for (jj = ir[i]; jj < end_rowR; ++jj) {
                    j = jr[jj];
                    // for each column in A
                    end_rowA = ia[j + 1];
                    for (k = ia[j]; k < end_rowA; ++k) {
                        if (index[ja[k]] == -2) {
                            index[ja[k]] = istart;
                            istart       = ja[k];
                            ++length;
                        }
                    }
                }
                // book-keeping [resetting length and setting iistart]
                // count = length;
                iistart = -1;
                // length = 0;
                //  use each column that would have resulted from R*A
                // for (j = 0; j < count; ++ j) {
                for (j = 0; j < length; ++j) {
                    jj        = istart;
                    istart    = index[istart];
                    index[jj] = -2;
                    // go across the row of P
                    end_row = ip[jj + 1];
                    for (k = ip[jj]; k < end_row; ++k) {
                        // pull out the appropriate columns of P
                        if (iindex[jp[k]] == -2) {
                            iindex[jp[k]] = iistart;
                            iistart       = jp[k];
                            //++length;
                        }
                    } // end for k
                }     // end for j
                // put the correct columns of p into the column list of the products
                end_row = iac[i1];
                for (j = iac[i]; j < end_row; ++j) {
                    // put the value in B->JA
                    jac[j] = iistart;
                    // set istart to the next value
                    iistart = iindex[iistart];
                    // set the iindex spot to 0
                    iindex[jac[j]] = -2;
                } // end j
            }
        }
        // Third loop: compute entries of R*A*P
        DBL* acj = (DBL*)faspxx_mem_calloc(iac[row], sizeof(DBL));
#if DEBUG_MODE > 1
        total_alloc_mem += iac[row] * sizeof(DBL);
#endif
        DBL* temps = (DBL*)faspxx_mem_calloc(minus_one_length_A, sizeof(DBL));
#if DEBUG_MODE > 1
        total_alloc_mem += minus_one_length_A * sizeof(DBL);
#endif

#ifdef _OPENMP
#pragma omp parallel for private(myid, index, FiveMyid, mybegin, myend, min_A, min_P,  \
                                     A_pos, P_pos, ic, FiveIc, BTindex, temp, i, i1,   \
                                     end_row, j, istart, length, end_rowR, jj,         \
                                     end_rowA, k)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            index    = index_array[myid];
            FiveMyid = myid * 5;
            mybegin  = icor_ysk[FiveMyid];
            if (myid == nthreads - 1) {
                myend = row;
            } else {
                myend = icor_ysk[FiveMyid + 5];
            }
            min_A = icor_ysk[FiveMyid + 2];
            min_P = icor_ysk[FiveMyid + 4];
            A_pos = 0;
            P_pos = 0;
            for (ic = myid - 1; ic >= 0; ic--) {
                FiveIc = ic * 5;
                A_pos += icor_ysk[FiveIc + 1];
                P_pos += icor_ysk[FiveIc + 3];
            }
            BTindex = BTindexs + P_pos - min_P;
            temp    = temps + A_pos - min_A;
            for (i = mybegin; i < myend; ++i) {
                i1 = i + 1;
                // each col of B
                end_row = iac[i1];
                for (j = iac[i]; j < end_row; ++j) {
                    BTindex[jac[j]] = j;
                }
                // reset istart and length at the beginning of each loop
                istart = -1;
                length = 0;
                // go across the rows in R
                end_rowR = ir[i1];
                for (jj = ir[i]; jj < end_rowR; ++jj) {
                    j = jr[jj];
                    // for each column in A
                    end_rowA = ia[j + 1];
                    for (k = ia[j]; k < end_rowA; ++k) {
                        if (index[ja[k]] == -2) {
                            index[ja[k]] = istart;
                            istart       = ja[k];
                            ++length;
                        }
                        temp[ja[k]] += rj[jj] * aj[k];
                    }
                }
                // book-keeping [resetting length and setting iistart]
                // use each column that would have resulted from R*A
                for (j = 0; j < length; ++j) {
                    jj        = istart;
                    istart    = index[istart];
                    index[jj] = -2;
                    // go across the row of P
                    end_row = ip[jj + 1];
                    for (k = ip[jj]; k < end_row; ++k) {
                        // pull out the appropriate columns of P
                        acj[BTindex[jp[k]]] += temp[jj] * pj[k];
                    }
                    temp[jj] = 0.0;
                }
            }
        }
        // setup coarse matrix B
        B->row = row;
        B->col = col;
        B->IA  = iac;
        B->JA  = jac;
        B->val = acj;
        B->nnz = B->IA[B->row] - B->IA[0];

        faspxx_mem_free(temps);
        temps = NULL;
        faspxx_mem_free(iindexs);
        iindexs = NULL;
        faspxx_mem_free(part_end);
        part_end = NULL;
        faspxx_mem_free(iindex_array);
        iindex_array = NULL;
        faspxx_mem_free(index_array);
        index_array = NULL;
    } else {
        faspxx_blas_dcsr_rap(R, A, P, B);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
