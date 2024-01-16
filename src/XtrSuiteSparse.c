/*! \file  XtrSuiteSparse.c
 *
 *  \brief Interface to SuiteSparse (or UMFpack) direct solvers
 *
 *  Reference for SuiteSparse:
 *  http://faculty.cse.tamu.edu/davis/suitesparse.html
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include "faspxx_functs.h"

#if WITH_SUITESPARSE
#include "umfpack.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT faspxx_solver_suitesparse (dCSRmat *ptrA, dvector *b, dvector *u,
 *                              const SHORT prtlvl)
 *
 * \brief Solve Au=b by SuiteSparse (or UMFpack)
 *
 * \param ptrA         Pointer to a dCSRmat matrix
 * \param b            Pointer to the dvector of right-hand side term
 * \param u            Pointer to the dvector of solution
 * \param prtlvl       Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_solver_suitesparse(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl)
{

#if WITH_SUITESPARSE

    const INT n   = ptrA->col;
    const INT m   = ptrA->row;
    const INT nnz = ptrA->nnz;

    INT*    Ap = ptrA->IA;
    INT*    Ai = ptrA->JA;
    double* Ax = ptrA->val;
    void *  Symbolic, *Numeric;
    INT     status = SUCCESS;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    umfpack_di_free_symbolic(&Symbolic);
    status =
        umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, u->val, b->val, Numeric, NULL, NULL);
    umfpack_di_free_numeric(&Numeric);

    if (prtlvl > PRINT_MIN) {
        faspxx_gettime(&end_time);
        double solve_time = end_time - start_time;
        printf("SuiteSparse (or UMFPACK) costs %f seconds.\n", solve_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;

#else

    printf("### ERROR: SuiteSparse (or UMFPACK) is not available!\n");
    return ERROR_SOLVER_EXIT;

#endif
}

#if WITH_SUITESPARSE
/**
 * \fn void* faspxx_suitesparse_factorize (dCSRmat *ptrA, const SHORT prtlvl)
 * \brief factorize A by SuiteSparse (or UMFpack)
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param Numeric   Pointer to the numerical factorization
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
void* faspxx_suitesparse_factorize(dCSRmat* ptrA, const SHORT prtlvl)
{
    const INT n   = ptrA->col;
    const INT m   = ptrA->row;
    const INT nnz = ptrA->nnz;

    INT*    Ap = ptrA->IA;
    INT*    Ai = ptrA->JA;
    double* Ax = ptrA->val;
    void*   Symbolic;
    void*   Numeric;
    INT     status = SUCCESS;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    umfpack_di_free_symbolic(&Symbolic);

    if (prtlvl > PRINT_MIN) {
        faspxx_gettime(&end_time);
        double factorize_time = end_time - start_time;
        printf("SuiteSparse (or UMFPACK) factorize costs %f seconds.\n",
               factorize_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return Numeric;
}

/**
 * \fn INT faspxx_suitesparse_solve (dCSRmat *ptrA, dvector *b, dvector *u,
 *                             void *Numeric, const SHORT prtlvl)
 * \brief Solve Au=b by SuiteSparse (or UMFpack), numerical factorization is given
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param b         Pointer to the dvector of right hand side term
 * \param u         Pointer to the dvector of dofs
 * \param Numeric   Pointer to the numerical factorization
 * \param prtlvl    Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_suitesparse_solve(dCSRmat* ptrA, dvector* b, dvector* u, void* Numeric,
                             const SHORT prtlvl)
{
    const INT n   = ptrA->col;
    const INT m   = ptrA->row;
    const INT nnz = ptrA->nnz;

    INT*    Ap     = ptrA->IA;
    INT*    Ai     = ptrA->JA;
    double* Ax     = ptrA->val;
    INT     status = SUCCESS;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    status =
        umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, u->val, b->val, Numeric, NULL, NULL);

    if (prtlvl > PRINT_NONE) {
        faspxx_gettime(&end_time);
        double solve_time = end_time - start_time;
        printf("SuiteSparse (or UMFPACK) costs %f seconds.\n", solve_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn INT faspxx_suitesparse_free (void *Numeric)
 * \brief Solve Au=b by SuiteSparse (or UMFpack)
 *
 * \param Numeric   Pointer to the numerical factorization
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_suitesparse_free(void* Numeric)
{
    INT status = SUCCESS;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    umfpack_di_free_numeric(&Numeric);

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
