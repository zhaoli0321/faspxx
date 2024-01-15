/*! \file  XtrMumps.c
 *
 *  \brief Interface to MUMPS direct solvers
 *
 *  Reference for MUMPS:
 *  http://mumps.enseeiht.fr/
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include "faspxx_functs.h"

#if WITH_MUMPS
#include "dmumps_c.h"
#endif

#define ICNTL(I) icntl[(I)-1] /**< macro s.t. indices match documentation */

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT faspxx_solver_mumps (dCSRmat *ptrA, dvector *b, dvector *u,
 *                            const SHORT prtlvl)
 *
 * \brief Solve Ax=b by MUMPS directly
 *
 * \param ptrA      Pointer to a dCSRmat matrix
 * \param b         Pointer to the dvector of right-hand side term
 * \param u         Pointer to the dvector of solution
 * \param prtlvl    Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 *
 */
INT faspxx_solver_mumps(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl)
{

#if WITH_MUMPS

    DMUMPS_STRUC_C id;

    const INT n  = ptrA->row;
    const INT nz = ptrA->nnz;
    INT*      IA = ptrA->IA;
    INT*      JA = ptrA->JA;
    DBL*      AA = ptrA->val;
    DBL*      b1 = b->val;
    DBL*      x  = u->val;

    INT* irn;
    INT* jcn;
    DBL* a;
    DBL* rhs;
    INT  i, j;
    INT  begin_row, end_row;

#if DEBUG_MODE
    printf("### DEBUG: faspxx_solver_mumps ... [Start]\n");
    printf("### DEBUG: nr=%d,  nnz=%d\n", n, nz);
#endif

    // First check the matrix format
    if (IA[0] != 0 && IA[0] != 1) {
        printf("### ERROR: Matrix format is wrong -- IA[0] = %d\n", IA[0]);
        return ERROR_SOLVER_EXIT;
    }

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    /* Define A and rhs */
    irn = (INT*)malloc(sizeof(INT) * nz);
    jcn = (INT*)malloc(sizeof(INT) * nz);
    a   = (DBL*)malloc(sizeof(DBL) * nz);
    rhs = (DBL*)malloc(sizeof(DBL) * n);

    if (IA[0] == 0) { // C-convention
        for (i = 0; i < n; i++) {
            begin_row = IA[i];
            end_row   = IA[i + 1];
            for (j = begin_row; j < end_row; j++) {
                irn[j] = i + 1;
                jcn[j] = JA[j] + 1;
                a[j]   = AA[j];
            }
        }
    } else { // For-convention
        for (i = 0; i < n; i++) {
            begin_row = IA[i] - 1;
            end_row   = IA[i + 1] - 1;
            for (j = begin_row; j < end_row; j++) {
                irn[j] = i + 1;
                jcn[j] = JA[j];
                a[j]   = AA[j];
            }
        }
    }

    /* Initialize a MUMPS instance. */
    id.job          = -1;
    id.par          = 1;
    id.sym          = 0;
    id.comm_fortran = 0;
    dmumps_c(&id);

    /* Define the problem on the host */
    id.n   = n;
    id.nz  = nz;
    id.irn = irn;
    id.jcn = jcn;
    id.a   = a;
    id.rhs = rhs;

    /* No outputs */
    id.ICNTL(1) = -1;
    id.ICNTL(2) = -1;
    id.ICNTL(3) = -1;
    id.ICNTL(4) = 0;

    /* Call the MUMPS package */
    for (i = 0; i < n; i++) rhs[i] = b1[i];

    id.job = 6;
    dmumps_c(&id);

    for (i = 0; i < n; i++) x[i] = id.rhs[i];

    id.job = -2;
    dmumps_c(&id); /* Terminate instance */

    free(irn);
    free(jcn);
    free(a);
    free(rhs);

    if (prtlvl > PRINT_MIN) {
        faspxx_gettime(&end_time);
        double solve_time = end_time - start_time;
        printf("MUMPS costs %f seconds.\n", solve_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: faspxx_solver_mumps ... [Finish]\n");
#endif
    return SUCCESS;
#else

    printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;

#endif
}

#if 0
/**
 * \fn INT faspxx_solver_mumps_steps (dCSRmat *ptrA, dvector *b, dvector *u,
 *                                  Mumps_data *mumps)
 *
 * \brief Solve Ax=b by MUMPS in three steps
 *
 * \param ptrA   Pointer to a dCSRmat matrix
 * \param b      Pointer to the dvector of right-hand side term
 * \param u      Pointer to the dvector of solution
 * \param mumps  Pointer to MUMPS data
 *
 * \author Li Zhao
 * \date   01/15/2024
 *
 */
INT faspxx_solver_mumps_steps(dCSRmat* ptrA, dvector* b, dvector* u, Mumps_data* mumps)
{
#if WITH_MUMPS

    DMUMPS_STRUC_C id;

    INT job = mumps->job;

    static INT job_stat = 0;
    INT        i, j;

    INT* irn;
    INT* jcn;
    DBL* a;
    DBL* rhs;

    switch (job) {

        case 1:
            {
#if DEBUG_MODE
                printf("### DEBUG: %s, step %d, job_stat = %d... [Start]\n",
                       __FUNCTION__, job, job_stat);
#endif
                INT       begin_row, end_row;
                const INT n  = ptrA->row;
                const INT nz = ptrA->nnz;
                INT*      IA = ptrA->IA;
                INT*      JA = ptrA->JA;
                DBL*      AA = ptrA->val;

                irn = id.irn = (INT*)malloc(sizeof(INT) * nz);
                jcn = id.jcn = (INT*)malloc(sizeof(INT) * nz);
                a = id.a = (DBL*)malloc(sizeof(DBL) * nz);
                rhs = id.rhs = (DBL*)malloc(sizeof(DBL) * n);

                // First check the matrix format
                if (IA[0] != 0 && IA[0] != 1) {
                    printf("### ERROR: Matrix format is wrong, IA[0] = %d!\n", IA[0]);
                    return ERROR_SOLVER_EXIT;
                }

                // Define A and rhs
                if (IA[0] == 0) { // C-convention
                    for (i = 0; i < n; i++) {
                        begin_row = IA[i];
                        end_row   = IA[i + 1];
                        for (j = begin_row; j < end_row; j++) {
                            irn[j] = i + 1;
                            jcn[j] = JA[j] + 1;
                            a[j]   = AA[j];
                        }
                    }
                } else { // For-convention
                    for (i = 0; i < n; i++) {
                        begin_row = IA[i] - 1;
                        end_row   = IA[i + 1] - 1;
                        for (j = begin_row; j < end_row; j++) {
                            irn[j] = i + 1;
                            jcn[j] = JA[j];
                            a[j]   = AA[j];
                        }
                    }
                }

                /* Initialize a MUMPS instance. */
                id.job          = -1;
                id.par          = 1;
                id.sym          = 0;
                id.comm_fortran = 0;
                dmumps_c(&id);
                /* Define the problem on the host */
                id.n   = n;
                id.nz  = nz;
                id.irn = irn;
                id.jcn = jcn;
                id.a   = a;
                id.rhs = rhs;

                /* No outputs */
                id.ICNTL(1) = -1;
                id.ICNTL(2) = -1;
                id.ICNTL(3) = -1;
                id.ICNTL(4) = 0;

                id.job = 4;
                dmumps_c(&id);
                job_stat = 1;

                mumps->id = id;

#if DEBUG_MODE
                printf("### DEBUG: %s, step %d, job_stat = %d... [Finish]\n",
                       __FUNCTION__, job, job_stat);
#endif
                break;
            }

        case 2:
            {
#if DEBUG_MODE
                printf("### DEBUG: %s, step %d, job_stat = %d... [Start]\n",
                       __FUNCTION__, job, job_stat);
#endif
                id = mumps->id;

                if (job_stat != 1)
                    printf("### ERROR: %s setup failed!\n", __FUNCTION__);

                /* Call the MUMPS package. */
                for (i = 0; i < id.n; i++) id.rhs[i] = b->val[i];

                id.job = 3;
                dmumps_c(&id);

                for (i = 0; i < id.n; i++) u->val[i] = id.rhs[i];

#if DEBUG_MODE
                printf("### DEBUG: %s, step %d, job_stat = %d... [Finish]\n",
                       __FUNCTION__, job, job_stat);
#endif
                break;
            }

        case 3:
            {
#if DEBUG_MODE
                printf("### DEBUG: %s, step %d, job_stat = %d... [Start]\n",
                       __FUNCTION__, job, job_stat);
#endif
                id = mumps->id;

                if (job_stat != 1)
                    printf("### ERROR: %s setup failed!\n", __FUNCTION__);

                free(id.irn);
                free(id.jcn);
                free(id.a);
                free(id.rhs);
                id.job = -2;
                dmumps_c(&id); /* Terminate instance */

#if DEBUG_MODE
                printf("### DEBUG: %s, step %d, job_stat = %d... [Finish]\n",
                       __FUNCTION__, job, job_stat);
#endif

                break;
            }

        default:
            printf("### ERROR: Parameter job = %d. Should be 1, 2, or 3!\n", job);
            return ERROR_SOLVER_EXIT;
    }

    return SUCCESS;

#else

    printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;

#endif
}
#endif

#if WITH_MUMPS
/**
 * \fn DMUMPS_STRUC_C faspxx_mumps_factorize (dCSRmat *ptrA, dvector *b, dvector *u,
 *                                          const SHORT prtlvl)
 * \brief Factorize A by MUMPS
 *
 * \param ptrA     Pointer to stiffness matrix of levelNum levels
 * \param b        Pointer to the dvector of right hand side term
 * \param u        Pointer to the dvector of dofs
 * \param prtlvl   output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
Mumps_data faspxx_mumps_factorize(dCSRmat* ptrA, dvector* b, dvector* u,
                                  const SHORT prtlvl)
{
    Mumps_data     mumps;
    DMUMPS_STRUC_C id;

    INT       i, j;
    const INT m  = ptrA->row;
    const INT n  = ptrA->col;
    const INT nz = ptrA->nnz;
    INT*      IA = ptrA->IA;
    INT*      JA = ptrA->JA;
    DBL*      AA = ptrA->val;

    INT* irn = id.irn = (INT*)malloc(sizeof(INT) * nz);
    INT* jcn = id.jcn = (INT*)malloc(sizeof(INT) * nz);
    DBL* a = id.a = (DBL*)malloc(sizeof(DBL) * nz);
    DBL* rhs = id.rhs = (DBL*)malloc(sizeof(DBL) * n);

    INT begin_row, end_row;

#if DEBUG_MODE
    printf("### DEBUG: %s ... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nz);
#endif

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    if (IA[0] == 0) { // C-convention
        for (i = 0; i < n; i++) {
            begin_row = IA[i];
            end_row   = IA[i + 1];
            for (j = begin_row; j < end_row; j++) {
                irn[j] = i + 1;
                jcn[j] = JA[j] + 1;
                a[j]   = AA[j];
            }
        }
    } else { // For-convention
        for (i = 0; i < n; i++) {
            begin_row = IA[i] - 1;
            end_row   = IA[i + 1] - 1;
            for (j = begin_row; j < end_row; j++) {
                irn[j] = i + 1;
                jcn[j] = JA[j];
                a[j]   = AA[j];
            }
        }
    }

    /* Initialize a MUMPS instance. */
    id.job          = -1;
    id.par          = 1;
    id.sym          = 0;
    id.comm_fortran = 0;
    dmumps_c(&id);
    /* Define the problem on the host */
    id.n   = n;
    id.nz  = nz;
    id.irn = irn;
    id.jcn = jcn;
    id.a   = a;
    id.rhs = rhs;

    /* No outputs */
    id.ICNTL(1) = -1;
    id.ICNTL(2) = -1;
    id.ICNTL(3) = -1;
    id.ICNTL(4) = 0;

    id.job = 4;
    dmumps_c(&id);

    if (prtlvl > PRINT_MIN) {
        faspxx_gettime(&end_time);
        double fac_time = end_time - start_time;
        printf("MUMPS factorize costs %f seconds.\n", fac_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ... [Finish]\n", __FUNCTION__);
#endif

    mumps.id = id;

    return mumps;
}

/**
 * \fn INT faspxx_mumps_solve (dCSRmat *ptrA, dvector *b, dvector *u,
 *                            Mumps_data mumps, const SHORT prtlvl)
 * \brief Solve A by MUMPS
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param b         Pointer to the dvector of right hand side term
 * \param u         Pointer to the dvector of dofs
 * \param mumps     Pointer to mumps data
 * \param prtlvl    Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_mumps_solve(dCSRmat* ptrA, dvector* b, dvector* u, Mumps_data mumps,
                       const SHORT prtlvl)
{
    INT i, j;
    INT status = SUCCESS;

    DMUMPS_STRUC_C id = mumps.id;

    const INT m  = ptrA->row;
    const INT n  = ptrA->row;
    const INT nz = ptrA->nnz;
    INT*      IA = ptrA->IA;
    INT*      JA = ptrA->JA;
    DBL*      AA = ptrA->val;

    INT* irn = id.irn;
    INT* jcn = id.jcn;
    DBL* a   = id.a;
    DBL* rhs = id.rhs;

#if DEBUG_MODE
    printf("### DEBUG: %s ... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nz);
#endif

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    DBL* b1 = b->val;
    DBL* x  = u->val;

    /* Call the MUMPS package. */
    for (i = 0; i < id.n; i++) rhs[i] = b1[i];

    id.job = 3;
    dmumps_c(&id);

    for (i = 0; i < id.n; i++) x[i] = id.rhs[i];

    if (prtlvl > PRINT_NONE) {
        faspxx_gettime(&end_time);
        double solve_time = end_time - start_time;
        printf("MUMPS costs %f seconds.\n", solve_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ... [Finish]\n", __FUNCTION__);
#endif
    return status;
}

/**
 * \fn INT faspxx_mumps_free (Mumps_data *mumps)
 *
 * \brief Free MUMPS memory
 *
 * \param mumps   Pointer to mumps data
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_mumps_free(Mumps_data* mumps)
{
    INT            status = SUCCESS;
    DMUMPS_STRUC_C id     = mumps->id;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    free(id.irn);
    free(id.jcn);
    free(id.a);
    free(id.rhs);

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}
#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
