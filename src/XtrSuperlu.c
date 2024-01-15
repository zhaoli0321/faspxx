/*! \file  XtrSuperlu.c
 *
 *  \brief Interface to SuperLU direct solvers
 *
 *  Reference for SuperLU:
 *  http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include "faspxx_functs.h"

#if WITH_SUPERLU
#include "slu_ddefs.h"

static void faspxx_superlu_factorize_internal(SuperLU_data* superlu_data);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT faspxx_solver_superlu (dCSRmat *ptrA, dvector *b, dvector *u,
 *                              const SHORT prtlvl)
 *
 * \brief Solve Au=b by SuperLU
 *
 * \param ptrA      Pointer to a dCSRmat matrix
 * \param b         Pointer to the dvector of right-hand side term
 * \param u         Pointer to the dvector of solution
 * \param prtlvl    Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 *
 * \note  Factorization and solution are combined together!!! Not efficient!!!
 */
INT faspxx_solver_superlu(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl)
{

#if WITH_SUPERLU

    SuperMatrix A, L, U, B;

    INT* perm_r; /* row permutations from partial pivoting */
    INT* perm_c; /* column permutation vector */
    INT  nrhs = 1, info, m = ptrA->row, n = ptrA->col, nnz = ptrA->nnz;

    // if (prtlvl > PRINT_NONE) printf("superlu: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);

    clock_t LU_start = clock();

    dCSRmat tempA = faspxx_dcsr_create(m, n, nnz);
    faspxx_dcsr_cp(ptrA, &tempA);

    dvector tempb = faspxx_dvec_create(n);
    faspxx_dvec_cp(b, &tempb);

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, tempA.val, tempA.JA, tempA.IA, SLU_NR, SLU_D,
                           SLU_GE);

    /* Create right-hand side B. */
    dCreate_Dense_Matrix(&B, m, nrhs, tempb.val, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = COLAMD; // MMD_AT_PLUS_A; MMD_ATA; NATURAL;

    /* Initialize the statistics variables. */
    SuperLUStat_t stat;
    StatInit(&stat);

    /* SuperLU */
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    DNformat* BB = (DNformat*)B.Store;
    u->val       = (double*)BB->nzval;
    u->row       = n;

    if (prtlvl > PRINT_MIN) {
        clock_t LU_end  = clock();
        double  LU_time = (double)(LU_end - LU_start) / (double)(CLOCKS_PER_SEC);
        printf("SuperLU totally costs %f seconds.\n", LU_time);
    }

    /* De-allocate storage */
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    return SUCCESS;

#else

    printf("### ERROR: SuperLU is not available!\n");
    return ERROR_SOLVER_EXIT;

#endif
}

#if WITH_SUPERLU
/**
 * \fn INT faspxx_superlu_factorize(dCSRmat* ptrA, dvector* b, SuperLU_data*
                                    superlu_data, const SHORT prtlvl)

 * \brief factorize A by SuperLU.
 *
 * \param ptrA          Pointer to matrix A
 * \param b             Pointer to rhs b
 * \param superlu_data  Pointer to numerical factorization data
 * \param prtlvl        Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_superlu_factorize(dCSRmat* ptrA, dvector* b, SuperLU_data* superlu_data,
                             const SHORT prtlvl)
{
    INT status = SUCCESS;
    DBL start_time, end_time;
    faspxx_gettime(&start_time);

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif

    SuperMatrix A, B;

    INT* perm_r; /* row permutations from partial pivoting */
    INT* perm_c; /* column permutation vector */
    INT  nrhs = 1, m = ptrA->row, n = ptrA->col, nnz = ptrA->nnz;

    // if (prtlvl > PRINT_NONE) printf("superlu: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);

    dCSRmat tempA = faspxx_dcsr_create(m, n, nnz);
    faspxx_dcsr_cp(ptrA, &tempA);

    dvector tempb = faspxx_dvec_create(n);
    faspxx_dvec_cp(b, &tempb);

    //! Note that matirx A is a pointer variable, rhs B also is a pointer variable
    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, tempA.val, tempA.JA, tempA.IA, SLU_NR, SLU_D,
                           SLU_GE);

    /* Create right-hand side B. */
    dCreate_Dense_Matrix(&B, m, nrhs, tempb.val, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = COLAMD; // MMD_AT_PLUS_A; MMD_ATA; NATURAL;

    /* Initialize the statistics variables. */
    SuperLUStat_t stat;
    StatInit(&stat);

    //! Pointer
    superlu_data->perm_r  = perm_r;
    superlu_data->perm_c  = perm_c;
    superlu_data->stat    = stat;
    superlu_data->options = options;
    superlu_data->A       = A;
    superlu_data->B       = B;

    faspxx_superlu_factorize_internal(superlu_data);

    if (prtlvl > PRINT_MIN) {
        faspxx_gettime(&end_time);
        double factorize_time = end_time - start_time;
        printf("SuperLU factoization costs %f seconds.\n", factorize_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn INT faspxx_superlu_solve (dCSRmat *ptrA, dvector *b, dvector *u,
 *                             Pardiso_data *pdata, const SHORT prtlvl)
 * \brief Solve Au=b by Intel MKL PARDISO, numerical factorization is given.
 *        Each row of A should be in ascending order w.r.t. column indices.
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param b         Pointer to the dvector of right hand side term
 * \param u         Pointer to the dvector of dofs
 * \param pdata     Pointer to the numerical factorization data
 * \param prtlvl    Output level
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_superlu_solve(dCSRmat* ptrA, dvector* b, dvector* u,
                         SuperLU_data* superlu_data, const SHORT prtlvl)
{
    INT            status = SUCCESS;
    INT*           perm_c = superlu_data->perm_c; // input/output
    INT*           perm_r = superlu_data->perm_r; // input/output
    SuperMatrix*   L      = &superlu_data->L;     // output
    SuperMatrix*   U      = &superlu_data->U;     // output
    SuperMatrix*   B      = &superlu_data->B;     // input/output
    SuperLUStat_t* stat   = &superlu_data->stat;  // output
    INT*           info   = &superlu_data->info;  // output
    trans_t        trans  = superlu_data->trans;  // output

    DBL start_time, end_time;
    faspxx_gettime(&start_time);

    /* Solve the system A*X=B, overwriting B with X. */
    dgstrs(trans, L, U, perm_c, perm_r, B, stat, info);

    if (*info != 0) {
        printf("### ERROR: Solution failed %d!\n", *info);
        exit(3);
    }

    DNformat* BB = (DNformat*)B->Store;
    u->val       = (double*)BB->nzval;
    u->row       = B->nrow;

    if (prtlvl > PRINT_MIN) {
        faspxx_gettime(&end_time);
        double solve_time = end_time - start_time;
        printf("SuperLU costs %f seconds.\n", solve_time);
    }

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn INT faspxx_superlu_free (Pardiso_data *pdata)
 * \brief Free internal solver memory for PARDISO
 *
 * \param  pdata  Pointer to the numerical factorization data
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
INT faspxx_superlu_free(SuperLU_data* superlu_data)
{
    INT status = SUCCESS;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    INT*           perm_c = superlu_data->perm_c;
    INT*           perm_r = superlu_data->perm_r;
    SuperMatrix*   A      = &superlu_data->A;
    SuperMatrix*   L      = &superlu_data->L;
    SuperMatrix*   U      = &superlu_data->U;
    SuperMatrix*   B      = &superlu_data->B;
    SuperLUStat_t* stat   = &superlu_data->stat;

    /* De-allocate storage */
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(A);
    Destroy_SuperMatrix_Store(B);
    Destroy_SuperNode_Matrix(L);
    Destroy_CompCol_Matrix(U);
    StatFree(stat);

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 *
 * \fn static void faspxx_superlu_factorize_internal(SuperLU_data* superlu_data)
 *
 * \brief using the LU factorization from DGSTRF. It performs
 * the following steps:
 *
 *   1. If A is stored column-wise (A->Stype = SLU_NC):
 *
 *      1.1. Permute the columns of A, forming A*Pc, where Pc
 *           is a permutation matrix. For more details of this step,
 *           see sp_preorder.c.
 *
 *      1.2. Factor A as Pr*A*Pc=L*U with the permutation Pr determined
 *           by Gaussian elimination with partial pivoting.
 *           L is unit lower triangular with offdiagonal entries
 *           bounded by 1 in magnitude, and U is upper triangular.
 *
 *      1.3. Solve the system of equations A*X=B using the factored
 *           form of A.
 *
 *   2. If A is stored row-wise (A->Stype = SLU_NR), apply the
 *      above algorithm to the transpose of A:
 *
 *      2.1. Permute columns of transpose(A) (rows of A),
 *           forming transpose(A)*Pc, where Pc is a permutation matrix.
 *           For more details of this step, see sp_preorder.c.
 *
 *      2.2. Factor A as Pr*transpose(A)*Pc=L*U with the permutation Pr
 *           determined by Gaussian elimination with partial pivoting.
 *           L is unit lower triangular with offdiagonal entries
 *           bounded by 1 in magnitude, and U is upper triangular.
 *
 *      2.3. Solve the system of equations A*X=B using the factored
 *           form of A.
 *
 *   See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * options (input) superlu_options_t*
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed and how the
 *         system will be solved.
 *
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of linear equations is A->nrow. Currently, the type of A can be:
 *         Stype = SLU_NC or SLU_NR; Dtype = SLU_D; Mtype = SLU_GE.
 *         In the future, more general A may be handled.
 *
 * perm_c  (input/output) INT*
 *         If A->Stype = SLU_NC, column permutation vector of size A->ncol
 *         which defines the permutation matrix Pc; perm_c[i] = j means
 *         column i of A is in position j in A*Pc.
 *         If A->Stype = SLU_NR, column permutation vector of size A->nrow
 *         which describes permutation of columns of transpose(A)
 *         (rows of A) as described above.
 *
 *         If options->ColPerm = MY_PERMC or options->Fact = SamePattern or
 *            options->Fact = SamePattern_SameRowPerm, it is an input argument.
 *            On exit, perm_c may be overwritten by the product of the input
 *            perm_c and a permutation that postorders the elimination tree
 *            of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
 *            is already in postorder.
 *         Otherwise, it is an output argument.
 *
 * perm_r  (input/output) INT*
 *         If A->Stype = SLU_NC, row permutation vector of size A->nrow,
 *         which defines the permutation matrix Pr, and is determined
 *         by partial pivoting.  perm_r[i] = j means row i of A is in
 *         position j in Pr*A.
 *         If A->Stype = SLU_NR, permutation vector of size A->ncol, which
 *         determines permutation of rows of transpose(A)
 *         (columns of A) as described above.
 *
 *         If options->RowPerm = MY_PERMR or
 *            options->Fact = SamePattern_SameRowPerm, perm_r is an
 *            input argument.
 *         otherwise it is an output argument.
 *
 * L       (output) SuperMatrix*
 *         The factor L from the factorization
 *             Pr*A*Pc=L*U              (if A->Stype = SLU_NC) or
 *             Pr*transpose(A)*Pc=L*U   (if A->Stype = SLU_NR).
 *         Uses compressed row subscripts storage for supernodes, i.e.,
 *         L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
 *
 * U       (output) SuperMatrix*
 *	   The factor U from the factorization
 *             Pr*A*Pc=L*U              (if A->Stype = SLU_NC) or
 *             Pr*transpose(A)*Pc=L*U   (if A->Stype = SLU_NR).
 *         Uses column-wise storage scheme, i.e., U has types:
 *         Stype = SLU_NC, Dtype = SLU_D, Mtype = SLU_TRU.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics on runtime and floating-point operation count.
 *        See util.h for the definition of 'SuperLUStat_t'.
 *
 * info    (output) INT*
 *	   = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *
 * \author Li Zhao
 * \date   01/15/2024
 */
static void faspxx_superlu_factorize_internal(SuperLU_data* superlu_data)
{

    superlu_options_t* options = &superlu_data->options; // input
    SuperMatrix*       A       = &superlu_data->A;       // input
    INT*               perm_c  = superlu_data->perm_c;   // input/output
    INT*               perm_r  = superlu_data->perm_r;   // input/output
    SuperMatrix*       L       = &superlu_data->L;       // output
    SuperMatrix*       U       = &superlu_data->U;       // output
    SuperMatrix*       B       = &superlu_data->B;       // input/output
    SuperLUStat_t*     stat    = &superlu_data->stat;    // output
    INT*               info    = &superlu_data->info;    // output
    trans_t            trans   = superlu_data->trans;    // output

    // local variables
    DNformat*    Bstore;
    SuperMatrix* AA; /* A in SLU_NC format used by the factorization routine.*/
    SuperMatrix  AC; /* Matrix postmultiplied by Pc */
    INT          lwork = 0, *etree, i;
    GlobalLU_t   Glu; /* Not needed on return. */

    /* Set default values for some parameters */
    INT panel_size; /* panel size */
    INT relax;      /* no of columns in a relaxed snodes */
    INT permc_spec;
    // trans_t trans = NOTRANS;
    trans = NOTRANS;
    double* utime;
    double  t; /* Temporary time */

    /* Test the input parameters ... */
    *info  = 0;
    Bstore = B->Store;
    if (options->Fact != DOFACT)
        *info = -1;
    else if (A->nrow != A->ncol || A->nrow < 0 ||
             (A->Stype != SLU_NC && A->Stype != SLU_NR) || A->Dtype != SLU_D ||
             A->Mtype != SLU_GE)
        *info = -2;
    else if (B->ncol < 0 || Bstore->lda < SUPERLU_MAX(0, A->nrow) ||
             B->Stype != SLU_DN || B->Dtype != SLU_D || B->Mtype != SLU_GE)
        *info = -7;
    if (*info != 0) {
        i = -(*info);
        input_error("faspxx_superlu_factorize_internal", &i);
        return;
    }

    utime = stat->utime;

    /* Convert A to SLU_NC format when necessary. */
    if (A->Stype == SLU_NR) {
        NRformat* Astore = A->Store;
        AA               = (SuperMatrix*)SUPERLU_MALLOC(sizeof(SuperMatrix));
        dCreate_CompCol_Matrix(AA, A->ncol, A->nrow, Astore->nnz, Astore->nzval,
                               Astore->colind, Astore->rowptr, SLU_NC, A->Dtype,
                               A->Mtype);
        trans = TRANS;
    } else {
        if (A->Stype == SLU_NC) AA = A;
    }

    t = SuperLU_timer_();
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = NATURAL:  natural ordering
     *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
     *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
     *   permc_spec = COLAMD:   approximate minimum degree column ordering
     *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
     */
    permc_spec = options->ColPerm;
    if (permc_spec != MY_PERMC && options->Fact == DOFACT)
        get_perm_c(permc_spec, AA, perm_c);
    utime[COLPERM] = SuperLU_timer_() - t;

    etree = intMalloc(A->ncol);

    t = SuperLU_timer_();
    sp_preorder(options, AA, perm_c, etree, &AC);
    utime[ETREE] = SuperLU_timer_() - t;

    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);

    /*printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n",
      relax, panel_size, sp_ienv(3), sp_ienv(4));*/
    t = SuperLU_timer_();
    /* Compute the LU factorization of A. */
    dgstrf(options, &AC, relax, panel_size, etree, NULL, lwork, perm_c, perm_r, L, U,
           &Glu, stat, info);
    utime[FACT] = SuperLU_timer_() - t;

    if (*info != 0) {
        printf("\n### ERROR: SuperLU factorization failed %d!\n", *info);
        exit(1);
    }

    SUPERLU_FREE(etree);
    Destroy_CompCol_Permuted(&AC);
    if (A->Stype == SLU_NR) {
        Destroy_SuperMatrix_Store(AA);
        SUPERLU_FREE(AA);
    }
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
