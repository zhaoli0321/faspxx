#include "faspxx.h"
#include "faspxx_functs.h"

#define FACTSOLCOMB FALSE
// TRUE:  Factorization and solution are combined together
// FALSE: Separate factorization and solution independently

int main()
{
    dCSRmat A;
    dvector b, x, u, r;
    SHORT   PrtLvl = PRINT_MORE;
    DBL     errL2 = LARGE_DBL, residual_norm2 = LARGE_DBL;
    DBL     CheckTOL = 1e-10;
    INT     status;

    printf("Direct Solver (Pardiso, Mumps, SuperLU, Umfpack) Test...\n");
#if FACTSOLCOMB
    printf("Factorization and solution are combined together!\n\n");
#else
    printf("Separate factorization and solution independently!\n\n");
#endif

    const char csrfilename[512] = "../data/fem_small.dat";
    faspxx_dcsr_read(csrfilename, &A);

    u = faspxx_dvec_create(A.row); // exact solution
    faspxx_dvec_set(u.row, &u, 1.0);
    printf("A.row: %d, A.col: %d, A.nnz: %d\n", A.row, A.col, A.nnz);

    b = faspxx_dvec_create(A.row);
    faspxx_blas_dcsr_mxv(&A, u.val, b.val);

    x = faspxx_dvec_create(A.row);
    faspxx_dvec_set(x.row, &x, 0.0);

    r = faspxx_dvec_create(A.row);

#if WITH_PARDISO // use Pardiso directly
    printf("\nPardiso Test...\n");
    faspxx_dcsr_sort(&A);

#if FACTSOLCOMB
    status = faspxx_solver_pardiso(&A, &b, &x, PrtLvl);
#else
    Pardiso_data pdata;
    status = faspxx_pardiso_factorize(&A, &pdata, PrtLvl);     // factorize
    status = faspxx_pardiso_solve(&A, &b, &x, &pdata, PrtLvl); // solve
    status = faspxx_pardiso_free(&pdata);                      // free
#endif

    //! Checking
    // ||x-u||
    errL2 = faspxx_dvec_norm2diff(&u, &x);
    // u = b-A*x
    faspxx_darray_cp(b.row, b.val, r.val);
    faspxx_blas_dcsr_aAxpy(-1.0, &A, x.val, r.val);
    // ||b-Ax||
    residual_norm2 = faspxx_blas_darray_norm2(r.row, r.val);
    printf("||x-u|| = %e, ||b-Ax|| = %e\n", errL2, residual_norm2);

    if (errL2 < CheckTOL && residual_norm2 < CheckTOL) printf("Pardiso Pass!\n");
#endif

#if WITH_MUMPS // use Mumps directly
    printf("\nMumps Test...\n");

#if FACTSOLCOMB
    status = faspxx_solver_mumps(&A, &b, &x, PrtLvl);
#else
    Mumps_data mumps_data;
    mumps_data = faspxx_mumps_factorize(&A, &b, &x, PrtLvl);
    status     = faspxx_mumps_solve(&A, &b, &x, mumps_data, PrtLvl);
    status     = faspxx_mumps_free(&mumps_data);
#endif

    //! Checking
    // ||x-u||
    errL2 = faspxx_dvec_norm2diff(&u, &x);
    // u = b-A*x
    faspxx_darray_cp(b.row, b.val, r.val);
    faspxx_blas_dcsr_aAxpy(-1.0, &A, x.val, r.val);
    // ||b-Ax||
    residual_norm2 = faspxx_blas_darray_norm2(r.row, r.val);
    printf("||x-u|| = %e, ||b-Ax|| = %e\n", errL2, residual_norm2);

    if (errL2 < CheckTOL && residual_norm2 < CheckTOL) printf("Mumps Pass!\n");
#endif

#if WITH_SUPERLU // use SuperLU directly
    printf("\nSuperLU Test...\n");

#if FACTSOLCOMB
    status = faspxx_solver_superlu(&A, &b, &x, PrtLvl);
#else
    SuperLU_data superlu_data;
    status = faspxx_superlu_factorize(&A, &b, &superlu_data, PrtLvl);
    status = faspxx_superlu_solve(&A, &b, &x, &superlu_data, PrtLvl);
    status = faspxx_superlu_free(&superlu_data);
#endif

    //! Checking
    // ||x-u||
    errL2 = faspxx_dvec_norm2diff(&u, &x);
    // u = b-A*x
    faspxx_darray_cp(b.row, b.val, r.val);
    faspxx_blas_dcsr_aAxpy(-1.0, &A, x.val, r.val);
    // ||b-Ax||
    residual_norm2 = faspxx_blas_darray_norm2(r.row, r.val);
    printf("||x-u|| = %e, ||b-Ax|| = %e\n", errL2, residual_norm2);

    if (errL2 < CheckTOL && residual_norm2 < CheckTOL) printf("SuperLU Pass!\n");
#endif

#if WITH_UMFPACK // use UMFPACK directly
    printf("\nUMFPACK Test...\n");
    dCSRmat A_tran;
    faspxx_dcsr_trans(&A, &A_tran);
    faspxx_dcsr_sort(&A_tran);
    faspxx_dcsr_cp(&A_tran, &A);
    faspxx_dcsr_free(&A_tran);

#if FACTSOLCOMB
    status = faspxx_solver_umfpack(&A, &b, &x, PrtLvl);
#else
    void* Numeric = faspxx_umfpack_factorize(&A, PrtLvl);              // factorize
    status        = faspxx_umfpack_solve(&A, &b, &x, Numeric, PrtLvl); // solve
    status        = faspxx_umfpack_free(Numeric);                      // free
#endif

    //! Checking
    // ||x-u||
    errL2 = faspxx_dvec_norm2diff(&u, &x);
    // u = b-A*x
    faspxx_darray_cp(b.row, b.val, r.val);
    faspxx_blas_dcsr_aAxpy(-1.0, &A, x.val, r.val);
    // ||b-Ax||
    residual_norm2 = faspxx_blas_darray_norm2(r.row, r.val);
    printf("||x-u|| = %e, ||b-Ax|| = %e\n", errL2, residual_norm2);

    if (errL2 < CheckTOL && residual_norm2 < CheckTOL) printf("UMFPACK Pass!\n");
#endif

    faspxx_dcsr_free(&A);
    faspxx_dvec_free(&b);
    faspxx_dvec_free(&x);
    faspxx_dvec_free(&u);
    faspxx_dvec_free(&r);

    return SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Li Zhao          Jan/15/2024       Create file                            */
/*----------------------------------------------------------------------------*/