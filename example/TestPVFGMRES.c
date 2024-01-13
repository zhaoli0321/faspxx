#include "faspxx.h"
#include "faspxx_functs.h"

int main()
{
    dCSRmat A;
    dvector b, x, u;
    DBL     tol      = 1e-8;
    INT     MaxIt    = 2000;
    SHORT   StopType = STOP_REL_RES;
    SHORT   PrtLvl   = PRINT_MORE;
    INT     Iter;
    SHORT   restart = 30;

    printf("PVFGMRES Test...\n");

    const char csrfilename[512] = "../data/fem_small.dat";
    faspxx_dcsr_read(csrfilename, &A);

    u = faspxx_dvec_create(A.row); // exact solution
    faspxx_dvec_set(u.row, &u, 1.0);
    printf("A.row: %d, A.col: %d, A.nnz: %d\n", A.row, A.col, A.nnz);

    b = faspxx_dvec_create(A.row);
    faspxx_blas_dcsr_mxv(&A, u.val, b.val);

    x = faspxx_dvec_create(A.row);
    faspxx_dvec_set(x.row, &x, 0.0);

    Iter = faspxx_solver_dcsr_pvfgmres(&A, &b, &x, NULL, tol, MaxIt, restart, StopType,
                                       PrtLvl);

    DBL errL2 = faspxx_dvec_norm2diff(&u, &x);

    // u = b-A*x
    faspxx_darray_cp(b.row, b.val, u.val);
    faspxx_blas_dcsr_aAxpy(-1.0, &A, x.val, u.val);
    DBL residual_norm2 = faspxx_blas_darray_norm2(u.row, u.val);
    printf("||x-u|| = %e, Iter = %d, ||b-Ax|| = %e\n", errL2, Iter, residual_norm2);

    faspxx_dcsr_free(&A);
    faspxx_dvec_free(&b);
    faspxx_dvec_free(&x);
    faspxx_dvec_free(&u);

    printf("Pass!\n");
    return SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Li Zhao          Jan/13/2024       Create file                            */
/*----------------------------------------------------------------------------*/