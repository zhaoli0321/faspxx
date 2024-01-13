#include "faspxx.h"
#include "faspxx_functs.h"

int main()
{
    dCSRmat A1, A2;
    dvector b;

    printf("Read Matrices and Vector Test...\n");

    const char csrfilename[512] = "../data/fdm_10X10.dat";
    faspxx_dcsr_read(csrfilename, &A1);

    const char coofilename[512] = "../data/nos7.mtx";
    faspxx_dmtx_read(coofilename, &A2);

    const char rhsfilename[512] = "../data/rhs_10X10.dat";
    faspxx_dvec_read(rhsfilename, &b);

    faspxx_dcsr_free(&A1);
    faspxx_dcsr_free(&A2);
    faspxx_dvec_free(&b);

    printf("Pass!\n");
    return SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Li Zhao          Jan/01/2024       Create file                            */
/*----------------------------------------------------------------------------*/