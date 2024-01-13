/*! \file  AuxVector.c
 *
 *  \brief Simple vector operations -- init, set, copy, etc
 *
 *  \note  This file contains Level-0 (Aux) functions. It requires:
 *         AuxThreads.c
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
 * \fn dvector faspxx_dvec_create (const INT m)
 *
 * \brief Create dvector data space of DBL type
 *
 * \param m    Number of rows
 *
 * \return u   The new dvector
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
dvector faspxx_dvec_create(const INT m)
{
    dvector u;

    u.row = m;
    u.val = (DBL*)faspxx_mem_calloc(m, sizeof(DBL));

    return u;
}

/**
 * \fn ivector faspxx_ivec_create (const INT m)
 *
 * \brief Create vector data space of INT type
 *
 * \param m   Number of rows
 *
 * \return u  The new ivector
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
ivector faspxx_ivec_create(const INT m)
{
    ivector u;

    u.row = m;
    u.val = (INT*)faspxx_mem_calloc(m, sizeof(INT));

    return u;
}

/**
 * \fn void faspxx_dvec_alloc (const INT m, dvector *u)
 *
 * \brief Create dvector data space of DBL type
 *
 * \param m    Number of rows
 * \param u    Pointer to dvector (OUTPUT)
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_alloc(const INT m, dvector* u)
{
    u->row = m;
    u->val = (DBL*)faspxx_mem_calloc(m, sizeof(DBL));

    return;
}

/**
 * \fn void faspxx_ivec_alloc (const INT m, ivector *u)
 *
 * \brief Create vector data space of INT type
 *
 * \param m   Number of rows
 * \param u   Pointer to ivector (OUTPUT)
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_ivec_alloc(const INT m, ivector* u)
{

    u->row = m;
    u->val = (INT*)faspxx_mem_calloc(m, sizeof(INT));

    return;
}

/**
 * \fn void faspxx_dvec_free (dvector *u)
 *
 * \brief Free vector data space of DBL type
 *
 * \param u   Pointer to dvector which needs to be deallocated
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_free(dvector* u)
{
    if (u == NULL) return;

    faspxx_mem_free(u->val);
    u->val = NULL;
    u->row = 0;
}

/**
 * \fn void faspxx_ivec_free (ivector *u)
 *
 * \brief Free vector data space of INT type
 *
 * \param u   Pointer to ivector which needs to be deallocated
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \note This function is same as faspxx_dvec_free except input type.
 */
void faspxx_ivec_free(ivector* u)
{
    if (u == NULL) return;

    faspxx_mem_free(u->val);
    u->val = NULL;
    u->row = 0;
}

/**
 * \fn void faspxx_dvec_rand (const INT n, dvector *x)
 *
 * \brief Generate random DBL vector in the range from 0 to 1
 *
 * \param n    Size of the vector
 * \param x    Pointer to dvector
 *
 * \note Sample usage:
 * \par
 *   dvector xapp;
 * \par
 *   faspxx_dvec_create(100,&xapp);
 * \par
 *   faspxx_dvec_rand(100,&xapp);
 * \par
 *   faspxx_dvec_print(100,&xapp);
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_rand(const INT n, dvector* x)
{
    const INT va = (DBL)0;
    const INT vb = (DBL)n;

    INT s = 1, i, j;

    srand(s);
    for (i = 0; i < n; ++i) {
        j         = 1 + (INT)(((DBL)n) * rand() / (RAND_MAX + 1.0));
        x->val[i] = (((DBL)j) - va) / (vb - va);
    }
    x->row = n;
}

/**
 * \fn void faspxx_dvec_set (INT n, dvector *x, const DBL val)
 *
 * \brief Initialize dvector x[i]=val for i=0:n-1
 *
 * \param n      Number of variables
 * \param x      Pointer to dvector
 * \param val    Initial value for the vector
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_set(INT n, dvector* x, const DBL val)
{
    INT  i;
    DBL* xpt = x->val;

    if (n > 0)
        x->row = n;
    else
        n = x->row;

#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend;
    INT nthreads = faspxx_get_num_threads();
#endif

    if (val == 0.0) {

#ifdef _OPENMP
        if (n > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend)
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
                memset(&xpt[mybegin], 0x0, sizeof(DBL) * (myend - mybegin));
            }
        } else {
#endif
            memset(xpt, 0x0, sizeof(DBL) * n);
#ifdef _OPENMP
        }
#endif

    }

    else {

#ifdef _OPENMP
        if (n > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend)
            for (myid = 0; myid < nthreads; myid++) {
                faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) xpt[i] = val;
            }
        } else {
#endif
            for (i = 0; i < n; ++i) xpt[i] = val;
#ifdef _OPENMP
        }
#endif
    }
}

/**
 * \fn void faspxx_ivec_set (INT n, ivector *u, const INT m)
 *
 * \brief Set ivector value to be m
 *
 * \param n    Number of variables
 * \param m    Integer value of ivector
 * \param u    Pointer to ivector (MODIFIED)
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_ivec_set(INT n, ivector* u, const INT m)
{
    SHORT nthreads = 1, use_openmp = FALSE;
    INT   i;

    if (n > 0)
        u->row = n;
    else
        n = u->row;

#ifdef _OPENMP
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
        INT mybegin, myend, myid;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) u->val[i] = m;
        }
    } else {
        for (i = 0; i < n; ++i) u->val[i] = m;
    }
}

/**
 * \fn void faspxx_dvec_cp (const dvector *x, dvector *y)
 *
 * \brief Copy dvector x to dvector y
 *
 * \param x  Pointer to dvector
 * \param y  Pointer to dvector (MODIFIED)
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_cp(const dvector* x, dvector* y)
{
    y->row = x->row;
    memcpy(y->val, x->val, x->row * sizeof(DBL));
}

/**
 * \fn DBL faspxx_dvec_maxdiff (const dvector *x, const dvector *y)
 *
 * \brief Maximal difference of two dvector x and y
 *
 * \param  x    Pointer to dvector
 * \param  y    Pointer to dvector
 *
 * \return      Maximal norm of x-y
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
DBL faspxx_dvec_maxdiff(const dvector* x, const dvector* y)
{
    const INT  length = x->row;
    const DBL *xpt = x->val, *ypt = y->val;

    SHORT use_openmp = FALSE;
    INT   i;
    DBL   Linf = 0.0, diffi = 0.0;

#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if (length > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP
        DBL temp = 0.;
#pragma omp parallel firstprivate(temp) private(myid, mybegin, myend, i, diffi)
        {
            myid = omp_get_thread_num();
            faspxx_get_start_end(myid, nthreads, length, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                if ((diffi = ABS(xpt[i] - ypt[i])) > temp) temp = diffi;
            }
#pragma omp critical
            if (temp > Linf) Linf = temp;
        }
#endif
    } else {
        for (i = 0; i < length; ++i) {
            if ((diffi = ABS(xpt[i] - ypt[i])) > Linf) Linf = diffi;
        }
    }

    return Linf;
}

/**
 * \fn void faspxx_dvec_symdiagscale (dvector *b, const dvector *diag)
 *
 * \brief Symmetric diagonal scaling D^{-1/2}b
 *
 * \param b       Pointer to dvector
 * \param diag    Pointer to dvector: the diagonal entries
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_symdiagscale(dvector* b, const dvector* diag)
{
    // information about dvector
    const INT n   = b->row;
    DBL*      val = b->val;

    // local variables
    SHORT use_openmp = FALSE;
    INT   i;

    if (diag->row != n) {
        printf("### ERROR: Sizes of diag = %d != dvector = %d!", diag->row, n);
        faspxx_chkerr(ERROR_MISC, __FUNCTION__);
    }

#ifdef _OPENMP
    INT mybegin, myend, myid, nthreads;
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend)
        for (myid = 0; myid < nthreads; myid++) {
            faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) val[i] = val[i] / sqrt(diag->val[i]);
        }
#endif
    } else {
        for (i = 0; i < n; ++i) val[i] = val[i] / sqrt(diag->val[i]);
    }

    return;
}

/**
 * \fn SHORT faspxx_dvec_isnan (const dvector *u)
 *
 * \brief Check a dvector whether there is NAN
 *
 * \param u    Pointer to dvector
 *
 * \return     Return TRUE if there is NAN
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
BOOL faspxx_dvec_isnan(const dvector* u)
{
    INT i;

    for (i = 0; i < u->row; i++) {
        if (ISNAN(u->val[i])) return TRUE;
    }

    return FALSE;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
