/*! \file  BlaArray.c
 *
 *  \brief BLAS1 operations for arrays
 *
 *  \note  This file contains Level-1 (Bla) functions. It requires:
 *         AuxThreads.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "faspxx.h"
#include "faspxx_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void faspxx_blas_darray_ax (const INT n, const DBL a, DBL *x)
 *
 * \brief x := a*x
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 *
 * \author Li Zhao
 * \date   01/13/2024
 *
 * \warning x is reused to store the resulting array!
 */
void faspxx_blas_darray_ax(const INT n, const DBL a, DBL* x)
{
    if (a == 1.0) return; // do nothing

    {
        SHORT use_openmp = FALSE;
        INT   i;

#ifdef _OPENMP
        INT myid, mybegin, myend, nthreads;
        if (n > OPENMP_HOLDS) {
            use_openmp = TRUE;
            nthreads   = faspxx_get_num_threads();
        }
#endif

        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i)
            {
                myid = omp_get_thread_num();
                faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) x[i] *= a;
            }
#endif
        } else {
            for (i = 0; i < n; ++i) x[i] *= a;
        }
    }
}

/**
 * \fn void faspxx_blas_darray_axpy (const INT n, const DBL a,
 *                                 const DBL *x, DBL *y)
 *
 * \brief y = a*x + y
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 * \param y    Pointer to y, reused to store the resulting array
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_darray_axpy(const INT n, const DBL a, const DBL* x, DBL* y)
{
    SHORT use_openmp = FALSE;
    INT   i;

#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (a == 1.0) {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) y[i] += x[i];
            }
#endif
        } else {
            for (i = 0; i < n; ++i) y[i] += x[i];
        }
    }

    else if (a == -1.0) {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) y[i] -= x[i];
            }
#endif
        } else {
            for (i = 0; i < n; ++i) y[i] -= x[i];
        }
    }

    else {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) y[i] += a * x[i];
            }
#endif
        } else {
            for (i = 0; i < n; ++i) y[i] += a * x[i];
        }
    }
}

/**
 * \fn void faspxx_blas_darray_axpyz (const INT n, const DBL a, const DBL *x,
 *                                  const DBL *y, DBL *z)
 *
 * \brief z = a*x + y
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 * \param y    Pointer to y
 * \param z    Pointer to z
 *
 * \author Chensong Zhang
 * \date   07/01/2009
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
void faspxx_blas_darray_axpyz(const INT n, const DBL a, const DBL* x, const DBL* y,
                              DBL* z)
{
    SHORT use_openmp = FALSE;
    INT   i;

#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) z[i] = a * x[i] + y[i];
        }
#endif
    } else {
        for (i = 0; i < n; ++i) z[i] = a * x[i] + y[i];
    }
}

/**
 * \fn void faspxx_blas_darray_axpby (const INT n, const DBL a, const DBL *x,
 *                                  const DBL b, DBL *y)
 *
 * \brief y = a*x + b*y
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 * \param b    Factor b
 * \param y    Pointer to y, reused to store the resulting array
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_blas_darray_axpby(const INT n, const DBL a, const DBL* x, const DBL b,
                              DBL* y)
{
    SHORT use_openmp = FALSE;
    INT   i;

#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) y[i] = a * x[i] + b * y[i];
        }
#endif
    } else {
        for (i = 0; i < n; ++i) y[i] = a * x[i] + b * y[i];
    }
}

/**
 * \fn DBL faspxx_blas_darray_norm1 (const INT n, const DBL *x)
 *
 * \brief L1 norm of array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     L1 norm of x
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
DBL faspxx_blas_darray_norm1(const INT n, const DBL* x)
{
    register DBL onenorm = 0.0;
    INT          i;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : onenorm) private(i)
#endif
    for (i = 0; i < n; ++i) onenorm += ABS(x[i]);

    return onenorm;
}

/**
 * \fn DBL faspxx_blas_darray_norm2 (const INT n, const DBL *x)
 *
 * \brief L2 norm of array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     L2 norm of x
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
DBL faspxx_blas_darray_norm2(const INT n, const DBL* x)
{
    register DBL twonorm = 0.0;
    INT          i;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : twonorm) private(i)
#endif
    for (i = 0; i < n; ++i) twonorm += x[i] * x[i];

    return sqrt(twonorm);
}

/**
 * \fn DBL faspxx_blas_darray_norminf (const INT n, const DBL *x)
 *
 * \brief Linf norm of array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     L_inf norm of x
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
DBL faspxx_blas_darray_norminf(const INT n, const DBL* x)
{
    SHORT        use_openmp = FALSE;
    register DBL infnorm    = 0.0;
    INT          i;

#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if (n > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads   = faspxx_get_num_threads();
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP
        DBL infnorm_loc = 0.0;
#pragma omp parallel firstprivate(infnorm_loc) private(myid, mybegin, myend, i)
        {
            myid = omp_get_thread_num();
            faspxx_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) infnorm_loc = MAX(infnorm_loc, ABS(x[i]));

            if (infnorm_loc > infnorm) {
#pragma omp critical
                infnorm = MAX(infnorm_loc, infnorm);
            }
        }
#endif
    } else {
        for (i = 0; i < n; ++i) infnorm = MAX(infnorm, ABS(x[i]));
    }

    return infnorm;
}

/**
 * \fn DBL faspxx_blas_darray_dotprod (const INT n, const DBL *x, const DBL *y)
 *
 * \brief Inner product of two arraies x and y
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 * \param y    Pointer to y
 *
 * \return     Inner product (x,y)
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
DBL faspxx_blas_darray_dotprod(const INT n, const DBL* x, const DBL* y)
{
    SHORT        use_openmp = FALSE;
    register DBL value      = 0.0;
    INT          i;

#ifdef _OPENMP
    if (n > OPENMP_HOLDS) use_openmp = TRUE;
#endif

    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : value) private(i)
#endif
        for (i = 0; i < n; ++i) value += x[i] * y[i];
    } else {
        for (i = 0; i < n; ++i) value += x[i] * y[i];
    }

    return value;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
