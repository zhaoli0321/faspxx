/*! \file  AuxTiming.c
 *
 *  \brief Timing subroutines
 *
 *  \note  This file contains Level-0 (Aux) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include <sys/time.h>

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void faspxx_gettime (REAL *time)
 *
 * \brief Get system time
 *
 * \author Li Zhao
 * \date   01/15/2024
 */

#if 0
void faspxx_gettime(DBL* time)
{
    if (time != NULL) {
#ifdef _OPENMP
        *time = omp_get_wtime();
#else
        *time = (DBL)clock() / CLOCKS_PER_SEC;
#endif
    }
}
#else

void faspxx_gettime(DBL* time)
{
    struct timeval tp;

    gettimeofday(&tp, NULL);
    *time = tp.tv_sec + tp.tv_usec / 1000000.0;
}
#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
