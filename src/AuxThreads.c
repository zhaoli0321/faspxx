/*! \file  AuxThreads.c
 *
 *  \brief Get and set number of threads and assign work load for each thread
 *
 *  \note  This file contains Level-0 (Aux) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

#ifdef _OPENMP

INT thread_ini_flag = 0;

/**
 * \fn     INT faspxx_get_num_threads ( void )
 *
 * \brief  Get the number of threads for thread related functions.
 *
 * \return The number of threads to run
 *
 * \author Li Zhao
 * \date   01/10/2024
 */
INT faspxx_get_num_threads(void)
{
    static INT nthreads;

    if (thread_ini_flag == 0) {
        nthreads = 1;
#pragma omp parallel
        nthreads = omp_get_num_threads();

        printf("\nFASPXX is running on %d thread(s).\n\n", nthreads);
        thread_ini_flag = 1;
    }

    return nthreads;
}

/**
 * \fn     INT faspxx_set_num_threads (const INT nthreads)
 *
 * \brief  Set the number of threads for thread related functions.
 *
 * \param  nthreads  Desirable number of threads
 *
 * \return The number of threads to run
 *
 * \author Li Zhao
 * \date   01/10/2024
 */
INT faspxx_set_num_threads(const INT nthreads)
{
    omp_set_num_threads(nthreads);

    return nthreads;
}

#endif

/**
 * \fn    void faspxx_get_start_end (const INT procid, const INT nprocs, const INT n,
 *                                 INT *start, INT *end)
 *
 * \brief Assign Load to each thread.
 *
 * \param procid Index of thread
 * \param nprocs Number of threads
 * \param n      Total workload
 * \param start  Pointer to the begin of each thread in total workload
 * \param end    Pointer to the end of each thread in total workload
 *
 * \author Li Zhao
 * \date   01/10/2024
 */
void faspxx_get_start_end(const INT procid, const INT nprocs, const INT n, INT* start,
                          INT* end)
{
    INT chunk_size = n / nprocs;
    INT mod        = n % nprocs;
    INT start_loc, end_loc;

    if (procid < mod) {
        end_loc   = chunk_size + 1;
        start_loc = end_loc * procid;
    } else {
        end_loc   = chunk_size;
        start_loc = end_loc * procid + mod;
    }
    end_loc = end_loc + start_loc;

    *start = start_loc;
    *end   = end_loc;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
