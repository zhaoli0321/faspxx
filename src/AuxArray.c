/*! \file  AuxArray.c
 *
 *  \brief Simple array operations -- init, set, copy, etc
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
 * \fn void faspxx_darray_set (const INT n, DBL *x, const DBL val)
 *
 * \brief Set initial value for an double array to be x=val
 *
 * \param n    Number of variables
 * \param x    Pointer to the vector
 * \param val  Initial value for the DBL array
 *
 * \author Li Zhao
 * \date   01/10/2024
 *
 */
void faspxx_darray_set(const INT n, DBL* x, const DBL val)
{
    if (val == 0.0) {
        memset(x, 0x0, sizeof(DBL) * n);
    } else {
        INT i;
        for (i = 0; i < n; ++i) x[i] = val;
    }
}

/**
 * \fn void faspxx_iarray_set (const INT n, INT *x, const INT val)
 *
 * \brief Set initial value for an interge array to be x=val
 *
 * \param n    Number of variables
 * \param x    Pointer to the vector
 * \param val  Initial value for the DBL array
 *
 * \author Li Zhao
 * \date   01/10/2024
 *
 */
void faspxx_iarray_set(const INT n, INT* x, const INT val)
{
    if (val == 0) {
        memset(x, 0, sizeof(INT) * n);
    } else {
        INT i;
        for (i = 0; i < n; ++i) x[i] = val;
    }
}

/**
 * \fn void faspxx_darray_cp (const INT n, const DBL *x, DBL *y)
 *
 * \brief Copy a double array to the other y=x
 *
 * \param n    Number of variables
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 *
 * \author Li Zhao
 * \date   01/10/2024
 */
void faspxx_darray_cp(const INT n, const DBL* x, DBL* y)
{
    memcpy(y, x, n * sizeof(DBL));
}

/**
 * \fn void faspxx_iarray_cp (const INT n, const INT *x, INT *y)
 *
 * \brief Copy an interge array to the other y=x
 *
 * \param n    Number of variables
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 *
 * \author Li Zhao
 * \date   01/10/2024
 */
void faspxx_iarray_cp(const INT n, const INT* x, INT* y)
{
    memcpy(y, x, n * sizeof(INT));
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
