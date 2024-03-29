/*! \file  AuxMessage.c
 *
 *  \brief Output some useful messages
 *
 *  \note  This file contains Level-0 (Aux) functions.
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
 * \fn void faspxx_itinfo (const INT ptrlvl, const INT stop_type, const INT iter,
 *                       const DBL relres, const DBL absres, const DBL factor)
 *
 * \brief Print out iteration information for iterative solvers
 *
 * \param ptrlvl     Level for output
 * \param stop_type  Type of stopping criteria
 * \param iter       Number of iterations
 * \param relres     Relative residual of different kinds
 * \param absres     Absolute residual of different kinds
 * \param factor     Contraction factor
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 *
 * Modified by Chensong Zhang on 03/28/2013: Output initial guess
 * Modified by Chensong Zhang on 04/05/2013: Fix a typo
 */
void faspxx_itinfo(const INT ptrlvl, const INT stop_type, const INT iter,
                   const DBL relres, const DBL absres, const DBL factor)
{
    if (ptrlvl >= PRINT_SOME) {

        if (iter > 0) {
            printf("%6d | %13.6e   | %13.6e  | %10.4f\n", iter, relres, absres, factor);
        } else { // iter = 0: initial guess
            printf("-----------------------------------------------------------\n");
            switch (stop_type) {
                case STOP_REL_RES:
                    printf(
                        "It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
                    break;
                case STOP_REL_PRECRES:
                    printf(
                        "It Num | ||r||_B/||b||_B |    ||r||_B     |  Conv. Factor\n");
                    break;
                case STOP_MOD_REL_RES:
                    printf(
                        "It Num |   ||r||/||x||   |     ||r||      |  Conv. Factor\n");
                    break;
            }
            printf("-----------------------------------------------------------\n");
            printf("%6d | %13.6e   | %13.6e  |     -.-- \n", iter, relres, absres);
        } // end if iter

    } // end if ptrlvl
}

/**
 * \fn void void faspxx_cputime (const char *message, const DBL cputime)
 *
 * \brief Print CPU walltime
 *
 * \param message   Some string to print out
 * \param cputime   Walltime since start to end
 *
 * \author Chensong Zhang
 * \date   04/10/2012
 */
void faspxx_cputime(const char* message, const DBL cputime)
{
    printf("%s costs %.4f seconds.\n", message, cputime);
}

/**
 * \fn void faspxx_message (const INT ptrlvl, const char *message)
 *
 * \brief Print output information if necessary
 *
 * \param ptrlvl   Level for output
 * \param message  Error message to print
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
void faspxx_message(const INT ptrlvl, const char* message)
{
    if (ptrlvl > PRINT_NONE) printf("%s", message);
}

/**
 * \fn void faspxx_chkerr (const SHORT status, const char *fctname)
 *
 * \brief Check error status and print out error messages before quit
 *
 * \param status   Error status
 * \param fctname  Function name where this routine is called
 *
 * \author Chensong Zhang
 * \date   01/10/2012
 */
void faspxx_chkerr(const SHORT status, const char* fctname)
{
    if (status >= 0) return; // No error found!!!

    switch (status) {
        // case ERROR_READ_FILE:
        //     printf("### ERROR: Cannot read file! [%s]\n", fctname);
        //     break;
        case ERROR_OPEN_FILE:
            printf("### ERROR: Cannot open file! [%s]\n", fctname);
            break;
        case ERROR_INPUT_FILE:
            printf("### ERROR: Unknown file format! [%s]\n", fctname);
            break;
        case ERROR_INPUT_PAR:
            printf("### ERROR: Unknown input argument! [%s]\n", fctname);
            break;
        // case ERROR_REGRESS:
        //     printf("### ERROR: Regression test failed! [%s]\n", fctname);
        //     break;
        case ERROR_ALLOC_MEM:
            printf("### ERROR: Cannot allocate memory! [%s]\n", fctname);
            break;
        // case ERROR_NUM_BLOCKS:
        //     printf("### ERROR: Unexpected number of blocks! [%s]\n", fctname);
        //     break;
        // case ERROR_DATA_STRUCTURE:
        //     printf("### ERROR: Wrong data structure! [%s]\n", fctname);
        //     break;
        // case ERROR_DATA_ZERODIAG:
        //     printf("### ERROR: Matrix has zero diagonal entries! [%s]\n", fctname);
        //     break;
        case ERROR_DUMMY_VAR:
            printf("### ERROR: Unknown input argument! [%s]\n", fctname);
            break;
        case ERROR_AMG_INTERP_TYPE:
            printf("### ERROR: Unknown AMG interpolation type! [%s]\n", fctname);
            break;
        case ERROR_AMG_COARSE_TYPE:
            printf("### ERROR: Unknown AMG coarsening type! [%s]\n", fctname);
            break;
        case ERROR_AMG_SMOOTH_TYPE:
            printf("### ERROR: Unknown AMG smoother type! [%s]\n", fctname);
            break;
        case ERROR_SOLVER_TYPE:
            printf("### ERROR: Unknown solver type! [%s]\n", fctname);
            break;
        // case ERROR_SOLVER_PRECTYPE:
        //     printf("### ERROR: Unknown preconditioner type! [%s]\n", fctname);
        //     break;
        case ERROR_SOLVER_STAG:
            printf("### ERROR: Solver stagnation! [%s]\n", fctname);
            break;
        case ERROR_SOLVER_SOLSTAG:
            printf("### ERROR: Solution close to zero! [%s]\n", fctname);
            break;
        case ERROR_SOLVER_TOLSMALL:
            printf("### ERROR: Convergence tolerance too small! [%s]\n", fctname);
            break;
        // case ERROR_SOLVER_ILUSETUP:
        //     printf("### ERROR: ILU setup failed! [%s]\n", fctname);
        //     break;
        case ERROR_SOLVER_MAXIT:
            printf("### ERROR: Max iteration number reached! [%s]\n", fctname);
            break;
        // case ERROR_SOLVER_EXIT:
        //     printf("### ERROR: Iterative solver failed! [%s]\n", fctname);
        //     break;
        // case ERROR_SOLVER_MISC:
        //     printf("### ERROR: Unknown solver runtime error! [%s]\n", fctname);
        //     break;
        case ERROR_MISC:
            printf("### ERROR: Miscellaneous error! [%s]\n", fctname);
            break;
        // case ERROR_QUAD_TYPE:
        //     printf("### ERROR: Unknown quadrature rules! [%s]\n", fctname);
        //     break;
        // case ERROR_QUAD_DIM:
        //     printf("### ERROR: Num of quad points not supported! [%s]\n", fctname);
        //     break;
        case ERROR_UNKNOWN:
            printf("### ERROR: Unknown error! [%s]\n", fctname);
            break;
        default:
            printf("### ERROR: Unknown error! [%s]\n", fctname);
            break;
    }

    exit(status);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
