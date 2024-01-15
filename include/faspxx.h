/** \file    Faspxx.hxx
 *  \brief   Public FASP++ header file for users
 *  \author  Li Zhao
 *  \date    Jan/10/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FASPXX_HEADER__ /*-- allow multiple inclusions --*/
#define __FASPXX_HEADER__ ///< indicate faspxx.hxx has been included before */

#ifdef _MSC_VER /*-- Define inline attribute for MSVC --*/
#define __inline__ __inline
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if WITH_MUMPS
#include "dmumps_c.h"
#endif

#if WITH_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

#if WITH_SUPERLU
#include "slu_ddefs.h"
#endif

/**
 * \brief FASP++ version information
 */
#define FASPXX_VERSION 0.1 ///< fasp++ version
#define DEBUG_MODE 0       ///< debug mode

/*----------------------------------------------------------------------------*/
/* Definition of logic type                                                   */
/*----------------------------------------------------------------------------*/
#define TRUE 1  ///< logic TRUE
#define FALSE 0 ///< logic FALSE

/*----------------------------------------------------------------------------*/
/* Definition of switch                                                       */
/*----------------------------------------------------------------------------*/
#define ON 1  ///< turn on certain parameter
#define OFF 0 ///< turn off certain parameter

#define DLMALLOC OFF  ///< use dlmalloc instead of standard malloc
#define NEDMALLOC OFF ///< use nedmalloc instead of standard malloc

/*----------------------------------------------------------------------------*/
/* Definition of data-type length                                             */
/*----------------------------------------------------------------------------*/

typedef int          INT;      ///< Regular integer numbers
typedef unsigned int USI;      ///< Unsigned integer numbers
typedef double       DBL;      ///< Double precision numbers
typedef int          BOOL;     ///< Boolean value
typedef short        SHORT;    ///< short integer type
typedef long         LONG;     ///< long integer type
typedef long long    LONGLONG; ///< long long integer type

/*----------------------------------------------------------------------------*/
/* Definition of OpenMP and declarations                                      */
/*----------------------------------------------------------------------------*/
#ifdef _OPENMP

#include "omp.h"

#define OPENMP_HOLDS 2000 ///< Smallest size for OpenMP version

#endif /* end if for _OPENMP */

/*----------------------------------------------------------------------------*/
/* Definition of constants for range, time units, and tolerance               */
/*----------------------------------------------------------------------------*/
#define SMALL_TOL 1e-14      ///< Small positive real for tolerance
#define SMALL_TOL2 1e-40     ///< Small positive real for tolerance
#define LARGE_DBL 1e+60      ///< Largest double number
#define SMALL_DBL -1e+60     ///< Smallest double number
#define CLOSE_ZERO 1e-20     ///< Tolerance for almost zero
#define CLOCK_USE_SEC 5000   ///< Show clock time in seconds
#define CLOCK_USE_MIN 200000 ///< Show clock time in minutes
#define STRLEN 256           ///< Length of strings

/*----------------------------------------------------------------------------*/
/* Definition of constants for solvers and preconditioners                    */
/*----------------------------------------------------------------------------*/
#define MAX_ITER_NUM 100000 ///< Maximal number of multigrid levels
#define KSM_CHK_RATIO 0.95  ///< Check ratio for Krylov space methods
#define MAX_STAG_NUM 20     ///< Maximal number of stagnation checks
#define PRT_STEP_NUM 20     ///< Print iteration info every N steps
#define MAX_MG_LEVEL 20     ///< Maximal number of multigrid levels
#define MAX_RESTART 20      ///< Maximal restarting number
#define STAG_RATIO 1e-4     ///< Stagnation tolerance = tol*STAGRATIO

/// Level of output
enum Output {
    PRINT_NONE = 0, ///< no output
    PRINT_MIN  = 2, ///< minimal output
    PRINT_SOME = 4, ///< some output
    PRINT_MORE = 6, ///< more output
    PRINT_MOST = 8, ///< maximal output
    PRINT_ALL  = 10 ///< all output
};

/// Solver types avaiable.
enum SOLType {
    SOLVER_CG       = 1,  ///< Conjugate Gradient
    SOLVER_BICGSTAB = 2,  ///< Bi-Conjugate Gradient Stabilized
    SOLVER_MINRES   = 3,  ///< Minimal Residual
    SOLVER_GMRES    = 4,  ///< Generalized Minimal Residual
    SOLVER_FGMRES   = 5,  ///< Flexible GMRES
    SOLVER_VFGMRES  = 6,  ///< Variable-restarting FGMRES
    SOLVER_JACOBI   = 11, ///< Jacobi method
    SOLVER_GS       = 12, ///< Gauss-Seidel method
    SOLVER_SGS      = 13, ///< Symmetrized Gauss-Seidel method
    SOLVER_SOR      = 14, ///< Successive over-relaxation method
    SOLVER_SSOR     = 15, ///< Symmetrized successive over-relaxation method
    SOLVER_MG       = 21, ///< Multigrid method
    SOLVER_FMG      = 22, ///< Full multigrid method
    SOLVER_UMFPACK  = 91, ///< Direct method from UMFPACK
    SOLVER_MUMPS    = 92, ///< Direct method from MUMPS
    SOLVER_SUPERLU  = 93, ///< Direct method from SUPERLU
    SOLVER_PARDISO  = 94  ///< Direct method from PARDISO
};

/// Definition of preconditioner type for iterative methods
enum PCType {
    PC_NULL    = 0, ///< with no preconditioner
    PC_DIAG    = 1, ///< with diagonal preconditioner
    PC_AMG     = 2, ///< with AMG preconditioner
    PC_FMG     = 3, ///< with full AMG preconditioner
    PC_ILU     = 4, ///< with ILU preconditioner
    PC_SCHWARZ = 5  ///< with Schwarz preconditioner
};

/// Solver stop criteria.
enum SOLStop {
    STOP_REL_RES     = 1, ///< Relative residual ||r||/||b||
    STOP_REL_PRECRES = 2, ///< Relative B-residual ||r||_B/||b||_B
    STOP_MOD_REL_RES = 3, ///< Modified relative residual ||r||/||x||
};

/// Return code definition
enum FaspRetCode {
    SUCCESS               = 0,   ///< Everything is fine
    ERROR_OPEN_FILE       = -10, ///< Failed to open a file
    ERROR_INPUT_FILE      = -11, ///< Wrong input file
    ERROR_INPUT_PAR       = -12, ///< Wrong input argument
    ERROR_VEC_SIZE        = -14, ///< Wrong vector size
    ERROR_MAT_SIZE        = -15, ///< Wrong matrix size
    ERROR_NONMATCH_SIZE   = -16, ///< Two sizes do not match
    ERROR_MAT_DATA        = -17, ///< Wrong matrix format
    ERROR_DIVIDE_ZERO     = -18, ///< Divided by zero!
    ERROR_MAT_ZERODIAG    = -19, ///< MAT has zero diagonal entries
    ERROR_ALLOC_MEM       = -20, ///< Failed to allocate memory
    ERROR_UNKNOWN         = -21, ///< Unknown error
    ERROR_MISC            = -22, ///< Other error
    ERROR_DUMMY_VAR       = -23, ///< Unknown function dummy variables
    ERROR_SOLVER_TYPE     = -30, ///< Unknown solver type
    ERROR_SOLVER_PCD_TYPE = -31, ///< Unknown preconditioner type
    ERROR_SOLVER_STAG     = -32, ///< Iterative solver stagnates
    ERROR_SOLVER_SOLSTAG  = -33, ///< Iterative solver's solution is too small
    ERROR_SOLVER_TOLSMALL = -34, ///< Iterative solver's tolerance is too small
    ERROR_SOLVER_EXIT     = -38, ///< Solver does not quit successfully
    ERROR_SOLVER_MAXIT    = -39, ///< Maximal iteration number reached
    ERROR_AMG_INTERP_TYPE = -40, ///< Unknown AMG interpolation type
    ERROR_AMG_SMOOTH_TYPE = -41, ///< Unknown AMG smoother type
    ERROR_AMG_COARSE_TYPE = -42, ///< Unknown AMG coarsening type
    ERROR_AMG_COARSEING   = -43, ///< AMG coarsening step failed to complete
    ERROR_AMG_SETUP       = -49, ///< AMG setup failed to complete
    ERROR_ILU_TYPE        = -50, ///< Unknown ILU method type
    ERROR_ILU_SETUP       = -59, ///< ILU setup failed to complete
    ERROR_SWZ_TYPE        = -60, ///< Unknown Schwarz method type
    ERROR_SWZ_SETUP       = -69, ///< Schwarz method setup failed to complete
    ERROR_DSOLVER_SETUP   = -91, ///< Fail to setup in direct solver
    ERROR_DSOLVER_SOLVE   = -92, ///< Fail to solve in direct solver
    PRINT_HELP            = 0,   ///< Print help message
};

/**
 * \brief Definition of max, min, abs
 */
#define MAX(a, b) (((a) > (b)) ? (a) : (b)) ///< bigger one in a and b
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) ///< smaller one in a and b
#define ABS(a) (((a) >= 0.0) ? (a) : -(a))  ///< absolute value of a

/**
 * \brief Definition of >, >=, <, <=, and isnan
 */
#define GT(a, b) (((a) > (b)) ? (TRUE) : (FALSE))  ///< is a > b?
#define GE(a, b) (((a) >= (b)) ? (TRUE) : (FALSE)) ///< is a >= b?
#define LS(a, b) (((a) < (b)) ? (TRUE) : (FALSE))  ///< is a < b?
#define LE(a, b) (((a) <= (b)) ? (TRUE) : (FALSE)) ///< is a <= b?
#define ISNAN(a) (((a) != (a)) ? (TRUE) : (FALSE)) ///< is a == NAN?

/*! Important Note:
 *-----------------------------------------------------------------------------------
 *  This class defines the basic MAT data structure and its basic operations. The
 *  CSRx data structure is an extension of the wellknown CSR sparse matrix format.
 *  The differences lie in the following two aspects:
 *
 *  1. Unlike the classical CSR format, the CSRx format requires the column indices
 *  in each row are in a rough ascending order (namely, the nonzero entries should be
 *  ordered as lower trig, diag, uppper trig);
 *  2. The CSRx format has a diagPtr array which stores the locations of the diagonal
 *  entries in each row.
 *
 *  Note that the CSRx format stores the diagonal entries even if they are zero.
 *  Furthermore, it is compatible with CSR subroutines!
 */

/**
 * \struct dCSRxmat
 * \brief  Sparse matrix of DBL type in CSRx format
 *
 * CSR Format (IA,JA,A) in DBL
 *
 * \note The starting index of A is 0.
 */
// typedef struct dCSRxmat {

//     //! row number of matrix A, m
//     INT row;

//     //! column of matrix A, n
//     INT col;

//     //! number of nonzero entries
//     INT nnz;

//     //! integer array of row pointers, the size is m+1
//     INT* IA;

//     //! integer array of column indexes, the size is nnz
//     INT* JA;

//     //! nonzero entries of A
//     DBL* val;

//     //! pointers to diagonal entries in values, the size is m
//     INT* diagPtr;

// } dCSRxmat; ///< Sparse matrix of DBL type in CSR format */

/**
 * \struct dCSRmat
 * \brief  Sparse matrix of DBL type in CSR format
 *
 * CSR Format (IA,JA,A) in DBL
 *
 * \note The starting index of A is 0.
 */
typedef struct dCSRmat {

    //! row number of matrix A, m
    INT row;

    //! column of matrix A, n
    INT col;

    //! number of nonzero entries
    INT nnz;

    //! integer array of row pointers, the size is m+1
    INT* IA;

    //! integer array of column indexes, the size is nnz
    INT* JA;

    //! nonzero entries of A
    DBL* val;

    //! pointers to diagonal entries in values, the size is m
    INT* diagPtr;

} dCSRmat; ///< Sparse matrix of DBL type in CSR format */

/**
 * \struct dCOOmat
 * \brief  Sparse matrix of DBL type in COO (IJ) format
 *
 * Coordinate Format (I,J,A)
 *
 * \note The starting index of A is 0.
 * \note Change I to rowind, J to colind. To avoid with complex.h confliction on I.
 */
typedef struct dCOOmat {

    //! row number of matrix A, m
    INT row;

    //! column of matrix A, n
    INT col;

    //! number of nonzero entries
    INT nnz;

    //! integer array of row indices, the size is nnz
    INT* rowind;

    //! integer array of column indices, the size is nnz
    INT* colind;

    //! nonzero entries of A
    DBL* val;

} dCOOmat; ///< Sparse matrix of DBL type in COO format */

/**
 * \struct dvector
 * \brief  Vector with n entries of DBL type
 */
typedef struct dvector {

    //! number of rows
    INT row;

    //! actual vector entries
    DBL* val;

} dvector; ///< Vector of DBL type

/**
 * \struct ivector
 * \brief  Vector with n entries of INT type
 */
typedef struct ivector {

    //! number of rows
    INT row;

    //! actual vector entries
    INT* val;

} ivector; ///< Vector of INT type

/**
 * \struct precond
 * \brief  Preconditioner data and action
 *
 * \note This is the preconditioner structure for preconditioned iterative methods.
 */
typedef struct {

    //! data for preconditioner, void pointer
    void* data;

    //! action for preconditioner, void function pointer
    void (*fct)(DBL*, DBL*, void*);

} precond; ///< Vector of INT type

/**
 * \struct input_param
 * \brief  Input parameters
 *
 * Input parameters, reading from disk file
 */
typedef struct {

    // output flags
    SHORT print_level; ///< print level
    SHORT output_type; ///< type of output stream

    // problem parameters
    char inifile[STRLEN]; ///< ini file name
    char workdir[STRLEN]; ///< working directory for data files
    INT  problem_num;     ///< problem number to solve

    // parameters for iterative solvers
    SHORT solver_type;    ///< type of iterative solvers
    SHORT precond_type;   ///< type of preconditioner for iterative solvers
    SHORT stop_type;      ///< type of stopping criteria for iterative solvers
    DBL   itsolver_tol;   ///< tolerance for iterative linear solver
    INT   itsolver_maxit; ///< maximal number of iterations for iterative solvers
    INT   restart;        ///< restart number used in GMRES

    // parameters for ILU
    SHORT ILU_type;    ///< ILU type for decomposition*/
    INT   ILU_lfil;    ///< level of fill-in
    DBL   ILU_droptol; ///< drop tolerance
    DBL   ILU_relax;   ///< scaling factor: add the sum of dropped entries to diagonal
    DBL   ILU_permtol; ///< permutation tolerance

#if 0
    // parameter for Schwarz
    INT SWZ_mmsize;    ///< maximal block size
    INT SWZ_maxlvl;    ///< maximal levels
    INT SWZ_type;      ///< type of Schwarz method
    INT SWZ_blksolver; ///< type of Schwarz block solver

    // parameters for AMG
    SHORT AMG_type;              ///< Type of AMG
    SHORT AMG_levels;            ///< maximal number of levels
    SHORT AMG_cycle_type;        ///< type of cycle
    SHORT AMG_smoother;          ///< type of smoother
    SHORT AMG_smooth_order;      ///< order for smoothers
    DBL   AMG_relaxation;        ///< over-relaxation parameter for SOR
    SHORT AMG_polynomial_degree; ///< degree of the polynomial smoother
    SHORT AMG_presmooth_iter;    ///< number of presmoothing
    SHORT AMG_postsmooth_iter;   ///< number of postsmoothing
    DBL   AMG_tol;               ///< tolerance for AMG if used as preconditioner
    INT   AMG_coarse_dof;        ///< max number of coarsest level DOF
    INT   AMG_maxit;          ///< number of iterations for AMG used as preconditioner
    SHORT AMG_ILU_levels;     ///< how many levels use ILU smoother
    SHORT AMG_coarse_solver;  ///< coarse solver type
    SHORT AMG_coarse_scaling; ///< switch of scaling of the coarse grid correction
    SHORT AMG_amli_degree;    ///< degree of the polynomial used by AMLI cycle
    SHORT
    AMG_nl_amli_krylov_type; ///< type of Krylov method used by nonlinear AMLI cycle
    INT AMG_SWZ_levels;      ///< number of levels use Schwarz smoother

    // parameters for classical AMG
    SHORT AMG_coarsening_type;      ///< coarsening type
    SHORT AMG_aggregation_type;     ///< aggregation type
    SHORT AMG_interpolation_type;   ///< interpolation type
    DBL   AMG_strong_threshold;     ///< strong threshold for coarsening
    DBL   AMG_truncation_threshold; ///< truncation factor for interpolation
    DBL   AMG_max_row_sum;          ///< maximal row sum
    INT   AMG_aggressive_level;     ///< number of levels use aggressive coarsening
    INT   AMG_aggressive_path; ///< number of paths to determine strongly coupled C-set
    INT   AMG_pair_number;     ///< number of pairs in matching algorithm
    DBL   AMG_quality_bound;   ///< threshold for pair wise aggregation

    //  parameters for smoothed aggregation AMG
    DBL AMG_strong_coupled;   ///< strong coupled threshold for aggregate
    INT AMG_max_aggregation;  ///< max size of each aggregate
    DBL AMG_tentative_smooth; ///< relaxation factor for smoothing the tentative
                              ///< prolongation
    SHORT AMG_smooth_filter; ///< use filter for smoothing the tentative prolongation or
                             ///< not
    SHORT AMG_smooth_restriction; ///< smoothing the restriction or not
#endif

} input_param; ///< Input parameters

/*---------------------------*/
/*--- Parameter structures --*/
/*---------------------------*/

/**
 * \struct ITS_param
 * \brief  Parameters for iterative solvers
 */
typedef struct {

    SHORT print_level;   ///< print level: 0--10
    SHORT itsolver_type; ///< solver type
    SHORT precond_type;  ///< preconditioner type
    SHORT stop_type;     ///< stopping criteria type
    INT   restart;       ///< number of steps for restarting: for GMRES etc
    INT   maxit;         ///< max number of iterations
    DBL   tol;           ///< convergence tolerance

} ITS_param; ///< Parameters for iterative solvers

/*-----------------------------------*/
/*--- structures of direct solver ---*/
/*-----------------------------------*/

/**
 * \struct Mumps_data
 * \brief  Data for MUMPS interface
 */
typedef struct {

#if WITH_MUMPS
    //! solver instance for MUMPS
    DMUMPS_STRUC_C id;
#endif

    //! work for MUMPS
    INT job;

} Mumps_data; /**< Data for MUMPS */

/**
 * \struct Pardiso_data
 * \brief  Data for Intel MKL PARDISO interface
 *
 */
typedef struct {

    //! Internal solver memory pointer
    void* pt[64];

#if WITH_PARDISO
    //! Pardiso control parameters
    MKL_INT iparm[64];

    //! Type of the matrix
    MKL_INT mtype;

    //! Maximum number of numerical factorizations
    MKL_INT maxfct;

    //! Indicate the actual matrix for the solution phase, 1 <= mnum <= maxfct
    MKL_INT mnum;

#endif

} Pardiso_data; /**< Data for PARDISO */

/**
 * \struct SuperLU_data
 * \brief  Data for SuperLU interface
 */
typedef struct {

    //! Row permutations from partial pivoting
    INT* perm_r;

    //! Column permutation vector
    INT* perm_c;

    //! Return variable
    INT info;

#if WITH_SUPERLU
    //! Matrices A, L, U, right-hand side B
    SuperMatrix A, L, U, B;

    //! The statistics variables, e.g., running time, flops, tiny pivots, iterative
    //! refinement steps, memory
    SuperLUStat_t stat;

    //! Set the default input options
    superlu_options_t options;

    //! Convert A to SLU_NC format when necessary
    trans_t trans;
#endif

} SuperLU_data; /**< Data for SuperLU */

#endif /* end if for __FASPXX_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Li Zhao           Sep/01/2019        Create file                          */
/*----------------------------------------------------------------------------*/