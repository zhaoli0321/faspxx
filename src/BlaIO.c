/*! \file  BlaIO.c
 *
 *  \brief Matrix/vector input/output subroutines
 *
 *  \note  Read, write or print a matrix or a vector in various formats
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "faspxx.h"
#include "faspxx_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/
void skip_comments(FILE* fp);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void faspxx_dcsr_read (const char *filename, dCSRmat *A)
 *
 * \brief Read A from matrix disk file in CSR format -- indices starting from 0 or 1
 *
 * \param filename  Char for matrix file name
 * \param A         Pointer to the CSR matrix
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_read(const char* filename, dCSRmat* A)
{
    int i, m, idata;
    DBL ddata;

    // Open input disk file
    FILE* fp = fopen(filename, "r");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: reading file %s...\n", __FUNCTION__, filename);

    skip_comments(fp); // skip the comments in the beginning --zcs 06/30/2020

    // Read CSR matrix
    if (fscanf(fp, "%d", &m) > 0)
        A->row = A->col = m;
    else {
        faspxx_chkerr(ERROR_INPUT_FILE, filename);
    }

    A->IA = (INT*)faspxx_mem_calloc(m + 1, sizeof(INT));
    for (i = 0; i <= m; ++i) {
        if (fscanf(fp, "%d", &idata) > 0)
            A->IA[i] = idata;
        else {
            faspxx_chkerr(ERROR_INPUT_FILE, filename);
        }
    }

    // If IA starts from 1, shift by -1
    if (A->IA[0] == 1)
        for (i = 0; i <= m; ++i) A->IA[i]--;

    INT nnz = A->IA[m] - A->IA[0];

    A->nnz = nnz;
    A->JA  = (INT*)faspxx_mem_calloc(nnz, sizeof(INT));
    A->val = (DBL*)faspxx_mem_calloc(nnz, sizeof(DBL));

    for (i = 0; i < nnz; ++i) {
        if (fscanf(fp, "%d", &idata) > 0)
            A->JA[i] = idata;
        else {
            faspxx_chkerr(ERROR_INPUT_FILE, filename);
        }
    }

    // If JA starts from 1, shift by -1
    if (A->JA[0] == 1)
        for (i = 0; i < nnz; ++i) A->JA[i]--;

    for (i = 0; i < nnz; ++i) {
        if (fscanf(fp, "%lf", &ddata) > 0)
            A->val[i] = ddata;
        else {
            faspxx_chkerr(ERROR_INPUT_FILE, filename);
        }
    }

    fclose(fp);
}

/**
 * \fn void faspxx_dcoo_read (const char *filename, dCSRmat *A)
 *
 * \brief Read A from matrix disk file in IJ format -- indices starting from 0
 *
 * \param filename  File name for matrix
 * \param A         Pointer to the CSR matrix
 *
 * \note File format:
 *   - nrow ncol nnz     % number of rows, number of columns, and nnz
 *   - i  j  a_ij        % i, j a_ij in each line
 *
 * \note After reading, it converts the matrix to dCSRmat format.
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcoo_read(const char* filename, dCSRmat* A)
{
    int i, j, k, m, n, nnz;
    DBL value;

    FILE* fp = fopen(filename, "r");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: reading file %s...\n", __FUNCTION__, filename);

    skip_comments(fp); // skip the comments in the beginning

    if (fscanf(fp, "%d %d %d", &m, &n, &nnz) <= 0) {
        faspxx_chkerr(ERROR_INPUT_FILE, filename);
    }

    dCOOmat Atmp = faspxx_dcoo_create(m, n, nnz);

    for (k = 0; k < nnz; k++) {
        if (fscanf(fp, "%d %d %le", &i, &j, &value) != EOF) {
            Atmp.rowind[k] = i; //! indices starting from 0
            Atmp.colind[k] = j; //! indices starting from 0
            Atmp.val[k]    = value;
        } else {
            faspxx_chkerr(ERROR_INPUT_FILE, filename);
        }
    }

    fclose(fp);

    faspxx_format_dcoo_dcsr(&Atmp, A);
    faspxx_dcoo_free(&Atmp);
}

/**
 * \fn void faspxx_dmtx_read (const char *filename, dCSRmat *A)
 *
 * \brief Read A from matrix disk file in MatrixMarket general format -- Indices start
 * from 1
 *
 * \param filename  File name for matrix
 * \param A         Pointer to the CSR matrix
 *
 * \note File format:
 *   This routine reads a MatrixMarket general matrix from a mtx file.
 *   And it converts the matrix to dCSRmat format. For details of mtx format,
 *   please refer to http://math.nist.gov/MatrixMarket/.
 *
 * \note Indices start from 1, NOT 0!!!
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dmtx_read(const char* filename, dCSRmat* A)
{
    int i, j, m, n, nnz;
    INT innz; // index of nonzeros
    DBL value;

    FILE* fp = fopen(filename, "r");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: reading file %s...\n", __FUNCTION__, filename);

    skip_comments(fp); // skip the comments in the beginning

    if (fscanf(fp, "%d %d %d", &m, &n, &nnz) <= 0) {
        faspxx_chkerr(ERROR_INPUT_FILE, filename);
    }

    dCOOmat Atmp = faspxx_dcoo_create(m, n, nnz);

    innz = 0;

    while (innz < nnz) {
        if (fscanf(fp, "%d %d %le", &i, &j, &value) != EOF) {
            Atmp.rowind[innz] = i - 1; //! indices starting from 1
            Atmp.colind[innz] = j - 1; //! indices starting from 1
            Atmp.val[innz]    = value;
            innz              = innz + 1;
        } else {
            faspxx_chkerr(ERROR_INPUT_FILE, filename);
        }
    }

    fclose(fp);

    faspxx_format_dcoo_dcsr(&Atmp, A);
    faspxx_dcoo_free(&Atmp);
}

/**
 * \fn void faspxx_dmtxsym_read (const char *filename, dCSRmat *A)
 *
 * \brief Read A from matrix disk file in MatrixMarket sym format
 *
 * \param filename  File name for matrix
 * \param A         Pointer to the CSR matrix
 *
 * \note File format:
 *   This routine reads a MatrixMarket symmetric matrix from a mtx file.
 *   And it converts the matrix to dCSRmat format. For details of mtx format,
 *   please refer to http://math.nist.gov/MatrixMarket/.
 *
 * \note Indices start from 1, NOT 0!!!
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dmtxsym_read(const char* filename, dCSRmat* A)
{
    int i, j, m, n, nnz;
    int innz; // index of nonzeros
    DBL value;

    FILE* fp = fopen(filename, "r");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: reading file %s...\n", __FUNCTION__, filename);

    skip_comments(fp); // skip the comments in the beginning --zcs 06/30/2020

    if (fscanf(fp, "%d %d %d", &m, &n, &nnz) <= 0) {
        faspxx_chkerr(ERROR_INPUT_FILE, filename);
    }

    nnz          = 2 * (nnz - m) + m; // adjust for sym problem
    dCOOmat Atmp = faspxx_dcoo_create(m, n, nnz);

    innz = 0;

    while (innz < nnz) {
        if (fscanf(fp, "%d %d %le", &i, &j, &value) != EOF) {

            if (i == j) {
                Atmp.rowind[innz] = i - 1;
                Atmp.colind[innz] = j - 1;
                Atmp.val[innz]    = value;
                innz              = innz + 1;
            } else {
                Atmp.rowind[innz]     = i - 1;
                Atmp.rowind[innz + 1] = j - 1;
                Atmp.colind[innz]     = j - 1;
                Atmp.colind[innz + 1] = i - 1;
                Atmp.val[innz]        = value;
                Atmp.val[innz + 1]    = value;
                innz                  = innz + 2;
            }

        } else {
            faspxx_chkerr(ERROR_INPUT_FILE, filename);
        }
    }

    fclose(fp);

    faspxx_format_dcoo_dcsr(&Atmp, A);
    faspxx_dcoo_free(&Atmp);
}

/**
 * \fn void faspxx_dvec_read (const char *filename, dvector *b)
 *
 * \brief Read b from a disk file in array format
 *
 * \param filename  File name for vector b
 * \param b         Pointer to the dvector b (output)
 *
 * \note File Format:
 *   - nrow
 *   - val_j, j=0:nrow-1
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_read(const char* filename, dvector* b)
{
    int    i, n;
    DBL    value;
    size_t status;

    FILE* fp = fopen(filename, "r");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: reading file %s...\n", __FUNCTION__, filename);

    skip_comments(fp); // skip the comments in the beginning --zcs 06/30/2020

    status = fscanf(fp, "%d", &n);

    faspxx_dvec_alloc(n, b);

    for (i = 0; i < n; ++i) {

        status    = fscanf(fp, "%le", &value);
        b->val[i] = value;

        if (value > LARGE_DBL) {
            faspxx_dvec_free(b);
            fclose(fp);

            printf("### ERROR: Wrong value = %lf!\n", value);
            faspxx_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
        }
    }

    fclose(fp);
    faspxx_chkerr(status, filename);
}

/**
 * \fn void faspxx_dcoo_write (const char *filename, dCSRmat *A)
 *
 * \brief Write a matrix to disk file in IJ format (coordinate format)
 *
 * \param A         pointer to the dCSRmat matrix
 * \param filename  char for vector file name
 *
 * \note
 *      The routine writes the specified DBL vector in COO format.
 *      Refer to the reading subroutine \ref faspxx_dcoo_read.
 *
 * \note File format:
 *   - The first line of the file gives the number of rows, the
 *   number of columns, and the number of nonzeros.
 *   - Then gives nonzero values in i j a(i,j) format.
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcoo_write(const char* filename, dCSRmat* A)
{
    const INT m = A->row, n = A->col;
    INT       i, j;

    FILE* fp = fopen(filename, "w");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: writing to file %s...\n", __FUNCTION__, filename);

    fprintf(fp, "%d  %d  %d\n", m, n, A->nnz);
    for (i = 0; i < m; ++i) {
        for (j = A->IA[i]; j < A->IA[i + 1]; j++)
            fprintf(fp, "%d  %d  %0.15e\n", i, A->JA[j], A->val[j]);
    }

    fclose(fp);
}

/**
 * \fn void faspxx_dvec_write (const char *filename, dvector *vec)
 *
 * \brief Write a dvector to disk file
 *
 * \param vec       Pointer to the dvector
 * \param filename  File name
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_write(const char* filename, dvector* vec)
{
    INT m = vec->row, i;

    FILE* fp = fopen(filename, "w");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: writing to file %s...\n", __FUNCTION__, filename);

    fprintf(fp, "%d\n", m);

    for (i = 0; i < m; ++i) fprintf(fp, "%0.15e\n", vec->val[i]);

    fclose(fp);
}

/**
 * \fn void faspxx_ivec_write (const char *filename, ivector *vec)
 *
 * \brief Write a ivector to disk file in coordinate format
 *
 * \param vec       Pointer to the dvector
 * \param filename  File name
 *
 * \note The routine writes the specified INT vector in IJ format.
 *   - The first line of the file is the length of the vector;
 *   - After that, each line gives index and value of the entries.
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_ivec_write(const char* filename, ivector* vec)
{
    INT m = vec->row, i;

    FILE* fp = fopen(filename, "w");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: writing to file %s...\n", __FUNCTION__, filename);

    // write number of nonzeros
    fprintf(fp, "%d\n", m);

    // write index and value each line
    for (i = 0; i < m; ++i) fprintf(fp, "%d %d\n", i, vec->val[i] + 1);

    fclose(fp);
}

/**
 * \fn void faspxx_dvec_print (const INT n, dvector *u)
 *
 * \brief Print first n entries of a vector of DBL type
 *
 * \param n   An interger (if n=0, then print all entries)
 * \param u   Pointer to a dvector
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dvec_print(const INT n, dvector* u)
{
    INT i;
    INT NumPrint = n;

    if (n <= 0) NumPrint = u->row; // print all

    for (i = 0; i < NumPrint; ++i) printf("vec_%d = %15.10E\n", i, u->val[i]);
}

/**
 * \fn void faspxx_ivec_print (const INT n, ivector *u)
 *
 * \brief Print first n entries of a vector of INT type
 *
 * \param n   An interger (if n=0, then print all entries)
 * \param u   Pointer to an ivector
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_ivec_print(const INT n, ivector* u)
{
    INT i;
    INT NumPrint = n;

    if (n <= 0) NumPrint = u->row; // print all

    for (i = 0; i < NumPrint; ++i) printf("vec_%d = %d\n", i, u->val[i]);
}

/**
 * \fn void faspxx_dcsr_print (const dCSRmat *A)
 *
 * \brief Print out a dCSRmat matrix in coordinate format
 *
 * \param A   Pointer to the dCSRmat matrix A
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_print(const dCSRmat* A)
{
    const INT m = A->row, n = A->col;
    INT       i, j;

    printf("nrow = %d, ncol = %d, nnz = %d\n", m, n, A->nnz);
    for (i = 0; i < m; ++i) {
        for (j = A->IA[i]; j < A->IA[i + 1]; j++)
            printf("A_(%d,%d) = %+.10E\n", i, A->JA[j], A->val[j]);
    }
}

/**
 * \fn void faspxx_dcoo_print (const dCOOmat *A)
 *
 * \brief Print out a dCOOmat matrix in coordinate format
 *
 * \param A   Pointer to the dCOOmat matrix A
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcoo_print(const dCOOmat* A)
{
    INT k;

    printf("nrow = %d, ncol = %d, nnz = %d\n", A->row, A->col, A->nnz);
    for (k = 0; k < A->nnz; k++) {
        printf("A_(%d,%d) = %+.10E\n", A->rowind[k], A->colind[k], A->val[k]);
    }
}

/**
 * \fn void faspxx_dcsr_write_coo (const char *filename, const dCSRmat *A)
 *
 * \brief Print out a dCSRmat matrix in coordinate format for matlab spy
 *
 * \param filename   Name of file to write to
 * \param A          Pointer to the dCSRmat matrix A
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
void faspxx_dcsr_write_coo(const char* filename, const dCSRmat* A)
{

    INT i, j;

#if DEBUG_MODE > PRINT_MIN
    printf("nrow = %d, ncol = %d, nnz = %d\n", A->row, A->col, A->nnz);
#endif

    FILE* fp = fopen(filename, "w");

    if (fp == NULL) faspxx_chkerr(ERROR_OPEN_FILE, filename);

    printf("%s: writing to file %s...\n", __FUNCTION__, filename);

    // write dimension of the block matrix
    fprintf(fp, "%% dimension of the block matrix and nonzeros %d  %d  %d\n", A->row,
            A->col, A->nnz);

    for (i = 0; i < A->row; i++) {
        for (j = A->IA[i]; j < A->IA[i + 1]; j++) {
            fprintf(fp, "%d %d %+.10E\n", i + 1, A->JA[j] + 1, A->val[j]);
        }
    }

    fclose(fp);
}

/**
 * \fn static inline void skip_comments (FILE * fp)
 *
 * \brief Skip the comments in the beginning of data file
 *
 * \param FILE        File handler
 *
 * \author Li Zhao
 * \date   01/13/2024
 */
static inline void skip_comments(FILE* fp)
{
    while (1) {
        char buffer[500];
        int  loc = ftell(fp);                // record the current position
        int  val = fscanf(fp, "%s", buffer); // read in a string
        if (val != 1 || val == EOF) {
            printf("### ERROR: Could not get any data!\n");
            exit(ERROR_INPUT_FILE);
        }
        if (buffer[0] == '%' || buffer[0] == '!' || buffer[0] == '/') {
            if (fscanf(fp, "%*[^\n]")) { /* skip rest of line and do nothing */
            };
            continue;
        } else {
            fseek(fp, loc, SEEK_SET); // back to the beginning of this line
            break;
        }
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
