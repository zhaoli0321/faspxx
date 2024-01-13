/*! \file  AuxMemory.c
 *
 *  \brief Memory allocation and deallocation subroutines
 *
 *  \note  This file contains Level-0 (Aux) functions.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

/*---------------------------------*/
/*-- Declare External Functions  --*/
/*---------------------------------*/

#include "faspxx.h"

#if DLMALLOC
#include "dlmalloc.h"
#elif NEDMALLOC
#include "nedmalloc.h"
#ifdef __cplusplus
extern "C" {
#endif
void* nedcalloc(size_t no, size_t size);
void* nedrealloc(void* mem, size_t size);
void  nedfree(void* mem);
#ifdef __cplusplus
}
#endif
#endif

#if DEBUG_MODE > 1
extern unsigned long total_alloc_mem;
extern unsigned long total_alloc_count;
#endif

/*---------------------------------*/
/*--      Global Variables       --*/
/*---------------------------------*/

const int Million = 1048576; //! 1M = 1024*1024

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void * faspxx_mem_calloc (const LONGLONG size, const int type)
 *
 * \brief Allocate, initiate, and check memory
 *
 * \param size    Number of memory blocks
 * \param type    Size of memory blocks
 *
 * \return        Void pointer to the allocated memory
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 *
 * Modified by Chensong Zhang on 07/30/2013: print warnings if failed
 */
void* faspxx_mem_calloc(const unsigned int size, const unsigned int type)
{
    const LONGLONG tsize = size * type;
    void*          mem   = NULL;

#if DEBUG_MODE > 1
    printf("### DEBUG: Trying to allocate %.3lfMB RAM!\n", (DBL)tsize / Million);
#endif

    if (tsize > 0) {

#if DLMALLOC
        mem = dlcalloc(size, type);
#elif NEDMALLOC
        mem = nedcalloc(size, type);
#else
        mem = calloc(size, type);
#endif

#if DEBUG_MODE > 1
        total_alloc_mem += tsize;
        total_alloc_count++;
#endif
    }

    if (mem == NULL) {
        printf("### WARNING: Trying to allocate %lldB RAM...\n", tsize);
        printf("### WARNING: Cannot allocate %.4fMB RAM!\n", (DBL)tsize / Million);
    }

    return mem;
}

/**
 * \fn void * faspxx_mem_realloc (void * oldmem, const LONGLONG tsize)
 *
 * \brief Reallocate, initiate, and check memory
 *
 * \param oldmem  Pointer to the existing mem block
 * \param tsize   Size of memory blocks
 *
 * \return        Void pointer to the reallocated memory
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 *
 * Modified by Chensong Zhang on 07/30/2013: print error if failed
 */
void* faspxx_mem_realloc(void* oldmem, const LONGLONG tsize)
{
    void* mem = NULL;

#if DEBUG_MODE > 1
    printf("### DEBUG: Trying to allocate %.3lfMB RAM!\n", (DBL)tsize / Million);
#endif

    if (tsize > 0) {

#if DLMALLOC
        mem = dlrealloc(oldmem, tsize);
#elif NEDMALLOC
        mem = nedrealloc(oldmem, tsize);
#else
        mem = realloc(oldmem, tsize);
#endif
    }

    if (mem == NULL) {
        printf("### WARNING: Trying to allocate %lldB RAM!\n", tsize);
        printf("### WARNING: Cannot allocate %.3lfMB RAM!\n", (DBL)tsize / Million);
    }

    return mem;
}

/**
 * \fn void faspxx_mem_free (void *mem)
 *
 * \brief Free up previous allocated memory body and set pointer to NULL
 *
 * \param mem   Pointer to the memory body need to be freed
 *
 * \author Chensong Zhang
 * \date   2010/12/24
 *
 * Modified on 2018/01/10 by Chensong: Add output when mem is NULL
 */
void faspxx_mem_free(void* mem)
{
    if (mem) {
#if DLMALLOC
        dlfree(mem);
#elif NEDMALLOC
        nedfree(mem);
#else
        free(mem);
#endif

#if DEBUG_MODE > 1
        total_alloc_count--;
#endif
    } else {
#if DEBUG_MODE > 1
        printf("### WARNING: Trying to free an empty pointer!\n");
#endif
    }
}

/**
 * \fn void faspxx_mem_usage ( void )
 *
 * \brief Show total allocated memory currently
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 */
void faspxx_mem_usage(void)
{
#if DEBUG_MODE > 1
    printf("### DEBUG: Number of alloc = %ld, allocated memory = %.3fMB.\n",
           total_alloc_count, (DBL)total_alloc_mem / Million);
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
