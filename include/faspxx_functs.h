/*! \file  faspxx_functs.h
 *
 *  \brief Function decoration for the FASP++ package
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2024--Present by the FASP++ team. All rights reserved.           
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  \warning DO NOT EDIT!!! This file is automatically generated!
 */ 

#include "faspxx.h" 
 

/*-------- In file: AuxArray.c --------*/

void faspxx_darray_set(const INT n, DBL* x, const DBL val);

void faspxx_iarray_set(const INT n, INT* x, const INT val);

void faspxx_darray_cp(const INT n, const DBL* x, DBL* y);

void faspxx_iarray_cp(const INT n, const INT* x, INT* y);


/*-------- In file: AuxMemory.c --------*/

void* faspxx_mem_calloc(const unsigned int size, const unsigned int type);

void* faspxx_mem_realloc(void* oldmem, const LONGLONG tsize);

void faspxx_mem_free(void* mem);

void faspxx_mem_usage(void);


/*-------- In file: AuxMessage.c --------*/

void faspxx_itinfo(const INT ptrlvl, const INT stop_type, const INT iter,
                   const DBL relres, const DBL absres, const DBL factor);

void faspxx_cputime(const char* message, const DBL cputime);

void faspxx_message(const INT ptrlvl, const char* message);

void faspxx_chkerr(const SHORT status, const char* fctname);


/*-------- In file: AuxSort.c --------*/

INT faspxx_aux_BiSearch(const INT nlist, const INT* list, const INT value);

INT faspxx_aux_unique(INT numbers[], const INT size);

void faspxx_aux_merge(INT numbers[], INT work[], INT left, INT mid, INT right);

void faspxx_aux_msort(INT numbers[], INT work[], INT left, INT right);

void faspxx_aux_iQuickSort(INT* a, INT left, INT right);

void faspxx_aux_dQuickSort(DBL* a, INT left, INT right);

void faspxx_aux_iQuickSortIndex(INT* a, INT left, INT right, INT* index);

void faspxx_aux_dQuickSortIndex(DBL* a, INT left, INT right, INT* index);


/*-------- In file: AuxThreads.c --------*/

INT faspxx_get_num_threads(void);

INT faspxx_set_num_threads(const INT nthreads);

void faspxx_get_start_end(const INT procid, const INT nprocs, const INT n, INT* start,
                          INT* end);


/*-------- In file: AuxTiming.c --------*/

void faspxx_gettime(DBL* time);

void faspxx_gettime(DBL* time);


/*-------- In file: AuxVector.c --------*/

dvector faspxx_dvec_create(const INT m);

ivector faspxx_ivec_create(const INT m);

void faspxx_dvec_alloc(const INT m, dvector* u);

void faspxx_ivec_alloc(const INT m, ivector* u);

void faspxx_dvec_free(dvector* u);

void faspxx_ivec_free(ivector* u);

void faspxx_dvec_rand(const INT n, dvector* x);

void faspxx_dvec_set(INT n, dvector* x, const DBL val);

void faspxx_ivec_set(INT n, ivector* u, const INT m);

void faspxx_dvec_cp(const dvector* x, dvector* y);

DBL faspxx_dvec_maxdiff(const dvector* x, const dvector* y);

DBL faspxx_dvec_norm2diff(const dvector* x, const dvector* y);

void faspxx_dvec_symdiagscale(dvector* b, const dvector* diag);

BOOL faspxx_dvec_isnan(const dvector* u);


/*-------- In file: BlaArray.c --------*/

void faspxx_blas_darray_ax(const INT n, const DBL a, DBL* x);

void faspxx_blas_darray_axpy(const INT n, const DBL a, const DBL* x, DBL* y);

void faspxx_blas_darray_axpyz(const INT n, const DBL a, const DBL* x, const DBL* y,
                              DBL* z);

void faspxx_blas_darray_axpby(const INT n, const DBL a, const DBL* x, const DBL b,
                              DBL* y);

DBL faspxx_blas_darray_norm1(const INT n, const DBL* x);

DBL faspxx_blas_darray_norm2(const INT n, const DBL* x);

DBL faspxx_blas_darray_norminf(const INT n, const DBL* x);

DBL faspxx_blas_darray_dotprod(const INT n, const DBL* x, const DBL* y);


/*-------- In file: BlaFormat.c --------*/

SHORT faspxx_format_dcoo_dcsr(const dCOOmat* A, dCSRmat* B);

SHORT faspxx_format_dcsr_dcoo(const dCSRmat* A, dCOOmat* B);


/*-------- In file: BlaIO.c --------*/

void faspxx_dcsr_read(const char* filename, dCSRmat* A);

void faspxx_dcoo_read(const char* filename, dCSRmat* A);

void faspxx_dmtx_read(const char* filename, dCSRmat* A);

void faspxx_dmtxsym_read(const char* filename, dCSRmat* A);

void faspxx_dvec_read(const char* filename, dvector* b);

void faspxx_dcoo_write(const char* filename, dCSRmat* A);

void faspxx_dvec_write(const char* filename, dvector* vec);

void faspxx_ivec_write(const char* filename, ivector* vec);

void faspxx_dvec_print(const INT n, dvector* u);

void faspxx_ivec_print(const INT n, ivector* u);

void faspxx_dcsr_print(const dCSRmat* A);

void faspxx_dcoo_print(const dCOOmat* A);

void faspxx_dcsr_write_coo(const char* filename, const dCSRmat* A);


/*-------- In file: BlaSparseCOO.c --------*/

dCOOmat faspxx_dcoo_create(const INT m, const INT n, const INT nnz);

void faspxx_dcoo_alloc(const INT m, const INT n, const INT nnz, dCOOmat* A);

void faspxx_dcoo_free(dCOOmat* A);


/*-------- In file: BlaSparseCSR.c --------*/

dCSRmat faspxx_dcsr_create(const INT m, const INT n, const INT nnz);

void faspxx_dcsr_alloc(const INT m, const INT n, const INT nnz, dCSRmat* A);

void faspxx_dcsr_free(dCSRmat* A);

void faspxx_dcsr_sort(dCSRmat* A);

void faspxx_dcsr_getdiag(INT n, const dCSRmat* A, dvector* diag);

void faspxx_dcsr_diagpref(dCSRmat* A);

void faspxx_dcsr_cp(const dCSRmat* A, dCSRmat* B);

INT faspxx_dcsr_trans(const dCSRmat* A, dCSRmat* AT);

void faspxx_dcsr_compress(const dCSRmat* A, dCSRmat* B, const DBL dtol);

SHORT faspxx_dcsr_compress_inplace(dCSRmat* A, const DBL dtol);


/*-------- In file: BlaSparseUtil.c --------*/

void faspxx_sparse_abybms_(INT* ia, INT* ja, INT* ib, INT* jb, INT* nap, INT* map,
                           INT* mbp, INT* ic, INT* jc);

void faspxx_sparse_abyb_(INT* ia, INT* ja, DBL* a, INT* ib, INT* jb, DBL* b, INT* nap,
                         INT* map, INT* mbp, INT* ic, INT* jc, DBL* c);

void faspxx_sparse_iit_(INT* ia, INT* ja, INT* na, INT* ma, INT* iat, INT* jat);

void faspxx_sparse_aat_(INT* ia, INT* ja, DBL* a, INT* na, INT* ma, INT* iat, INT* jat,
                        DBL* at);

void faspxx_sparse_aplbms_(INT* ia, INT* ja, INT* ib, INT* jb, INT* nab, INT* mab,
                           INT* ic, INT* jc);

void faspxx_sparse_aplusb_(INT* ia, INT* ja, DBL* a, INT* ib, INT* jb, DBL* b, INT* nab,
                           INT* mab, INT* ic, INT* jc, DBL* c);

void faspxx_sparse_rapms_(INT* ir, INT* jr, INT* ia, INT* ja, INT* ip, INT* jp,
                          INT* nin, INT* ncin, INT* iac, INT* jac, INT* maxrout);

void faspxx_sparse_wtams_(INT* jw, INT* ia, INT* ja, INT* nwp, INT* map, INT* jv,
                          INT* nvp, INT* icp);

void faspxx_sparse_wta_(INT* jw, DBL* w, INT* ia, INT* ja, DBL* a, INT* nwp, INT* map,
                        INT* jv, DBL* v, INT* nvp);

void faspxx_sparse_ytxbig_(INT* jy, DBL* y, INT* nyp, DBL* x, DBL* s);

void faspxx_sparse_ytx_(INT* jy, DBL* y, INT* jx, DBL* x, INT* nyp, INT* nxp, INT* icp,
                        DBL* s);

void faspxx_sparse_rapcmp_(INT* ir, INT* jr, DBL* r, INT* ia, INT* ja, DBL* a, INT* ipt,
                           INT* jpt, DBL* pt, INT* nin, INT* ncin, INT* iac, INT* jac,
                           DBL* ac, INT* idummy);

ivector faspxx_sparse_mis(dCSRmat* A);


/*-------- In file: BlaSpmvCSR.c --------*/

SHORT faspxx_blas_dcsr_add(const dCSRmat* A, const DBL alpha, const dCSRmat* B,
                           const DBL beta, dCSRmat* C);

void faspxx_blas_dcsr_axm(dCSRmat* A, const DBL alpha);

void faspxx_blas_dcsr_mxv(const dCSRmat* A, const DBL* x, DBL* y);

void faspxx_blas_dcsr_mxv_agg(const dCSRmat* A, const DBL* x, DBL* y);

void faspxx_blas_dcsr_aAxpy(const DBL alpha, const dCSRmat* A, const DBL* x, DBL* y);

void faspxx_blas_dcsr_aAxpy_agg(const DBL alpha, const dCSRmat* A, const DBL* x, DBL* y);

DBL faspxx_blas_dcsr_vmv(const dCSRmat* A, const DBL* x, const DBL* y);

void faspxx_blas_dcsr_mxm(const dCSRmat* A, const dCSRmat* B, dCSRmat* C);

void faspxx_blas_dcsr_rap(const dCSRmat* R, const dCSRmat* A, const dCSRmat* P,
                          dCSRmat* RAP);

void faspxx_blas_dcsr_rap_agg(const dCSRmat* R, const dCSRmat* A, const dCSRmat* P,
                              dCSRmat* RAP);

void faspxx_blas_dcsr_rap_agg1(const dCSRmat* R, const dCSRmat* A, const dCSRmat* P,
                               dCSRmat* B);

void faspxx_blas_dcsr_ptap(const dCSRmat* Pt, const dCSRmat* A, const dCSRmat* P,
                           dCSRmat* Ac);

dCSRmat faspxx_blas_dcsr_rap2(INT* ir, INT* jr, DBL* r, INT* ia, INT* ja, DBL* a,
                              INT* ipt, INT* jpt, DBL* pt, INT n, INT nc, INT* maxrpout,
                              INT* ipin, INT* jpin);

void faspxx_blas_dcsr_rap4(dCSRmat* R, dCSRmat* A, dCSRmat* P, dCSRmat* B,
                           INT* icor_ysk);


/*-------- In file: KryPCG.c --------*/

INT faspxx_solver_dcsr_pcg(dCSRmat* A, dvector* b, dvector* u, precond* pc,
                           const DBL tol, const INT MaxIt, const SHORT StopType,
                           const SHORT PrtLvl);


/*-------- In file: KryPVFGMRES.c --------*/

INT faspxx_solver_dcsr_pvfgmres(dCSRmat* A, dvector* b, dvector* x, precond* pc,
                                const DBL tol, const INT MaxIt, const SHORT restart,
                                const SHORT StopType, const SHORT PrtLvl);


/*-------- In file: XtrMumps.c --------*/

INT faspxx_solver_mumps(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl);

INT faspxx_solver_mumps_steps(dCSRmat* ptrA, dvector* b, dvector* u, Mumps_data* mumps);

Mumps_data faspxx_mumps_factorize(dCSRmat* ptrA, dvector* b, dvector* u,
                                  const SHORT prtlvl);

INT faspxx_mumps_solve(dCSRmat* ptrA, dvector* b, dvector* u, Mumps_data mumps,
                       const SHORT prtlvl);

INT faspxx_mumps_free(Mumps_data* mumps);


/*-------- In file: XtrPardiso.c --------*/

INT faspxx_solver_pardiso(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl);

INT faspxx_pardiso_factorize(dCSRmat* ptrA, Pardiso_data* pdata, const SHORT prtlvl);

INT faspxx_pardiso_solve(dCSRmat* ptrA, dvector* b, dvector* u, Pardiso_data* pdata,
                         const SHORT prtlvl);

INT faspxx_pardiso_free(Pardiso_data* pdata);


/*-------- In file: XtrSuperlu.c --------*/

INT faspxx_solver_superlu(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl);

INT faspxx_superlu_factorize(dCSRmat* ptrA, dvector* b, SuperLU_data* superlu_data,
                             const SHORT prtlvl);

INT faspxx_superlu_solve(dCSRmat* ptrA, dvector* b, dvector* u,
                         SuperLU_data* superlu_data, const SHORT prtlvl);

INT faspxx_superlu_free(SuperLU_data* superlu_data);


/*-------- In file: XtrUmfpack.c --------*/

INT faspxx_solver_umfpack(dCSRmat* ptrA, dvector* b, dvector* u, const SHORT prtlvl);

void* faspxx_umfpack_factorize(dCSRmat* ptrA, const SHORT prtlvl);

INT faspxx_umfpack_solve(dCSRmat* ptrA, dvector* b, dvector* u, void* Numeric,
                         const SHORT prtlvl);

INT faspxx_umfpack_free(void* Numeric);

 
/* End of faspxx_functs.h */
