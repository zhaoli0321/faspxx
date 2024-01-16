########################################################################
# Fast Auxiliary Space Preconditioners (FASP++)
# 
# This is the Makefile for FASP++ test! 
#
########################################################################

#==============================================================================#
# User compilers                                                               #
# FASP++ has been tested with many different compilers                         #
#==============================================================================#
CC  = icc 
CPP = icpc 
FC  = ifort -nofor_main
AR  = ar ruc

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-O2

#==============================================================================#
# OpenMP Support                                                               #
# Uncomment the following line to turn on the OpenMP support                   #
#==============================================================================#
# BOPT+=-fopenmp

########################################################################
# Root directory for FASPXX package
########################################################################
FASPXXDIR = ../faspxx
INCLUDE = -I$(FASPXXDIR)/include -I./src 
FASXXPLIB = $(FASPXXDIR)/lib/libfaspxx.a 


########################################################################
# External libraries for FASPXX package 
########################################################################
# We recommend using Pardiso.
# (1) If you wish to use Pardiso as a direct solver, uncomment the following
MKL_DIR = /opt/intel/mkl
MKL_INC = -I $(MKL_DIR)/include
MKL_LIB = -L $(MKL_DIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


# (2) If you wish to use MUMPS as a direct solver, uncomment the following
# MUMPS_DIR = /home/spring/zhaoli/faspTest/MUMPS_4.10.0.ifort # The root directory of [MUMPS_DIR] is specified by the user
# MUMPS_INC = -I$(MUMPS_DIR)/libseq -I$(MUMPS_DIR)/include
# MUMPS_LIB = -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -lpord -L$(MUMPS_DIR)/libseq -lmpiseq -lblas


# (3) If you wish to use SuperLU as a direct solver, uncomment the following
# SUPERLU_DIR = /home/spring/zhaoli/faspTest/superlu-install # The root directory of [SUPERLU_DIR] is specified by the user
# SUPERLU_INC = -I$(SUPERLU_DIR)/include
# SUPERLU_LIB = -L$(SUPERLU_DIR)/lib -lsuperlu


# (4) If you want to use SUITESPARSE (or UMFPACK) as a direct solver, uncomment the following
# SUITESPARSE_DIR=/home/spring/zhaoli/faspTest/SuiteSparse # The root directory of [SUITESPARSE_DIR] is specified by the user
# SUITESPARSE_INC=-I$(SUITESPARSE_DIR)/include -I$(SUITESPARSE_DIR)/AMD/Include -I$(SUITESPARSE_DIR)/UMFPACK/Include 
# SUITESPARSE_LIB=-L$(SUITESPARSE_DIR)/lib -lsuitesparseconfig -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd -L/opt/local/lib -lmetis


#==============================================================================#
# User preprocessing definitions                                               #
# e.g., using PARDISO, set -DWITH_PARDISO=1, otherwise -DWITH_PARDISO    =0    #
# e.g., using MUMPS,   set -DWITH_MUMPS  =1, otherwise -DWITH_MUMPS      =0    #
# e.g., using SUPERLU, set -DWITH_SUPERLU=1, otherwise -DWITH_SUPERLU    =0    #
# e.g., using UMFPACK, set -DWITH_UMFPACK=1, otherwise -DWITH_SUITESPARSE=0    #
#==============================================================================#
CDEFS=
CDEFS+=-DWITH_PARDISO=1 -DWITH_MUMPS=0 -DWITH_SUPERLU=0 -DWITH_SUITESPARSE=0 

COPTS=$(BOPT)
CINCLUDES=$(INCLUDE)
CINCLUDES+=$(MKL_INC) $(MUMPS_INC) $(SUPERLU_INC) $(SUITESPARSE_INC) 
CFLAGS=$(CDEFS) $(COPTS) $(CINCLUDES)

FOPTS=$(BOPT)
FDEFS=$(CDEFS)
FINCLUDES=$(CINCLUDES)
FFLAGS=$(FDEFS) $(FOPTS) $(FINCLUDES)

########################################################################
# Link options
########################################################################
LINKOPTS=$(BOPT)

LIBS=$(TESTLIB) $(FASPLIB)
LIBS+=$(MKL_LIB) $(MUMPS_LIB) $(SUPERLU_LIB) $(SUITESPARSE_LIB)

CLFLAGS=-lstdc++ $(LINKOPTS) $(LIBS)
FLFLAGS=-lm $(LINKOPTS) $(LIBS)

CSRCDIR = ./src
FSRCDIR = ./src
TESTLIB = ./lib/libfaspxx.a

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .cc .cpp .for .f .f77 .f90 .f95
#
FSRC := $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.for))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f77))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f90))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f95))
CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))
CSRC += $(foreach dir,$(EXTRDIR),$(wildcard $(EXTRDIR)/*.c))
#
OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSF += $(patsubst %.f77,%.o,$(FSRC))
OBJSF += $(patsubst %.f90,%.o,$(FSRC))
OBJSF += $(patsubst %.f95,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))
#
.for.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f90.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F90 object $@'
	@$(AR) $(TESTLIB) $@
#
.f95.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F95 object $@'
	@$(AR) $(TESTLIB) $@
#
.c.o:
	@$(CC) -c $< -o $@ $(CFLAGS)
	@echo 'Building C object $@'
	@$(AR) $(TESTLIB) $@
#
.cpp.o:
	@$(CPP) -c $< -o $@ $(CFLAGS)
	@echo 'Building CPP object $@'
	@$(AR) $(TESTLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(TESTLIB)

########################################################################
# Link
########################################################################

all: $(ALLPROG) TestReadData TestPCG TestPVFGMRES TestDirectSolver

Default:
	regression

headers: 
	cat $(CSRCDIR)/*.c \
	| awk -v name="faspxx_functs.h" -f ./util/mkheaders.awk > ./include/faspxx_functs.h

$(TESTLIB): $(OBJSC) $(OBJSF)
	@ranlib $(TESTLIB)
	@echo 'Generating library $@'

lib: $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)
	@echo 'Generating library $@'

########################################################################
# Some test problems
########################################################################

TestReadData:
	@$(CC) $(CFLAGS) -c example/TestReadData.c -o example/TestReadData.o
	@$(FC) $(LOPT) example/TestReadData.o $(FLFLAGS) -o example/TestReadData.ex
	@echo 'Building executable $@'

TestPCG:
	@$(CC) $(CFLAGS) -c example/TestPCG.c -o example/TestPCG.o
	@$(FC) $(LOPT) example/TestPCG.o $(FLFLAGS) -o example/TestPCG.ex
	@echo 'Building executable $@'

TestPVFGMRES:
	@$(CC) $(CFLAGS) -c example/TestPVFGMRES.c -o example/TestPVFGMRES.o
	@$(FC) $(LOPT) example/TestPVFGMRES.o $(FLFLAGS) -o example/TestPVFGMRES.ex
	@echo 'Building executable $@'

TestDirectSolver:
	@$(CC) $(CFLAGS) -c example/TestDirectSolver.c -o example/TestDirectSolver.o
	@$(FC) $(LOPT) example/TestDirectSolver.o $(FLFLAGS) -o example/TestDirectSolver.ex
	@echo 'Building executable $@'

########################################################################
# Clean up
########################################################################

.PHONY : clean distclean help

clean:
	@rm -f $(CSRCDIR)/*.o
	@rm -f $(FSRCDIR)/*.o
	@rm -f example/*.o

distclean:
	@make clean
	@rm -f lib/*.a
	@rm -f *~ *.ex *.out
	@rm -f $(CSRCDIR)/*~
	@rm -f $(FSRCDIR)/*~
	@rm -f example/*~
	@rm -f example/*.ex

help:
	@echo "======================================================"
	@echo " Fast Auxiliary Space Preconditioners (FASP++)        "
	@echo "======================================================"
	@echo " "
	@echo " make            : build all exe files "
	@echo " make headers    : build the header file automatically"
	@echo " make clean      : clean all obj files "
	@echo " make distclean  : clean all obj, exe, bak, out files "
	@echo " make help       : show this screen "
	@echo " "
