#
# Makefile for Eve
#
# Usage
#    make
# Variables can be set as command line input or in the --- SWITCHES --- section below
#
# make targets are: all, clean, cleanclean
# 
# 
# Netcdf is:           netcdf=netcdf3, netcdf4 or ""
# Static is:           static=static, shared, dynamic
#
# Written Matthias Cuntz & Juliane Mai, UFZ Leipzig, Germany, Aug. 2011 - matthias.cuntz@ufz.de

export
SHELL = /bin/bash

# --- SWITCHES -------------------------------------------------------
MAKEPATH = . # where are the make files (. is current directory, .. is parent directory)
SRCPATH  = . # where are the source files; use fortran_test to run test directory
#SRCPATH  = fortran_test
PROGPATH = . # where shall be the executable

PROGNAME = RNG # Name of executable

# Switches
# Possible releases: release, debug
release  = release
# Possible netcdf versions: netcdf3, netcdf4
netcdf   = 
# Possible linking: static, shared, dynamic (last two are equal)
static   = shared

# --- CHECKS ---------------------------------------------------
# check input
ifneq ($(MAKECMDGOALS),$(findstring $(MAKECMDGOALS),clean cleanclean))
    ifeq (,$(findstring $(release),debug release))
        $(error Error: release '$(release)' not found; must be in 'debug release')
    endif
    ifneq ($(netcdf),)
        ifeq (,$(findstring $(netcdf),netcdf3 netcdf4))
            $(error Error: netcdf '$(netcdf)' not found; must be in 'netcdf3 netcdf4')
        endif
    endif
    ifeq (,$(findstring $(static),static shared dynamic))
        $(static Error: static '$(static)' not found; must be in 'static shared dynamic')
    endif
endif

OBJPATH = $(strip $(SRCPATH))/.$(strip $(release))

# --- DEFAULTS ---------------------------------------------------
# These variables will be used to compile
F90      =
F90FLAGS =
DEFINES  =
INCLUDES =
# and link, and therefore set below
LD       =
LDFLAGS  =
LIBS     =

# --- COMPILER ---------------------------------------------------
# v12
# INTEL = /usr/local/intel/composerxe-2011.4.191
# INTELLIB = $(INTEL)/compiler/lib/intel64
# v11
INTEL    = /opt/intel/Compiler/11.1/075
INTELBIN = $(INTEL)/bin/intel64
INTELLIB = /usr/local/intel/11.1.075

F90 = $(INTELBIN)/ifort
ifeq ($(release),debug)
    F90FLAGS = -check all -warn all -g -debug -traceback -fp-stack-check -O0 -debug
else
    # -vec-report1 to see vectorized loops; -vec-report2 to see also non-vectorized loops
    F90FLAGS  = -O3 -vec-report0 -override-limits
endif
F90FLAGS += -cpp -fp-model precise -openmp -m64 -module $(OBJPATH)
LDFLAGS  += -openmp
DEFINES  += -DINTEL

ifeq ($(static),static)
    LIBS += -static-intel -Bstatic -Wl,--start-group
else
    LIBS += -Bdynamic
endif
LIBS += -L$(INTELLIB) -limf -lm -lsvml
ifneq ($(static),static)
    LIBS += -lintlc
endif
RPATH += -Wl,-rpath,$(INTELLIB)

# --- IMSL ---------------------------------------------------
IMSL    = /global/local/imsl/imsl/fnl700/rdhin111e64
IMSLINC = $(IMSL)/include
IMSLLIB = $(IMSL)/lib

LIBS += -z muldefs
ifneq ($(static),static)
    LIBS += -i_dynamic
endif
LIBS += -L$(IMSLLIB) -limsl -limslscalar -limsllapack_imsl -limslblas_imsl -limsls_err -limslmpistub -limslsuperlu -limslhpc_l
#LIBS  += -L$(IMSLLIB) -limsl -limslscalar -limsllapack_vendor -limslblas_vendor -limsls_err -limslmpistub -limslsuperlu -limslhpc_l
RPATH += -Wl,-rpath,$(IMSLLIB)

# --- NETCDF ---------------------------------------------------
ifneq ($(netcdf),)
    ifeq ($(netcdf),$(findstring $(netcdf),netcdf3 netcdf4))
        NCDIR =
        ifeq ($(netcdf),netcdf3)
            #NCDIR = /usr/local/netcdf/3.6.3_intel_12.0.4
            NCDIR = /usr/local/netcdf/3.6.3_intel11.1.075
        else
            #NCDIR = /usr/local/netcdf/4.1.1_intel_12.0.4
            NCDIR = /usr/local/netcdf/4.1.1_intel11.1.075
        endif
        NCINC = $(strip $(NCDIR))/include
        NCLIB = $(strip $(NCDIR))/lib

        INCLUDES += -I$(NCINC)
        LIBS     += -L$(NCLIB) -lnetcdf -lnetcdff
        RPATH    += -Wl,-rpath,$(NCLIB)

        # libraries for netcdf4, ignored for netcdf3
        ifeq ($(netcdf),netcdf4)
            SZLIB     = /usr/local/szip/2.1/lib
            HDF5LIB   = /usr/local/hdf5/1.8.6/lib
            LIBS     += -lz -L$(SZLIB) -lsz -L$(HDF5ZLIB) -lhdf5 -lhdf5_hl
            RPATH    += -Wl,-rpath,$(SZLIB) -Wl,-rpath,$(HDF5LIB)
        endif
    endif
endif

# --- MKL ---------------------------------------------------
MKL       = /usr/local/intel/composerxe-2011.4.191/mkl
MKLINC    = $(MKL)/include
MKLLIB    = $(MKL)/lib/intel64
INCLUDES += -I$(MKLINC)
LIBS     += -L$(MKLLIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential # -lmkl_intel_thread
RPATH    += -Wl,-rpath,$(MKLLIB)

# --- LAPACK ---------------------------------------------------
# LAPACK    = /usr
# LAPACKLIB = $(LAPACK)/lib64
# LIBS     += -L$(LAPACKLIB) -lblas -llapack
# RPATH    += -Wl,-rpath,$(LAPACKLIB)

# --- FINISH SETUP ---------------------------------------------------
ifeq ($(release),debug)
    DEFINES += -DDEBUG
endif

ifeq ($(static),static)
    LIBS += -Wl,--end-group
else
    LIBS += $(RPATH)
endif

# Progs
PROG          = $(strip $(PROGPATH))/$(strip $(PROGNAME))
MAKEPROG      = $(strip $(MAKEPATH))/Makefile2
MAKEDEPSPROG  = $(strip $(MAKEPATH))/makedeps.pl

# --- TARGETS ---------------------------------------------------
LD       := $(F90)
# A vars contain source dir informations
ASRCS    := $(wildcard $(strip $(SRCPATH))/*.f90)
SRCS     := $(notdir $(ASRCS))
AOBJS    := $(SRCS:.f90=.o)
# main driver routines can be excluded and added at targets in $(MAKEPROG)
EXCL     := 
FAOBJS   := $(filter-out $(EXCL), $(AOBJS))
OBJS     := $(addprefix $(OBJPATH)/, $(FAOBJS))

.SUFFIXES: .f90 .o

# targets for executables
all: makedeps makedirs
	$(MAKE) -f $(MAKEPROG)

# helper targets
makedeps:
	rm -f make.deps
	rm -f tmp.gf3.*
	$(MAKEDEPSPROG) $(OBJPATH) $(SRCPATH)

makedirs:
	if [ ! -d $(OBJPATH) ] ; then mkdir $(OBJPATH) ; fi

clean:
	rm -f $(OBJPATH)/*.o $(OBJPATH)/*.mod make.deps $(PROG)

cleanclean: clean
	rm -rf $(strip $(SRCPATH))/.release $(strip $(SRCPATH))/.debug $(PROG).dSYM
