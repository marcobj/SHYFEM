
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#------------------------------------------------------------
#
# This is the Rules.make file for shyfem
#
#------------------------------------------------------------

#------------------------------------------------------------
# This file defines various parameters to be used
# during compilation. Please customize the first part
# of this file to your needs.
#------------------------------------------------------------

##############################################
# User defined parameters and flags
##############################################

##############################################
# Compiler profile
##############################################
#
# You can either compile with maximum speed
# or with checks enabled, depending on your
# application.
# If in doubt, please leave as it is.
#
##############################################

COMPILER_PROFILE = NORMAL
#COMPILER_PROFILE = CHECK
#COMPILER_PROFILE = SPEED

##############################################
# Compiler
##############################################
#
# Please choose a compiler. Compiler options are
# usually correct. You might check below.
#
# Available options for the Fortran compiler are:
#
# GNU_G77		->	g77
# GNU_GFORTRAN		->	gfortran
# INTEL			->	ifort
# PORTLAND		->	pgf90
# IBM			->	xlf
# PGI			->	nvfortran
#
# Available options for the C compiler are:
#
# GNU_GCC		->	gcc
# INTEL			->	icc
# IBM			->	xlc
# PGI			->	nvc
#
##############################################

#FORTRAN_COMPILER = GNU_G77
FORTRAN_COMPILER = GNU_GFORTRAN
#FORTRAN_COMPILER = INTEL
#FORTRAN_COMPILER = PORTLAND
#FORTRAN_COMPILER = IBM
#FORTRAN_COMPILER = PGI

C_COMPILER = GNU_GCC
#C_COMPILER = INTEL
#C_COMPILER = IBM
#C_COMPILER = PGI

##############################################
# Parallel compilation
##############################################
#
# For some compilers you can specify
# parallel execution of some parts of the
# code. This can be specified here. Please
# switch back to serial execution if you are
# in doubt of the results.
#
# If running in parallel you get a segmentation fault
# then you might have to run one of the following
# commands on the command line, before you run
# the model:
#
#	ulimit -a		# gives you the system settings
#	ulimit -s 32000		# sets stack to 32MB
#	ulimit -s unlimited	# sets stack to unlimited
#
# For the Intel compiler you may also try inserting a
# command similar to the following in your .bashrc file:
#
#       export KMP_STACKSIZE=32M
#
# There are two ways of parallelizing. One is OMP
# and the other is MPI. Please note that MPI is
# highly experimental and not recommended to be used
# at this stage.
#
##############################################

PARALLEL_OMP = false
#PARALLEL_OMP = true

PARALLEL_MPI = NONE
#PARALLEL_MPI = NODE
#PARALLEL_MPI = ELEM

##############################################
# Partition library for domain decomposition
##############################################
#
# Here you specify the external module to be used
# for the partition of the grid. The software
# should be downloaded and installed separately.
#
# There are different options for the software:
#
#  - METIS: http://glaros.dtc.umn.edu/gkhome/views/metis
#  - ...
#
# The variable PARTSDIR indicates the directory
# where the library and its include files can be found.
# Please leave out the final lib specification.
# This is mandatory only if the library has been
# installed in a non-standard place.
#
##############################################

PARTS = NONE
#PARTS = METIS
#PARTSDIR = /usr/local
#PARTSDIR = $(HOME)/lib/metis

##############################################
# Solver for matrix solution
##############################################
#
# Here you have to specify what solver you want
# to use for the solution of the system matrix.
# You have a choice between Gaussian elimination,
# iterative solution (Sparskit), the Pardiso
# solver, and the Paralution solver.
# The gaussian elimination is the most robust.
# Sparskit is very fast on big matrices.
# To use the Pardiso solver you must have the 
# libraries installed. It is generally faster
# then the default method. Please note that
# in order to use Pardiso you must also use
# the INTEL compiler. The option to use the
# Pardiso solver is considered experimental.
# Finally you can use the Paralution solver,
# which allows you to use the GPU. 
# This solver is still in a testing phase,
# it should be faster with large matrices
# and with a good GPU.
# Moreover, it can reduce the CPU usage.
# If you are unsure what solver to use, leave
# the default which is SPARSKIT.
#
##############################################

#SOLVER = GAUSS
SOLVER = SPARSKIT
#SOLVER = PARDISO
#SOLVER = PARALUTION
#SOLVER = PETSC
#SOLVER = PETSC_AmgX

##############################################
#
# PETSC and PETSC_AmgX solvers
#
##############################################

# PETSC_DIR it the path to the PETSc installation folder, it is 
# needed for both the PETSc and the PETSc_AmgX solvers
PETSC_DIR =

# The next 4 paths must be filled in for the PETSc_AmgX solver only.

# AMGX_C_WRAPPER_DIR is the path to the amgx-c-wrapper folder
# (https://github.com/tobiashuste/amgx-c-wrapper/ : 
# amgx-c-wrapper is a C interface for the C++ AmgXWrapper) 
AMGX_C_WRAPPER_DIR = ../amgx-c-wrapper/amgx-c-wrapper

# AMGX_WRAPPER_DIR is the path to the C++ AmgXWrapper folder 
# (https://github.com/barbagroup/AmgXWrapper : a wrapper that 
# simplifies the usage of AmgX when using AmgX together with PETSc)
AMGX_WRAPPER_DIR =

# AMGX_DIR is the path to the installation folder of NVIDIA/AmgX, 
# a GPU accelerated core solver library (https://github.com/NVIDIA/AMGX)
AMGX_DIR =

# CUDA_DIR is the path to the installation folder of CUDA where the following 
# libraries must be installed: cusolver cusparse cublas cuda cudart
CUDA_DIR =

##############################################
#
# Paralution solver
#
# In order to use the paralution solver
# please see the README file in fempara,
# set SOLVER = PARALUTION, set the GPU variable,
# and define the PARADIR directory below.
#
# this feature is still experimental - no support
#
##############################################

#PARADIR = $(HOME)/my_paralution

GPU=NONE
#GPU=OpenCL
#GPU=CUDA
#GPU=MIC

##############################################
# NetCDF library
##############################################
#
# If you want output in NetCDF format and have
# the library installed, then you can specify
# it here. Normally the place where the netcdf
# files reside can be found automatically. However,
# if the libraries are in some non standard 
# place the directory where the netcdf files
# (include and libraries) reside must also be 
# indicated.
#
# You can normally find the directory by one
# of the following commands:
#   ldconfig -p | grep libnetcdff
#   whereis libnetcdff
#   locate libnetcdff.a
# Do not include the final /lib part of the directory.
#
##############################################

NETCDF = false
#NETCDF = true
#NETCDFDIR =

##############################################
# GOTM library
##############################################
#
# This software comes with a version of the
# GOTM turbulence model. It is needed if you
# want to run the model in 3D mode, unless you
# use constant vertical viscosity and diffusivity.
# If you have problems compiling the GOTM
# library, you can set GOTM to false. This will use
# an older version of the program without the
# need of the external library. However, this 
# option is not recommended.
#
##############################################

#GOTM = false
GOTM = true

##############################################
# Ecological models
##############################################
#
# The model also comes with code for some
# ecological models. You can activiate this
# code here. Choices are between EUTRO,
# ERSEM, AQUABC, and BFM.
# The BFM model is still experimental.
#
##############################################

ECOLOGICAL = NONE
#ECOLOGICAL = EUTRO
#ECOLOGICAL = ERSEM
#ECOLOGICAL = AQUABC
#ECOLOGICAL = BFM

##############################################
#
# BFM model - in order to use the BFM model
# please see the README file in fembfm
# and set the BFMDIR directory below.
#
# this feature is still experimental - no support
#
##############################################

#BFMDIR = /gpfs/work/OGS18_PRACE_P_0/SHYFEM_BFM/bfm
#BFMDIR = /home/georg/appl/donata/bfm/bfmv5
#BFMDIR = /home/georg/appl/donata/bfm/BiogeochemicalFluxModel-5.1.0
#BFMDIR = $(HOME)/BFM

##############################################
# Experimental features
##############################################

FLUID_MUD = false
#FLUID_MUD = true

##############################################
# end of user defined parameters and flags
##############################################

#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------
# Normally it should not be necessary to change anything beyond here.
#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------

##############################################
# DEFINE VERSION
##############################################

RULES_MAKE_VERSION = 1.6
DISTRIBUTION_TYPE = experimental

##############################################
# DEFINE DIRECTORIES
##############################################

DEFDIR  = $(HOME)
DIRLIB  = $(FEMDIR)/femlib
MODDIR  = 
MODDIR  = $(DIRLIB)/mod

LIBX = -L/usr/X11R6/lib -L/usr/X11/lib -L/usr/lib/X11  -lX11

##############################################
# check compatibility of options
##############################################

RULES_MAKE_PARAMETERS = RULES_MAKE_OK
RULES_MAKE_MESSAGE = ""

ifeq ($(FORTRAN_COMPILER),GNU_G77)
  ifeq ($(GOTM),true)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "g77 compiler and GOTM=true are incompatible"
  endif
  ifeq ($(PARALLEL_OMP),true)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "g77 and PARALLEL_OMP=true are incompatible"
  endif
  ifneq ($(PARALLEL_MPI),NONE)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "g77 and PARALLEL_MPI/=NONE are incompatible"
  endif
endif

ifneq ($(FORTRAN_COMPILER),INTEL)
  ifeq ($(SOLVER),PARDISO)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "Pardiso solver needs Intel compiler"
  endif
endif

ifneq ($(SOLVER),PARALUTION)
  ifneq ($(GPU),NONE)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "Use GPU=NONE without PARALUTION solver"
  endif
endif

ifeq ($(SOLVER),PARALUTION)
  ifneq ($(PARALLEL_OMP),true)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "Paralution solver needs PARALLEL_OMP=true"
  endif
  ifeq ($(PARADIR),)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "PARALUTION solver needs PARADIR directory"
  endif
endif

ifeq ($(C_COMPILER),INTEL)
  ifneq ($(FORTRAN_COMPILER),INTEL)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "INTEL C works only with INTEL Fortran compiler"
  endif
endif

ifeq ($(ECOLOGICAL),BFM)
  ifeq ($(BFMDIR),)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "BFM model needs BFMDIR directory"
  endif
  ifeq ($(NETCDF),false)
    RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
    RULES_MAKE_MESSAGE = "BFM model needs NETCDF support"
  endif
endif

##############################################
# some utilities
##############################################

# do "make print-VARIABLE" to see value of $VARIABLE

print-% : ; @echo $* = $($*)

# can also use following construct:
#
# $(warning this is a warning)
# $(warning var=$(VAR))

#------------------------------------------------------------

##############################################
# COMPILER OPTIONS
##############################################

##############################################
# General Compiler options
##############################################

# if unsure please leave defaults
#
# PROFILE      insert profiling instructions
# DEBUG        insert debug information and run time checks
# OPTIMIZE     optimize program for speed
# WARNING      generate compiler warnings for unusual constructs
# BOUNDS       generate bounds check during run

CPROF = false

ifeq ($(COMPILER_PROFILE),NORMAL)
  CPROF = true
  PROFILE = false
  DEBUG = true
  OPTIMIZE = MEDIUM
  WARNING = true
  BOUNDS = false
endif

ifeq ($(COMPILER_PROFILE),CHECK)
  CPROF = true
  PROFILE = true
  DEBUG = true
  OPTIMIZE = NONE
  WARNING = true
  BOUNDS = true
endif

ifeq ($(COMPILER_PROFILE),SPEED)
  CPROF = true
  PROFILE = false
  DEBUG = false
  OPTIMIZE = HIGH
  WARNING = false
  BOUNDS = false
endif

ifeq ($(CPROF),false)
  RULES_MAKE_PARAMETERS = RULES_MAKE_PARAMETER_ERROR
  RULES_MAKE_MESSAGE = "COMPILER_PROFILE must be one of NORMAL,CHECK,SPEED"
  $(warning COMPILER_PROFILE=$(COMPILER_PROFILE))
endif

##############################################
#
# GNU compiler (g77, f77 or gfortran)
#
##############################################
#
# -Wall			warnings
# -pedantic		a lot of warnings
# -fautomatic		forget local variables
# -static		save local variables
# -no-automatic		save local variables
# -O			optimization
# -g			debug code
# -r8			double precision for old compiler
# -fdefault-real-8	double precision for new compiler
# -p			use profiling
#
# problems with stacksize (segmentation fault)
#	ulimit -a
#	ulimit -s 32000		(32MB)
#	ulimit -s unlimited
#
##############################################

# determines major version for compilers

GMV := $(shell $(FEMBIN)/cmv.sh -quiet gfortran)
IMV := $(shell $(FEMBIN)/cmv.sh -quiet intel)
GMV_LE_4  := $(shell [ $(GMV) -le 4 ] && echo true || echo false )
IMV_LE_14 := $(shell [ $(IMV) -le 14 ] && echo true || echo false )

MVDEBUG := true
MVDEBUG := false
ifeq ($(MVDEBUG),true)
  $(info gfortran major version = $(GMV) )
  $(info gfortran major version <= 4: $(GMV_LE_4) )
  $(info intel major version = $(IMV) )
  $(info intel major version <= 14: $(IMV_LE_14) )
endif

# next solves incompatibility of option -Wtabs between version 4 and higher

WTABS = -Wno-tabs
ifeq ($(GMV_LE_4),true)
  WTABS = -Wtabs
endif
ifeq ($(MVDEBUG),true)
  $(info WTABS = $(WTABS) )
endif

FGNU_GENERAL = -cpp
ifdef MODDIR
  FGNU_GENERAL = -cpp -J$(MODDIR)
endif

FGNU_PROFILE = 
ifeq ($(PROFILE),true)
  FGNU_PROFILE = -p
endif

FGNU_WARNING = 
ifeq ($(WARNING),true)
  FGNU_WARNING = -Wall -pedantic
  FGNU_WARNING = -Wall $(WTABS) -Wno-unused -Wno-uninitialized
  FGNU_WARNING = -Wall $(WTABS) -Wno-unused
  FGNU_WARNING = -Wall $(WTABS) -Wno-unused \
			-Wno-conversion -Wno-unused-dummy-argument \
			-Wno-zerotrip
endif

FGNU_BOUNDS = 
ifeq ($(BOUNDS),true)
  FGNU_BOUNDS = -fbounds-check
  FGNU_BOUNDS = -fcheck=all
endif

FGNU_NOOPT = 
ifeq ($(DEBUG),true)
  TRAP_LIST = zero,invalid,overflow,underflow,denormal
  TRAP_LIST = zero
  TRAP_LIST = zero,invalid,overflow,denormal
  TRAP_LIST = zero,invalid,overflow
  FGNU_NOOPT = -g
  #FGNU_NOOPT = -g -fbacktrace -ffpe-trap=$(TRAP_LIST)
  FGNU_NOOPT = -g -fbacktrace -ffpe-trap=$(TRAP_LIST) $(FGNU_BOUNDS)
endif

FGNU_OPT   = -O
ifeq ($(OPTIMIZE),HIGH)
  FGNU_OPT   = -O3
endif
ifeq ($(OPTIMIZE),NONE)
  FGNU_OPT   = 
endif

FGNU_OMP   =
ifeq ($(PARALLEL_OMP),true)
  FGNU_OMP   =  -fopenmp
endif

#----------------------------------

ifeq ($(FORTRAN_COMPILER),GNU_G77)
  FGNU		= g77
  FGNU95	= g95
  F77		= $(FGNU)
  F95		= $(FGNU95)
  LINKER	= $(F77)
  LFLAGS	= $(FGNU_OPT) $(FGNU_PROFILE) $(FGNU_OMP)
  FFLAGS	= $(LFLAGS) $(FGNU_NOOPT) $(FGNU_WARNING)
  FFLAG_SPECIAL	= $(LFLAGS) $(FGNU_WARNING)
  FINFOFLAGS	= --version
endif

ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  FGNU		= gfortran
  FGNU95	= gfortran
  ifneq ($(PARALLEL_MPI),NONE)
    FGNU        = mpif90
    FGNU95      = mpif90
  endif
  F77		= $(FGNU)
  F95		= $(FGNU95)
  LINKER	= $(F77)
  LFLAGS	= $(FGNU_OPT) $(FGNU_PROFILE) $(FGNU_OMP)
  FFLAGS	= $(LFLAGS) $(FGNU_NOOPT) $(FGNU_WARNING) $(FGNU_GENERAL)
  FFLAG_SPECIAL	= $(LFLAGS) $(FGNU_WARNING) $(FGNU_GENERAL)
  FINFOFLAGS	= --version
endif

##############################################
#
# PGI compiler (nvfortran)
#
##############################################
#
# for download see: https://developer.nvidia.com/nvidia-hpc-sdk-download
#
##############################################

FPGI_GENERAL = 
ifdef MODDIR
  FPGI_GENERAL = -module $(MODDIR)
endif

FPGI_OMP   =
ifeq ($(PARALLEL_OMP),true)
  FPGI_OMP   = -mp
endif

FPGI_BOUNDS = 
ifeq ($(BOUNDS),true)
  FPGI_BOUNDS = -check uninit 
  FPGI_BOUNDS = -Mbounds -Mchkptr -Mchkstk
endif

FPGI_PROFILE = 
ifeq ($(PROFILE),true)
  FPGI_PROFILE = -Mprof
endif

FPGI_NOOPT = -cpp
ifeq ($(DEBUG),true)
  FPGI_NOOPT = -g -traceback -Ktrap=fp -cpp
endif

FPGI_OPT   = -O
ifeq ($(OPTIMIZE),HIGH)
  FPGI_OPT   = -O3
endif
ifeq ($(OPTIMIZE),NONE)
  FPGI_OPT   = 
endif

FGNU_OMP   =
ifeq ($(PARALLEL_OMP),true)
  FGNU_OMP   =  -fopenmp
endif

FPGI_WARNING =

ifeq ($(FORTRAN_COMPILER),PGI)
  FPGI		= nvfortran
  F77		= $(FPGI)
  F95		= nvfortran
  LINKER	= $(FPGI)
  LFLAGS	= $(FPGI_OPT) $(FPGI_PROFILE) $(FPGI_OMP) $(FPGI_BOUNDS)
  FFLAGS	= $(LFLAGS) $(FPGI_NOOPT) $(FPGI_WARNING) $(FPGI_GENERAL)
  FFLAG_SPECIAL	= $(FFLAGS)
  FINFOFLAGS	= --version
endif
 
##############################################
#
# IBM compiler (xlf)
#
##############################################
#
# if you use xlf95  "-qnosave" is a default option
# xlf_r is thread safe
# all the compiler options are included in FIBM_OMP
# set PARALLEL_OMP = TRUE
##############################################

FIBM_PROFILE = 
ifeq ($(PROFILE),true)
  FIBM_PROFILE = 
endif

FIBM_WARNING = 
ifeq ($(WARNING),true)
  FIBM_NOOPT = 
endif

FIBM_NOOPT = 
ifeq ($(DEBUG),true)
  FIBM_NOOPT =
endif

FIBM_OPT   = -O
ifeq ($(OPTIMIZE),HIGH)
  FIBM_OPT   = -O3
endif
ifeq ($(OPTIMIZE),NONE)
  FIBM_OPT   = 
endif

FIBM_OMP   =
ifeq ($(PARALLEL_OMP),true)
     FIBM_OMP    = -qsmp=omp -qnosave -q64 -qmaxmem=-1 -NS32648 -qextname -qsource -qcache=auto -qstrict -O3 -qarch=pwr6 -qtune=pwr6
endif

#----------------------------------

ifeq ($(FORTRAN_COMPILER),IBM)
  FIBM		= xlf_r
  F77		= $(FIBM)
  F95		= xlf_r
  LINKER	= $(FIBM)
  FFLAGS	= $(FIBM_OMP)
  FFLAG_SPECIAL	= $(FFLAGS)
  LFLAGS	= $(FIBM_OMP) -qmixed  -b64 -bbigtoc -bnoquiet -lpmapi -lessl -lmass -lmassvp4
endif
 
##############################################
#
# Portland compiler
#
##############################################

FPG_PROFILE = 
ifeq ($(PROFILE),true)
  FPG_PROFILE = -Mprof=func
endif

FPG_WARNING = 
ifeq ($(WARNING),true)
  FPG_WARNING =
endif

FPG_NOOPT = -cpp
ifeq ($(DEBUG),true)
  FPG_NOOPT = -g -cpp
endif

FPG_OPT   = -O
ifeq ($(OPTIMIZE),HIGH)
  FPG_OPT   = -O3
endif
ifeq ($(OPTIMIZE),NONE)
  FPG_OPT   = 
endif

FPG_OMP   =
ifeq ($(PARALLEL_OMP),true)
  FPG_OMP   =  -mp
endif

#----------------------------------

ifeq ($(FORTRAN_COMPILER),PORTLAND)
  FPG		= pgf90
  FPG95		= pgf90
  F77		= $(FPG)
  F95		= $(FPG95)
  LINKER	= $(F77)
  LFLAGS	= $(FPG_OPT) $(FPG_PROFILE) $(FPG_OMP)
  FFLAGS	= $(LFLAGS) $(FPG_NOOPT) $(FPG_WARNING)
  FFLAG_SPECIAL	= $(FFLAGS)
  FINFOFLAGS	= -v
endif

##############################################
#
# INTEL compiler
#
##############################################
#
# -w    warnings	no warnings
# -CU   run time exception
# -dn   level n of diagnostics
# -O    optimization (O2 O3)
# -i8   integer 8 byte (-i2 -i4 -i8)
# -r8   real 8 byte    (-r4 -r8 -r16), also -autodouble
#
# -check none
# -check all
# -check bounds		(run time exception for out of bounds arrays)
# -check uninit		(run time exception for uninitialized values)
# -check pointer	(run time exception for zero pointers)
#
# -implicitnone
# -debug
# -openmp
# -openmp-profile
# -openmp-report	(1 2)
#
# -warn interfaces -gen-interfaces
#
# problems with stacksize (segmentation fault)
#	export KMP_STACKSIZE=32M
#
#############################################

# ERSEM FLAGS -------------------------------------

REAL_4B = real\(4\)
DEFINES += -DREAL_4B=$(REAL_4B)
DEFINES += -DFORTRAN95 
DEFINES += -DPRODUCTION -static
FINTEL_ERSEM = $(DEFINES) 

#-------------------------------------------------

FINTEL_GENERAL = -fpp
ifdef MODDIR
  FINTEL_GENERAL = -fpp -module $(MODDIR)
endif

FINTEL_PROFILE = 
ifeq ($(PROFILE),true)
  FINTEL_PROFILE = -p
endif

FINTEL_WARNING = 
ifeq ($(WARNING),true)
  FINTEL_WARNING =
  FINTEL_WARNING = -w
  FINTEL_WARNING = -warn interfaces,nouncalled -gen-interfaces
endif

FINTEL_BOUNDS = 
ifeq ($(BOUNDS),true)
  FINTEL_BOUNDS = -check uninit -check bounds -check pointer
endif

FINTEL_NOOPT = -g -traceback
ifeq ($(DEBUG),true)
  FINTEL_TRAP = -fp-trap-all=common
  FINTEL_TRAP = -ftrapuv -debug all -fpe0
  FINTEL_NOOPT = -xP
  FINTEL_NOOPT = -CU -d1
  FINTEL_NOOPT = -CU -d5
  FINTEL_NOOPT = -g -traceback -O0
  FINTEL_NOOPT = -g -traceback
  FINTEL_NOOPT = -g -traceback $(FINTEL_BOUNDS) $(FINTEL_TRAP)
endif

# FINTEL_OPT   = -O -g -Mprof=time
# FINTEL_OPT   = -O3 -g -axSSE4.2 #-mcmodel=medium -shared-intel
# FINTEL_OPT   = -O3 -g -axAVX -mcmodel=medium -shared-intel
# FINTEL_OPT   = -O -g -fp-model precise -no-prec-div

FINTEL_OPT   = -O -mcmodel=large
FINTEL_OPT   = -O 
ifeq ($(OPTIMIZE),HIGH)
  FINTEL_OPT   = -O3
  FINTEL_OPT   = -O3 -xhost
  FINTEL_OPT   = -O2 -xhost
  #FINTEL_OPT   = -O3 -mcmodel=medium
  #FINTEL_OPT   = -O3 -mcmodel=large
endif
ifeq ($(OPTIMIZE),NONE)
  FINTEL_OPT   = 
endif

FINTEL_OMP   =
ifeq ($(PARALLEL_OMP),true)
  FINTEL_OMP   = -threads -qopenmp
  FINTEL_OMP   = -qopenmp
  ifeq ($(IMV_LE_14),true)
    FINTEL_OMP   = -openmp
  endif
  ifeq ($(MVDEBUG),true)
    $(info FINTEL_OMP = $(FINTEL_OMP) )
  endif
endif

ifeq ($(FORTRAN_COMPILER),INTEL)
  FINTEL	= ifort
  ifneq ($(PARALLEL_MPI),NONE)
    FINTEL      = mpiifort
  endif
  F77		= $(FINTEL)
  F95     	= $(F77)
  LINKER	= $(F77)
  LFLAGS	= $(FINTEL_OPT) $(FINTEL_PROFILE) $(FINTEL_OMP)
  FFLAGS	= $(LFLAGS) $(FINTEL_NOOPT) $(FINTEL_WARNING) $(FINTEL_GENERAL)
  FFLAG_SPECIAL	= $(FFLAGS)
  FINFOFLAGS	= -v
endif

##############################################
#
# C compiler
#
##############################################

ifeq ($(C_COMPILER),GNU_GCC)
  CC     = gcc
  CFLAGS = -O -Wall -pedantic
  CFLAGS = -O -Wall -pedantic -std=gnu99  #no warnings for c++ style comments
  LCFLAGS = -O 
  CINFOFLAGS = --version
endif

ifeq ($(C_COMPILER),PGI)
  CC     = nvc
  CFLAGS = -O -Wall -pedantic
  CFLAGS = -O -Wall -pedantic -std=gnu99  #no warnings for c++ style comments
  CFLAGS = -O -Wall
  LCFLAGS = -O 
  CINFOFLAGS = --version
endif

ifeq ($(C_COMPILER),INTEL)
  CC     = icc
  CFLAGS = -O -g -traceback -check-uninit
  CFLAGS = -O -g -traceback
  LCFLAGS = -O 
  CINFOFLAGS = -v
endif

ifeq ($(C_COMPILER),IBM)
  CC     = xlc
#  CFLAGS = -O -traceback -check-uninit
  CFLAGS = -O
  LCFLAGS = -O 
  CINFOFLAGS = -v
endif

##############################################
#
# old stuff - do not bother
#
##############################################

##############################################
#
# profiling: use -p for linking ; example for ht
#
#   f77 -p -c subany.f ... 	(compiles with profiling)
#   f77 -p -o ht ...     	(links with profiling)
#   ht                   	(creates mon.out)
#   prof ht              	(elaborates mon.out)
#   gprof ht              	(elaborates mon.out)
#	-b	brief
#	-p	flat profile
#
##############################################


##############################################
#
# lahey
#
#	help for options :
#		f77l3
#		386link
#		up l32 /help
#	remove kernel :
#		os386 /remove
#	compiler switches :
#		/R	remember local variables
#		/H	hardcopy listing
#		/S	create sold file .sld
#		/L	line-number traceback table
#		/O	output options
#
# compiler (all versions)
# 
# FLC             = f77l3
# FLFLAGS         = /nO /nS /L /nR /nH
# 
# linker lahey 3.01
# 
# LLINKER         = up l32
# LDFLAGS         = /nocommonwarn /stack:5004
# LDMAP           = nul
# LLIBS		=
# LDEXTRA		= ,$@,$(LDMAP),$(LLIBS) $(LDFLAGS);
# EXE		= exp
# 
# linker lahey 5.20
# 
# LLINKER		= 386link
# LDFLAGS		= -stack 2000000 -fullsym
# LDEXTRA		= $(LDFLAGS)
# EXE    		= exe
#
##############################################

##############################################
# Private stuff
##############################################

#MAKEGENERAL = $(FEMDIR)/general/Makefile.general

