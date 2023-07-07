################################################################################
# 
#     Makefile for Eta Coordinate Post
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         TARGET   - name of the executable
#         FC       - name of Fortran compiler
#         CPP      - name of CPP
#         ARCH     - architecture
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         MEM      - user data area and stack size
#         MAP      - load map
#         W3LIB    - w3lib
#         BACIO    - bacio lib
#         ESSL     - ESSL library
#         MASS     - MASS library
#         SEARCH   - library search location
#
#################################################################################
#
# Define the name of the executable
#
include ../../configure/make.inc_pos
TARGET = ../../../exe/etapost_new.x
#
# CPP, Compiler, and Linker Options
#
# FC       = f90
CPP      = /lib/cpp -P
CPPFLAGS = 
OPTS     = $(FFLAGS)
W3LIB    = ../../../../libraries/w3lib
#BACIO    = ../../../libraries/bacio
# MPILIB   = ../../../libraries/dummyMPI/libmpi.a
#MPILIB = -lmpi
#
# Assemble Options
#
# FFLAGS   = $(OPTS) 
FFLAGST  = $(OPTS) 
LDFLAGS  =  $(FFLAGS)
#GSMLIBS     = $(W3LIB) $(BACIO) $(MPILIB) 
LIBS     = $(W3LIB) $(MPILIB)
#
#
# Non-threaded object files
#
# OBJS=	CALCAPE.o    MISCLN.o    NGMFLD2.o
#
# Threaded object files
#
OBJST=	BLOSFC2.o    BNDLYR2.o   BOUND.o BOUNDL.o CALCAPE.o  CALDRG.o  \
	CALDWP2.o    CALEKM.o    CALEKM2.o  CALHEL.o    CALLCL.o  \
	CALMCVG.o    CALMXW.o    CALPOT2.o  CALPW.o     CALRCH.o  \
	CALRH2.o     CALSTRM.o   CALTAU.o   CALTHTE.o   CALVIS.o  \
	CALVOR.o     CALWXT.o    CALCWI.o   CALQUV.o    CALHPBL.o CETLIH4.o  CHR2INT.o   CLDRAD.o  \
	CLMAX.o      COEFI.o     COEFJ.o    COLLECT.o   DEWPOINT.o  DIR.o     \
	DIST.o       E2OUT.o     E2TLL.o    EFILL.o     EFILT.o     ETA2P.o   \
	ETAFLD2.o  ETAPOST.o EUVGUV.o EXCH.o EXCH2.o EXTEND.o EXTRACT.o \
	FDLVL.o  FILLH.o  FILLV.o  FILTER.o FIXED.o     FRZLVL.o \
	FRZLVL2.o    GD2EG.o     GD2EGK.o   GD2TLL.o    GENBIL.o  \
	GETGDEF.o \
	GENEUV.o     GENLL.o      GET_BITS.o GRIBIT.o    H2V.o     \
	INITPOST.o   INTERP3.o    LFMFLD.o   LL2EG.o     MINMAX.o  \
	MISCLN.o     MIXLEN.o     MPI_FIRST.o MPI_LAST.o NETAL.o NGMFLD2.o NGMSLP2.o \
	NOPACK.o     OTLIFT2.o    OTLFT2.o   OUTPUT.o  OUT_MASKS.o PARA_RANGE.o P2FILT.o  \
	P2FLTV.o     POSTDATA.o   PROCESS.o  RDRST2D.o  \
	READCNTRL2.o SCLFLD.o     STAT.o    SURFCE2.o  TRNSEG.o    TRNSGD.o  \
	TRPAUS.o     TTBLEX.o     WDIR.o WETBULB.o P2FILTF.o P2FLTVF.o 
#
# Includes
#
INCLUDES= parm.tbl parmeta parmout parmsoil cuparm
#
# Common Blocks
#
COMMS=  ACMCLD.comm BITMAP.comm DYNAMD.comm IOUNIT.comm OMGAOT.comm PHYS2.comm \
	ACMCLH.comm BOCO.comm E2PFLG.comm LLGRDS.comm OPTIONS.comm PVRBLS.comm \
	ACMPRE.comm CLDWTR.comm EGRID.comm LOOPS.comm OUTFIL.comm RQSTFLD.comm \
	ACMRDL.comm CNVCLD.comm EXTRA.comm MAPOT.comm OUTGRD.comm SOIL.comm \
	ACMRDS.comm CTLBLK.comm GRBDAT.comm MAPOT1.comm PHYS.comm SRFDSP.comm \
	ACMSFC.comm DYNAM.comm INDX.comm MASKS.comm PHYS1.comm VRBLS.comm

DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES:	.F .f .o

.F.f:
	$(CPP) $(CPPFLAGS) $< > $*.f

# $(TARGET):	$(OBJS) $(OBJST)
$(TARGET):	$(OBJST)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(LIBS)

# $(OBJST):	$(DEPS)
.f.o:
	$(FC) $(FFLAGST) -c $<

clean:	
	/bin/rm -f  $(TARGET) *.lst *.o *.i
#
