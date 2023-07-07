################################################################################
# 
#     Makefile for Eta Coordinate Model - Dummy MPI Version
#                                          Tuccillo Dec 1999
# 	
#	Updated for internal quilting version of the model
#					   Pyle July 2000
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
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         INCLUDE  - search location for include files
#         W3LIB    - w3lib
#         MPI      - dummy MPI library
#         SEARCH   - library search location
#
#################################################################################

.SUFFIXES: .f90 .o
# INC = ../include_mpi
INC=./
include ../configure/make.inc

#
# Define the name of the executable
#
TARGET = ../../exe/etafcst.x

COMPILE.f = $(FC) -c $(FFLAGS) -I$(INC)
PREPROC.F = $(CPP) -P $(CFLAGS)

CPP      = /lib/cpp -P
CPPFLAGS = 
OPTS= -O1
SMP      = 
TRAPS    =
PROFILE  =
DEBUG=
INCLUDE = ./
W3LIB= ../../../libraries/w3lib
# MPI=../../libraries/dummyMPI/libmpi.a
# SEARCH   = -L../../libraries/dummyMPI
#
# Assemble Options
#
# FFLAGS should be defined in make.inc
# FFLAGS   = $(OPTS) $(LIST) $(TRAPS) $(PROFILE) $(DEBUG) $(INCLUDE)
FFLAGST  = $(FFLAGS) $(SMP)
LDFLAGS  = $(SMP) $(PROFILE) $(FFLAGS)
#LIBS     = $(MPI) $(W3LIB) -lmpi
LIBS     = $(MPI) $(W3LIB)
#
# Non-threaded object files
#
OBJS=	MODULE_F77KINDS.o MODULE_PARMEXCH.o MODULE_CUINIT.o MODULE_TIMMING.o MODULE_TIMCHK.o MODULE_PARMETA.o  MODULE_MPPCOM.o MODULE_GLB_TABLE.o MODULE_TEMPCOM.o MODULE_TOPO.o MODULE_MAPPINGS.o \
        MODULE_EXCH_BUF_INTEGER.o MODULE_EXCH_BUF_REAL.o MODULE_EXCHM.o MODULE_RDPARM.o MODULE_BANDTA.o MODULE_BDWIDE.o MODULE_BDCOMB.o MODULE_RNDDTA.o MODULE_TEMPV.o MODULE_SLOPES.o MODULE_PRFHLD.o \
	MODULE_ACMCLD.o MODULE_ACMPRE.o MODULE_ACMSFC.o MODULE_CNVCLD.o MODULE_PARMSOIL.o MODULE_ACMRDL.o MODULE_ACMRDS.o MODULE_ACMCLH.o  MODULE_BOCO.o MODULE_NHYDRO.o MODULE_CLDWTR.o MODULE_LOOPS.o \
	 MODULE_CONTIN.o MODULE_INDX.o MODULE_PVRBLS.o MODULE_PARM_TBL.o MODULE_CUPARM.o MODULE_MASKS.o MODULE_SOIL.o MODULE_VRBLS.o MODULE_INPUT.o MODULE_CTLBLK.o MODULE_MAPOT.o  MODULE_PHYS.o MODULE_DYNAM.o \
	 MODULE_RD1TIM.o MODULE_SWRSAV.o MODULE_TBLTMP.o MODULE_DIUCON.o MODULE_SSALB.o MODULE_RDFSAV.o MODULE_SAVMEM.o MODULE_SCRTCH.o MODULE_CMICRO_STATS.o MODULE_CMY600.o MODULE_CMICRO_CONS.o MODULE_IACCR_TABLES.o \
	 MODULE_IMASS_TABLES.o MODULE_IVENT_TABLES.o MODULE_IRATE_TABLES.o MODULE_IRIME_TABLES.o MODULE_RACCR_TABLES.o MODULE_SDENS_TABLES.o MODULE_RMASS_TABLES.o MODULE_RRATE_TABLES.o MODULE_RVELR_TABLES.o MODULE_RVENT_TABLES.o \
	 MODULE_COMPVS0.o MODULE_COMPVS.o MODULE_CMICRO_START.o MODULE_ASTSAV.o MODULE_MOMENTO.o MODULE_RITE.o MODULE_ABCI.o MODULE_UPDATE_FLDS.o  MODULE_Z0_DATA.o \
	MODULE_HCON.o MODULE_PHYCON.o BOCOH.o BOCOHF.o BOCOV.o MODULE_PARMBUF.o MODULE_BUFFER.o MODULE_NSOILTYPE.o MODULE_O3DATA.o MODULE_SEASO3.o MODULE_OUTFIL.o CHKOUT.o MODULE_C_TADJ.o CLTEND.o \
	MODULE_Z0EFFT.o MODULE_PARA.o CLO89.o COLLECT.o  MODULE_CO2BD2.o MODULE_CO2BD3.o MODULE_CO2BD4.o MODULE_CO2BD5.o MODULE_CO2DTA.o CONRAD.o CONRAD2.o MODULE_PPTASM.o CUCNVC.o CUCNVC_SHALLOW.o DDAMP.o \
	DIFCOF.o DIGFLT.o DIST.o DSTRB.o MODULE_TABCOM.o E1E290.o E290.o \
	E2SPEC.o E3V88.o EBU.o EPS.o FILT25.o FST88.o \
	GFDLRD.o GOSSIP.o GRADFS.o \
	GSMCOLUMN.o GSMCONST.o MODULE_C_FRACN.o GSMDRIVE.o GSCOND.o \
	HADZ.o HCONST.o \
	IDSTRB.o INIT.o INITS.o LOC2GLB.o LWR88.o \
	MIXLEN.o MPI_FIRST.o MPPINIT.o O3CLIM.o O3INT.o \
	OZON2D.o PARA_RANGE.o PDNEW.o PDTEDT.o \
	PGCOR.o PRECPD.o PRODQ2.o \
	QUILT.o  RADFS.o  RDTEMP.o READ_NHB.o  \
	READ_RESTRT.o READ_RESTRT2.o SETUP_SERVERS.o \
	SFCDIF.o SFLX.o SGETMO.o SLP.o SOLARD.o SPA88.o MODULE_OPTIONS.o SURFCE.o SWR93.o \
	TABLE.o TTBLEX.o UPDATE.o VADZ.o VDIFH.o VDIFQ.o VDIFV.o \
	VTADV.o VTADVF.o ZENITH.o ZERO2.o ZERO3.o ZERO3_T.o TWR.o  \
	SLPSIG.o SLPSIGSPLINE.o VWR.o \
	READ_OUT2RESTRT.o OUT2RESTRT.o \
	
#
# Threaded object files
#
OBJST=  DIVHOA.o DIVHOAST.o DIVHOASTQL.o HDIFF.o HDIFFS.o HZADV.o HZADV_LM1.o HZADVS.o HZADV2.o RADTN.o TURBL.o

#
# Includes
#
INCLUDES=  MPP.h MODULE_PARMETA.f90 EXCHM.h MODULE_PARMBUF.f90 MODULE_PARMEXCH.f90
#
# Common Blocks
#
COMMS=  MODULE_BUFFER.f90 COMM_ACMCLD.f90 COMM_CTLBLK.f90 COMM_MAPOT.f90  COMM_MAPOT1.f90  COMM_PHYS.f90   COMM_SEASO3.f90 \
        COMM_ACMCLH.f90 COMM_BOCO.f90   COMM_DYNAM.f90  COMM_MASKS.f90   COMM_PHYS1.f90  COMM_SOIL.f90 \
        COMM_ACMPRE.f90 COMM_CLDWTR.f90 COMM_DYNAMD.f90 COMM_MSLP.f90    COMM_PHYS2.f90  COMM_TABCOM.f90 \
        COMM_ACMRDL.f90 COMM_CNVCLD.f90 COMM_HCON.f90   COMM_NHYDRO.f90  COMM_PRFHLD.f90 COMM_TEMPV.f90 \
        COMM_ACMRDS.f90 COMM_CO2DTA.f90 COMM_INDX.f90   COMM_OPTIONS.f90 COMM_PVRBLS.f90 COMM_VRBLS.f90 \
        COMM_ACMSFC.f90 COMM_CONTIN.f90 COMM_LOOPS.f90  COMM_OUTFIL.f90  COMM_RNDDTA.f90 COMM_Z0EFFT.f90 \
        COMM_PARA.f90

DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES: .F .f .o

.f90.o:
#	echo "USING CPP STUFF"
#	$(PREPROC.F) $< > $*.f
#	mv $(<:.F=.i)  $(<:.F=.f)
	echo "Compiling $(<:.f90=.f90)"
	$(COMPILE.f) $(<:.f90=.f90)
#	rm $(<:.F=.f)

$(TARGET):	$(OBJS) $(OBJST)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(OBJST) $(LIBS)
#	$(FC) $(LDFLAGS) $(LIBS) -o $@ $(OBJS) $(OBJST)

# $(OBJS):	$(DEPS)
#	$(COMPILE.f) $<

# $(OBJST):	$(DEPS)
#	$(FC) $(FFLAGST) -c $<

clean:	
	/bin/rm -f  *.i *.lst *.o *.mod

exe: $(TARGET)
