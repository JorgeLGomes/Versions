      SUBROUTINE GRIBIT(IFLD,ILVL,GRID,IMOUT,JMOUT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    GRIBIT      POST FIELDS IN GRIB1
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-06-18       
C     
C ABSTRACT:
C     THIS ROUTINE POSTS THE DATA IN THE PASSED ARRAY GRID
C     TO THE OUTPUT FILE IN GRIB1 FORMAT.
C     
C PROGRAM HISTORY LOG:
C   93-06-18  RUSS TREADON
C   93-11-23  RUSS TREADON - REMOVED CODE GENERATING GRIB INDEX FILE.
C   98-07-17  MIKE BALDWIN - REMOVED LABL84, NOW USING ID
C     
C USAGE:    CALL GRIBIT(IFLD,ILVL,GRID,IMOUT,JMOUT)
C   INPUT ARGUMENT LIST:
C     IFLD     - FIELD ID TAG.
C     ILVL     - INTEGER TAG FOR LEVEL OF FIELD.
C     GRID     - FIELD TO BE POSTED IN GRIB.
C     IMOUT    - FIRST DIMENSION OF OUTPUT GRID.
C     JMOUT    - SECOND DIMENSION OF OUTPUT GRID.
C
C   OUTPUT ARGUMENT LIST: 
C     
C   OUTPUT FILES:
C     STDOUT    - RUN TIME STANDARD OUT.
C     LUNOUT+1  - UNIT TO RECEIVE GRIB1 DATA.
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C     GETENV   - CRAY SUBROUTINE TO GET VALUE OF ENVIRONMENT VARIABLE.
C     MINMAX   - DETERMINES MIN/MAX VALUES IN AN ARRAY.
C     WRYTE    - WRITE DATA OUT BY BYTES.
C     GET_BITS   - COMPUTE NUMBER OF BITS 
C     VARIOUS W3LIB ROUTINES
C     LIBRARY:
C       COMMON   - CTLBLK
C                  MAPOT
C                  RQSTFLD
C                  IOUNIT
C                  OUTGRD
C                  GRBDAT
C                  AVBLFLDS
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C
C
C     INCLUDE GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
C
      INCLUDE "parmeta"
      INCLUDE "parmout"
      PARAMETER (LP1=LM+1,D01=0.01,D50=0.5E0)
      PARAMETER (IMT=2*IM-1,JMT=JM,IMJMT=IMT*JMT)
C     
C     GRIB1 PARAMETERS.
C        MNBIT  = MINIMUM NUMBER OF BITS TO USE IN PACKING.
C        MXBIT  = MAXIMUM NUMBER OF BITS TO USE IN PACKING.
C        LENPDS = LENGTH OF GRIB1 PDS.
C        LENGDS = LENGTH OF GRIB1 GDS.
C     
      PARAMETER (MNBIT=0,MXBIT=16,LENPDS=28,LENGDS=32)
C     
C     DECLARE VARIABLES.
C     
      LOGICAL RUN,FIRST,RESRT,SIGMA,OLDRD,STRD
      LOGICAL NORTH
      CHARACTER*1  KBUF(30+LENPDS+LENGDS+IMOUT*JMOUT*(MXBIT+2)/8)
      CHARACTER*1  IFLAG
      CHARACTER*4  RESTHR,BLANK
      CHARACTER*6  CRUN,PROJ
      CHARACTER*7  DESCR2
      CHARACTER*28 PDS
Cmp      CHARACTER*50 ENVAR
      CHARACTER*90 ENVAR
      CHARACTER*110 FNAME,OPATH
      CHARACTER*90 CMD
      INTEGER IBDSFL(9)
Cmp      INTEGER IGRD(IMOUT,JMOUT),IGDS(18),IBMASK(IMOUT,JMOUT)
      INTEGER IGRD(IMOUT,JMOUT),IBMASK(IMOUT,JMOUT)
Cmp
      REAL GRID(IMOUT,JMOUT),GRIDOT(IMOUT,JMOUT)
C     
C     THE BELOW VARIABLE ARE ONLY NEEDED FOR THE CALL TO W3FI63.
      REAL DATAFLD(IMOUT,JMOUT)
      INTEGER KGDS(20),KPTR(16)
      LOGICAL KBMS(IMOUT,JMOUT)
C     
C     INCLUDE COMMON BLOCKS.
      INCLUDE "CTLBLK.comm"
      INCLUDE "MAPOT.comm"
      INCLUDE "RQSTFLD.comm"
      INCLUDE "IOUNIT.comm"
      INCLUDE "OUTGRD.comm"
      INCLUDE "BITMAP.comm"
      INCLUDE "GRBDAT.comm"
      INCLUDE "OUTFIL.comm"
Cmp
        INCLUDE "GDS.com"
Cmp
C     
C     SET DEFAULT GRIB1 PARAMETERS.  
C     PARAMETERS MNBIT, MXBIT, IBX, AND NBIT ARE USED 
C     IN THE CALL TO GET_BITS.
C        IBX    = DESIRED BINARY PRECISION.
C        NBIT   = NUMBER OF BITS TO USE IN PACKING DATA.
C     
      DATA IBX,NBIT / 0, 12 /
      DATA BLANK /'    '/
      SAVE OPATH
C
C*****************************************************************************
C     START GRIBIT HERE.
C
C     SET NUMBER OF OUTPUT GRID POINTS.
      IJOUT = IMOUT*JMOUT
C     
C     PREPARE GRIB PDS
C     
C     SET ARRAY ID VALUES TO GENERATE GRIB1 PDS.  
C        ID(1)  = NUMBER OF BYTES IN PRODUCT DEFINITION SECTION (PDS)
C        ID(2)  = PARAMETER TABLE VERSION NUMBER
C        ID(3)  = IDENTIFICATION OF ORIGINATING CENTER
C        ID(4)  = MODEL IDENTIFICATION (ALLOCATED BY ORIGINATING CENTER)
C        ID(5)  = GRID IDENTIFICATION
C        ID(6)  = 0 IF NO GDS SECTION, 1 IF GDS SECTION IS INCLUDED
C        ID(7)  = 0 IF NO BMS SECTION, 1 IF BMS SECTION IS INCLUDED
C        ID(8)  = INDICATOR OF PARAMETER AND UNITS (TABLE 2)
C        ID(9)  = INDICATOR OF TYPE OF LEVEL       (TABLE 3)
C        ID(10) = VALUE 1 OF LEVEL (=0 FOR 1-100,102,103,105,107,
C          109,111,113,115,117,119,125,160,200,201 LEVEL IS IN ID WORD 11)
C        ID(11) = VALUE 2 OF LEVEL
C        ID(12) = YEAR OF CENTURY
C        ID(13) = MONTH OF YEAR
C        ID(14) = DAY OF MONTH
C        ID(15) = HOUR OF DAY
C        ID(16) = MINUTE OF HOUR   (IN MOST CASES SET TO 0)
C        ID(17) = FCST TIME UNIT
C        ID(18) = P1 PERIOD OF TIME
C        ID(19) = P2 PERIOD OF TIME
C        ID(20) = TIME RANGE INDICATOR
C        ID(21) = NUMBER INCLUDED IN AVERAGE
C        ID(22) = NUMBER MISSING FROM AVERAGES
C        ID(23) = CENTURY  (20, CHANGE TO 21 ON JAN. 1, 2001)
C        ID(24) = RESERVED - SET TO 0
C        ID(25) = SCALING POWER OF 10
C
      IF (IOUTYP.EQ.3.OR.IOUTYP.EQ.5) THEN
C     
C        PREPARE DATE PART OF GRIB PDS RECORD.
C         IFHR       = NTSD/TSPH+D50
         IFHR       = ITAG
	 ICENT      = (IDAT(3)-1)/100 + 1
         IYY        = IDAT(3) - (ICENT-1)*100
         IMM        = IDAT(1)
         IDD        = IDAT(2)
         AYEAR0     = IYY
         AMNTH0     = IMM
         ADAY0      = IDD
         AGMT0      = IHRST
         ID(01)     = 28
         ID(02)     = 2
         ID(03)     = 7
         ID(12)     = IYY
         ID(13)     = IMM
         ID(14)     = IDD
         ID(15)     = IHRST
         ID(16)     = 0

C	ID(17)=1 (hourly time increment.  limits to 256 hours)
C	ID(17)=10 (3 hour increment.  limits to 768 hours)
C	
C	below options are listed, but currently would have problems
C	beyond 999 hours due to filename constraints.  
C
C	ID(17)=11 (6 hour increment.  limits to 1536 hours)
C	ID(17)=12 (12 hour increment.  allows 3072 hours)
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         ID(17)     = 1
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	IF (ID(17) .eq. 1) IFACT=1
	IF (ID(17) .eq. 10) IFACT=3
	IF (ID(17) .eq. 11) IFACT=6
	IF (ID(17) .eq. 12) IFACT=12
C
C    ASSUMING ID(18-20), (P1, P2, TIME RANGE INDICATOR) 
C    ARE PASSED IN CORRECTLY IF NOT AN INSTANTANEOUS FIELD
C   
         IF (ID(20).EQ.0) THEN
          ID(18)     = IFHR/IFACT
          ID(19)     = 0
         ENDIF

!tst
	IF (ID(20).EQ.3 .OR. ID(20).EQ.4) THEN
	  ID(18)     = ID(18)/IFACT
          ID(19)     = ID(19)/IFACT
	ENDIF
!endtst

         ID(21)     = 0
         ID(22)     = 0
         ID(23)     = ICENT
         ID(24)     = 0
C
C     
C        SET OUTPUT GRID TYPE.  WE ASSUME KGYTPE HOLDS THE GRIB
C        ID FOR THE OUTPUT GRID.  
C
         KGTYP = KGTYPE
C     
C        SET GRID TYPE ID(5)
C        GENERATING PROGRAM ID(4)
C
         IJOUT      = IMOUT*JMOUT
         ID(4) = IMDLTY
         ID(5) = KGTYP
C
C        ID(6) =0 IF NO GDS SECTION, =1 IF GDS INCLUDED, 
C                 ALWAYS INCLUDE GDS
C
         ID(6) = 1
C     
C        SET DATA TYPE ID(8) AND SURFACE ID(9).
C
C     DON'T SET PARAMETER IF PRECIP TYPE, SINCE THERE ARE
C     4 PARAMETER NUMBERS FOR THE SAME IFLD
C
         IF (ID(8).LT.140.OR.ID(8).GT.143) ID(8) = IQ(IDENT(IFLD))
         IF (ID(9).EQ.0) ID(9) = IS(IDENT(IFLD))
C     
C        SET VALUE OF LEVEL IF ON PRESSURE OR ETA SURFACE.
C        OTHERWISE, WE ASSUME ID(10) AND (11) ARE SET 
C        APPROPRIATELY PRIOR TO ENTERING GRIBIT.
C     
         IF (ID(9).EQ.100)  THEN
            ISVALUE = NINT(SPL(ILVL)*D01)
            ID(10) = 0
            ID(11) = ISVALUE
         ELSEIF (ID(9).EQ.119) THEN
            ISVALUE = NINT(AETA(ILVL)*10000.)
C
C   TKE IS ON THE ETA INTERFACE AT THE BOTTOM OF THE LAYER ILVL
C
            IF (ID(8).EQ.158) ISVALUE = NINT(ETA(ILVL+1)*10000.)
            ID(10) = 0
            ID(11) = ISVALUE
         ELSEIF (ID(9) .EQ. 109) THEN
            ISVALUE = ILVL
            ID(10) = 0
            ID(11) = ISVALUE
         ENDIF
C     
C     END OF GRIB PDS LABEL PREPARATION.
C
      ENDIF

C     
C     SET DECIMAL SCALING (IDECI) FROM LIST IN INCLUDE FILE 
C     RQSTFLD.  A CALL TO GET_BITS WILL COMPUTE THE NUMBER OF
C     BITS NECESSARY TO PACK THE DATA BASED ON THE RANGE OF 
C     THE FIELD.  THE FIELD IS SCALED TO THIS PRECISION AND
C     RETURNED FOR PACKING BY THE GRIB PACKER.
C     
      DO JJ = 1,JMOUT
       DO II = 1,IMOUT
        IBMASK(II,JJ)=IBMAP(II,JJ)
       ENDDO
      ENDDO
      IBM = 0
      IBITM = 0
      SGDG  = DEC(IFLD)
!$omp  parallel do
      DO J=1,JMOUT
      DO I=1,IMOUT
        GRIDOT(I,J)=GRID(I,J)
      ENDDO
      ENDDO
C
      DO J=1,JMOUT
      DO I=1,IMOUT
        IBITM=IBITM+IBMASK(I,J)
      ENDDO
      ENDDO
C
C        ID(7) =0 IF NO BMS SECTION, =1 IF BMS INCLUDED
C
      IF (IBITM.EQ.IJOUT) THEN
        ID(7) = 0
        IBM = 0
      ELSE
        ID(7) = 1
        IBM = 1
      ENDIF
      CALL GET_BITS(IBM,SGDG,IJOUT,IBMASK,GRID,
     &                IDECI,GRIDOT,GMIN,GMAX,NBIT)
C
C        ID(25) = SCALING POWER OF 10
C
      ID(25) = IDECI
C     
C     GENERATE COMPLETE GRIB1 MESSAGE USING W3FI72.
C        ITYPE  = 0 SPECIFIES REAL DATA TO BE PACKED.
C        IGRD   = DUMMY ARRAY FOR INTEGER DATA.
C        IBITL  = NBIT TELLS W3FI72 TO PACK DATA USING NBIT BITS.
C        IPFLAG = 0 IS PDS INFORMATION IN USER ARRAY ID.
C                 1 IS PDS (GENERATED ABOVE BY W3FP12).
C        ID     = (DUMMY) ARRAY FOR USER DEFINED PDS.
C        IGFLAG = 0 TELLS W3FI72 TO MAKE GDS USING IGRID.
C                 1 IS GDS GENERATED BY USER IN ARRAY IGDS
C        IGRID  = GRIB1 GRID TYPE (TABLE B OF ON388).
C        IGDS   = ARRAY FOR USER DEFINED GDS.
C        ICOMP  = 0 FOR EARTH ORIENTED WINDS,
C                 1 FOR GRID ORIENTED WINDS.
C        IBFLAG = 0 TELLS W3FI72 TO MAKE BIT MAP FROM USER
C                 SUPPLIED DATA.
C        IBMASK = ARRAY CONTAINING USER DEFINED BIT MAP.
C        IBLEN  = LENGTH OF ARRAY IBMASK.
C        IBDSFL = ARRAY CONTAINING TABLE 11 (ON388) FLAG INFORMATION.
C        NPTS   = LENGTH OF ARRAY GRID OR IGRD.  MUST AGREE WITH IBLEN.
C     
C     INTIALIZE VARIABLES.
      ITYPE  = 0
!$omp  parallel do
      DO J=1,JMOUT
      DO I=1,IMOUT
        IGRD(I,J)=0
      ENDDO
      ENDDO
C
      IBITL  = MIN(NBIT,MXBIT)
C
      IPFLAG = 0
C
      IGFLAG = 0
Cwrkst      IGRID  = ID(5)
      IF (IGRID.EQ.26) IGRID=6
Cwrkst      DO 20 K = 1,18
Cwrkst        IGDS(K) = 0
Cwrkst 20   CONTINUE
      ICOMP  = 1
      IF (INDEX(PROJ,'LOLA').NE.0) ICOMP = 0
      IBFLAG = 0
      IBLEN  = IJOUT
      DO 30 K = 1,9
         IBDSFL(K) = 0
 30   CONTINUE

Cmp     this is where things need to be defined
Cmp
Cmp     want IGRID=255 (user defined type)
Cmp     IGDS needs to have the w3fi71 style GDS (18 elements)
Cmp     also NEED IGFLAG=1
Cmp
Cmp     what to do with bitmap IBFLAG/IBLEN?
        IGFLAG=1
        IGRID=255
Cmp

C
!	write(6,*) 'ID= ', ID
      CALL W3FI72(ITYPE,GRIDOT,IGRD,IBITL,
     X            IPFLAG,ID,PDS,
     X            IGFLAG,IGRID,IGDS,ICOMP,
     X            IBFLAG,IBMASK,IBLEN,
     X            IBDSFL,
     X            NPTS,KBUF,ITOT,IER)
C     
C     EXPLICITLY SET BYTE 12 OF KBUF (BYTE 4 OF THE PDS)
C     TO 2.  THIS WILL REFER ALL QUANTITIES TO PARAMETER
C     TABLE VERSION 2 OF WHICH TABLE VERSION 1 IS A SUBSET.
C     THIS IS NEEDED BECAUSE THE W3 ROUTINES HARDWIRE THIS
C     VALUE TO 1 YET SOME OF THE OUTPUT VARIABLES ARE ONLY 
C     DEFINED IN VERSION 2 OF THE PARAMETER TABLE.
C
      KBUF(12)=CHAR(2)
C
      IF (IER.NE.0) THEN
         WRITE(STDOUT,1040) IER,FIELD(IFLD)
 1040    FORMAT('GRIBIT:  ***W3FI72 ERROR IER=',I1,
     X        ' FOR ',A20)
         WRITE(STDOUT,*)'GRIBIT:  DID NOT POST THIS FIELD'
         RETURN
      ENDIF
C     
C     ON FIRST ENTRY MAKE OUTPUT DIRECTORY.  SET SWITCH (RITEHD)
C     TO FALSE FOR SUBSEQUENT ENTRIES.
      IF ( ((IOUTYP.EQ.3).AND.RITEHD) .OR.
     X     ((IOUTYP.EQ.5).AND.RITEHD) .OR.
     X     ((IOUTYP.EQ.4).AND.RITE2 ) ) THEN
C
C        PUT FORECAST HOUR INTO DIR PREFIX FOR GRIB FILE.
C         IHR = NTSD/TSPH + 0.5
          IHR = ITAG
C     
C        GET FULL PATH FOR OUTPUT FILE FROM ENVIRONMENT VARIABLE
C        COMSP WHICH IS SET IN THE SCRIPT RUNNING THE MODEL.
C     
C        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
         ENVAR = ' '
         RESTHR = ' '
         CALL GETENV('COMSP',ENVAR)
         CALL GETENV('tmmark',RESTHR)
         KDAT = INDEX(DATSET,' ') -1
         IF (KDAT.LE.0) KDAT = LEN(DATSET)
         KENV = INDEX(ENVAR,' ') -1
         IF (KENV.LE.0) KENV = LEN(ENVAR)
         KTHR = INDEX(RESTHR,' ') -1
         IF (KTHR.LE.0) KTHR = LEN(RESTHR)
       IF (IOUTYP.EQ.5) THEN
         WRITE(DESCR2,1010) IHR
 1010    FORMAT('f',I2.2)
         IF (ENVAR(1:4).EQ.BLANK.AND.RESTHR(1:4).EQ.BLANK) THEN
          OPATH = DATSET(1:KDAT) // '/' // DESCR2(1:3) // '/'
         ELSEIF (ENVAR(1:4).NE.BLANK.AND.RESTHR(1:4).EQ.BLANK) THEN
          OPATH = ENVAR(1:KENV) // DATSET(1:KDAT) // '/' 
     &              // DESCR2(1:3) // '/'
         ELSEIF (ENVAR(1:4).EQ.BLANK.AND.RESTHR(1:4).NE.BLANK) THEN
          OPATH = DATSET(1:KDAT) // '/' // DESCR2(1:3) // '.' //
     &              RESTHR(1:KTHR) // '/'
         ELSE
          OPATH = ENVAR(1:KENV) // DATSET(1:KDAT) // '/' 
     &              // DESCR2(1:3) // '.' // RESTHR(1:KTHR) // '/'
         ENDIF
C
         WRITE(STDOUT,*)'GRIBIT:  DIRECTORY ',OPATH,
     X        ' CREATED FOR GRIB DATA '
       ELSE
C     
C        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
         IF (ENVAR(1:4).EQ.BLANK.AND.RESTHR(1:4).EQ.BLANK) THEN
          WRITE(DESCR2,1011) IHR
 1011     FORMAT('.GrbF',I6.6)
          FNAME = DATSET(1:KDAT) // DESCR2
         ELSEIF(ENVAR(1:4).EQ.BLANK.AND.RESTHR(1:4).NE.BLANK) THEN
          WRITE(DESCR2,1012) IHR
          FNAME = DATSET(1:KDAT) // DESCR2(1:2)  //'.'// RESTHR
         ELSE
          WRITE(DESCR2,1012) IHR
Cmp 1012     FORMAT(I2.2)
 1012     FORMAT(I6.6)
          FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // DESCR2(1:6)
     &              //'.'// RESTHR
         ENDIF
C
C        ASSIGN AND OPEN UNIT FOR GRIB DATA FILE.
         CLOSE(LUNOUT+1)
C        CALL ASNUNIT(LUNOUT+1,'-s unblocked',IER)
         CALL BAOPEN(LUNOUT+1,FNAME,IER)
         IF (IER.NE.0) WRITE(STDOUT,*)
     X        'GRIBIT:  BAOPEN ERROR FOR GRIB DATA ',
     X        'FILE.  IER=',IER
         WRITE(STDOUT,*)'GRIBIT:  OPENED ',LUNOUT+1,
     X        ' FOR GRIB DATA ',FNAME
       ENDIF
C     
C        SET OPEN-UNIT FLAGS TO FALSE.
         RITEHD = .FALSE.
         RITE2  = .FALSE.
      ENDIF
C     
C     WRITE GRIB1 MESSAGE TO OUTPUT FILE.
      CALL WRYTE(LUNOUT+1,ITOT,KBUF)
C     
C     WRITE DIAGNOSTIC MESSAGE.
C        ID(8)  = INDICATOR OF PARAMETER AND UNITS (TABLE 2)
C        ID(9)  = INDICATOR OF TYPE OF LEVEL       (TABLE 3)
C        ID(10) = VALUE 1 OF LEVEL  (0 FOR 1-100,102,103,105,107
C              111,160   LEVEL IS IN ID WORD 11)
C        ID(11) = VALUE 2 OF LEVEL
      WRITE(STDOUT,1050) ID(8),FIELD(IFLD),ID(9),ID(10),ID(11)
 1050 FORMAT('GRIBIT:  ',I3,1X,A20,1X,I3,1X,I5,1X,I5)
C     
C     END OF ROUTINE.
C     
      RETURN
      END
