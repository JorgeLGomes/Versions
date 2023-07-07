    SUBROUTINE CONRAD2(NFILE)
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE CONRAD2 
!>
!> SUBROUTINE: CONRAD2 - READ AND UPDATE CO2 TRANSMISSION DATA FROM UNIT(NFILE) FOR NEW VERTICAL 
!>                       COORDINATE TESTS THESE ARRAYS USED TO BE IN BLOCK DATA.     
!> PROGRAMMER: K.CAMPANA
!> ORG: ?????
!> DATE: 90-03-??
!>
!> ABSTRACT:
!> CO2 DATA TABLES FOR USER'S VERTICAL COORDINATE.
!>
!> THE FOLLOWING COMMON BLOCKS CONTAIN PRETABULATED CO2 TRANSMISSION FUNCTIONS, EVALUATED USING THE
!> METHODS OF FELS AND SCHWARZKOPF (1981) AND SCHWARZKOPF AND FELS (1985).
!> THE 2-DIMENSIONAL ARRAYS ARE CO2 TRANSMISSION FUNCTIONS AND THEIR DERIVATIVES FROM 109-LEVEL 
!> LINE-BY-LINE CALCULATIONS MADE USING THE 1982 MCCLATCHY TAPE (12511 LINES), CONSOLIDATED, 
!> INTERPOLATED TO THE NMC MRF VERTICAL COORDINATTE, AND RE-CONSOLIDATED TO A 200 CM-1 BANDWIDTH. 
!> THE INTERPOLATION METHOD IS DESCRIBED IN SCHWARZKOPF AND FELS (J.G.R., 1985).
!>
!> THE 1-DIM ARRAYS ARE CO2 TRANSMISSION FUNCTIONS AND THEIR DERIVATIVES FOR TAU(I,I+1), I=1, L, 
!> WHERE THE VALUES ARE NOT OBTAINED BY QUADRATURE,BUT ARE THE ACTUAL TRANSMISSIVITIES, ETC, 
!> BETWEEN A PAIR OF PRESSURES.
!> THESE USED ONLY FOR NEARBY LAYER CALCULATIONS INCLUDING QH2O.
!>
!> THE WEIGHTING FUNCTION GTEMP = P(K) ** 0.2 * (1. + P(K) / 30000.) ** 0.8 / 1013250.,
!> WHERE P(K) = PRESSURE, NMC MRF(NEW) L18 DATA LEVELS FOR PSTAR = 1013250.
!>
!> STEMP IS US STANDARD ATMOSPHERES, 1976, AT DATA PRESSURE LEVELS USING NMC MRF SIGMAS, 
!> WHERE PSTAR = 1013.25 MB (PTZ PROGRAM)
!>
!> MODULE CO2BD3 CONTAINS CO2 TRANSMISSION FUNCTIONS AND TEMPERATURE AND PRESSURE DERIVATIVES FOR 
!> THE 560-800 CM-1 BAND. 
!> ALSO INCLUDED ARE THE STANDARD TEMPERATURES AND THE WEIGHTING FUNCTION. THESE DATA ARE IN BLOCK
!> DATA BD3:
!>
!> CO251  - TRANSMISSION FCTNS FOR T0 (STD. PROFILE) WITH P(SFC)=1013.25 MB
!> CO258  - TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)WITH P(SFC)= 810 MB
!> CDT51  - FIRST TEMPERATURE DERIVATIVE OF CO251
!> CDT58  - FIRST TEMPERATURE DERIVATIVE OF CO258
!> C2D51  - SECOND TEMPERATURE DERIVATIVE OF CO251
!> C2D58  - SECOND TEMPERATURE DERIVATIVE OF CO251
!> CO2M51 - TRANSMISSION FCTNS FOR T0 FOR ADJACENT PRESSURE LEVELS, WITH NO PRESSURE QUADRATURE. 
!>          USED FOR NEARBY LAYER COMPUTATIONS. P(SFC)=1013.25 MB
!> CO2M58 - SAME AS CO2M51,WITH P(SFC)= 810 MB
!> CDTM51 - FIRST TEMPERATURE DERIVATIVE OF CO2M51
!> CDTM58 - FIRST TEMPERATURE DERIVATIVE OF CO2M58
!> C2DM51 - SECOND TEMPERATURE DERIVATIVE OF CO2M51
!> C2DM58 - SECOND TEMPERATURE DERIVATIVE OF CO2M58
!> STEMP  - STANDARD TEMPERATURES FOR MODEL PRESSURE LEVEL STRUCTURE WITH P(SFC)=1013.25 MB
!> GTEMP  - WEIGHTING FUNCTION FOR MODEL PRESSURE LEVEL STRUCTURE WITH P(SFC)=1013.25 MB.
!>
!> THE FOLLOWING ARE STILL IN BLOCK DATA
!>
!> B0 - TEMP. COEFFICIENT USED FOR CO2 TRANS. FCTN. CORRECTION FOR T(K). (SEE REF. 4 AND BD3)
!> B1 - TEMP. COEFFICIENT, USED ALONG WITH B0
!> B2 - TEMP. COEFFICIENT, USED ALONG WITH B0
!> B3 - TEMP. COEFFICIENT, USED ALONG WITH B0
!>
!> PROGRAM HISTORY LOG:
!> 00-01-20  K.CAMPANA - ORIGINATOR
!> 18-01-15  LUCCI     - MODERNIZATION OF THE CODE, INCLUDING:
!>                       * F77 TO F90/F95
!>                       * INDENTATION & UNIFORMIZATION CODE
!>                       * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                       * DOCUMENTATION WITH DOXYGEN
!>                       * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NFILE - 
!>
!> OUTPUT ARGUMENT LIST:
!> NONE 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE 
!>
!> USE MODULES: CTLBLK
!>              CO2DTA
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : GRADFS
!>
!> CALLS      : MPI_BCAST
!>              TABLE
!>-------------------------------------------------------------------------------------------------- 
    USE CTLBLK
    USE CO2DTA
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
    USE UPDATE_FLDS
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IX   = 2 * IM - 1
    INTEGER(KIND=I4KIND), PARAMETER :: KX   = LM
    INTEGER(KIND=I4KIND), PARAMETER :: KP   = KX + 1  
    INTEGER(KIND=I4KIND), PARAMETER :: LP12 = LP1 * LP1
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND), DIMENSION(LP1, 2)                                                     ::&
    & SGTMP
!
    REAL   (KIND=R4KIND), DIMENSION(L, 6)                                                       ::&
    & CO21D
!
    REAL   (KIND=R4KIND), DIMENSION(LP1, LP1, 6)                                                ::&
    & CO22D
!
    REAL   (KIND=R4KIND), DIMENSION(LP1, 6)                                                     ::&
    & CO21D3  , CO21D7
!
    REAL   (KIND=R4KIND), DIMENSION(LP12)                                                       ::&
    & DATA2   
!
    INTEGER(KIND=I4KIND), DIMENSION(3)                                                          ::&
    & RSZE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & RSIZE   , KK      , I       , IRTN    , N       , I1      , I2      , K       , J
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & NFILE
!
    INTEGER (KIND=I4KIND)                                                                       ::&
    & IYR     , IMO     , IDY     , IUTC    , YR      , MON     , DAY     , UTC
!
    INTEGER (KIND=I4KIND)                                                                       ::&
    & NREC1   , NREC2   , KKREC1  , KKREC2
!---------------------------------------------------
! THE ABOVE NOT USED IN CURRENT VERSION OF RADIATION
!
! BEGIN HERE TO GET CONSTANTS FOR RADIATION PACKAGE
!--------------------------------------------------- 
    IYR   = IDAT(3)
    IMO   = IDAT(1)
    IDY   = IDAT(2)
    IUTC = NTSD*DT/3600

    CALL GETDATE(IYR, IMO, IDY, IUTC, YR, MON, DAY, UTC)

	IF(MYPE.EQ.0) then
	print*,"CONRAD2", YR, MON, DAY, UTC, NTSD
	ENDIF
	
	IF (NTSD.eq.1) THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 for ntsd=1'
    	  NREC1=0
	  NREC2=0    
        ENDIF

	IF (yr-iyr.eq.3.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
          IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 3yr'
	  NREC1=2
	  NREC2=6
        ENDIF	  	
	
	IF (yr-iyr.eq.6.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 6yr'	
	  NREC1=4
	  NREC2=12	
	ENDIF
	
	IF (yr-iyr.eq.9.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 9yr'	
	  NREC1=6
	  NREC2=18	
	ENDIF

	IF (yr-iyr.eq.12.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 12yr'	
	  NREC1=8
	  NREC2=24	
	ENDIF

	IF (yr-iyr.eq.15.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 15yr'
	  NREC1=10
	  NREC2=30	
	ENDIF

	IF (yr-iyr.eq.18.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 18yr'	
	  NREC1=12
	  NREC2=36	
	ENDIF

	IF (yr-iyr.eq.21.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 21yr'	
	  NREC1=14
	  NREC2=42	
	ENDIF

	IF (yr-iyr.eq.24.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 24yr'	
	  NREC1=16
	  NREC2=48	
	ENDIF

	IF (yr-iyr.eq.27.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 27yr'	
	  NREC1=18
	  NREC2=54	
	ENDIF

	IF (yr-iyr.eq.30.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 30yr'	
	  NREC1=20
	  NREC2=60	
	ENDIF

	IF (yr-iyr.eq.33.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 33yr'
	  NREC1=22
	  NREC2=66	
	ENDIF

	IF (yr-iyr.eq.36.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 36yr'	
	  NREC1=24
	  NREC2=72	
	ENDIF

	IF (yr-iyr.eq.39.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN	
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 39yr'	
	  NREC1=26
	  NREC2=78	
	ENDIF

	IF (yr-iyr.eq.42.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN	
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 42yr'	
	  NREC1=26
	  NREC2=84	
	ENDIF
	
	IF (yr-iyr.eq.45.and.mon.eq.1.and.day.eq.1.and.utc.eq.0)  THEN	
	  IF (MYPE .eq. 0) print*,' Update CO2 file record in CONRAD2 45yr'	
	  NREC1=30
	  NREC2=90	
	ENDIF

!    REWIND NFILE  !Chou 20080207 - no rewind allows to update the CO2 fctns.
!
!---------------------------------------------------------------------- 
! READ IN PRE-COMPUTED CO2 TRANSMISSION DATA AND CONVERT TO CYBER WORDS
!---------------------------------------------------------------------- 
    RSZE(1) = LP1
    RSZE(2) = L
    RSZE(3) = LP1 * LP1
!
!----------------------------------------------------------------------
      OPEN(UNIT=67,FILE='co2_SGTMP_cluster.bin',FORM='unformatted',STATUS='unknown',access='direct',RECL=4*LP1)
     
      RSIZE = RSZE(1)
      DO 10 KK=1,2
      KKREC1=KK+NREC1
      READ(NFILE,REC=KKREC1)(SGTMP(I,KK),I=1,RSIZE)
!      IF(MYPE.EQ.0)READ(NFILE)(SGTMP(I,KK),I=1,RSIZE)
!      CALL MPI_BCAST(SGTMP(1,KK),RSIZE,MPI_REAL,0,
!     1               MPI_COMM_COMP,IRTN)
     
   10 CONTINUE
      CLOSE(67)
!
!----------------------------------------------------------------------
      OPEN(UNIT=67,FILE='co2_CO21D_cluster.bin',FORM='unformatted',STATUS='unknown',access='direct',RECL=4*L)
     
      RSIZE = RSZE(2)
      DO 15 KK=1,6
      KKREC2=KK+NREC2
      READ(NFILE,REC=KKREC2)(CO21D(I,KK),I=1,RSIZE)
!      IF(MYPE.EQ.0)READ(NFILE)(CO21D(I,KK),I=1,RSIZE)
!      CALL MPI_BCAST(CO21D(1,KK),RSIZE,MPI_REAL,0,
!     1               MPI_COMM_COMP,IRTN)
     
   15 CONTINUE
      CLOSE(67)
!
!----------------------------------------------------------------------
      OPEN(UNIT=67,FILE='co2_DATA_cluster.bin',FORM='unformatted',STATUS='unknown',access='direct',RECL=4*LP12)
     
      RSIZE = RSZE(3)
      DO 20 KK=1,6
      KKREC2=KK+NREC2
      READ(NFILE,REC=KKREC2)(DATA2(I),I=1,RSIZE)
!      IF(MYPE.EQ.0)READ(NFILE)(DATA2(I),I=1,RSIZE)
!      CALL MPI_BCAST(DATA2(1),RSIZE,MPI_REAL,0,
!     1               MPI_COMM_COMP,IRTN)

      N=0
      DO 5673 I1=1,LP1
      DO 5673 I2=1,LP1
      N=N+1
      CO22D(I1,I2,KK)=DATA2(N)
 5673 CONTINUE
   20 CONTINUE
      CLOSE(67)
!
!----------------------------------------------------------------------
      OPEN(UNIT=67,FILE='co2_CO21D3_cluster.bin',FORM='unformatted',STATUS='unknown',access='direct',RECL=4*LP1)
     
      RSIZE = RSZE(1)
      DO 25 KK=1,6
      KKREC2=KK+NREC2
      READ(NFILE,REC=KKREC2)(CO21D3(I,KK),I=1,RSIZE)
!      IF(MYPE.EQ.0)READ(NFILE)(CO21D3(I,KK),I=1,RSIZE)
!      CALL MPI_BCAST(CO21D3(1,KK),RSIZE,MPI_REAL,0,
!     1               MPI_COMM_COMP,IRTN)
	 
   25 CONTINUE
      CLOSE(67)
!
!----------------------------------------------------------------------
      OPEN(UNIT=67,FILE='co2_CO21D7_cluster.bin',FORM='unformatted',STATUS='unknown',access='direct',RECL=4*LP1)
     
      DO 30 KK=1,6
      KKREC2=KK+NREC2
      READ(NFILE,REC=KKREC2)(CO21D7(I,KK),I=1,RSIZE)
!      IF(MYPE.EQ.0)READ(NFILE)(CO21D7(I,KK),I=1,RSIZE)
!      CALL MPI_BCAST(CO21D7(1,KK),RSIZE,MPI_REAL,0,
!     1               MPI_COMM_COMP,IRTN)
      
   30 CONTINUE
      CLOSE(67)
!----------------------------------------------------------------------
!
!    REWIND NFILE  !Chou 20080207 - no rewind allows to update the CO2 fctns.
!
    DO 35 K=1,LP1
        STEMP(K) = SGTMP(K,1)
        GTEMP(K) = SGTMP(K,2)
    35  END DO
!
    DO 40 K=1,L
        CDTM51(K) = CO21D(K,1)
        CO2M51(K) = CO21D(K,2)
        C2DM51(K) = CO21D(K,3)
        CDTM58(K) = CO21D(K,4)
        CO2M58(K) = CO21D(K,5)
        C2DM58(K) = CO21D(K,6)
    40  END DO
!
    DO 45 J=1,LP1
        DO 45 I=1,LP1
            CDT51(I,J) = CO22D(I,J,1)
            CO251(I,J) = CO22D(I,J,2)
            C2D51(I,J) = CO22D(I,J,3)
            CDT58(I,J) = CO22D(I,J,4)
            CO258(I,J) = CO22D(I,J,5)
            C2D58(I,J) = CO22D(I,J,6)
    45  END DO
!   
    DO 50 K=1,LP1
        CDT31(K) = CO21D3(K,1)
        CO231(K) = CO21D3(K,2)
        C2D31(K) = CO21D3(K,3)
        CDT38(K) = CO21D3(K,4)
        CO238(K) = CO21D3(K,5)
        C2D38(K) = CO21D3(K,6)
    50  END DO
!
    DO 55 K=1,LP1
        CDT71(K) = CO21D7(K,1)
        CO271(K) = CO21D7(K,2)
        C2D71(K) = CO21D7(K,3)
        CDT78(K) = CO21D7(K,4)
        CO278(K) = CO21D7(K,5)
        C2D78(K) = CO21D7(K,6)
    55  END DO
!
    IF (MYPE == 0) PRINT 66,NFILE
66  FORMAT(1H ,'----READ CO2 TRANSMISSION FUNCTIONS FROM UNIT ',I2)
!-------------------------------
! DEFINE TABLES FOR LW RADIATION
!-------------------------------
    CALL TABLE
!
    RETURN
    END SUBROUTINE CONRAD2
