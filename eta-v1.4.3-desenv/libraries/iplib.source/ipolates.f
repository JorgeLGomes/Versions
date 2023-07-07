C-----------------------------------------------------------------------
      SUBROUTINE IPOLATES(IP,IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                    NO,RLAT,RLON,IBO,LO,GO,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  IPOLATES   IREDELL'S POLATE FOR SCALAR FIELDS
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
C
C ABSTRACT: THIS SUBPROGRAM INTERPOLATES SCALAR FIELDS
C           FROM ANY GRID TO ANY GRID (JOE IRWIN'S DREAM).
C           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
C           THE FOLLOWING INTERPOLATION METHODS ARE POSSIBLE:
C             (IP=0) BILINEAR
C             (IP=1) BICUBIC
C             (IP=2) NEIGHBOR
C             (IP=3) BUDGET
C             (IP=4) SPECTRAL
C             (IP=6) NEIGHBOR-BUDGET
C           SOME OF THESE METHODS HAVE INTERPOLATION OPTIONS AND/OR
C           RESTRICTIONS ON THE INPUT OR OUTPUT GRIDS, BOTH OF WHICH
C           ARE DOCUMENTED MORE FULLY IN THEIR RESPECTIVE SUBPROGRAMS.
C           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
C           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
C           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
C             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
C             (KGDS(1)=001) MERCATOR CYLINDRICAL
C             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
C             (KGDS(1)=004) GAUSSIAN CYLINDRICAL
C             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
C             (KGDS(1)=201) ROTATED EQUIDISTANT CYLINDRICAL
C             (KGDS(1)=202) ROTATED EQUIDISTANT CYLINDRICAL
C           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
C           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
C           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
C           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
C           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
C           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
C           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
C           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
C           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
C           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
C        
C PROGRAM HISTORY LOG:
C   96-04-10  IREDELL
C
C USAGE:    CALL IPOLATES(IP,IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
C    &                    NO,RLAT,RLON,IBO,LO,GO,IRET)
C
C   INPUT ARGUMENT LIST:
C     IP       - INTEGER INTERPOLATION METHOD
C                (IP=0 FOR BILINEAR;
C                 IP=1 FOR BICUBIC;
C                 IP=2 FOR NEIGHBOR;
C                 IP=3 FOR BUDGET;
C                 IP=4 FOR SPECTRAL;
C                 IP=6 FOR NEIGHBOR-BUDGET)
C     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
C                (IP=0: (NO OPTIONS)
C                 IP=1: CONSTRAINT OPTION
C                 IP=2: (NO OPTIONS)
C                 IP=3: NUMBER IN RADIUS, RADIUS WEIGHTS ...
C                 IP=4: SPECTRAL SHAPE, SPECTRAL TRUNCATION
C                 IP=6: NUMBER IN RADIUS, RADIUS WEIGHTS ...)
C     KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
C     KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
C     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
C                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
C     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
C                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
C     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
C     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
C     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF RESPECTIVE IBI(K)=1)
C     GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
C     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
C     RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
C     RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
C
C   OUTPUT ARGUMENT LIST:
C     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)>=0)
C     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)>=0)
C     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)>=0)
C     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
C     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
C     GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
C     IRET     - INTEGER RETURN CODE
C                0    SUCCESSFUL INTERPOLATION
C                1    UNRECOGNIZED INTERPOLATION METHOD
C                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
C                3    UNRECOGNIZED OUTPUT GRID
C                1X   INVALID BICUBIC METHOD PARAMETERS
C                3X   INVALID BUDGET METHOD PARAMETERS
C                4X   INVALID SPECTRAL METHOD PARAMETERS
C
C SUBPROGRAMS CALLED:
C   POLATES0     INTERPOLATE SCALAR FIELDS (BILINEAR)
C   POLATES1     INTERPOLATE SCALAR FIELDS (BICUBIC)
C   POLATES2     INTERPOLATE SCALAR FIELDS (NEIGHBOR)
C   POLATES3     INTERPOLATE SCALAR FIELDS (BUDGET)
C   POLATES4     INTERPOLATE SCALAR FIELDS (SPECTRAL)
C
C REMARKS: EXAMPLES DEMONSTRATING RELATIVE CPU COSTS.
C   THIS EXAMPLE IS INTERPOLATING 12 LEVELS OF TEMPERATURES
C   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
C   TO THE 93 X 68 HAWAIIAN MERCATOR GRID (NCEP GRID 204).
C   THE EXAMPLE TIMES ARE FOR THE C90.  AS A REFERENCE, THE CP TIME
C   FOR UNPACKING THE GLOBAL 12 TEMPERATURE FIELDS IS 0.04 SECONDS.
C
C   METHOD      IP  IPOPT          CP SECONDS
C   --------    --  -------------  ----------
C   BILINEAR    0                   0.03
C   BICUBIC     1   0               0.07
C   BICUBIC     1   1               0.07
C   NEIGHBOR    2                   0.01
C   BUDGET      3   -1,-1           0.48
C   SPECTRAL    4   0,40            0.22
C   SPECTRAL    4   1,40            0.24
C   SPECTRAL    4   0,-1            0.42
C   N-BUDGET    6   -1,-1           0.15
C
C   THE SPECTRAL INTERPOLATION IS FAST FOR THE MERCATOR GRID.
C   HOWEVER, FOR SOME GRIDS THE SPECTRAL INTERPOLATION IS SLOW.
C   THE FOLLOWING EXAMPLE IS INTERPOLATING 12 LEVELS OF TEMPERATURES
C   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
C   TO THE 93 X 65 CONUS LAMBERT CONFORMAL GRID (NCEP GRID 211).
C
C   METHOD      IP  IPOPT          CP SECONDS
C   --------    --  -------------  ----------
C   BILINEAR    0                   0.03
C   BICUBIC     1   0               0.07
C   BICUBIC     1   1               0.07
C   NEIGHBOR    2                   0.01
C   BUDGET      3   -1,-1           0.51
C   SPECTRAL    4   0,40            3.94
C   SPECTRAL    4   1,40            5.02
C   SPECTRAL    4   0,-1           11.36
C   N-BUDGET    6   -1,-1           0.18
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      INTEGER IPOPT(20)
      INTEGER KGDSI(200),KGDSO(200)
      INTEGER IBI(KM),IBO(KM)
	LOGICAL*1 LI(MI,KM),LO(MO,KM)
      REAL GI(MI,KM),GO(MO,KM)
      REAL RLAT(MO),RLON(MO)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  BILINEAR INTERPOLATION
      IF(IP.EQ.0) THEN
        CALL POLATES0(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                NO,RLAT,RLON,IBO,LO,GO,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  BICUBIC INTERPOLATION
      ELSEIF(IP.EQ.1) THEN
        CALL POLATES1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                NO,RLAT,RLON,IBO,LO,GO,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  NEIGHBOR INTERPOLATION
      ELSEIF(IP.EQ.2) THEN
        CALL POLATES2(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                NO,RLAT,RLON,IBO,LO,GO,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  BUDGET INTERPOLATION
      ELSEIF(IP.EQ.3) THEN
        CALL POLATES3(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                NO,RLAT,RLON,IBO,LO,GO,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SPECTRAL INTERPOLATION
      ELSEIF(IP.EQ.4) THEN
Cmp        CALL POLATES4(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
Cmp     &                NO,RLAT,RLON,IBO,LO,GO,IRET)
	write(6,*) 'spectral option removed...try something else!'
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  NEIGHBOR-BUDGET INTERPOLATION
      ELSEIF(IP.EQ.6) THEN
        CALL POLATES6(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                NO,RLAT,RLON,IBO,LO,GO,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  UNRECOGNIZED INTERPOLATION METHOD
      ELSE
        IF(KGDSO(1).GE.0) NO=0
CMIC$ DO ALL AUTOSCOPE
        DO K=1,KM
          IBO(K)=1
          DO N=1,NO
            LO(N,K)=.FALSE.
            GO(N,K)=0.
          ENDDO
        ENDDO
        IRET=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
