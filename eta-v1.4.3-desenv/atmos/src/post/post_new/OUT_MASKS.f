      SUBROUTINE OUT_MASKS(EGFUL,GDOUT,IMOT,JMOT)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    OUT_MASKS   
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-23
C     
C ABSTRACT: 
C    
C   .     
C     
C PROGRAM HISTORY LOG:
C   ??-??-??  DAVID PLUMMER - SUBROUTINE INTERP
C   92-12-23  RUSS TREADON - BROKE INTERP INTO SEVERAL PIECES
C                            ONE OF WHICH BECAME THIS
C                            INTERPOLATION ROUTINE.
C   95-05-03  MIKE BALDWIN - ADDED BIT MAP ARRAY
C     
C USAGE:    CALL INTERP3(EGFUL,GDOUT,IMOT,JMOT)
C   INPUT ARGUMENT LIST:
C     EGFUL    - DATA ON FILLED E-GRID.
C     IMOT     - FIRST DIMENSION OF OUTPUT GRID.
C     JMOT     - SECOND DIMENSION OF OUTPUT GRID.
C
C   OUTPUT ARGUMENT LIST: 
C     GDOUT    - DATA INTERPOLATED TO OUTPUT GRID.
C     
C   OUTPUT FILES:
C     STDOUT  - RUN TIME STANDARD OUT.
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       EXTEND   - FILLS MISSING VALUES ON OUTPUT GRID WITH 
C                  EITHER AN EXTENSION OF BORDER VALUES OR
C                  THE FIELD MEAN.
C     LIBRARY:
C       COMMON   - OPTIONS
C                  LLGRDS
C                  IOUNIT
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY Y-MP
C$$$  
C     
C     
C     INCLUDE DECLARED GRID DIMENSIONS.
C-----------------------------------------------------------------
      INCLUDE "parmeta"
      INCLUDE "parmout"
C-----------------------------------------------------------------
      PARAMETER (IMT=2*IM-1,JMT=JM)
C     
C     DECLARE VARIABLES.
C-----------------------------------------------------------------
      INTEGER IMOT,JMOT
      REAL EGFUL(IMT,JMT),GDOUT(IMOT,JMOT)
      DOUBLE PRECISION SUM
C-----------------------------------------------------------------
C     
C     INCLUDE COMMON.
      INCLUDE "OPTIONS.comm"
      INCLUDE "LLGRDS.comm"
      INCLUDE "BITMAP.comm"
      INCLUDE "IOUNIT.comm"
C
C     SET TOLERANCE LIMITS.
C     
      DATA MXPASS,SPVC,SMALL /2,1.E20,1.E-4/
C     
C
      GDOUT(1:IMOT,1:JMOT)=SPVC
C
!$omp  parallel do
!$omp& private(m,n)
      doout120: DO J = 1,JMOT
      doin120: DO  I = 1,IMOT
        IF (IWGT(I,J).EQ.1) THEN
          M = IEGRD(I,J)
          N = JEGRD(I,J) 
          GDOUT(I,J) = EGFUL(M,N) 
        ENDIF
       END DO doin120
       END DO doout120
!$omp end parallel do
C
C
      RETURN
      END
