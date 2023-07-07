      SUBROUTINE CALPW(PW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    CALPW       COMPUTES 
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-24       
C     
C ABSTRACT:  
C     THIS ROUTINE COMPUTES PRECIPITABLE WATER IN A COLUMN
C     EXTENDING FROM THE FIRST ATMOSPHERIC ETA LAYER TO THE
C     MODEL TOP.  THE DEFINITION USED IS
C                                 TOP
C            PRECIPITABLE WATER = SUM (Q+CLDW) DP*HTM/G
C                                 BOT
C     WHERE,
C        BOT IS THE FIRST ETA LAYER,
C        TOP IS THE MODEL TOP,
C        Q IS THE SPECIFIC HUMIDITY (KG/KG) IN THE LAYER
C        CLDW IS THE CLOUD WATER (KG/KG) IN THE LAYER
C        DP (Pa) IS THE LAYER THICKNESS.
C        HTM IS THE HEIGHT MASK AT THAT LAYER (=0 IF BELOW GROUND)
C        G IS THE GRAVITATIONAL CONSTANT
C     
C PROGRAM HISTORY LOG:
C   92-12-24  RUSS TREADON
C   96-03-04  MIKE BALDWIN - ADD CLOUD WATER AND SPEED UP CODE
C   98-06-15  T BLACK      - CONVERSION FROM 1-D TO 2-D
C   00-01-04  JIM TUCCILLO - MPI VERSION                 
C     
C USAGE:    CALL CALPW(PW)
C   INPUT ARGUMENT LIST:
C     PW       - ARRAY OF PRECIPITABLE WATER.
C
C   OUTPUT ARGUMENT LIST: 
C     NONE     
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       NONE
C     LIBRARY:
C       COMMON   - LOOPS
C                  EXTRA
C                  VRBLS
C                  MASKS
C                  CLDWTR
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C     
C     INCLUDE/SET PARAMETERS.
      INCLUDE "parmeta"
      INCLUDE "params"
C     
C     SET DENSITY OF WATER AT 1 ATMOSPHERE PRESSURE, 0C.
C     UNITS ARE KG/M**3.
      PARAMETER (RHOWAT=1.E3)
C     
C     DECLARE VARIABLES.
C     
      INTEGER LLMH
      REAL ALPM,DZ,PM,PWSUM,RHOAIR
      REAL PW(IM,JM)
C     
C     INCLUDE COMMON BLOCKS.
      INCLUDE "VRBLS.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "LOOPS.comm"
      INCLUDE "MASKS.comm"
      INCLUDE "CLDWTR.comm"
      INCLUDE "CTLBLK.comm"
C
C***************************************************************
C     START CALPW HERE.
C
C     INITIALIZE PW TO 0.    
C     
      PW = 0.
C     
C     OUTER LOOP OVER VERTICAL DIMENSION.
C     INNER LOOP OVER HORIZONTAL GRID.
C     
!$omp  parallel do
!$omp& private(dp)
      DO L = 1,LM
        DO J=JSTA,JEND
        DO I=1,IM
          DP   =PINT(I,J,L+1)-PINT(I,J,L)
          PW(I,J)=PW(I,J)+(Q(I,J,L)+CWM(I,J,L))*DP*GI*HTM(I,J,L)
        ENDDO
        ENDDO
      ENDDO
C     
C        AT ONE TIME THE SUM WAS DIVIDED BY THE DENSITY OF
C        WATER.  (TO GET PW IN M) THIS IS NO LONGER DONE.
C        PW(I,J) = PWSUM/RHOWAT
C     
C     END OF ROUTINE.
C     
      RETURN
      END
