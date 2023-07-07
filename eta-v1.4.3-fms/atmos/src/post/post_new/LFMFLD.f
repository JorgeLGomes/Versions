      SUBROUTINE LFMFLD(RH3310,RH6610,RH3366,PW3310)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    LFMFLD      COMPUTES LAYER MEAN LFM FIELDS
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22
C     
C ABSTRACT:
C     THIS ROUTINE COMPUTES THREE LAYER MEAN RELATIVE HUMIDITIES
C     AND A PRECIPITABLE WATER FIELD FROM ETA LEVEL DATA.  THE
C     COMPUTED FIELDS ARE INTENDED TO MIMIC SIMILAR FIELDS COM-
C     PUTED BY THE LFM.  THE ALGORITHM USED HERE IS FAIRLY PRI-
C     MATIVE.  IN EACH COLUMN ABOVE A MASS POINT ON THE ETA GRID
C     WE SET THE FOLLOWING TARGET PRESSURES:
C         SIGMA LAYER 1.00 PRESSURE:  SURFACE PRESSURE
C         SIGMA LAYER 0.66 PRESSURE:  0.50 * SURFACE PRESSURE
C         SIGMA LAYER 0.33 PRESSURE:  0.4356 * SURFACE PRESSURE
C     GIVEN THESE PRESSURES A SURFACE UP SUMMATION IS MADE OF 
C     RELATIVE HUMIDITY AND/OR PRECIPITABLE WATER BETWEEN THESE
C     TARGET PRESSURES.  EACH TERM IN THE SUMMATION IS WEIGHTED
C     BY THE THICKNESS OF THE ETA LAYER.  THE FINAL LAYER MEAN
C     IS THIS SUM NORMALIZED BY THE TOTAL DEPTH OF THE LAYER.  
C     THERE IS, OBVIOUSLY, NO NORMALIZATION FOR PRECIPITABLE WATER.
C
C     
C PROGRAM HISTORY LOG:
C   92-12-22  RUSS TREADON
C   93-07-27  RUSS TREADON - MODIFIED SUMMATION LIMITS FROM
C                            0.66*PSFC TO 0.75*PSFC AND 0.33*PSFC 
C                            TO 0.50*PSFC, WHERE PSFC IS THE
C                            SURFACES PRESSURE.  THE REASON FOR
C                            THIS CHANGE WAS RECOGNITION THAT IN 
C                            THE LFM 0.33 AND 0.66 WERE MEASURED
C                            FROM THE SURFACE TO THE TROPOPAUSE,
C                            NOT THE TOP OF THE MODEL.
C   93-09-13  RUSS TREADON - RH CALCULATIONS WERE MADE INTERNAL
C                            TO THE ROUTINE.
C   96-03-04  MIKE BALDWIN - CHANGE PW CALC TO INCLUDE CLD WTR 
C   98-06-16  T BLACK      - CONVERSION FROM 1-D TO 2-D
C   98-08-17  MIKE BALDWIN - COMPUTE RH OVER ICE
C   98-12-22  MIKE BALDWIN - BACK OUT RH OVER ICE
C   00-01-04  JIM TUCCILLO - MPI VERSION
C     
C     
C USAGE:    CALL LFMFLD(RH3310,RH6610,RH3366,PW3310)
C   INPUT ARGUMENT LIST:
C     NONE
C
C   OUTPUT ARGUMENT LIST: 
C     RH3310   - SIGMA LAYER 0.33-1.00 MEAN RELATIVE HUMIDITY.
C     RH6610   - SIGMA LAYER 0.66-1.00 MEAN RELATIVE HUMIDITY.
C     RH3366   - SIGMA LAYER 0.33-0.66 MEAN RELATIVE HUMIDITY.
C     PW3310   - SIGMA LAYER 0.33-1.00 PRECIPITABLE WATER.
C     
C   OUTPUT FILES:
C     NONE
C     
C   LIBRARY:
C     COMMON   - VRBLS
C                MAPOT
C                EXTRA
C                LOOPS
C                OPTIONS
C
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C
C     
C     INCLUDE PARAMETERS.
      INCLUDE "parmeta"
      INCLUDE "params"
C
      PARAMETER (RHOWAT=1.E3)
C     
C     DECLARE VARIABLES.
C     
      REAL ALPM, DZ, ES, PM, PWSUM, QM, QS, TM
      REAL RH3310(IM,JM), RH6610(IM,JM), RH3366(IM,JM)
      REAL PW3310(IM,JM), IW(IM,JM,LM)
C     
C     INCLUDE COMMON BLOCKS.
      INCLUDE "VRBLS.comm"
      INCLUDE "MAPOT.comm"
      INCLUDE "CLDWTR.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "LOOPS.comm"
      INCLUDE "OPTIONS.comm"
      INCLUDE "CTLBLK.comm"
C
C***********************************************************************
C     START LFMFLD HERE
C     
C
C    COMPUTE IW
C
        CLIMIT =1.0E-20
        IW=0.
C
        DO L=2,LM 
        DO J=JSTA,JEND
        DO I=1,IM
           IF(CWM(I,J,L).GT.CLIMIT) THEN
             IF(T(I,J,L).LT.258.15)THEN
               IW(I,J,L)=1.
             ELSEIF(T(I,J,L).GE.273.15)THEN
               IW(I,J,L)=0.
             ELSE
               IF(IW(I,J,L-1).EQ.1.0)IW(I,J,L)=1.
             ENDIF
           ELSE
             IW(I,J,L)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
C
C     LOOP OVER HORIZONTAL GRID.
C     
      doout30: DO J=JSTA,JEND
      doin30: DO I=1,IM
C     
C        ZERO VARIABLES.
         RH3310(I,J) = D00
         PW3310(I,J) = D00
         RH6610(I,J) = D00
         RH3366(I,J) = D00
         Z3310     = D00
         Z6610     = D00
         Z3366     = D00
C     
C        SET BOUNDS FOR PRESSURES AND SURFACE L.
         P10  = PD(I,J) + PT
         P66  = 0.75*P10
         P33  = 0.50*P10
         LLMH = LMH(I,J)
C     
C        ACCULMULATE RELATIVE HUMIDITIES AND PRECIPITABLE WATER.
C
         do10: DO L = LLMH,1,-1
C     
C           GET P, Z, T, AND Q AT MIDPOINT OF ETA LAYER.
            ALPM = D50*(ALPINT(I,J,L)+ALPINT(I,J,L+1))
            DZ   = ZINT(I,J,L)-ZINT(I,J,L+1)
            DP   = PINT(I,J,L+1)-PINT(I,J,L)
            PM   = EXP(ALPM)
            TM   = T(I,J,L)
            QM   = Q(I,J,L)
            QM   = AMAX1(QM,D00)
C
            TMT0=TM-273.16
            TMT15=AMIN1(TMT0,-15.)
            AI=0.008855
            BI=1.
            IF(TMT0.LT.-20.)THEN
              AI=0.007225
              BI=0.9674
            ENDIF
            QW=PQ0/PM*EXP(A2*(TM-A3)/(TM-A4))
            QI=QW*(BI+AI*AMIN1(TMT0,0.))
            QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
            IF(TMT0.LT.-15.)THEN
               QS=QI
            ELSEIF(TMT0.GE.0.)THEN
               QS=QINT
            ELSE
               IF(IW(I,J,L).GT.0.0) THEN
                 QS=QI
               ELSE
                 QS=QINT
               ENDIF
            ENDIF
CMEB 12/22/98 SWITCH TO RH VS WATER NO MATTER WHAT
C             DELETE THIS LINE TO SWITCH BACK TO RH VS ICE
            QS=QW
CMEB 12/22/98 SWITCH TO RH VS WATER NO MATTER WHAT

C
            RH   = QM/QS
            IF (RH.GT.H1) THEN
               RH = H1
               QM = RH*QS
            ENDIF
            IF (RH.LT.D01) THEN
               RH = D01
               QM = RH*QS
            ENDIF
C
C           JUMP OUT OF THIS LOOP IF WE ARE ABOVE THE HIGHEST TARGET PRESSURE.
            IF (PM.LE.P33) EXIT do10 
C     
C           0.66-1.00 RELATIVE HUMIDITY.
            IF ((PM.LE.P10).AND.(PM.GE.P66)) THEN
               Z6610     = Z6610 + DZ
               RH6610(I,J) = RH6610(I,J) + RH*DZ
            ENDIF
C     
C           0.33-1.00 RELATIVE HUMIDITY AND PRECIPITABLE WATER.
            IF ((PM.LE.P10).AND.(PM.GE.P33)) THEN
               Z3310      = Z3310 + DZ
               RH3310(I,J)= RH3310(I,J)+RH*DZ
               PW3310(I,J)= PW3310(I,J)+(Q(I,J,L)+CWM(I,J,L))*DP*GI
            ENDIF
C     
C           0.33-0.66 RELATIVE HUMIDITY.
            IF ((PM.LE.P66).AND.(PM.GE.P33)) THEN
               Z3366     = Z3366 + DZ
               RH3366(I,J) = RH3366(I,J) + RH*DZ
            ENDIF
C
         END DO do10
C     
C        NORMALIZE TO GET MEAN RELATIVE HUMIDITIES.  AT
C        ONE TIME WE DIVIDED PRECIPITABLE WATER BY DENSITY
C        TO GET THE EQUIVALENT WATER DEPTH IN METERS.  NO MORE.
         IF (Z6610.GT.D00) THEN
            RH6610(I,J) = RH6610(I,J)/Z6610
         ELSE
            RH6610(I,J) = SPVAL
         ENDIF
C     
         IF (Z3310.GT.D00) THEN
            RH3310(I,J) = RH3310(I,J)/Z3310
         ELSE
            RH3310(I,J) = SPVAL
         ENDIF
C     
         IF (Z3366.GT.D00) THEN
            RH3366(I,J) = RH3366(I,J)/Z3366
         ELSE
            RH3366(I,J) = SPVAL
         ENDIF
C
CWAS     PW3310(I,J) = PW3310(I,J)/RHOWAT
C     
      END DO doin30
      END DO doout30
C     
C     
C     END OF ROUTINE.
C     
      RETURN
      END
