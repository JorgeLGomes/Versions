      SUBROUTINE CALCWI(CLW,CLI)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C SUBPROGRAM:    CALCWI     COMPUTES vert integr  CLOUD WATER & CLOUD ICE
C SUBPROGRAM:               BASED ON CALQUV and CALCAPE
C
C ABSTRACT:
C     GIVEN SPECIFIC HUMIDITY, Q, and Cloud water/ice AND THE U-V WIND COMPONENTS
C     THIS ROUTINE CALCULATES  integral (qu) dp/g ;ntegral (qv) dp/g
C     WHERE,
C
C     THIS REQUIRES WINDS BE AT VELOCITY POINTS.
C     SPECIFIC HUMIDITY IS INTERPOLATED TO VELOCITY POINTS.
C
C USAGE:    CALL CALCWI(CLW,CLI)
C   INPUT ARGUMENT LIST:
C
C   OUTPUT ARGUMENT LIST:
C     CLW     -  INTEGRAL OF CLOUD WATER
C     CLI     -  INTEGRAL OF CLOUD ICE
C
C   OUTPUT FILES:
C     NONE
C
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       NONE
C     LIBRARY:
C       COMMON   - MASKS
C                  DYNAM
C                  OPTIONS
C                  INDX
C
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN 90
C     MACHINE : CRAY C-90
C$$$
C
C
C
C     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
C
      INCLUDE "parmeta"
      INCLUDE "params"
C
C     DECLARE VARIABLES.
C
      REAL CLW(IM,JM), CLI(IM,JM), QI(IM,JM)
      REAL IW(IM,JM,LM)
C
C     DECLARE COMMONS.
      INCLUDE "DYNAM.comm"
      INCLUDE "MASKS.comm"
      INCLUDE "OPTIONS.comm"
      INCLUDE "INDX.comm"
      INCLUDE "CTLBLK.comm"
      INCLUDE "CLDWTR.comm"
      INCLUDE "VRBLS.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "LOOPS.comm"
C
C
C***************************************************************************
C     START CALCWI  HERE.
C
C     DEFINE IW switch between water and ice
C
        CLIMIT =1.0E-20

        DO J=JSTA,JEND
        DO I=1,IM
          IW(I,J,1)=0.
        ENDDO
        ENDDO
C
      DO 125 L=2,LM
      DO J=JSTA,JEND
      DO I=1,IM
        LML=LM-LMH(I,J)
        HH=HTM(I,J,L)*HBM2(I,J)
        TKL=T(I,J,L)
        QKL=Q(I,J,L)
        CWMKL=CWM(I,J,L)
        TMT0=(TKL-273.15)*HH
        TMT15=AMIN1(TMT0,-15.)*HH
        PP=PDSL(I,J)*AETA(L)+PT
        QW=HH*PQ0/PP*EXP(HH*A2*(TKL-A3)/(TKL-A4))
        QI(I,J)=QW*(1.+0.01*AMIN1(TMT0,0.))
C
        U00KL=U00(I,J)+UL(L+LML)*(0.95-U00(I,J))*UTIM
        IF(TMT0.LT.-15.0)THEN
          FIQ=QKL-U00KL*QI(I,J)
          IF(FIQ.GT.D00.OR.CWMKL.GT.CLIMIT) THEN
            IW(I,J,L)=1.
          ELSE
            IW(I,J,L)=0.
          ENDIF
        ENDIF
        IF(TMT0.GE.0.0)IW(I,J,L)=0.
        IF(TMT0.LT.0.0.AND.TMT0.GE.-15.0)THEN
          IW(I,J,L)=0.
          IF(IW(I,J,L-1).EQ.1.0.AND.CWMKL.GT.CLIMIT)IW(I,J,L)=1.
        ENDIF
      ENDDO
      ENDDO
  125 CONTINUE

!$omp  parallel do
      DO J=JSTA,JEND
      DO I=1,IM
        CLW(I,J) = D00
        CLI(I,J) = D00
      ENDDO
      ENDDO

      DO L=1,LM

!$omp  parallel do
      DO J=JSTA_M,JEND_M
      DO I=2,IM-1
      DP     = -(PINT(I,J,L) -PINT(I,J,L+1))
      CLW(I,J)= CLW(I,J)+ (CWM(I,J,L)*(1-IW(I,J,L))*HTM(I,J,L))*DP*GI
      CLI(I,J)= CLI(I,J)+ (CWM(I,J,L)*(IW(I,J,L))  *HTM(I,J,L))*DP*GI

CHOU checking       WRITE(1098,*)"CalCWI:", i,j,l, CLW(i,j), CLI(i,j),
CHOU checking     1                                cwm(i,j,l),IW(i,j,l), htm(i,j,l)
      ENDDO
      ENDDO
C
      ENDDO    ! L LOOP
C
C     END OF ROUTINE.
C
      RETURN
      END
