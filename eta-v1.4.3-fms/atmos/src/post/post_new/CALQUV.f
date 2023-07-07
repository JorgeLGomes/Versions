      SUBROUTINE CALQUV(QUint,QVint)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C SUBPROGRAM:    CALQUV     COMPUTES vert integrtd MOISTURE TRANSPORT COMPONENTS
C SUBPROGRAM:               BASED ON CALMCVG
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-01-22       d
C
C ABSTRACT:
C     GIVEN SPECIFIC HUMIDITY, Q, and Cloud water/ice AND THE U-V WIND COMPONENTS
C     THIS ROUTINE CALCULATES  integral (qu) dp/g ;ntegral (qv) dp/g
C     WHERE,
C
C     THIS REQUIRES WINDS BE AT VELOCITY POINTS.
C     SPECIFIC HUMIDITY IS INTERPOLATED TO VELOCITY POINTS.
C
C USAGE:    CALL CALQUV(QUint, QVint)
C   INPUT ARGUMENT LIST:
C
C   OUTPUT ARGUMENT LIST:
C     QUint     -  INTEGRAL OF MOISTURE TRANSPORT IN ZONAL DIRECTION
C     QVint     -  INTEGRAL OF MOISTURE TRANSPORT IN MERIDIONAL DIRECTION
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
      REAL QUint(IM,JM), QVint(IM,JM)
      REAL U1D(IM,JM), V1D(IM,JM)
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
C
C
C***************************************************************************
C     START CALQUV  HERE.
C
C     INITIALIZE MOISTURE TRANSPORT ARRAY.  LOAD TEMPORARY WIND ARRAYS.
C
!$omp  parallel do
      DO J=JSTA,JEND
      DO I=1,IM
        QUINT(I,J) = D00
        QVINT(I,J) = D00
      ENDDO
      ENDDO

      DO L=1,LM

      DO J=JSTA,JEND
      DO I=1,IM
        U1D(I,J)  = U(I,J,L)
        V1D(I,J)  = V(I,J,L)
      ENDDO
      ENDDO
C
      CALL EXCH(U1D)
      CALL EXCH(V1D)
C
!$omp  parallel do
      DO J=JSTA_M,JEND_M
      DO I=2,IM-1
      DMEAN  =  (VTM(I,J-1,L) +VTM(I+IHW(J),J,L)+VTM(I+IHE(J),J,L)
     1          +VTM(I,J+1,L))
      IF (DMEAN .LT. 1) THEN   ! inside topography
         DMEAN  =  0
      ELSE
         DMEAN  =  1./(VTM(I,J-1,L) +VTM(I+IHW(J),J,L)+VTM(I+IHE(J),J,L)
     1                +VTM(I,J+1,L))
      ENDIF
Chou DP is the DPi at the wind point, averaged from the mass points
      DP     = -(PINT(I,J,L) -PINT(I,J,L+1))
C
C     CALC  U and V at MASS POINTS times DP/G
C
      UMDP     = (U1D(I,J-1)     *VTM(I,J-1,L)      +
     1            U1D(I+IHW(J),J)*VTM(I+IHW(J),J,L) +
     2            U1D(I+IHE(J),J)*VTM(I+IHE(J),J,L) +
     1            U1D(I,J+1)     *VTM(I,J+1,L)     )* DMEAN*DP*GI

      VMDP     = (V1D(I,J-1)     *VTM(I,J-1,L)      +
     1            V1D(I+IHW(J),J)*VTM(I+IHW(J),J,L) +
     2            V1D(I+IHE(J),J)*VTM(I+IHE(J),J,L) +
     1            V1D(I,J+1)     *VTM(I,J+1,L)     )* DMEAN* DP*GI
C
C     VERTICAL INTEGRAL  (U* q , V* q) * DP/G
C
      QUINT(I,J)= QUINT(I,J)+ UMDP*(Q(I,J,L)+CWM(I,J,L))*HTM(I,J,L)
      QVINT(I,J)= QVINT(I,J)+ VMDP*(Q(I,J,L)+CWM(I,J,L))*HTM(I,J,L)

CHOU checking       WRITE(1099,*)"alquv:", i,j,l, quint(i,j), qvint(i,j), u1D(i,j),
CHOU checking     2                 v1d(i,j),q(i,j,l), cwm(i,j,l),vtm(i,j,l)
      ENDDO
      ENDDO
C
      ENDDO    ! L LOOP
C
C     END OF ROUTINE.
C
      RETURN
      END
