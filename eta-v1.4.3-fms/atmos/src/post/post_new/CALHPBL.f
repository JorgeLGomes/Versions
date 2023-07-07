      SUBROUTINE CALHPBL(HPBL)
C
C     CALCULATES THE  Planetary boundary layer height consistent w/ etafcst code
C
C
C     INPUT:
C     ------
C
C     Q2   (IM,JM,LM)  - TWICE THE TURBULENT ENERGY FIELD
C     ZINT (IM,JM,LP1) - ETA INTERFACES HEIGHT FIELD
C     LMH  (IM,JM)     - TOPOGRAPHY INDEXES ARRAY
C     IM,JM            - ARRAY DIMENSIONS FOR HORIZONTAL GRIDS
C                        IN A VECTORIZED FORM
C     LM               - ARRAY DIMENSION FOR VERTICAL GRIDS
C     LP1              = LM+1
C
C
C     OUTPUT:
C     -------
C
C     HPBL(IM,JM)      - PBL Height    !Chou added
C
C
C     RELEVANT CONSTANTS:
C     -------------------
C
      INCLUDE "parmeta"
      INCLUDE "params"

C     MINIMAL VALUE OF TURBULENT ENERGY:
      PARAMETER (EPSQ2F = 0.12 * 1.01)    !in etafcst/MIXLEN.f90,EPSQ2=0.12, FH=1.01

      INCLUDE "CTLBLK.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "LOOPS.comm"
      INCLUDE "PVRBLS.comm"
C
      REAL HPBL(IM,JM)
C
      DO  J=1,JM
      DO  I=1,IM
     
       LMHK = LMH(I,J) 
       LPBL = LMHK
!
!teste JOrge - 28-01-21 - T4-NEW
       NEPSQ2=0
       DO 100 L=LMHK, LMHK/2, -1   !Daniela 
          IF ((Q2(I,J,L) < EPSQ2F).and.(NEPSQ2<2)) THEN
             NEPSQ2=NEPSQ2+1
             CYCLE             ! raise L either ecreasing or increasing tke with height
          END IF
          IF (Q2(I,J,L) >= EPSQ2F) THEN
             LPBL=L
             NEPSQ2=0
             CYCLE
          ELSEIF (NEPSQ2>=2) THEN
             EXIT   ! raise L either in decreasing or in increasing tke with height
             END IF
 100   END DO
       IF (ZINT(I,J,LMHK+1). lt. 2500.) LPBL = LPBL - 1  !raise pbl thickness 
!
!----------------------
! THE HEIGHT OF THE PBL
!----------------------
       HPBL(I,J) = ZINT(I,J,LPBL) - ZINT(I,J,LMHK+1)
       HPBL(I,J) = AMAX1(300., HPBL(I,J))
      END DO
      END DO

      RETURN
      END
