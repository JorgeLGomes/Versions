      SUBROUTINE CALHEL(UST,VST,HELI)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    CALHEL       COMPUTES STORM RELATIVE HELICITY
C   PRGRMMR: BALDWIN         ORG: W/NP2      DATE: 94-08-22       
C     
C ABSTRACT:
C     THIS ROUTINE COMPUTES ESTIMATED STORM MOTION AND
C     STORM-RELATIVE ENVIRONMENTAL HELICITY.  
C     (DAVIES-JONES ET AL 1990) THE ALGORITHM PROCEEDS AS 
C     FOLLOWS.
C     
C     AT EACH WIND POINT MOVE UP VERTICALLY FROM THE LMV(K)-TH
C     ETA LAYER.  USE ALL LAYER WINDS WHOSE PRESSURES ARE 
C     FOUND BETWEEN 150 MB ABOVE THE SURFACE AND 300 MB LEVEL
C     TO COMPUTE A MASS WEIGHTED MEAN WIND.  THIS IS USED 
C     TO ESTIMATE THE STORM MOTION.  STORM MOTION IS SPECIFIED
C     AS 30 DEG TO THE RIGHT OF THE MEAN WIND AND 75% OF
C     THE MEAN WIND SPEED IF THE MEAN WIND IS LESS THAN 15 M/S.
C     IF THE MEAN WIND IS MORE THAN 15 M/S, THE STORM MOTION
C     IS 20 DEG TO THE RIGHT OF THE MEAN WIND AND 80% OF
C     ITS SPEED.  THIS STORM MOTION IS USED TO COMPUTE THE
C     STORM-RELATIVE HELICITY.  MOVING UP VERTICALLY FROM
C     THE LMV(K)-TH LAYER TO THE HIGHEST LAYER BELOW 3 KM
C     ABOVE GROUND, SUM UP THE FOLLOWING:
C        HELI = (V-VMEAN) * DU/DZ - (U-UMEAN) * DV/DZ
C     
C     
C PROGRAM HISTORY LOG:
C   94-08-22  MICHAEL BALDWIN
C   97-03-27  MICHAEL BALDWIN - SPEED UP CODE
C   98-06-15  T BLACK         - CONVERSION FROM 1-D TO 2-D
C   00-01-04  JIM TUCCILLO    - MPI VERSION               
C     
C USAGE:    CALHEL(UST,VST,HELI)
C   INPUT ARGUMENT LIST:
C     NONE     
C
C   OUTPUT ARGUMENT LIST: 
C     UST      - ESTIMATED U COMPONENT (M/S) OF STORM MOTION.
C     VST      - ESTIMATED V COMPONENT (M/S) OF STORM MOTION.
C     HELI     - STORM-RELATIVE HELICITY (M**2/S**2)
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C
C     LIBRARY:
C       COMMON   - VRBLS
C                  LOOPS
C                  PHYS 
C                  EXTRA
C                  MASKS
C                  OPTIONS
C                  INDX
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN 90
C     MACHINE : CRAY C-90
C$$$  
C     
C     
C     INCLUDE PARAMETERS.
      INCLUDE "parmeta"
      INCLUDE "params"
      INCLUDE "parm.tbl"
      PARAMETER (PI=3.141592654)
      PARAMETER (P150=15000.0,P300=30000.0,S15=15.0)
      PARAMETER (D3000=3000.0,PI6=0.5235987756,PI9=0.34906585)
C     
C     DECLARE VARIABLES
C     
      REAL UST(IM,JM),VST(IM,JM),HELI(IM,JM),ETOT(IM,JM)
C     
C     INCLUDE COMMON BLOCKS.
      INCLUDE "VRBLS.comm"
      INCLUDE "LOOPS.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "MASKS.comm"
      INCLUDE "PHYS.comm"
      INCLUDE "OPTIONS.comm"
      INCLUDE "INDX.comm"
      INCLUDE "CTLBLK.comm"
C     
C****************************************************************
C     START CALHEL HERE
C     
C     INITIALIZE ARRAYS.
C     
!$omp  parallel do
      DO J=JSTA,JEND
      DO I=1,IM
         UST(I,J)    = 0.0
         VST(I,J)    = 0.0
         HELI(I,J)   = 0.0
         ETOT(I,J)   = 0.0
      ENDDO
      ENDDO
C     
C     LOOP OVER HORIZONTAL GRID.
C 
      CALL EXCH(RES)
      CALL EXCH(PD)
      DO L = 1,LP1
        CALL EXCH(ZINT(1,1,L))
      END DO
C
!$omp  parallel do
!$omp& private(htsfc,ie,iw,pdslvk,pkl,psfck)
      DO L = 1,LM
        DO J=JSTA_M,JEND_M
        DO I=2,IM-1
          IE=I+IVE(J)
          IW=I+IVW(J)
          PDSLVK=(PD(IW,J)*RES(IW,J)+PD(IE,J)*RES(IE,J)+
     1           PD(I,J+1)*RES(I,J+1)+PD(I,J-1)*RES(I,J-1))*0.25
          PSFCK=AETA(LMV(I,J))*PDSLVK+PT
c         HTSFC=0.25*(ZINT(IW,J,LMV(I,J)+1)+ZINT(IE,J,LMV(I,J)+1)+
c    1                ZINT(I,J+1,LMV(I,J)+1)+ZINT(I,J-1,LMV(I,J)+1))
C     
C     COMPUTE MASS WEIGHTED MEAN WIND FROM 150 MB ABOVE
C      SURFACE TO 300 MB LAYER
C
         PKL = AETA(L)*PDSLVK+PT
         IF (PKL.LT.PSFCK-P150.AND.PKL.GT.P300) THEN
               UST(I,J) = UST(I,J) + U(I,J,L) * DETA(L)
               VST(I,J) = VST(I,J) + V(I,J,L) * DETA(L)
               ETOT(I,J) = ETOT(I,J) + DETA(L)
         ENDIF
        ENDDO
        ENDDO
      ENDDO


!$omp  parallel do
!$omp& private(rot,stspd,umean,unew,vmean,vnew)
      DO J=JSTA_M,JEND_M
      DO I=2,IM-1
         IF (ETOT(I,J).GT.0.0) THEN
           UMEAN = UST(I,J) / ETOT(I,J)
           VMEAN = VST(I,J) / ETOT(I,J)
           STSPD = SQRT(UMEAN*UMEAN+VMEAN*VMEAN)
C
C      COMPUTE STORM MOTION VECTOR
C
           IF (STSPD.LE.S15) THEN
            ROT=PI6
            UNEW=UMEAN*COS(ROT)+VMEAN*SIN(ROT)
            VNEW=VMEAN*COS(ROT)-UMEAN*SIN(ROT)
            UMEAN=0.75*UNEW
            VMEAN=0.75*VNEW
           ELSE
            ROT=PI9
            UNEW=UMEAN*COS(ROT)+VMEAN*SIN(ROT)
            VNEW=VMEAN*COS(ROT)-UMEAN*SIN(ROT)
            UMEAN=0.80*UNEW
            VMEAN=0.80*VNEW
           ENDIF
            UST(I,J)=UMEAN
            VST(I,J)=VMEAN
         ELSE
            UST(I,J)=0.0
            VST(I,J)=0.0
         ENDIF
      ENDDO
      ENDDO
C
C
C       COMPUTE STORM-RELATIVE HELICITY
C
!$omp  parallel do
!$omp& private(du1,du2,dv1,dv2,dz,dz1,dz2,dzabv,ie,iw,z1,z2,z3)
      DO L = 2,LM-1
        DO J=JSTA_M,JEND_M
        DO I=2,IM-1
          IW=I+IVW(J)
          IE=I+IVE(J)
          Z2=0.125*(ZINT(IW,J,L)+ZINT(IW,J,L+1)+
     &              ZINT(IE,J,L)+ZINT(IE,J,L+1)+
     &              ZINT(I,J+1,L)+ZINT(I,J+1,L+1)+
     &              ZINT(I,J-1,L)+ZINT(I,J-1,L+1))
          HTSFC=0.25*(ZINT(IW,J,LMV(I,J)+1)+ZINT(IE,J,LMV(I,J)+1)+
     1                ZINT(I,J+1,LMV(I,J)+1)+ZINT(I,J-1,LMV(I,J)+1))

          DZABV=Z2-HTSFC
C
          IF(DZABV.LT.D3000.AND.L.LE.LMV(I,J))THEN
            Z1=0.125*(ZINT(IW,J,L+1)+ZINT(IW,J,L+2)+
     &                ZINT(IE,J,L+1)+ZINT(IE,J,L+2)+
     &                ZINT(I,J+1,L+1)+ZINT(I,J+1,L+2)+
     &                ZINT(I,J-1,L+1)+ZINT(I,J-1,L+2))
            Z3=0.125*(ZINT(IW,J,L-1)+ZINT(IW,J,L)+
     &                ZINT(IE,J,L-1)+ZINT(IE,J,L)+
     &                ZINT(I,J+1,L-1)+ZINT(I,J+1,L)+
     &                ZINT(I,J-1,L-1)+ZINT(I,J-1,L))
            DZ=0.25*((ZINT(IW,J,L)+ZINT(IE,J,L)+
     &                ZINT(I,J-1,L)+ZINT(I,J+1,L))-
     &               (ZINT(IW,J,L+1)+ZINT(IE,J,L+1)+
     &                ZINT(I,J-1,L+1)+ZINT(I,J+1,L+1)))
            DZ1=Z1-Z2
            DZ2=Z2-Z3
            DU1=U(I,J,L+1)-U(I,J,L)
            DU2=U(I,J,L)-U(I,J,L-1)
            DV1=V(I,J,L+1)-V(I,J,L)
            DV2=V(I,J,L)-V(I,J,L-1)
            HELI(I,J)=((V(I,J,L)-VST(I,J))*
     1                (DZ2*(DU1/DZ1)+DZ1*(DU2/DZ2))
     2                -(U(I,J,L)-UST(I,J))*
     3                (DZ2*(DV1/DZ1)+DZ1*(DV2/DZ2)))
     4                *DZ/(DZ1+DZ2)+HELI(I,J) 
           ENDIF
        ENDDO
        ENDDO
      ENDDO
C
C
C     END OF ROUTINE.
C
      RETURN
      END
