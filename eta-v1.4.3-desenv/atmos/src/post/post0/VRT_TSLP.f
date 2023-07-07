      SUBROUTINE VRT_TSLP(TSLPB,TSLPB2,PMSLPB,PMSLPB2,PSLPB,PSLPB2
     1,                   IMI,JMI,LMI0,LMIMX,KBIMX
     2,                   NBCEX,PTOP,IUDETI)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  VRT_TSLP    INTERP TEMPERATURE IN VERTICAL
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 00-11-10
C
C  ABSTRACT:  VERTICALLY INTERPOLATES TEMEPRATURE FROM A
C             COARSER RESOLUTION DISTRIBUTION TO A FINER ONE.
C
C PROGRAM HISTORY LOG:
C    00-11-10  T BLACK - ORIGINATOR
C
C USAGE:  CALL PTETAE FROM SUBROUTINE BCEX
C
C   INPUT ARGUMENT LIST:
C      TSLPB - THE TEMPERATURE ON THE OUTER BOUNDARY ROW OF THE
C              NEST WITH THE PARENT VERTICAL DISTRIBUTION
C     TSLPB2 - THE TEMPERATURE ON THE SECOND BOUNDARY ROW OF THE
C              NEST WITH THE PARENT VERTICAL DISTRIBUTION
C     PMSLPB - THE MIDLAYER PRESSURE ON THE OUTER BOUNDARY ROW OF THE
C              NEST WITH THE PARENT VERTICAL DISTRIBUTION
C    PMSLPB2 - THE MIDLAYER PRESSURE ON THE SECOND BOUNDARY ROW OF THE
C              NEST WITH THE PARENT VERTICAL DISTRIBUTION
C      PSLPB - THE PARENT SLP ON THE OUTER BOUNDARY ROW OF THE NEST 
C     PSLPB2 - THE PARENT SLP ON THE SECOND BOUNDARY ROW OF THE NEST 
C      KBIMX - THE MAXIMUM LENGTH OF ANY NEST BOUNDARY'S HORIZONTAL
C              EXTENT
C      LMIMX - THE MAXIMUM NUMBER OF MODEL LEVELS FOR ANY NEST
C      NBCEX - THE NUMBER OF NESTS
C        IMI - THE ARRAY OF IM's FOR ALL NESTS
C        JMI - THE ARRAY OF JM's FOR ALL NESTS
C       LMI0 - THE ARRAY OF LM's FOR ALL NESTS
C     IUDETI - THE BASE UNIT NUMBER OF THE NEST DETA's
C
C   OUTPUT ARGUMENT LIST:
C      TSLPB - THE TEMPERATURE ON THE OUTER BOUNDARY ROW OF THE
C              NEST WITH THE NEST VERTICAL DISTRIBUTION
C     TSLPB2 - THE TEMPERATURE ON THE SECOND BOUNDARY ROW OF THE
C              NEST WITH THE NEST VERTICAL DISTRIBUTION
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:  NONE
C
C----------------------------------------------------------------------
      INCLUDE "parmeta"
C----------------------------------------------------------------------
                             P A R A M E T E R
     1 (LMOP=LM+1,LMOM=LM-1)
C
                             P A R A M E T E R
     1 (RD=287.04,G=9.80,RG=1./G)
C--------------------------------------------------------------------
                             P A R A M E T E R
     & (K15=SELECTED_REAL_KIND(15))
C
                             R E A L
     & (KIND=K15) B,C,ALPET,ALPETK
     &,           ALP(KBIMX,LM+1),ALPETA(KBIMX),ALPSQ(KBIMX,LM+1)
C--------------------------------------------------------------------
                             R E A L
     & TSLPB(KBIMX,LMIMX,NBCEX),TSLPB2(KBIMX,LMIMX,NBCEX)
     &,PMSLPB(KBIMX,LM,NBCEX),PMSLPB2(KBIMX,LM,NBCEX)
     &,PSLPB(KBIMX,NBCEX),PSLPB2(KBIMX,NBCEX)
C--------------------------------------------------------------------
                             R E A L
     1 DETAI(LMIMX)
     2,ALPMIDO(KBIMX,LM),ALPMIDI(KBIMX,LMIMX)
     3,TMIDI(KBIMX,LMIMX)
C-----------------------------------------------------------------------
                               I N T E G E R
     1 IMI(9),JMI(9),LMI0(9)
C-----------------------------------------------------------------------
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C-----------------------------------------------------------------------
C
      IUDETI0=IUDETI
C
C-------------------------------------------------------------
C-------------------------------------------------------------
C***
C***  THE GRAND LOOP OVER EACH OF THE NESTS
C***
C-------------------------------------------------------------
C-------------------------------------------------------------
      DO 1000 NB=1,NBCEX
C-------------------------------------------------------------
      KBI=2*IMI(NB)+JMI(NB)-3
      KBI2=KBI-4
      LMI=LMI0(NB)
C-------------------------------------------------------------
C***
C***  COMPUTE THE LOG PRESSURES AT THE PARENT MIDLAYERS
C***
      DO L=1,LM
        DO K=1,KBI
          ALPMIDO(K,L)=ALOG(PMSLPB(K,L,NB))
        ENDDO
      ENDDO
C***
C***  READ THE DETAs FOR THE NEST AND COMPUTE THE LOG PRESSURE
C***  AT THE NEST MIDLAYERS
C***
CCCC  IUDETI0=IUDETI0+1
      REWIND IUDETI0
      READ(IUDETI0)(DETAI(L),L=1,LMI)
C
      DETAU=0.
      AETAI=0.
C
      DO L=1,LMI
        AETAI=0.5*(DETAI(L)+DETAU)+AETAI
C
        DO K=1,KBI 
          ALPMIDI(K,L)=ALOG(AETAI*(PSLPB(K,NB)-PTOP)+PTOP)
        ENDDO
C
        DETAU=DETAI(L)
      ENDDO
C
C----------------------------------------------------------------------
C***  INTERPOLATE T LINEAR IN LN(P) ASSUMING THAT BOTH 
C***  VERTICAL DISTRIBUTIONS SHARE THE SAME TOP PRESSURE
C----------------------------------------------------------------------
C
      DO K=1,KBI
C
        LPARENT=2
        TDIF=TSLPB(K,2,NB)-TSLPB(K,1,NB)
        ALPDIF=ALPMIDO(K,2)-ALPMIDO(K,1)
        ALPBOT=ALPMIDO(K,LM)
C
        DO LNEST=1,LMI
          ALPNEST=ALPMIDI(K,LNEST)
C
          IF(ALPNEST.GT.ALPMIDO(K,LPARENT).AND.
     1       ALPNEST.LE.ALPBOT)THEN
            LPARENT=LPARENT+1
            TDIF=TSLPB(K,LPARENT,NB)-TSLPB(K,LPARENT-1,NB)
            ALPDIF=ALPMIDO(K,LPARENT)-ALPMIDO(K,LPARENT-1)
          ENDIF
C
          FRAC=(ALPMIDI(K,LNEST)-ALPMIDO(K,LPARENT-1))/ALPDIF
          TMIDI(K,LNEST)=TSLPB(K,LPARENT-1,NB)+FRAC*TDIF
        ENDDO
C
      ENDDO
C
C***  REPLACE THE PARENT VALUES IN TSLPB WITH THE NEST VALUES
C
      DO L=1,LMI
        DO K=1,KBI
          TSLPB(K,L,NB)=TMIDI(K,L)
        ENDDO
      ENDDO
C----------------------------------------------------------------------
C***
C***  NOW DO THE SECOND BOUNDARY ROW
C***
C----------------------------------------------------------------------
C
C***  COMPUTE THE LOG PRESSURES AT THE PARENT MIDLAYERS
C
      DO L=1,LM
        DO K=1,KBI2
          ALPMIDO(K,L)=ALOG(PMSLPB2(K,L,NB))
        ENDDO
      ENDDO
C
      DETAU=0.
      AETAI=0.
C
C***  AND AT THE NEST MIDLAYERS
C
      DO L=1,LMI
        AETAI=0.5*(DETAI(L)+DETAU)+AETAI
C
        DO K=1,KBI2
          ALPMIDI(K,L)=ALOG(AETAI*(PSLPB2(K,NB)-PTOP)+PTOP)
        ENDDO
C
        DETAU=DETAI(L)
      ENDDO
C
C----------------------------------------------------------------------
C***  INTERPOLATE T LINEAR IN LN(P) ASSUMING THAT BOTH 
C***  VERTICAL DISTRIBUTIONS SHARE THE SAME TOP PRESSURE
C----------------------------------------------------------------------
C
      DO K=1,KBI2
C
        LPARENT=2
        TDIF=TSLPB2(K,2,NB)-TSLPB2(K,1,NB)
        ALPDIF=ALPMIDO(K,2)-ALPMIDO(K,1)
        ALPBOT=ALPMIDO(K,LM)
C
        DO LNEST=1,LMI
          ALPNEST=ALPMIDI(K,LNEST)
C
          IF(ALPNEST.GT.ALPMIDO(K,LPARENT).AND.
     1       ALPNEST.LE.ALPBOT)THEN
            LPARENT=LPARENT+1
            TDIF=TSLPB2(K,LPARENT,NB)-TSLPB2(K,LPARENT-1,NB)
            ALPDIF=ALPMIDO(K,LPARENT)-ALPMIDO(K,LPARENT-1)
          ENDIF
C
          FRAC=(ALPNEST-ALPMIDO(K,LPARENT-1))/ALPDIF
          TMIDI(K,LNEST)=TSLPB2(K,LPARENT-1,NB)+FRAC*TDIF
        ENDDO
C
      ENDDO
C
C***  REPLACE THE PARENT VALUES IN TSLPB2 WITH THE NEST VALUES
C
      DO L=1,LMI
        DO K=1,KBI2
          TSLPB2(K,L,NB)=TMIDI(K,L)
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
 1000 CONTINUE
C
      RETURN
      END

