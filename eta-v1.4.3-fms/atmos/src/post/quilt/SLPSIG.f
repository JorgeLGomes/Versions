      SUBROUTINE SLPSIG(PD,FIS,SM,TSIG,QSIG,CWMSIG,HBM2
     +,            U00,SPL,LSL,UL,DETA,PT,PSLP)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  SLPSIG      SLP REDUCTION
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-04-22
C
C ABSTRACT:  THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE
C            REDUCTION USING EITHER THE MESINGER RELAXATION
C            METHOD OR THE STANDARD NCEP REDUCTION FOR
C            SIGMA COORDINATES.  A BY-PRODUCT IS THE
C            SET OF VALUES FOR THE UNDERGROUND TEMPERATURES
C            ON THE SPECIFIED PRESSURE LEVELS
C
C PROGRAM HISTORY LOG:
C   99-09-23  T BLACK - REWRITTEN FROM ROUTINE SLP (ETA
C                       COORDINATES)
C   00-08-17  H CHUANG - MODIFIED THE ROUTINE TO BE INCLUDED
C                        IN THE QUILT INSTEAD OF POST.  
C
C USAGE:  CALL SLPSIG FROM SUBROUITNE ETA2P
C
C   INPUT ARGUMENT LIST:
C     PD   - SFC PRESSURE MINUS PTOP
C     FIS  - SURFACE GEOPOTENTIAL
C     T    - TEMPERATURE 
C     Q    - SPECIFIC HUMIDITY
C     FI   - GEOPOTENTIAL
C     PT   - TOP PRESSURE OF DOMAIN
C
C   OUTPUT ARGUMENT LIST:
C     PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C             NONE
C
C-----------------------------------------------------------------------
      INCLUDE "parmeta"
      INCLUDE "PARA.comm"
      INCLUDE "mpif.h"
      INCLUDE "mpp.h"
      PARAMETER (LMP1=LM+1)
C-----------------------------------------------------------------------
                            P A R A M E T E R
     & (LP1=LSM+1,IMJM=IM*JM-JM/2,IM_JM=IM*JM,JM2=JM-2)
                            P A R A M E T E R
     & (NFILL=8,NRLX1=500,NRLX2=100,KSLPD=1
     &, OVERRC=1.5,AD05=OVERRC*0.05,CFT0=OVERRC-1.
     &, RD=287.04,ROG=RD/9.8,PQ0=379.90516,A2=17.2693882
     &, A3=273.16,A4=35.86,GAMMA=6.5E-3,RGAMOG=GAMMA*ROG
     &, H1M12=1.E-12)
C-----------------------------------------------------------------------
                            R E A L
     & PD(IM,MY_JSD:MY_JED),FIS(IM,MY_JSD:MY_JED)
     &,SM(IM,MY_JSD:MY_JED),HBM2(IM,MY_JSD:MY_JED)
     &,HTM2D(IM,MY_JSD:MY_JED),TG(IM,MY_JSD:MY_JED)
     &,U00(IM,MY_JSD:MY_JED),HTM(IM,MY_JSD:MY_JED,LSM)
     &,T(IM,MY_JSD:MY_JED,LSL),Q(IM,MY_JSD:MY_JED,LSL)
     &,FI(IM,MY_JSD:MY_JED,LSL),TSIG(IM,MY_JSD:MY_JED,LM)
     &,QSIG(IM,MY_JSD:MY_JED,LM),CWMSIG(IM,MY_JSD:MY_JED,LM)
     &,ALPINT(IM,MY_JSD:MY_JED,LMP1),ZINT(IM,MY_JSD:MY_JED,LMP1)
     &,PINT(IM,MY_JSD:MY_JED,LMP1),IW(IM,MY_JSD:MY_JED,LM)
     &,IWU,IWL
                            R E A L
     & PSLP(IM,MY_JSD:MY_JED),TTV(IM,MY_JSD:MY_JED)
     &,SLPX(IM,MY_JSD:MY_JED)
     &,P1(IM,MY_JSD:MY_JED)
c     &,T00B(IM,MY_JSD:MY_JED),T00A(IM,MY_JSD:MY_JED)
                            R E A L
     & SPL(LSM),SPLI(LSM)
                            R E A L
     & DETA(LM),RDETA(LM),AETA(LM),F4Q2(LM),ETA(LMP1),DFL(LMP1)
     &,UL(2*LM)
C-----------------------------------------------------------------------
      real tempt(7,7)
C-----------------------------------------------------------------------
                            I N T E G E R
     & KMNTM(LSM),IMNT(IMJM,LSM),JMNT(IMJM,LSM)
     &,LMH(IM,MY_JSD:MY_JED),NL1X(IM,MY_JSD:MY_JED)
     &,IHOLD(IM_JM),JHOLD(IM_JM)
                            I N T E G E R
     & IHE(JM),IHW(JM),IVE(JM)
     &,IVW(JM)
C-----------------------------------------------------------------------
                            L O G I C A L
     & SIGMA,STDRD,DONE(IM,MY_JSD:MY_JED)
C-----------------------------------------------------------------------
      STDRD=.FALSE.
C-----------------------------------------------------------------------
C***
C***  CALCULATE THE I-INDEX EAST-WEST INCREMENTS
C***
      print*,'in slpsig'
      DO J=1,JM
        IHE(J)=MOD(J+1,2)
        IHW(J)=IHE(J)-1
        IVE(J)=MOD(J,2)
        IVW(J)=IVE(J)-1
      ENDDO
C-----------------------------------------------------------------------
C***
C***  INITIALIZE ARRAYS.  LOAD SLP ARRAY WITH SURFACE PRESSURE.
C***
CC!$omp parallel do 
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        PSLP(I,J)=PD(I,J)+PT
        TTV(I,J)=0.
        LMH(I,J)=0
      ENDDO
      ENDDO
C      
C*** COMPUTE ETA AND AETA 
      ETA(1)=0.0
      DO L=2,LM+1
       ETA(L)=ETA(L-1)+DETA(L-1)
      END DO 
      DO L=1,LM
       AETA(L)=0.5*(ETA(L)+ETA(L+1))
      ENDDO
C
C***  CALCULATE SEA LEVEL PRESSURE FOR PROFILES (AND POSSIBLY
C***  FOR POSTING BY POST PROCESSOR).
C
C***  "STDRD" REFERS TO THE "STANDARD" SLP REDUCTION SCHEME.
C
      IF(STDRD)GO TO 400
C--------------------------------------------------------------------
C***
C***  CREATE A 3-D "HEIGHT MASK" FOR THE SPECIFIED PRESSURE LEVELS
C***  (1 => ABOVE GROUND) AND A 2-D INDICATOR ARRAY THAT SAYS 
C***  WHICH PRESSURE LEVEL IS THE LOWEST ONE ABOVE THE GROUND
C***
      DO 100 L=1,LSM
      SPLL=SPL(L)
C       
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        PSFC=PD(I,J)+PT
        PCHK=PSFC
        IF(NFILL.GT.0)THEN
          DO LL=1,NFILL
            PCHK=PCHK-DETA(LM+1-LL)*PD(I,J)
          ENDDO
        ENDIF
        IF(SM(I,J).GT.0.5.AND.FIS(I,J).LT.10.)PCHK=PSFC
C
c       IF(SPLL.LT.PSFC)THEN
        IF(SPLL.LT.PCHK)THEN
          HTM(I,J,L)=1.
        ELSE
          HTM(I,J,L)=0.
          IF(L.GT.1.AND.HTM(I,J,L-1).GT.0.5)LMH(I,J)=L-1
        ENDIF
C
        IF(L.EQ.LSM.AND.HTM(I,J,L).GT.0.5)LMH(I,J)=LSM
      ENDDO
      ENDDO
C
  100 CONTINUE
C--------------------------------------------------------------------
C  INTERPOLATE T AND Q FROM SIGMA TO PRESSURE LEVELS
C
C---------------------------------------------------------------------
C***
C***  VALUES OF T ON THE OUTPUT PRESSURE LEVELS ABOVE GROUND MUST BE
C***  KNOWN BEFORE THEY CAN BE FILLED IN BELOW GROUND IN THE ROUTINE 
C***  WHICH COMPUTES THE MESINGER SEA LEVEL PRESSURE.  THEREFORE DO ALL
C***  INTERPOLATION ABOVE GROUND NOW.
C***
C
C DEFINE ALPINT AND PINT
      DO L=1,LMP1
      DO J=JSTA_I,JEND_I
      DO I=1,IM
       PINT(I,J,L)=PD(I,J)*ETA(L)+PT
       ALPINT(I,J,L)=LOG(PD(I,J)*ETA(L)+PT)
      END DO
      END DO
      END DO
C
C COMPUTE ZINT 
CC!$omp  parallel do
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        ZINT(I,J,LMP1)=FIS(I,J)/9.8
      ENDDO
      ENDDO
C
C     COMPUTE VALUES FROM THE SURFACE UP.
      DO L=LM,1,-1
CC!$omp  parallel do
        DO J=JSTA_I,JEND_I
        DO I=1,IM
          ZINT(I,J,L)=ZINT(I,J,L+1)+TSIG(I,J,L)*(QSIG(I,J,L)*.608+1.)
     1	        *RD*(ALPINT(I,J,L+1)-ALPINT(I,J,L))/9.8	  
        ENDDO
        ENDDO
      END DO
c      
      if (i.eq.10.and.j.eq.jm-10)then
      print*,'writing sample input height, T'
c      do j=JSTA_I,JEND_I
c       j=10  
c       i=86
       do l=1,lmp1
        print*,'i,j,l,P,T,ZINT = ',i,j,l,pint(i,j,l),tsig(i,j,l)
     +	,zint(i,j,l)
       end do
c      end do
      end if      
C	         
C  SET UP UTIM FOR THIS TIME STEP
C
        UTIM=1.
        CLIMIT=1.E-20
C
        do75: DO L=1,LM
        IF(L.EQ.1)THEN
CC!$omp  parallel do
          DO J=JSTA_I,JEND_I
          DO I=1,IM
            IW(I,J,L)=0.
          ENDDO
          ENDDO
          CYCLE do75
        ENDIF
C
CC!$omp  parallel do
CC!$omp& private(cwmkl,fiq,hh,lml,pp,qi,qkl,qw,tkl,tmt0,tmt15,u00kl)
        doout70: DO J=JSTA_I,JEND_I
        doin70: DO I=1,IM
        LML=LM-LM   
        HH=1*HBM2(I,J)
        TKL=TSIG(I,J,L)
        QKL=QSIG(I,J,L)
        CWMKL=CWMSIG(I,J,L)
        TMT0=(TKL-273.16)*HH
        TMT15=AMIN1(TMT0,-15.)*HH    
        PP=PD(I,J)*AETA(L)+PT
        QW=HH*PQ0/PP*EXP(HH*A2*(TKL-A3)/(TKL-A4))
        QI=QW*(1.+0.01*AMIN1(TMT0,0.))     
        U00KL=U00(I,J)+UL(L+LML)*(0.95-U00(I,J))*UTIM
C
        IF(TMT0.LT.-15.0)THEN
          FIQ=QKL-U00KL*QI
          IF(FIQ.GT.0..OR.CWMKL.GT.CLIMIT) THEN
            IW(I,J,L)=1.
          ELSE
            IW(I,J,L)=0.
          ENDIF
        ENDIF
C
        IF(TMT0.GE.0.0)IW(I,J,L)=0.
        IF(TMT0.LT.0.0.AND.TMT0.GE.-15.0)THEN
          IW(I,J,L)=0.
          IF(IW(I,J,L-1).EQ.1..AND.CWMKL.GT.CLIMIT)IW(I,J,L)=1.
        ENDIF
C
        END DO doin70
        END DO doout70
        END DO do75
C-----------------------------------------------------------------
        DO 228 LP=1,LSL
	ALSL=LOG(SPL(LP))
        NHOLD=0
C
        doout125: DO J=JSTA_I,JEND_I
        doin125: DO I=1,IM
C
        T(I,J,LP)=-1.E6
        Q(I,J,LP)=-1.E6
        FI(I,J,LP)=-1.E6
C
C***  LOCATE VERTICAL INDEX OF MODEL INTERFACE JUST BELOW
C***  THE PRESSURE LEVEL TO WHICH WE ARE INTERPOLATING.
C
        do115: DO L=2,LM+1
c       IF((ALPINT(I,J,L)).LT.ALSL)GO TO 115
        IF(ALPINT(I,J,L).GE.ALSL)THEN
          NL1X(I,J)=L
          NHOLD=NHOLD+1
          IHOLD(NHOLD)=I
          JHOLD(NHOLD)=J
          CYCLE  doin125
        ELSEIF(ALPINT(I,J,LMP1).LT.ALSL)THEN
          NL1X(I,J)=LMP1
          NHOLD=NHOLD+1
          IHOLD(NHOLD)=I
          JHOLD(NHOLD)=J
          CYCLE doin125
        ENDIF
        END DO do115
C
        END DO doin125
        END DO doout125
        IF(NHOLD.EQ.0)GO TO 228
C
        TRF=2.*ALSL
C***
C***  NL1X WAS SET BASED ON INTERFACE PRESSURE VALUES.  HOWEVER, WE 
C***  ARE READY TO INTERPOLATE VALUES FROM MIDLAYER LOCATIONS SO 
C***  ADJUST NL1X UP ONE LAYER IF NEEDED SO THAT THE INDEX OF THE
C***  MIDLAYER LOCATION IS IMMEDIATELY BELOW THE PRESSURE LEVEL
C***  TO WHICH INTERPOLATION IS BEING DONE.
C***
CC!$omp  parallel do
CC!$omp& private(ahf,ahfo,ahfq,ahfq2,ahfqc,ahfqi,ai,b,bi,bom,
CC!$omp&         bq,bq2,bqc,bqc_2,bqi_2,bqi,fac,gmiw,gmiw_2,iwl,iwu,
CC!$omp&         lmp1,pl,pnl1,pu,q2a,q2b,qabv,qi,qint,ql,
CC!$omp&         qsat,qu,qw,rhu,tabv,tblo,tl,tmt0,tmt15,tu,
CC!$omp&         tvrabv,tvrblo,tvrl,tvru,zl,zu)
        DO 220 NN=1,NHOLD
        I=IHOLD(NN)
        J=JHOLD(NN)
        DIFU=ALSL-ALPINT(I,J,NL1X(I,J)-1)
        DIFL=ALPINT(I,J,NL1X(I,J))-ALSL
        IF(DIFU.LT.DIFL)NL1X(I,J)=NL1X(I,J)-1
        PNL1=PINT(I,J,NL1X(I,J))
C---------------------------------------------------------------------
C***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
C***  HUMIDITY, CLOUD WATER/ICE, AND TKE.
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
C---------------------------------------------------------------------
C
        IF(NL1X(I,J).EQ.1)THEN
          PU=PINT(I,J,2)
          ZU=ZINT(I,J,2)
          TU=0.5*(TSIG(I,J,1)+TSIG(I,J,2))
          QU=0.5*(QSIG(I,J,1)+QSIG(I,J,2))
C
          IWU=0.5*(IW(I,J,1)+IW(I,J,2))
          TMT0=TU-273.16
          TMT15=AMIN1(TMT0,-15.)
          AI=0.008855
          BI=1.
          IF(TMT0.LT.-20.)THEN
            AI=0.007225
            BI=0.9674
          ENDIF
          QW=PQ0/PU
     1      *EXP(A2*(TU-A3)/(TU-A4))
          QI=QW*(BI+AI*AMIN1(TMT0,0.))
          QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
          IF(TMT0.LT.-15.)THEN
            QSAT=QI
          ELSEIF(TMT0.GE.0.)THEN
            QSAT=QINT
          ELSE
            IF(IWU.GT.0.)THEN
              QSAT=QI
            ELSE
              QSAT=QINT
            ENDIF
          ENDIF
C
C***  USE RH FOR LIQUID WATER NO MATTER WHAT.
C***  DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
          QSAT=QW
C
          RHU=QU/QSAT
C
          IF(RHU.GT.1.)THEN
            RHU=1.
            QU =RHU*QSAT
          ENDIF
C
          IF(RHU.LT.0.01)THEN
            RHU=0.01
            QU =RHU*QSAT
          ENDIF
C
          TVRU=TU*(1.+0.608*QU)
          TVRABV=TVRU*(SPL(LP)/PU)**RGAMOG
          TABV=TVRABV/(1.+0.608*QU)
C     
          TMT0=TABV-273.16
          TMT15=AMIN1(TMT0,-15.)
          AI=0.008855
          BI=1.
          IF(TMT0.LT.-20.)THEN
            AI=0.007225
            BI=0.9674
          ENDIF
          QW=PQ0/SPL(LP)
     1      *EXP(A2*(TABV-A3)/(TABV-A4))
          QI=QW*(BI+AI*AMIN1(TMT0,0.))
          QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
          IF(TMT0.LT.-15.)THEN
            QSAT=QI
          ELSEIF(TMT0.GE.0.)THEN
            QSAT=QINT
          ELSE
            IF(IWU.GT.0.)THEN
              QSAT=QI
            ELSE
              QSAT=QINT
            ENDIF
          ENDIF
C
C***  USE RH FOR LIQUID WATER NO MATTER WHAT.
C***  DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
          QSAT=QW
C
          QABV =RHU*QSAT
          QABV =AMAX1(1.E-12,QABV)
          B    =TABV
          BQ   =QABV
          GMIW =IW(I,J,1) 
          AHF  =0.
          AHFQ =0.
          FAC  =0.
C
        ELSEIF(NL1X(I,J).EQ.LMP1)THEN 
C---------------------------------------------------------------------
C***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
C---------------------------------------------------------------------
C
          PL=PINT(I,J,LM-1)
          ZL=ZINT(I,J,LM-1)
          TL=0.5*(TSIG(I,J,LM-2)+TSIG(I,J,LM-1))
          QL=0.5*(QSIG(I,J,LM-2)+QSIG(I,J,LM-1))
C     
          IWL=0.5*(IW(I,J,LM-2)+IW(I,J,LM-1))

          TMT0=TL-273.16
          TMT15=AMIN1(TMT0,-15.)
          AI=0.008855
          BI=1.
          IF(TMT0.LT.-20.)THEN
            AI=0.007225
            BI=0.9674
          ENDIF
          QW=PQ0/PL
     1      *EXP(A2*(TL-A3)/(TL-A4))
          QI=QW*(BI+AI*AMIN1(TMT0,0.))
          QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
          IF(TMT0.LT.-15.)THEN
            QSAT=QI
          ELSEIF(TMT0.GE.0.)THEN
            QSAT=QINT
          ELSE
            IF(IWL.GT.0.)THEN
              QSAT=QI
            ELSE
              QSAT=QINT
            ENDIF
          ENDIF
C
C***  USE RH FOR LIQUID WATER NO MATTER WHAT.
C***  DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
          QSAT=QW
C
          RHL=QL/QSAT
C
          IF(RHL.GT.1.)THEN
            RHL=1.
            QL =RHL*QSAT
          ENDIF
C
          IF(RHL.LT.0.01)THEN
            RHL=0.01
            QL =RHL*QSAT
          ENDIF
C
          TVRL  =TL*(1.+0.608*QL)
          TVRBLO=TVRL*(SPL(LP)/PL)**RGAMOG
          TBLO  =TVRBLO/(1.+0.608*QL)
C     
          TMT0=TBLO-273.16
          TMT15=AMIN1(TMT0,-15.)
          AI=0.008855
          BI=1.
          IF(TMT0.LT.-20.)THEN
            AI=0.007225
            BI=0.9674
          ENDIF
          QW=PQ0/SPL(LP)
     1      *EXP(A2*(TBLO-A3)/(TBLO-A4))
          QI=QW*(BI+AI*AMIN1(TMT0,0.))
          QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
          IF(TMT0.LT.-15.)THEN
            QSAT=QI
          ELSEIF(TMT0.GE.0.)THEN
            QSAT=QINT
          ELSE
            IF(IWL.GT.0.) THEN
              QSAT=QI
            ELSE
              QSAT=QINT
            ENDIF
          ENDIF
C
C***  USE RH FOR LIQUID WATER NO MATTER WHAT.
C***  DELETE THE FOLLOWING LINE TO SWITCH BACK TO RH FOR ICE
          QSAT=QW
C
          QBLO =RHL*QSAT
          QBLO =AMAX1(1.E-12,QBLO)
          B    =TBLO
          BQ   =QBLO
          AHF  =0.
          AHFQ =0.
          FAC  =0.
C
        ELSE
C---------------------------------------------------------------------
C***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
C---------------------------------------------------------------------
C
          B     =TSIG(I,J,NL1X(I,J))
          BQ    =QSIG(I,J,NL1X(I,J))
          FAC  =2.*ALOG(PT+PD(I,J)*AETA(NL1X(I,J)))
          AHF  =(B-TSIG(I,J,NL1X(I,J)-1))/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
          AHFQ =(BQ-QSIG(I,J,NL1X(I,J)-1))/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
        ENDIF
C
        T(I,J,LP)=B+AHF*(TRF-FAC)
        Q(I,J,LP)=BQ+AHFQ*(TRF-FAC)
        Q(I,J,LP)=AMAX1(Q(I,J,LP),H1M12)
        FI(I,J,LP)=(PNL1-SPL(LP))/(SPL(LP)+PNL1)
     1       *((ALSL+ALPINT(I,J,NL1X(I,J))-FAC)*AHF+B)*RD*2.
     2       +ZINT(I,J,NL1X(I,J))*9.8
c
        if(i.eq.86.and.j.eq.10)then
	 print*,'LP,NL1X,T,Q,FI= ',lp,NL1X(I,J),T(I,J,LP)
     1,  Q(I,J,LP),FI(I,J,LP)   	
	END IF
  220   CONTINUE
C
c 225   CONTINUE
C     
  228   CONTINUE
C
C***  THE END OF VERTICAL INTERPOLATION FROM SIGMA TO PRESSURE
c put the 1000mb temp in the tg
c      do j=JSTA_I,JEND_I
c      do i=1,im      
c       t00b(i,j)=t(i,j,lsm)
c      end do
c      end do
CCC
C***
C***  WE REACH THIS LINE IF WE WANT THE MESINGER ETA SLP REDUCTION
C***  BASED ON RELAXATION TEMPERATURES.  THE FIRST STEP IS TO
C***  FIND THE HIGHEST LAYER CONTAINING MOUNTAINS.
C***
      DO 230 L=LSM,1,-1
C
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        IF(HTM(I,J,L).LT.0.5)GO TO 230
      ENDDO
      ENDDO
C
      LHMNT=L+1
      GO TO 235
  230 CONTINUE
C
  235 CONTINUE
      CALL MPI_ALLREDUCE
     &  (LHMNT,LXXX,1,MPI_INTEGER,MPI_MIN,MPI_COMM_COMP,IERR)
      LHMNT = LXXX
      IF(LHMNT.EQ.LP1)THEN
        GO TO 430
      ENDIF
C***
C***  NOW GATHER THE ADDRESSES OF ALL THE UNDERGROUND POINTS.
C***
CC!$omp parallel do private(kmn,kount)
      DO 250 L=LHMNT,LSM
      KMN=0
      KMNTM(L)=0
      KOUNT=0
      doout240: DO J=JSTA_IM2,JEND_IM2
      doin240: DO I=2,IM-1
      KOUNT=KOUNT+1
      IMNT(KOUNT,L)=0
      JMNT(KOUNT,L)=0
      IF(HTM(I,J,L).GT.0.5) CYCLE doin240
      KMN=KMN+1
      IMNT(KMN,L)=I
      JMNT(KMN,L)=J
      END DO doin240
      END DO doout240
      KMNTM(L)=KMN
  250 CONTINUE
C
C***  AS THE FIRST GUESS, SET THE UNDERGROUND TEMPERATURES EQUAL
C***  TO 0.0C.
C
c     IF(NTSD.EQ.1)THEN
        KMM=KMNTM(LSM)
CC!$omp parallel do private(i,j,lmap1),shared(t)
        doout260: DO KM=1,KMM
        I=IMNT(KM,LSM)
        J=JMNT(KM,LSM)
        LMAP1=LMH(I,J)+1
        doin260: DO L=LMAP1,LSM
        T(I,J,L)=273.15
        END DO doin260
        END DO doout260
c     ENDIF
C
C***  CREATE A TEMPORARY TV ARRAY, AND FOLLOW BY SEQUENTIAL
C***  OVERRELAXATION, DOING NRLX PASSES.
C
c     IF(NTSD.EQ.1)THEN
        NRLX=NRLX1
c     ELSE
c       NRLX=NRLX2
c     ENDIF
C
CC!$omp parallel do private(i,j,tinit,ttv)
      DO 300 L=LHMNT,LSM
C
      doout270: DO J=JSTA_I,JEND_I
      doin270: DO I=1,IM
      TTV(I,J)=T(I,J,L)
      HTM2D(I,J)=HTM(I,J,L)
      END DO doin270
      END DO doout270
C
C***  FOR GRID BOXES NEXT TO MOUNTAINS, COMPUTE TV TO USE AS
C***  BOUNDARY CONDITIONS FOR THE RELAXATION UNDERGROUND
C

      CALL UPDATE(HTM2D)
      DO J=JSTA_IM2,JEND_IM2
      DO I=2,IM-1
c        print*,'I,J,L,HTM2D= ',I,J,L,IHW(J),IHE(J)
c     1,HTM2D(I,J),HTM2D(I+IHW(J),J-1)	
c     2,HTM2D(I+IHE(J),J-1),HTM2D(I+IHW(J),J+1)
c     3,HTM2D(I+IHE(J),J+1),HTM2D(I-1     ,J)
c     4,HTM2D(I+1     ,J),HTM2D(I       ,J-2)
c     5,HTM2D(I       ,J+2)	
        IF(HTM2D(I,J).GT.0.5.AND.
     1     HTM2D(I+IHW(J),J-1)*HTM2D(I+IHE(J),J-1)
     2    *HTM2D(I+IHW(J),J+1)*HTM2D(I+IHE(J),J+1)
     3    *HTM2D(I-1     ,J)*HTM2D(I+1     ,J)
     4    *HTM2D(I,J-2)*HTM2D(I,J+2).LT.0.5)THEN
          TTV(I,J)=T(I,J,L)*(1.+0.608*Q(I,J,L))
c        if(i.eq.86.and.j.eq.10)print*,'l,TV before relaxation='
c     1, l, TTV(I,J)	
        ENDIF
      ENDDO
      ENDDO
      print*,'convert t to ttv'
C
      KMM=KMNTM(L)
C
      print*,'calling update',l
      DO 285 N=1,NRLX
      CALL UPDATE(TTV)
      if(n.eq.500)print*,'after calling update',l,n
      if(n.eq.500)print*,'l,kmm=',l,kmm
      DO 280 KM=1,KMM
      I=IMNT(KM,L)
      J=JMNT(KM,L)
      TINIT=TTV(I,J)
      TTV(I,J)=AD05*(4.*(TTV(I+IHW(J),J-1)+TTV(I+IHE(J),J-1)
     1                  +TTV(I+IHW(J),J+1)+TTV(I+IHE(J),J+1))
     2                  +TTV(I-1,J)       +TTV(I+1,J)
     3                  +TTV(I,J-2)       +TTV(I,J+2))
     4                  -CFT0*TTV(I,J)
      if(kmm.eq.1.and.n.eq.500)then
c       print*,'l,neighbor points' ,l,I+IHW(J),J-1
c     1,I+IHE(J),J-1,I+IHW(J),J+1,I+IHE(J),J+1 
       print*,'l,n,ttv neighbor',l,n
     1,TTV(I+IHW(J),J-1),TTV(I+IHE(J),J-1),TTV(I+IHW(J),J+1)
     2,TTV(I+IHE(J),J+1),TTV(I-1,J),TTV(I+1,J),TTV(I,J-2)
     3,TTV(I,J+2),TTV(I,J)       
      end if
  280 CONTINUE
C
  285 CONTINUE
c      print*,'after relaxation at level', l
c     1, TTV(I,J) 
C
      DO 290 KM=1,KMM
      I=IMNT(KM,L)
      J=JMNT(KM,L)
      T(I,J,L)=TTV(I,J)
  290 CONTINUE
  300 CONTINUE
c      do j=JSTA_I,JEND_I
c      do i=1,im
c       t00a(i,j)=t(i,j,lsm)
c      end do
c      end do
C----------------------------------------------------------------
C***
C***  CALCULATE THE SEA LEVEL PRESSURE AS PER THE NEW SCHEME.
C***  INTEGRATE THE HYDROSTATIC EQUATION DOWNWARD FROM THE
C***  GROUND THROUGH EACH OUTPUT PRESSURE LEVEL (WHERE TV
C***  IS NOW KNOWN) TO FIND GZ AT THE NEXT MIDPOINT BETWEEN
C***  PRESSURE LEVELS.  WHEN GZ=0 IS REACHED, SOLVE FOR THE
C***  PRESSURE.
C***
C
C***  COUNT THE POINTS WHERE SLP IS DONE BELOW EACH OUTPUT LEVEL
C
      KOUNT=0
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        P1(I,J)=SPL(LMH(I,J))
        DONE(I,J)=.FALSE.
        IF(FIS(I,J).LT.10.)THEN
          PSLP(I,J)=PD(I,J)+PT
          DONE(I,J)=.TRUE.
          KOUNT=KOUNT+1
        ENDIF
      ENDDO
      ENDDO
C
      KMM=KMNTM(LSM)
CC!$omp parallel do private(gz1,gz2,i,j,lmap1,p1,p2),shared(pslp)
      DO 320 KM=1,KMM
      I=IMNT(KM,LSM)
      J=JMNT(KM,LSM)
      LMHIJ=LMH(I,J)
      GZ1=FI(I,J,LMHIJ)
      P1(I,J)=SPL(LMHIJ)
C
      LMAP1=LMHIJ+1
      DO L=LMAP1,LSM
        P2=SPL(L)
        TLYR=0.5*(T(I,J,L)+T(I,J,L-1))
        GZ2=GZ1+RD*TLYR*ALOG(P1(I,J)/P2)
        FI(I,J,L)=GZ2
        IF(GZ2.LE.0.)THEN
          PSLP(I,J)=P1(I,J)/EXP(-GZ1/(RD*T(I,J,L-1)))
          DONE(I,J)=.TRUE.
          KOUNT=KOUNT+1
          GO TO 320
        ENDIF
        P1(I,J)=P2
        GZ1=GZ2
      ENDDO
  320 CONTINUE
C
C***  WHEN SEA LEVEL IS BELOW THE LOWEST OUTPUT PRESSURE LEVEL,
C***  SOLVE THE HYDROSTATIC EQUATION BY CHOOSING A TEMPERATURE
C***  AT THE MIDPOINT OF THE LAYER BETWEEN THAT LOWEST PRESSURE
C***  LEVEL AND THE GROUND BY EXTRAPOLATING DOWNWARD FROM T ON
C***  THE LOWEST PRESSURE LEVEL USING THE DT/DFI BETWEEN THE
C***  LOWEST PRESSURE LEVEL AND THE ONE ABOVE IT.
C
      TOTAL=(IM-2)*(JM-4)
C
      DO 340 LP=LSM,1,-1
      IF(KOUNT.EQ.TOTAL)GO TO 350
      doout330: DO J=JSTA_I,JEND_I
      doin330: DO I=1,IM
      IF(FI(I,J,LP).LT.0..OR.DONE(I,J)) CYCLE doin330
      SLOPE=(T(I,J,LP)-T(I,J,LP-1))/(FI(I,J,LP)-FI(I,J,LP-1))
c     SLOPE=-6.6E-4
      TLYR=T(I,J,LP)-0.5*FI(I,J,LP)*SLOPE
      PSLP(I,J)=P1(I,J)/EXP(-FI(I,J,LP)/(RD*TLYR))
      DONE(I,J)=.TRUE.
      KOUNT=KOUNT+1
      END DO doin330
      END DO doout330
  340 CONTINUE
C
  350 CONTINUE
      print*,'done with mesinger slp reduction'
c      if (me.eq.2)then
c      print*,'writing sample pslp, t, fi'
cc      do j=JSTA_I,JEND_I
c       j=220
c       i=86
c       print*,'i,j,pslp = ',i,j,pslp(i,j)
c       do l=1,lsm
c        print*,'l,SPL,T,FI = ',l,spl(l),t(i,j,l),fi(i,j,l)
c       end do
cc      end do
c      end if 	
C--------------------------------------------------------------------
C     SKIP THE STANDARD SCHEME.
C--------------------------------------------------------------------
      GO TO 430
C--------------------------------------------------------------------
C***
C***  IF YOU WANT THE "STANDARD" ETA/SIGMA REDUCTION
C***  THIS IS WHERE IT IS DONE.
C***
  400 CONTINUE
C
      doout410: DO J=JSTA_I,JEND_I
      doin410: DO I=1,IM
      IF(FIS(I,J).GE.1.)THEN
        LMA=LM
        ALPP1=ALOG(PD(I,J)+PT)
        SLOP=0.0065*ROG*TSIG(I,J,LM)
        IF(SLOP.LT.0.50)THEN
          SLPP=ALPP1+FIS(I,J)/(RD*TSIG(I,J,LMA))
        ELSE
          TTT=-(ALOG(PD(I,J)+PT)+ALPP1)
     1         *SLOP*0.50+TSIG(I,J,LMA)
          SLPP=(-TTT+SQRT(TTT*TTT+2.*SLOP*
     1          (FIS(I,J)/RD+
     2          (TTT+0.50*SLOP*ALPP1)*ALPP1)))/SLOP
        ENDIF
        PSLP(I,J)=EXP(SLPP)
      ENDIF
      END DO doin410
      END DO doout410
  410 CONTINUE
C
C****************************************************************
C     AT THIS POINT WE HAVE A SEA LEVEL PRESSURE FIELD BY
C     EITHER METHOD.  5-POINT AVERAGE THE FIELD ON THE E-GRID.
C****************************************************************
C
  430 CONTINUE
C
C***  EXTRAPOLATE VALUES TO THE OUTER 2 ROWS
C
      print*,'ME in slpsig = ',me
      IF(ME.EQ.0)THEN
      IF((JEND_I-JSTA_I).LT.5)THEN
      DO J=1,2
      IEND=IM-1-MOD(J+1,2)
      DO I=2,IEND
        PSLP(I,J)=PSLP(I,3)
      ENDDO
      ENDDO
      ELSE
      DO J=1,2
      IEND=IM-1-MOD(J+1,2)
      DO I=2,IEND
        PSLP(I,J)=1.5*PSLP(I,J+2)-0.5*PSLP(I,J+4)
      ENDDO
      ENDDO
      END IF
      END IF
C
      IF(ME.EQ.(NUM_PROCS-1))THEN
      IF((JEND_I-JSTA_I).LT.5)THEN
      DO J=JM-1,JM
      IEND=IM-1-MOD(J+1,2)
      DO I=2,IEND
        PSLP(I,J)=PSLP(I,IM-2)
      ENDDO
      ENDDO
      ELSE
      DO J=JM-1,JM
      IEND=IM-1-MOD(J+1,2)
      DO I=2,IEND
        PSLP(I,J)=1.5*PSLP(I,J-2)-0.5*PSLP(I,J-4)
      ENDDO
      ENDDO
      END IF
      END IF
C
      DO J=JSTA_I,JEND_I
        PSLP(1,J)=1.5*PSLP(2,J)-0.5*PSLP(3,J)
      ENDDO
      DO J=JSTA_I,JEND_I
        I=IM-MOD(J+1,2)
        PSLP(I,J)=1.5*PSLP(I-1,J)-0.5*PSLP(I-2,J)
      ENDDO
c      print*,'pslp before smoothing'
c      do j=JSTA_I,JEND_I
c      do i=1,im
c      print*,'i,j,pslp= ',i,j,pslp(i,j)
c      end do
c      end do
C
CC!$omp parallel do 
      SLPX(1:IM,JSTA_I:JEND_I)=PSLP(1:IM,JSTA_I:JEND_I)
C
      DO 480 KS=1,KSLPD
C
      CALL UPDATE(PSLP)
CC!$omp parallel do private(ihh2)
      doout460: DO J=JSTA_IM2,JEND_IM2
      IHH2=IM-1-MOD(J+1,2)
      doin460: DO I=2,IHH2
C
C***  EXTRA AVERAGING UNDER MOUNTAINS TAKEN OUT, FM, MARCH 96
C
      SLPX(I,J)=0.125*(PSLP(I+IHW(J),J-1)+PSLP(I+IHE(J),J-1)
     1                +PSLP(I+IHW(J),J+1)+PSLP(I+IHE(J),J+1)
     2                +4.*PSLP(I,J))
      END DO doin460
      END DO doout460
  460 CONTINUE
C
CC!$omp parallel do
c      print*,'pslp after smoothing'
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        PSLP(I,J)=SLPX(I,J)
c	print*,'i,j,pslp= ',i,j,pslp(i,j)
      ENDDO
      ENDDO
      print*,'sample pslp',pslp(im/2,JSTA_I),pslp(im/2,Jend_I)
C
  480 CONTINUE
C
c      DO L=LHMNT,LSM
c        KMN=KMNTM(L)
c        DO KM=1,KMM
c          I=IMNT(KM,L)
c          J=JMNT(KM,L)
c          T(I,J,L)=T(I,J,L)/(1.+0.608*Q(I,J,L))
c        ENDDO
c      ENDDO
C----------------------------------------------------------------
      RETURN
      END
