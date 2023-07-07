      SUBROUTINE SIG2P(IMOUT,JMOUT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    SIG2P       VERT INTRP OF SIGMA TO PRESSURE
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-09-23       
C     
C ABSTRACT:
C     FOR MOST APPLICATIONS THIS ROUTINE IS THE WORKHORSE
C     OF THE POST PROCESSOR.  IN A NUTSHELL IT INTERPOLATES
C     DATA FROM SIGMA TO PRESSURE SURFACES.  IT ORIGINATED
C     FROM THE VERTICAL INTERPOLATION CODE IN THE OLD ETA
C     POST PROCESSOR SUBROUTINE OUTMAP AND IS A REVISION
C     OF SUBROUTINE ETA2P.
C   .     
C     
C PROGRAM HISTORY LOG:
C   99-09-23  T BLACK       - REWRITTEN FROM ETA2P
C     
C USAGE:    CALL SIG2P(IMOUT,JMOUT)
C   INPUT ARGUMENT LIST:
C     IMOUT    - FIRST DIMENSION OF OUTPUT GRID
C     JMOUT    - SECOND DIMENSION OF OUTPUT GRID
C
C   OUTPUT ARGUMENT LIST: 
C     NONE       
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       SLPSIG   - COMPUTE MESINGER SEA LEVEL PRESSURE
C       SCLFLD   - SCALE ARRAY ELEMENTS BY CONSTANT.
C       E2OUT    - E-GRID TO OUTPUT GRID INTERPOLATION/SMOOTHING.
C       OUTPUT   - POST ARRAY TO OUTPUT FILE.
C       CALPOT2  - COMPUTE POTENTIAL TEMPERATURE.
C       CALRH2   - COMPUTE RELATIVE HUMIDITY.
C       CALDWP2  - COMPUTE DEWPOINT TEMPERATURE.
C       BOUND    - BOUND ARRAY ELEMENTS BETWEEN LOWER AND UPPER LIMITS.
C       CALMCVG  - COMPUTE MOISTURE CONVERGENCE.
C       CALVOR   - COMPUTE ABSOLUTE VORTICITY.
C       CALSTRM  - COMPUTE GEOSTROPHIC STREAMFUNCTION.
C
C     LIBRARY:
C       COMMON   - OMGAOT
C                  LOOPS
C                  MASKS
C                  MAPOT
C                  VRBLS
C                  PVRBLS
C                  RQSTFLD
C                  EXTRA
C                  CLDWTR
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN 90
C     MACHINE : IBM SP
C$$$  
C
C
C     
C     INCLUDE ETA MODEL DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
C     GAMMA AND RGAMOG ARE USED IN THE EXTRAPOLATION OF VIRTUAL
C     TEMPERATURES BEYOND THE UPPER OF LOWER LIMITS OF ETA DATA.
C     
      INCLUDE "parmeta"
      INCLUDE "parmout"
      INCLUDE "params"
C
      PARAMETER (IM_JM=IM*JM,LMP1=LM+1)
      PARAMETER (GAMMA=6.5E-3,RGAMOG=RD*GAMMA/G)
C     
C     DECLARE VARIABLES.
C     
      LOGICAL RUN,FIRST,RESTRT,SIGMA,OLDRD,STDRD
      LOGICAL IOOMG,IOALL
      REAL OSL(IM,JM),USL(IM,JM),VSL(IM,JM)
      REAL PRSI(IM,JM),QCSL(IM,JM),Q2SL(IM,JM)
      REAL ICE(IM,JM),IW(IM,JM,LM),IWU,IWL

      REAL EGRID1(IM,JM),EGRID2(IM,JM)
      REAL GRID1(IMOUT,JMOUT),GRID2(IMOUT,JMOUT)
C
      REAL TPRS(IM,JM,LSL),QPRS(IM,JM,LSL),FPRS(IM,JM,LSL)
C
      INTEGER IHOLD(IM_JM),JHOLD(IM_JM)
C
C     INCLUDE COMMON BLOCKS.
      INCLUDE "CTLBLK.comm"
      INCLUDE "OMGAOT.comm"
      INCLUDE "LOOPS.comm"
      INCLUDE "MASKS.comm"
      INCLUDE "MAPOT.comm"
      INCLUDE "VRBLS.comm"
      INCLUDE "PVRBLS.comm"
      INCLUDE "RQSTFLD.comm"
      INCLUDE "EXTRA.comm"
      INCLUDE "CLDWTR.comm"
      INCLUDE "E2PFLG.comm"
C
      COMMON/JIMA/NL1X(IM,JM),ALPETUX(IM,JM),ALPET2X(IM,JM)
C
C      DATA SIGFLG/.FALSE./
C     
C******************************************************************************
C
C     START SIG2P. 
C     
C     SET TOTAL NUMBER OF POINTS ON OUTPUT GRID.
C
C---------------------------------------------------------------
C
C     *** PART I ***
C
C     VERTICAL INTERPOLATION OF EVERYTHING ELSE.  EXECUTE ONLY
C     IF THERE'S SOMETHING WE WANT.
C
      IF((IGET(012).GT.0).OR.(IGET(013).GT.0).OR.
     X   (IGET(014).GT.0).OR.(IGET(015).GT.0).OR.
     X   (IGET(016).GT.0).OR.(IGET(017).GT.0).OR.
     X   (IGET(018).GT.0).OR.(IGET(019).GT.0).OR.
     X   (IGET(020).GT.0).OR.(IGET(030).GT.0).OR.
     X   (IGET(021).GT.0).OR.(IGET(022).GT.0).OR.
     X   (IGET(153).GT.0).OR.(IGET(166).GT.0))THEN
C
C  SET UP UTIM FOR THIS TIME STEP
C
        UTIM=1.
        CLIMIT=1.E-20
C
        do75: DO L=1,LM
        IF(L.EQ.1)THEN
!$omp  parallel do
          DO J=JSTA,JEND
          DO I=1,IM
            IW(I,J,L)=0.
          ENDDO
          ENDDO
          CYCLE  do75
        ENDIF
C
!$omp  parallel do
!$omp& private(cwmkl,fiq,hh,lml,pp,qi,qkl,qw,tkl,tmt0,tmt15,u00kl)
        doout70: DO J=JSTA,JEND
        doin70: DO I=1,IM
        LML=LM-LMH(I,J)
        HH=HTM(I,J,L)*HBM2(I,J)
        TKL=T(I,J,L)
        QKL=Q(I,J,L)
        CWMKL=CWM(I,J,L)
        TMT0=(TKL-273.16)*HH
        TMT15=AMIN1(TMT0,-15.)*HH    
        PP=PDSL(I,J)*AETA(L)+PT
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
C
C---------------------------------------------------------------------
C***
C***  BECAUSE SIGMA LAYERS DO NOT GO UNDERGROUND, SUBROUTINE BLOSFC2
C***  COULD DO NO WORK TO FILL IN VALUES BELOW THE GROUND SURFACE.
C***  VALUES OF T ON THE OUTPUT PRESSURE LEVELS ABOVE GROUND MUST BE
C***  KNOWN BEFORE THEY CAN BE FILLED IN BELOW GROUND IN THE ROUTINE 
C***  WHICH COMPUTES THE MESINGER SEA LEVEL PRESSURE.  THEREFORE DO ALL
C***  INTERPOLATION ABOVE GROUND NOW.
C***
C
        do310: DO LP=1,LSL
        NHOLD=0
C
        doout125: DO J=JSTA,JEND
        doin125: DO I=1,IM
C
        TSL(I,J)=-1.E6
        QSL(I,J)=-1.E6
        FSL(I,J)=-1.E6
C
C***  LOCATE VERTICAL INDEX OF MODEL INTERFACE JUST BELOW
C***  THE PRESSURE LEVEL TO WHICH WE ARE INTERPOLATING.
C
        do115: DO L=2,LMP1
c       IF((ALPINT(I,J,L)).LT.ALSL(LP))GO TO 115
        IF(ALPINT(I,J,L).GE.ALSL(LP))THEN
          NL1X(I,J)=L
          NHOLD=NHOLD+1
          IHOLD(NHOLD)=I
          JHOLD(NHOLD)=J
          CYCLE doin125
        ELSEIF(ALPINT(I,J,LMP1).LT.ALSL(LP))THEN
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
        IF(NHOLD.EQ.0)CYCLE do310
C
        TRF=2.*ALSL(LP)
C***
C***  NL1X WAS SET BASED ON INTERFACE PRESSURE VALUES.  HOWEVER, WE 
C***  ARE READY TO INTERPOLATE VALUES FROM MIDLAYER LOCATIONS SO 
C***  ADJUST NL1X UP ONE LAYER IF NEEDED SO THAT THE INDEX OF THE
C***  MIDLAYER LOCATION IS IMMEDIATELY BELOW THE PRESSURE LEVEL
C***  TO WHICH INTERPOLATION IS BEING DONE.
C***
!$omp  parallel do
!$omp& private(ahf,ahfo,ahfq,ahfq2,ahfqc,ahfqi,ai,b,bi,bom,
!$omp&         bq,bq2,bqc,bqc_2,bqi_2,bqi,fac,gmiw,gmiw_2,iwl,iwu,
!$omp&         lmp1,pl,pnl1,pu,q2a,q2b,qabv,qi,qint,ql,
!$omp&         qsat,qu,qw,rhu,tabv,tblo,tl,tmt0,tmt15,tu,
!$omp&         tvrabv,tvrblo,tvrl,tvru,zl,zu)
        do220: DO NN=1,NHOLD
        I=IHOLD(NN)
        J=JHOLD(NN)
        DIFU=ALSL(LP)-ALPINT(I,J,NL1X(I,J)-1)
        DIFL=ALPINT(I,J,NL1X(I,J))-ALSL(LP)
        IF(DIFU.LT.DIFL)NL1X(I,J)=NL1X(I,J)-1
        PNL1=PINT(I,J,NL1X(I,J))
C---------------------------------------------------------------------
C***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
C***  HUMIDITY, CLOUD WATER/ICE, AND TKE.
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
C***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
C---------------------------------------------------------------------
C
        IF(NL1X(I,J).EQ.1)THEN
          PU=PINT(I,J,2)
          ZU=ZINT(I,J,2)
          TU=0.5*(T(I,J,1)+T(I,J,2))
          QU=0.5*(Q(I,J,1)+Q(I,J,2))
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
          BOM  =OMGA(I,J,1)
          GMIW =IW(I,J,1) 
          BQC  =(1.-GMIW)*CWM(I,J,1)
          BQI  =GMIW*CWM(I,J,1)
          Q2A  =0.5*(Q2(I,J,1)+Q2(I,J,2))
          BQ2  =Q2A
          AHF  =0.
          AHFQ =0.
          AHFO =0.
          AHFQC=0.
          AHFQI=0.
          AHFQ2=0.
          FAC  =0.
C
        ELSEIF(NL1X(I,J).EQ.LMP1)THEN 
C---------------------------------------------------------------------
C***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
C---------------------------------------------------------------------
C
          PL=PINT(I,J,LM-1)
          ZL=ZINT(I,J,LM-1)
          TL=0.5*(T(I,J,LM-2)+T(I,J,LM-1))
          QL=0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
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
          BOM  =OMGA(I,J,LM)
          GMIW =IW(I,J,LMH(I,J))
          BQC  =(1.-GMIW)*CWM(I,J,LMH(I,J))
          BQI  =GMIW*CWM(I,J,LMH(I,J)) 
          Q2A  =0.5*(Q2(I,J,LMH(I,J)-1)+Q2(I,J,LMH(I,J)))
          BQ2  =Q2A
          AHF  =0.
          AHFQ =0.
          AHFO =0.
          AHFQC=0.
          AHFQI=0.
          AHFQ2=0.
          FAC  =0.
C
        ELSE
C---------------------------------------------------------------------
C***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
C---------------------------------------------------------------------
C
          B     =T(I,J,NL1X(I,J))
          BQ    =Q(I,J,NL1X(I,J))
          BOM   =OMGA(I,J,NL1X(I,J))
          GMIW  =IW(I,J,NL1X(I,J))
          BQC   =(1.-GMIW)*CWM(I,J,NL1X(I,J))
          BQI   =GMIW*CWM(I,J,NL1X(I,J)) 
          GMIW_2=IW(I,J,NL1X(I,J)-1)
          BQC_2 =(1.-GMIW)*CWM(I,J,NL1X(I,J)-1)
          BQI_2 =GMIW*CWM(I,J,NL1X(I,J)-1)
          Q2B   =0.5*(Q2(I,J,NL1X(I,J)-1)+Q2(I,J,NL1X(I,J)))
C
          IF(NL1X(I,J).GT.2)THEN
            Q2A=0.5*(Q2(I,J,NL1X(I,J)-2)+Q2(I,J,NL1X(I,J)-1))
          ELSE
            Q2A=Q2B
          ENDIF
C
          BQ2=Q2B*HTM(I,J,NL1X(I,J))
          FAC  =2.*ALOG(PT+PDSL(I,J)*AETA(NL1X(I,J)))
          AHF  =(B-T(I,J,NL1X(I,J)-1))/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
          AHFQ =(BQ-Q(I,J,NL1X(I,J)-1))/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
          AHFO =(BOM-OMGA(I,J,NL1X(I,J)-1))/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
          AHFQC=(BQC-BQC_2)/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
          AHFQI=(BQI-BQI_2)/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
          AHFQ2=(BQ2-Q2A*HTM(I,J,NL1X(I,J)-1))/
     1          (ALPINT(I,J,NL1X(I,J)+1)-ALPINT(I,J,NL1X(I,J)-1))
        ENDIF
C
        TSL(I,J)=B+AHF*(TRF-FAC)
        QSL(I,J)=BQ+AHFQ*(TRF-FAC)
        QSL(I,J)=AMAX1(QSL(I,J),H1M12)
        OSL(I,J)=BOM+AHFO*(TRF-FAC)
        QCSL(I,J)=BQC+AHFQC*(TRF-FAC)
        QCSL(I,J)=AMAX1(QCSL(I,J),H1M12)
        ICE(I,J)=BQI+AHFQI*(TRF-FAC)
        ICE(I,J)=AMAX1(ICE(I,J),1.E-12)
        Q2SL(I,J)=BQ2+AHFQ2*(TRF-FAC)
        Q2SL(I,J)=AMAX1(Q2SL(I,J),0.)
        FSL(I,J)=(PNL1-SPL(LP))/(SPL(LP)+PNL1)
     1       *((ALSL(LP)+ALPINT(I,J,NL1X(I,J))-FAC)*AHF+B)*R*2.
     2       +ZINT(I,J,NL1X(I,J))*G
        END DO do220
C
c 225   CONTINUE
C---------------------------------------------------------------------
C        LOAD GEOPOTENTIAL AND TEMPERATURE INTO STANDARD LEVEL 
C        ARRAYS FOR THE NEXT PASS.
C---------------------------------------------------------------------
C     
C***     SAVE 500MB TEMPERATURE FOR LIFTED INDEX.
C     
        IF((NINT(SPL(LP)).EQ.50000).AND.
     1    ((IGET(030).GT.0).OR.(IGET(031).GT.0).OR.
     2                         (IGET(075).GT.0)))THEN
!$omp  parallel do private(i,j)
          DO NN=1,NHOLD
            I=IHOLD(NN)
            J=JHOLD(NN)
            T500(I,J)=TSL(I,J)
          ENDDO
        ENDIF
C     
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C***  INTERPOLATE WIND COMPONENTS FROM ETA TO PRESSURE
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C     
        IF((IGET(018).GT.0).OR.(IGET(019).GT.0).OR.
     1     (IGET(021).GT.0).OR.(IGET(085).GT.0))THEN
C
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
            USL(I,J)=0.
            VSL(I,J)=0.
          ENDDO
          ENDDO
C
!!$omp  parallel do
!!$omp& private(alpet2,alpetl,alpetu,lmb,petal,petau)
!!$omp& shared (alpet2x,alpetux,nl1x)
          doout280: DO J=JSTA,JEND
          doin1280: DO  I=1,IM
CNOTE 
CNOTE         29 JANUARY 1993, RUSS TREADON.
CNOTE          - AS FOR THE OTHER FIELDS WE INTERPOLATE ONLY
CNOTE            BETWEEN THE FAL AND THE MODEL TOP.  BELOW 
CNOTE            SURFACE VALUES ARE FAL VALUES.
C
          LMVIJ=LMV(I,J)
C
cccccccc  DO 280 LV=2,LMVIJ
          doin2280: DO LV=2,LMVIJ+1
          PETAL=PT+PDVP1(I,J)*ETA(LV)
          PETAU=PT+PDVP1(I,J)*ETA(LV-1)
          ALPETL=ALOG(PETAL)
          ALPETU=ALOG(PETAU)
          ALPET2=SQRT(0.5*(ALPETL*ALPETL+ALPETU*ALPETU))
C     
C***  SEARCH FOR HIGHEST MID-LAYER MODEL SURFACE
C***  THAT IS BELOW THE GIVEN STANDARD PRESSURE LEVEL.
C
          IF(ALSL(LP).LT.ALPET2)THEN
            NL1X(I,J)=LV-1
            ALPETUX(I,J)=ALPETU
            ALPET2X(I,J)=ALPET2
            GO TO 281
          ENDIF
          END DO doin2280
          END DO doin1280
          END DO doout280
C
          NL1X(I,J)=LMVIJ+1
          ALPETUX(I,J)=ALPETU
          ALPET2X(I,J)=ALPET2
  281     CONTINUE
C     
C---------------------------------------------------------------------
C***      BELOW GROUND USE WIND IN LOWEST LAYER
C---------------------------------------------------------------------
C
!$omp  parallel do
!$omp& private(alpet1,alpetl,alpetu,fact,lmvij,petau)
          doout290: DO J=JSTA,JEND
          doin290: DO I=1,IM
          LMVIJ=LMV(I,J)
          IF(NL1X(I,J).GT.LMVIJ)THEN
            USL(I,J)=U(I,J,LMVIJ)
            VSL(I,J)=V(I,J,LMVIJ)
C     
C---------------------------------------------------------------------
C***  IF REQUESTED PRESSURE LEVEL IS NOT BELOW THE LOCAL GROUND
C***  THEN WE HAVE TWO POSSIBILITIES.  IF THE REQUESTED PRESSURE
C***  LEVEL IS BETWEEN THE LOCAL SURFACE PRESSURE AND TOP OF
C***  MODEL PRESSURE, VERTICALLY INTERPOLATE BETWEEN NEAREST
C***  BOUNDING ETA LEVELS TO GET THE WIND COMPONENTS.  IF THE   
C***  REQUESTED PRESSURE LEVEL IS ABOVE THE MODEL TOP, USE
C***  CONSTANT EXTRAPOLATION OF TOP ETA LAYER (L=1) WINDS.
C---------------------------------------------------------------------
C 
          ELSE
            IF(NL1X(I,J).GT.1)THEN
              ALPETL=ALPETUX(I,J)
              PETAU=PT+PDVP1(I,J)*ETA(NL1X(I,J)-1)
              ALPETU=ALOG(PETAU)
              ALPET1=SQRT(0.5*(ALPETL*ALPETL+ALPETU*ALPETU))
              FACT=(ALPET2X(I,J)-ALSL(LP))/(ALPET2X(I,J)-ALPET1)
              USL(I,J)=U(I,J,NL1X(I,J))
     1               +(U(I,J,NL1X(I,J)-1)-U(I,J,NL1X(I,J)))*FACT
              VSL(I,J)=V(I,J,NL1X(I,J))+(V(I,J,NL1X(I,J)-1)
     1                -V(I,J,NL1X(I,J)))*FACT
            ELSE
              USL(I,J)=U(I,J,NL1X(I,J))
              VSL(I,J)=V(I,J,NL1X(I,J))
            ENDIF
C     
C---------------------------------------------------------------------
C***  ALPET2 IS MID-LAYER ETA SURFACE JUST BELOW STANDARD PRESSURE
C***  LEVEL AND ALPET1 IS MID-LAYER ETA SURFACE JUST ABOVE.
C***  NOTE THAT IF THE STANDARD PRESSURE SURFACE IS SUBMERGED, THEN
C***  ALPET2 AND ALPET1 ARE THE LOWEST AND 2ND LOWEST MID-LAYER
C***  ETA SURFACES ABOVE THE TOPOGRAPHY (WITH OLDRD=.TRUE., ZJ).
C---------------------------------------------------------------------
C
          ENDIF
          END DO doin290
          END DO doout290
C
        ENDIF
C
C---------------------------------------------------------------------
C
C***  FILL THE 3-D-IN-PRESSURE ARRAYS FOR THE SLP REDUCTION
C
        DO J=JSTA,JEND
        DO I=1,IM
          TPRS(I,J,LP)=TSL(I,J)
          QPRS(I,J,LP)=QSL(I,J)
          FPRS(I,J,LP)=FSL(I,J)
        ENDDO
        ENDDO
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C        *** PART II ***
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C        INTERPOLATE/OUTPUT SELECTED FIELDS.
C
C---------------------------------------------------------------------
C     
C***  SPECIFIC HUMIDITY.
C
        IF(IGET(016).GT.0)THEN
          IF(LVLS(LP,IGET(016)).GT.0)THEN
            CALL E2OUT(016,000,QSL,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1M12,H99999,IMOUT,JMOUT)
            ID(1:25)=0
            print*,'calling output humidity from sig2p'
            CALL OUTPUT(IOUTYP,IGET(016),LP,GRID1,IMOUT,JMOUT)
            if(lp.eq.1)print*,'sample SH IOUTYP',IOUTYP,IGET(016)
          ENDIF
        ENDIF
C     
C***  OMEGA
C
        IF(IGET(020).GT.0)THEN
          IF(LVLS(LP,IGET(020)).GT.0)THEN
            CALL E2OUT(020,000,OSL,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            print*,'calling output omeg from sig2p'
            CALL OUTPUT(IOUTYP,IGET(020),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  MOISTURE CONVERGENCE
C
        IF(IGET(085).GT.0)THEN
          IF(LVLS(LP,IGET(085)).GT.0)THEN
            CALL CALMCVG(QSL,USL,VSL,-1,EGRID1)
            CALL E2OUT(085,000,EGRID1,EGRID2,
     1                 GRID1,GRID2,IMOUT,JMOUT)
C
C     CONVERT TO DIVERGENCE FOR GRIB UNITS
C
            CALL SCLFLD(GRID1,-1.0,IMOUT,JMOUT)
            ID(1:25)=0
            print*,'calling output moisture from sig2p'
            CALL OUTPUT(IOUTYP,IGET(085),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  U AND/OR V WIND
C
        IF(IGET(018).GT.0.OR.IGET(019).GT.0)THEN
          IF(LVLS(LP,IGET(018)).GT.0.OR.LVLS(LP,IGET(019)).GT.0)THEN
            CALL E2OUT(018,019,USL,VSL,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            IF(IGET(018).GT.0)THEN
              CALL OUTPUT(IOUTYP,IGET(018),LP,GRID1,IMOUT,JMOUT)
            ENDIF
            ID(1:25)=0
            IF(IGET(019).GT.0) 
     1       CALL OUTPUT(IOUTYP,IGET(019),LP,GRID2,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  ABSOLUTE VORTICITY
C
         IF (IGET(021).GT.0) THEN
          IF (LVLS(LP,IGET(021)).GT.0) THEN
            CALL CALVOR(USL,VSL,EGRID1)
            CALL E2OUT(021,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(021),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C     
C        GEOSTROPHIC STREAMFUNCTION.
         IF (IGET(086).GT.0) THEN
          IF (LVLS(LP,IGET(086)).GT.0) THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=FSL(I,J)*GI
            ENDDO
            ENDDO
            CALL CALSTRM(EGRID2,EGRID1)
            CALL E2OUT(086,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(086),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C     
C***  TURBULENT KINETIC ENERGY
C
         IF (IGET(022).GT.0) THEN
          IF (LVLS(LP,IGET(022)).GT.0) THEN
            CALL E2OUT(022,000,Q2SL,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(022),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C     
C***  TOTAL CLOUD WATER
C
         IF (IGET(153).GT.0) THEN
          IF (LVLS(LP,IGET(153)).GT.0) THEN
            CALL E2OUT(153,000,QCSL,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1M12,H99999,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(153),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF
C
C***  TOTAL CLOUD ICE 
C
         IF (IGET(166).GT.0) THEN
          IF (LVLS(LP,IGET(166)).GT.0) THEN
            CALL E2OUT(166,000,ICE,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1M12,H99999,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(166),LP,GRID1,IMOUT,JMOUT)
          ENDIF
         ENDIF

C     
C***  END OF MAIN VERTICAL LOOP
C     
        END DO do310
C---------------------------------------------------------------------
C***
C***  NOW CALL THE SEA LEVEL PRESSURE ROUTINE WHICH WILL ALSO
C***  GENERATE THE UNDERGROUND TEMPERATURES.  IF IT HAS ALREADY
C***  BEEN CALLED IN A PREVIOUS PASS THROUGH SIG2P, THEN WE
C***  DO NOT NEED TO CALL IT AGAIN.
C***
C
C        IF(.NOT.SIGFLG)THEN
C          CALL SLPSIG(PD,FIS,SM,TPRS,QPRS,FPRS,IM,JM,SPL,LSL
C     1,               DETA,PT,PSLP)
C          SIGFLG=.TRUE.
C        ENDIF
C
C***  OUTPUT SEA LEVEL PRESSURE IF REQUESTED.
C***  FIRST, MESINGER'S SEA LEVEL PRESSURE.
C
cc        IF(IGET(023).GT.0)THEN
cc          CALL E2OUT(023,000,PSLP,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
cc          ID(1:25)=0
cc          CALL OUTPUT(IOUTYP,IGET(023),LVLS(1,IGET(023)),
cc     1                GRID1,IMOUT,JMOUT)
cc          print*,'sample PSLP IOUTYP',IOUTYP,IGET(023)
cc     1    ,LVLS(1,IGET(023)) 
cc        ENDIF
C
C***  SECOND, STANDARD NGM SEA LEVEL PRESSURE.
C
        IF(IGET(105).GT.0)THEN
          CALL NGMSLP2
          CALL E2OUT(105,000,SLP,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
          ID(1:25)=0
          CALL OUTPUT(IOUTYP,IGET(105),LVLS(1,IGET(105)),
     1        GRID1,IMOUT,JMOUT)
        ENDIF
C
C---------------------------------------------------------------------
C***  OUTPUT THE TEMPERATURES AND QUANTITIES DERIVED FROM IT
C---------------------------------------------------------------------
C
        do400: DO LP=1,LSL
C     
C***  TEMPERATURE
C
        IF(IGET(013).GT.0) THEN
          IF(LVLS(LP,IGET(013)).GT.0)THEN
             CALL E2OUT(013,000,TPRS(1,1,LP),EGRID2,GRID1,GRID2
     1,                 IMOUT,JMOUT)
             ID(1:25)=0
             CALL OUTPUT(IOUTYP,IGET(013),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  POTENTIAL TEMPERATURE.
C
        IF(IGET(014).GT.0)THEN
          IF(LVLS(LP,IGET(014)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=SPL(LP)
            ENDDO
            ENDDO
C
            CALL CALPOT2(EGRID2,TPRS(1,1,LP),EGRID1,IM,JM)
            CALL E2OUT(014,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(014),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  RELATIVE HUMIDITY.
C
        IF(IGET(017).GT.0)THEN
          IF(LVLS(LP,IGET(017)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=SPL(LP)
            ENDDO
            ENDDO
C
            CALL CALRH2(EGRID2,TPRS(1,1,LP),QPRS(1,1,LP),ICE,EGRID1
     1,                 IM,JM)
            CALL E2OUT(017,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            CALL SCLFLD(GRID1,H100,IMOUT,JMOUT)
            CALL BOUND(GRID1,H1,H100,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(017),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C***  DEWPOINT TEMPERATURE.
C
        IF(IGET(015).GT.0)THEN
          IF(LVLS(LP,IGET(015)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID2(I,J)=SPL(LP)
            ENDDO
            ENDDO
C
            CALL CALDWP2(EGRID2,QPRS(1,1,LP),EGRID1,TPRS(1,1,LP))
            CALL E2OUT(015,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(015),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C     
C---------------------------------------------------------------------
C***  CALCULATE 1000MB GEOPOTENTIALS CONSISTENT WITH SLP OBTAINED 
C***  FROM THE MESINGER OR NWS SHUELL SLP REDUCTION.
C---------------------------------------------------------------------
C     
C***  FROM MESINGER SLP
C
        IF(IGET(023).GT.0.AND.LP.EQ.LSM)THEN
          ALPTH=ALOG(1.E5)
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
            PSLPIJ=PSLP(I,J)
            ALPSL=ALOG(PSLPIJ)
            PSFC=PD(I,J)+PT
            IEND=IM-MOD(J+1,2)
            IF(ABS(PSLPIJ-PSFC).LT.5.E2.OR.
     1        ((I.EQ.1.OR.I.LE.IEND.OR.J.LE.2.OR.J.GE.JM-1)
     2         .AND.SM(I,J).GT.0.5))THEN
              FPRS(I,J,LP)=R*TPRS(I,J,LSL)*(ALPSL-ALPTH)
            ELSE
              FPRS(I,J,LP)=FIS(I,J)/(ALPSL-ALOG(PD(I,J)+PT))*
     1                              (ALPSL-ALPTH)
            ENDIF
            Z1000(I,J)=FPRS(I,J,LP)*GI
          ENDDO
          ENDDO
C     
C***  FROM NWS SHUELL SLP. NGMSLP2 COMPUTES 1000MB GEOPOTENTIAL.
C
        ELSEIF(IGET(023).LE.0.AND.LP.EQ.LSM)THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
            FPRS(I,J,LP)=Z1000(I,J)*G
          ENDDO
          ENDDO
        ENDIF
C
C***  OUTPUT GEOPOTENTIAL (SCALE BY GI)
C
        IF(IGET(012).GT.0)THEN
          IF(LVLS(LP,IGET(012)).GT.0)THEN
!$omp  parallel do
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID1(I,J)=FPRS(I,J,LP)*GI
            ENDDO
            ENDDO
C
            CALL E2OUT(012,000,EGRID1,EGRID2,GRID1,GRID2,IMOUT,JMOUT)
            ID(1:25)=0
            CALL OUTPUT(IOUTYP,IGET(012),LP,GRID1,IMOUT,JMOUT)
          ENDIF
        ENDIF
C
        END DO do400
C
        IOALL=.TRUE.
C
C***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
C
      ENDIF
C     
C     END OF ROUTINE.
C
      RETURN
      END
