!
      PROGRAM POST0
!   
!     ##################################################################
!     #                                                                #    
!     #    LAST UPDATE: 2011-06-09 BY ISRAEL BORGES                    #
!     #                               LUIS THIAGO LUCCI                #
!     #                               JORGE GOMES                      # 
!     #                                                                #
!     #    MAIN PROGRAM DOCUMENTATION BLOCK                            #
!     #                                                                #
!     #    MAIN PROGRAM:  POST0       BCEX FOR NESTS, & PROFILES       #
!     #    PRGRMMR: BLACK         ORG: W/NP22     DATE: 99-04-21       #
!     #                                                                #
!     #    ABSTRACT:THIS PROGRAM PRODUCES THE RAW PROFILE FILES        #
!     #             AND EXTRACTS THE BOUNDARIES FOR THE NESTED         #
!     #             DOMAINS.  THESE JOBS WERE PREVIOUSLY DONE          #
!     #             INSIDE THE ETA MODEL FORECAST CODE IN              #
!     #             SUBROUTINE CHKOUT.                                 #
!     #                                                                #
!     #    PROGRAM HISTORY LOG:                                        #
!     #    99-04-21  T BLACK  - EXTRACTED THE RELEVANT PARTS OF CHKOUT #
!     #    11-06-09  I BORGES - DYNAMIC SCHEME TO OPEN RESTRT FILES    #
!     #              T LUCCI    DECLARATIONS CONVERTED TO FORTRAN90    #
!     #              J GOMES    REMOVED REDUNDANT CALLS OF SUBROUTINES #
!     #                                                                #
!     #    USAGE:  MAIN PROGRAM                                        #
!     #                                                                #
!     #    INPUT ARGUMENT LIST: NONE                                   #
!     #                                                                #
!     #    OUTPUT ARGUMENT LIST: NONE                                  #
!     #                                                                #
!     #    INPUT FILES: NONE                                           #
!     #                                                                #
!     #    OUTPUT FILES:                                               #
!     #                                                                #
!     #    SUBPROGRAMS CALLED: NONE                                    #
!     #                                                                #
!     #    EXIT STATES:                                                #
!     #    COND =   0 - NORMAL EXIT                                    #
!     #                                                                #
!     #    ATTRIBUTES:                                                 #
!     #    LANGUAGE: FORTRAN 90                                        #
!     #    MACHINE : IBM SP                                            #
!     #                                                                #
!     ##################################################################
!
!     ##################################################################     
!     #                                                                #
!     # PARMS FOR HOURLY PROFILER OUTPUT                               #
!     # NSTAT - MAX NUMBER OF STATIONS                                 #     
!     # NWORDM - DIMENSION OF OUTPUT ARRAY, MUST BE LARGE ENOUGH       #
!     # TO HOLD ALL VARIABLES                                          #
!     # (MAX NO MULTI-LAYER VARIABLES*LM + NO OF SINGLE LAYER VARS)    #
!     # LCL1ML - NUMBER OF MULTI-LAYER VARIABLES OUTPUT FOR CLASS 1    #
!     #	LCL1SL - NUMBER OF SINGLE LAYER VARIABLES OUTPUT FOR CLASS 1   #
!     #                                                                #
!     ##################################################################
!
      IMPLICIT NONE
!
!     ##################################################################
!     #                                                                #
!     #                   IMPLICIT NONE VARIABLES                      #
!     #                                                                #     
!     ##################################################################      
!
      INTEGER ::                                                       &
      IM,JM,LM,LSM,NSOIL,NROOT,ITAG,KOUNT,NINC,LCLAS1,IER,NUMSTA,N,NHB,&
      NFCST,NBC,LIST,L,NTSPH,ISLP,TMAX,JUMPMAX,NK,KPATH,KRST,LRSTRT,   &
      JUMP,LREC0,LREC,NKRST,IHRST,NTSD,NS,NHEAT,NPHS,NCNVC,NPREC,NRDSW,&
      NRDLW,NSRFC,IFHR,LML,IWKL,IHR,IYR,IMNTH,IDAY,IFCST,IH,JH,LMHK,   &
      NWORD2,NWORD3,NWORD4,NWORD5,NWORD6,NWORD7,NWORD8,NWORD9,NWORD10, &
      NWORD11,NWORD12,NWORD13,NWORD14,NWORD15,ISTAT,LVL,LV,LL,NN,NLEN, &
      NL,NREC      
!
      REAL  ::                                                         &
      DT,DY,CPGFV,EN,ENT,R,PT,TDDAMP,F4D,F4Q,EF4T,R1,PT1,TSPH,WBD,     &
      SBD,TLM0D,TPH0D,ACUTIM,ARATIM,APHTIM,UTIM,US,CCLIMIT,CLIMIT,     &
      HH,TKL,QKL,CWMKL,TMT0,TMT15,AI,BI,PP,QW,QI,QINT,U00KL,FIQ,       &
      FIW,QC,RQKL,ARG,TIME,RESET0,RESET1,                              &
      SINPH0,COSPH0,APEL,RTSPH,RTSCU,RTSRA,DLM,XX,YY,TLON,             &
      ALPHA,SINALP,COSALP,UT,VT,STAPRX,STACRX,PSFCEVP,PPOTEVP,PSFCSHX, &
      PSFCSUB,PSNOPCX,PRSWIN,PRSWOUT,PRLWIN,PRLWOUT,PRLWTOA,PRSWTOA,   &
      PACSNOW,PACSNOM,PSSROFF,PBGROFF	                   
!
!     ##################################################################
!
      INCLUDE "parmeta"
      INCLUDE "parmsoil"
!
!     ##################################################################
!
      INTEGER, PARAMETER ::                                            &                                            
      NSTAT=2000,LCL1ML=15,LCL1SL=50,NWORDM=(LCL1ML+1)*LM+2*LCL1SL,    &
      LRECPR=4*(8+9+LCL1ML*LM+LCL1SL),ITB=76,JTB=134
!
!     ##################################################################
!
      REAL, PARAMETER ::                                               &
      A2=17.2693882,A3=273.16,A4=35.86,PQ0=379.90516,DTR=1.74532925E-2,&
      GI=1./9.8,RD=287.04,CP=1004.6,CAPA=RD/CP
!
!     ##################################################################
!
      INTEGER, DIMENSION(NSTAT) ::                                     &
      IDSTN,IHINDX,JHINDX,IVINDX,JVINDX   
!      
      INTEGER, DIMENSION(IM,JM) ::                                     &
      IDUM,LMH
!
      INTEGER ::                                                       &
      IW(NSTAT,LM),IDAT(3),STATV(13) ! IDAT0(3)     	
!
!     ##################################################################
!
      REAL, DIMENSION(NSTAT) ::                                        &
      STNLAT,STNLON,RES,FIS,THS,HBOT,CFRACL,CFRACM,CFRACH,SNO,         &
      SOILTB,SFCEXC,SMSTAV,SMSTOT,Z0,CZEN,CZMEAN,U00,SR,ACPREC,        &
      CUPREC,ACSNOW,ACSNOM,SSROFF,BGROFF,SFCSHX,SFCLHX,SUBSHX,         &
      SNOPCX,ASWIN,ASWOUT,ASWTOA,ALWIN,ALWOUT,ALWTOA,TSHLTR,QSHLTR,    &
      TH10,Q10,U10,V10,TLMIN,TLMAX,CMC,VEGFRC,POTFLX,PSLP,PDSL1,       &
      EGRID2,SM,SICE,HBM2,FACTR,STATPR,STACPR,STAEVP,STAPOT,STASHX,    &
      STASUB,STAPCX,STASWI,STASWO,STALWI,STALWO,STALWT,STASWT,STASNM,  &
      STASRF,STABRF,STASNO,ACPREC0,CUPREC0,SFCLHX0,POTFLX0,SFCSHX0,    &
      SUBSHX0,SNOPCX0,ASWIN0,ASWOUT0,ALWIN0,ALWOUT0,ALWTOA0,ASWTOA0,   &
      ACSNOW0,ACSNOM0,SSROFF0,BGROFF0	
!
!     ##################################################################
!     
      REAL, DIMENSION(LM) ::                                           &
      DETA,RDETA,AETA,STADHC,STADHR
!
      REAL, DIMENSION(NWORDM) ::                                       &
      PRODAT,FPACK
!
      REAL, DIMENSION (IM,JM) ::                                       &
      PD,PDS
!
      REAL, DIMENSION(NSTAT,NSOIL) ::                                  &
      SMC,STC,SH2O
!
      REAL, DIMENSION(NSTAT,LM) ::                                     &
      T,Q,U,V,Q2,OMGALF,CWM,TRAIN,TCUCN,RSWTT,RLWTT,CCR,RTOP,HTM,      &
      TCUCN0,TRAIN0
!
      REAL, DIMENSION(LM,NSTAT) ::                                     &
      DHCNVC,DHRAIN
!
      REAL ::                                                          &
      UL(2*LM),PTBL(ITB,JTB),TTBL(JTB,ITB)
!
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: DUM
!
!     ##################################################################
!
      LOGICAL ::                                                       &
      RUN,RESTRT,FILEFOUND
!
!     ##################################################################
!
      CHARACTER ::                                                     &
      RSTFIL*80,NRSTFIL*80,RESTHR*4,LABEL*32,CISTAT*8,CIDSTN(NSTAT)*8, &
      FNAME*80,NFNAME*80,ENVAR*80,BLANK*4
!
!     ##################################################################
!
      DATA BLANK/'    '/
      DATA NHB/12/,LRSTRT/22/,LCLAS1/76/
!
!     ##################################################################
!
      REAL ::                                                          &
      TSTART,TEND,SPL(LSM),TPREC,THEAT,TCLOD,TRDSW,TRDLW,TSRFC
!
      INTEGER ::                                                       &
      TCP,NMAP,TSHDE(999999),NRADSH,NRADLH,NTDDMP
!
      LOGICAL ::                                                       &
      SINGLRST,SUBPOST,NEST,HYDRO,SPLINE
! 
      NAMELIST/FCSTDATA/                                               &
      TSTART,TEND,TCP,RESTRT,SINGLRST,SUBPOST,NMAP,TSHDE,SPL,NPHS,     &
      NCNVC,NRADSH,NRADLH,NTDDMP,TPREC,THEAT,TCLOD,TRDSW,TRDLW,TSRFC,  &
      NEST,HYDRO,SPLINE
!
!     ##################################################################
!     #                                                                #
!     # READ NAMELIST FCSTDATA WHICH CONTROLS TIMESTEPS,               #
!     # ACCUMULATION PERIODS, AND STANDARD OUTPUT                      #
!     #                                                                #
!     ##################################################################
!
      OPEN(11,FILE='fcstdata.meso',STATUS='OLD')
      REWIND 11
      READ(11,FCSTDATA)
      CLOSE(11)
!
!     ##################################################################     
!     #                                                                #
!     # READ IN THE FIRST FORECAST HOUR THAT THIS JOB WILL WORK ON     #
!     # PLUS THE NUMBER OF OUTPUT TIMES TO PROCESS AND THE TIME        #     
!     # INCREMENT IN HOURS BETWEEN PROCESS TIMES                       #
!     #                                                                #
!     ##################################################################
!
      READ(5,*) ITAG,KOUNT,NINC
!      
      OPEN(UNIT=LCLAS1,ACCESS='DIRECT',RECL=LRECPR,IOSTAT=IER,         &
           FILE='profilm.c1.t00s')
!
!     ##################################################################
!     #                                                                #
!     #        READ IN THE INFORMATION FILE ABOUT THE SOUNDINGS        #
!     #                                                                #     
!     ##################################################################      
!    
      REWIND 19
!
      READ(19)NUMSTA,IDSTN,STNLAT,STNLON,                              &
      IHINDX,JHINDX,IVINDX,JVINDX,CIDSTN
!
!      WRITE(6,20)NUMSTA
!   20 FORMAT('INIT:  NUMBER OF PROFILE STATIONS ',I5)
!
!      WRITE(6,30)(IDSTN(N),STNLAT(N)/DTR,STNLON(N)/DTR,                &
!      IHINDX(N),JHINDX(N),IVINDX(N),JVINDX(N),                         &
!      CIDSTN(N),N=1,NUMSTA)
!   30 FORMAT(2X,I6,2F8.2,4I8,4X,A8)
!
!     ##################################################################
!
      ALLOCATE(DUM(IM,JM,8))
!
!     ##################################################################
!     #                                                                #
!     #            READ QUANTITIES NEEDED FROM THE NHB FILE            #
!     #                                                                #     
!     ##################################################################      
! 
      OPEN(UNIT=NHB,FILE='cnst.file',ACCESS='SEQUENTIAL',              &
           FORM='UNFORMATTED')
!
      REWIND NHB
!
      READ(NHB) NFCST,NBC,LIST,DT
!
!      WRITE(6,*) 'NFCST,NBC,LIST,DT ', NFCST,NBC,LIST,DT
!
      READ(NHB)LMH
      READ(NHB)
      READ(NHB) DUM(:,:,1) 
!
      DO N=1,NUMSTA
       HBM2(N)=DUM(IHINDX(N),JHINDX(N),1)
      ENDDO
!
      READ(NHB)
      READ(NHB)
      READ(NHB) DUM(:,:,1) 
!
      DO N=1,NUMSTA
       SM(N)=DUM(IHINDX(N),JHINDX(N),1)
      ENDDO
!
      READ(NHB) DUM(:,:,1) 
!
      DO N=1,NUMSTA
       SICE(N)=DUM(IHINDX(N),JHINDX(N),1)
      ENDDO
!
      DO L=1,LM
       READ(NHB) DUM(:,:,1)
      DO N=1,NUMSTA
       HTM(N,L)=DUM(IHINDX(N),JHINDX(N),1)
      ENDDO     
      ENDDO
!
      DO L=1,LM
       READ(NHB)
      ENDDO
!
      READ(NHB)DY,CPGFV,EN,ENT,R,PT,TDDAMP,                            &
      F4D,F4Q,EF4T,DETA,RDETA,AETA
!
      NTSPH=INT(3600./DT+0.50)
!
      DO L=1,23
       READ(NHB)
      ENDDO
!
      READ(NHB) PTBL,TTBL,R1,PT1,TSPH,WBD,SBD,TLM0D,TPH0D
      READ(NHB)
      READ(NHB)
      READ(NHB)
!
      READ(NHB) DUM(:,:,1) 
!      
      DO N=1,NUMSTA
       VEGFRC(N)=DUM(IHINDX(N),JHINDX(N),1)
      ENDDO
!
!     ##################################################################
!
      DEALLOCATE(DUM)
!
!     ##################################################################
!     #                                                                #
!     #                 LOOP OVER ALL THE OUTPUT TIMES                 #
!     #                                                                #     
!     ##################################################################      
! 
      ISLP=5  ! SLEEP TIME PER ATTEMPT (SECONDS)
      TMAX=20 ! MAXIMUM SLEEP TIME PER FILE (MINUTES)
      JUMPMAX=TMAX*60/ISLP 
!
      DO 100 NK=1,KOUNT,NINC
!
       ALLOCATE(DUM(IM,JM,8))
!
!     ##################################################################
!     #                                                                #
!     #         GENERATE THE NAME OF THE CURRENT RESTRT FILE           #
!     #                                                                #     
!     ##################################################################      
! 
       ENVAR=' '
       CALL GETENV("RSTFNL",ENVAR)
       CALL GETENV("tmmark",RESTHR)
       KPATH=INDEX(ENVAR,' ') -1
       IF(KPATH <= 0) KPATH=LEN(ENVAR)
!
!     ##################################################################
!
       IF(RESTHR == '    ') THEN
       WRITE(RSTFIL,50)ITAG
       WRITE(NRSTFIL,50)ITAG+NINC
   50  FORMAT('restrt',I6.6)
       ELSE
       WRITE(RSTFIL,55)ITAG,RESTHR
       WRITE(NRSTFIL,55)ITAG+NINC,RESTHR
   55  FORMAT('restrt',I6.6,'.',a4)
       ENDIF
!
       KRST=INDEX(RSTFIL,' ') -1
       NKRST=INDEX(NRSTFIL,' ') -1
       IF(KRST <= 0) KRST=LEN(RSTFIL)
       IF(NKRST <= 0) NKRST=LEN(NRSTFIL)
!
       CLOSE(LRSTRT)
!
       IF(ENVAR(1:4) == BLANK) THEN
       FNAME=RSTFIL
       NFNAME=NRSTFIL
       ELSE
       FNAME=ENVAR(1:KPATH)//RSTFIL(1:KRST)
       NFNAME=ENVAR(1:KPATH)//NRSTFIL(1:NKRST)
       ENDIF
!
!     ##################################################################
!
       JUMP=0
!
!     ##################################################################
!    
       IF(NK/=KOUNT) THEN
 110   INQUIRE(FILE=NFNAME,EXIST=FILEFOUND)
       IF(.NOT.FILEFOUND) THEN
 115   JUMP=JUMP+1
!
       IF(JUMP>JUMPMAX) THEN
       WRITE(6,*)'FROM POST0.f90:'
       WRITE(6,*)'EXCEEDED THE MAXIMUM NUMBER OF ATTEMPTS TO OPEN FILE: ',FNAME
       WRITE(6,*)
       STOP
       ENDIF
!
       WRITE(6,*)'WAITING FOR FILE: ',FNAME
       CALL SLEEP(ISLP)      
       GOTO 110 
       ENDIF
       ENDIF
!
!     ##################################################################
!
       IF(NK==KOUNT-NINC) THEN 
       CALL STAT(FNAME,STATV)
       LREC0=STATV(8)
       ENDIF
!
!     ##################################################################
!
       IF(NK==KOUNT) THEN
 120   CALL STAT(FNAME,STATV)
       LREC=STATV(8)
       IF(LREC<LREC0) THEN
       JUMP=JUMP+1
!
       IF(JUMP>JUMPMAX) THEN
       WRITE(6,*)'FROM POST0.f90:'
       WRITE(6,*)'EXCEEDED THE MAXIMUM NUMBER OF ATTEMPTS TO OPEN FILE: ',FNAME
       WRITE(6,*)
       STOP
       ENDIF
!
       WRITE(6,*)'WAITING FOR FILE: ',FNAME
       CALL SLEEP(ISLP)      
       GOTO 120 
       ENDIF
       ENDIF
!
!     ##################################################################
!
       OPEN(UNIT=LRSTRT,FILE=FNAME,STATUS='OLD',FORM='UNFORMATTED')
!
!     ##################################################################
!     #                                                                #
!     #          READ QUANTITIES NEEDED FROM THE RESTRT FILE           #
!     #                                                                #     
!     ##################################################################      
!
       READ(LRSTRT)RUN,IDAT,IHRST,NTSD
!       WRITE(6,*) 'RUN,IDAT,IHRST,NTSD:', RUN,IDAT,IHRST,NTSD
       READ(LRSTRT)
!
       DO L=1,LM
        READ(LRSTRT) DUM(:,:,1) 
       DO N=1,NUMSTA
        OMGALF(N,L)=DUM(IHINDX(N),JHINDX(N),1)
       ENDDO
       ENDDO
!
       READ(LRSTRT)
       READ(LRSTRT)PD,DUM(:,:,1:2) 
!
       DO N=1,NUMSTA
        RES(N)=DUM(IHINDX(N),JHINDX(N),1)
        FIS(N)=DUM(IHINDX(N),JHINDX(N),2)
       ENDDO           
!
       READ(LRSTRT)
!
       DO L=1,LM
        READ(LRSTRT) DUM(:,:,1) 
        READ(LRSTRT) DUM(:,:,2) 
        READ(LRSTRT) DUM(:,:,3) 
        READ(LRSTRT) DUM(:,:,4)
        READ(LRSTRT) DUM(:,:,5) 
        READ(LRSTRT)
        READ(LRSTRT) DUM(:,:,6) 
        READ(LRSTRT) DUM(:,:,7) 
        READ(LRSTRT) DUM(:,:,8) 
!
       DO N=1,NUMSTA
        T(N,L)=DUM(IHINDX(N),JHINDX(N),1)
        Q(N,L)=DUM(IHINDX(N),JHINDX(N),2)
        U(N,L)=DUM(IVINDX(N),JVINDX(N),3)
        V(N,L)=DUM(IVINDX(N),JVINDX(N),4)
        Q2(N,L)=DUM(IHINDX(N),JHINDX(N),5)
        CWM(N,L)=DUM(IHINDX(N),JHINDX(N),6)
        TRAIN(N,L)=DUM(IHINDX(N),JHINDX(N),7)
        TCUCN(N,L)=DUM(IHINDX(N),JHINDX(N),8)
       ENDDO
!
       ENDDO
!
       READ(LRSTRT)RUN,IDAT,IHRST,NTSD,LABEL,DUM(:,:,1:6) 
!
       DO N=1,NUMSTA
        Z0(N)=DUM(IHINDX(N),JHINDX(N),4)
        CZEN(N)=DUM(IHINDX(N),JHINDX(N),6)
       ENDDO      
!
       READ(LRSTRT) DUM(:,:,1:7) 
!
       DO N=1,NUMSTA
        THS(N)=DUM(IHINDX(N),JHINDX(N),2)
        HBOT(N)=DUM(IHINDX(N),JHINDX(N),6)
        CFRACL(N)=DUM(IHINDX(N),JHINDX(N),7)
       ENDDO    
!
       READ(LRSTRT) DUM(:,:,1:7) 
!
       DO N=1,NUMSTA
        CFRACM(N)=DUM(IHINDX(N),JHINDX(N),7)
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1:7) 
!
       DO N=1,NUMSTA
        SNO(N)=DUM(IHINDX(N),JHINDX(N),1)
        PSLP(N)=DUM(IHINDX(N),JHINDX(N),5)
        CFRACH(N)=DUM(IHINDX(N),JHINDX(N),7)
       ENDDO            
!
       READ(LRSTRT) DUM(:,:,1:4) 
!
       DO N=1,NUMSTA
        SOILTB(N)=DUM(IHINDX(N),JHINDX(N),1)
        SFCEXC(N)=DUM(IHINDX(N),JHINDX(N),2)
        SMSTAV(N)=DUM(IHINDX(N),JHINDX(N),3)
        SMSTOT(N)=DUM(IHINDX(N),JHINDX(N),4)
       ENDDO             
!
       READ(LRSTRT) DUM(:,:,1:3) 
!
       DO N=1,NUMSTA
        CZMEAN(N)=DUM(IHINDX(N),JHINDX(N),3)
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1) ,UL,IDUM, DUM(:,:,2)
!
       DO N=1,NUMSTA
        U00(N)=DUM(IHINDX(N),JHINDX(N),1)
        SR(N)=DUM(IHINDX(N),JHINDX(N),2)
       ENDDO      
!
       READ(LRSTRT)RUN,IDAT,IHRST,NTSD,LABEL, DUM(:,:,1:4) 
!
       DO N=1,NUMSTA
        ACPREC(N)=DUM(IHINDX(N),JHINDX(N),2)
        CUPREC(N)=DUM(IHINDX(N),JHINDX(N),4)
       ENDDO
!
       READ(LRSTRT)
       READ(LRSTRT) DUM(:,:,1:4)
!
       DO N=1,NUMSTA
        ACSNOW(N)=DUM(IHINDX(N),JHINDX(N),1)
        ACSNOM(N)=DUM(IHINDX(N),JHINDX(N),2)
        SSROFF(N)=DUM(IHINDX(N),JHINDX(N),3)
        BGROFF(N)=DUM(IHINDX(N),JHINDX(N),4)
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1:4 )
!
       DO N=1,NUMSTA
        SFCSHX(N)=DUM(IHINDX(N),JHINDX(N),1)
        SFCLHX(N)=DUM(IHINDX(N),JHINDX(N),2)
        SUBSHX(N)=DUM(IHINDX(N),JHINDX(N),3)
        SNOPCX(N)=DUM(IHINDX(N),JHINDX(N),4)
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1:6)
!
       DO N=1,NUMSTA
        ASWIN(N)=DUM(IHINDX(N),JHINDX(N),1)
        ASWOUT(N)=DUM(IHINDX(N),JHINDX(N),2)
        ASWTOA(N)=DUM(IHINDX(N),JHINDX(N),3)
        ALWIN(N)=DUM(IHINDX(N),JHINDX(N),4)
        ALWOUT(N)=DUM(IHINDX(N),JHINDX(N),5)
        ALWTOA(N)=DUM(IHINDX(N),JHINDX(N),6)
       ENDDO         
!
       READ(LRSTRT)
       READ(LRSTRT) DUM(:,:,1:6) 
!
       DO N=1,NUMSTA
        TH10(N)=DUM(IHINDX(N),JHINDX(N),1)
        Q10(N)=DUM(IHINDX(N),JHINDX(N),2)
        U10(N)=DUM(IHINDX(N),JHINDX(N),3)
        V10(N)=DUM(IHINDX(N),JHINDX(N),4)
        TSHLTR(N)=DUM(IHINDX(N),JHINDX(N),5)
        QSHLTR(N)=DUM(IHINDX(N),JHINDX(N),6)
       ENDDO            
!
       READ(LRSTRT) DUM(:,:,1:NSOIL) 
!
       DO NS=1,NSOIL
       DO N=1,NUMSTA
        SMC(N,NS)=DUM(IHINDX(N),JHINDX(N),NS)
       ENDDO
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1) 
!
       DO N=1,NUMSTA
        CMC(N)=DUM(IHINDX(N),JHINDX(N),1)
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1:NSOIL) 
       DO NS=1,NSOIL
       DO N=1,NUMSTA
        STC(N,NS)=DUM(IHINDX(N),JHINDX(N),NS)
       ENDDO     
       ENDDO
!
       READ(LRSTRT) DUM(:,:,1:NSOIL) 
       DO NS=1,NSOIL
       DO N=1,NUMSTA
        SH2O(N,NS)=DUM(IHINDX(N),JHINDX(N),NS)
       ENDDO     
       ENDDO
!
       READ(LRSTRT)
!
!     ##################################################################
!
       READ(LRSTRT) DUM(:,:,1:3),                                      & 
       ACUTIM,ARATIM,APHTIM,                                           &
       NHEAT,NPHS,NCNVC,NPREC,NRDSW,NRDLW,NSRFC,                       &
       TPH0D,TLM0D,RESTRT
!
       DO N=1,NUMSTA
        POTFLX(N)=DUM(IHINDX(N),JHINDX(N),1)
        TLMIN(N)=DUM(IHINDX(N),JHINDX(N),2)
        TLMAX(N)=DUM(IHINDX(N),JHINDX(N),3)
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #             READ RADIATIVE TEMPERATURE TENDENCIES              #
!     #                                                                #     
!     ##################################################################      
!             
       DO L=1,LM
        READ(LRSTRT) DUM(:,:,1) 
        READ(LRSTRT) DUM(:,:,2) 
       DO N=1,NUMSTA
        RSWTT(N,L)=DUM(IHINDX(N),JHINDX(N),1)
        RLWTT(N,L)=DUM(IHINDX(N),JHINDX(N),2)
       ENDDO      
       ENDDO
!
       CLOSE(LRSTRT)
!
!     ##################################################################
!
       DEALLOCATE(DUM)
!
!     ##################################################################
!     #                                                                #
!     #                        THE FORECAST HOUR                       #
!     #                                                                #     
!     ##################################################################      
!  
       IFHR=NTSD/NTSPH
!
!     ##################################################################
!     #                                                                #
!     #                    ALL THE DATA IS NOW IN.                     #
!     #    CALCULATE CLOUD FRACTION AND CLOUD WATER/ICE ID NUMBER      #
!     #                                                                #     
!     ##################################################################      
!
       UTIM=1.
       US=1.
       CCLIMIT=1.E-3
       CLIMIT =1.E-20
!
!     ##################################################################
!

       DO N=1,NUMSTA
        IW(N,1)=0
        CCR(N,1)=0.
        PDSL1(N)=PD(IHINDX(N),JHINDX(N))*RES(N)
       ENDDO       
!
!     ##################################################################
!     #                                                                #
!     #                         QW, QI AND QINT                        #
!     #                                                                #     
!     ##################################################################      
!
       DO L=2,LM
!
       DO N=1,NUMSTA
        LML=LM-LMH(IHINDX(N),JHINDX(N))
        HH=HTM(N,L)*HBM2(N)
        TKL=T(N,L)
        QKL=Q(N,L)
        CWMKL=CWM(N,L)
        TMT0=(TKL-273.16)*HH
        TMT15=AMIN1(TMT0,-15.)*HH
        AI=0.008855
        BI=1.
!
        IF(TMT0 < -20.) THEN
        AI=0.007225
        BI=0.9674
        ENDIF
!
        PP=PDSL1(N)*AETA(L)+PT
        QW=HH*PQ0/PP*EXP(HH*A2*(TKL-A3)/(TKL-A4))
        QI=QW*(BI+AI*AMIN1(TMT0,0.))
        QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
        IF(TMT0<=-40.) QINT=QI
!
!     ##################################################################
!     #                                                                #
!     #                     ICE-WATER ID NUMBER IW                     #
!     #                                                                #     
!     ##################################################################      
!
        U00KL=U00(N)+UL(L+LML)*(0.95-U00(N))*UTIM
        IF(TMT0 < -15.) THEN
        FIQ=QKL-U00KL*QI
!
        IF(FIQ > 0. .OR. CWMKL > CLIMIT) THEN
        IW(N,L)=1
        ELSE
        IW(N,L)=0
        ENDIF
        ENDIF
!
        IF(TMT0 >= 0.) THEN
        IW(N,L)=0
        ENDIF
!
        IF(TMT0 < 0. .AND. TMT0 >= -15.) THEN
        IW(N,L)=0
        IF(IW(N,L-1)==1.AND.CWMKL>CLIMIT) IW(N,L)=1
        ENDIF
!
        IWKL=IW(N,L)
!
!     ##################################################################
!     #                                                                #
!     #               THE SATUATION SPECIFIC HUMIDITY                  #
!     #                                                                #     
!     ##################################################################      
!
        FIW=FLOAT(IWKL)
        QC=(1.-FIW)*QINT+FIW*QI
!
!     ##################################################################
!     #                                                                #
!     #                      THE RELATIVE HUMIDITY                     #
!     #                                                                #     
!     ##################################################################      
!
        IF(QC <= 0.) THEN
        RQKL=0.
        ELSE
        RQKL=QKL/QC
        ENDIF
!
!     ##################################################################
!     #                                                                #
!     #                     CLOUD COVER RATIO CCR                      #
!     #                                                                #     
!     ##################################################################      
!
        IF(RQKL >= 0.9999) THEN
        CCR(N,L)=AMIN1(US,RQKL)
        ELSE
        ARG=-1000.*CWMKL/(US-RQKL)
        ARG=AMAX1(ARG,-25.)
        CCR(N,L)= RQKL*(1.-EXP(ARG))
        ENDIF
!
!     ##################################################################
!
       ENDDO
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #                BEGIN THE PROFILE POSTING CODE                  #
!     #                                                                #
!     #      USE ZERO IN ACCUMULATION ARRAYS AT APPROPRIATE TIMES      #
!     #                                                                #     
!     ##################################################################      
!
       IF(RESTRT.OR.NTSD<2) THEN
       STATPR=0.
       STACPR=0.
       STAEVP=0.
       STAPOT=0.
       STASHX=0.
       STASUB=0.
       STAPCX=0.
       STASWI=0.
       STASWO=0.
       STALWI=0.
       STALWO=0.
       STALWT=0.
       STASWT=0.
       STASNM=0.
       STASNO=0.
       STASRF=0.
       STABRF=0.
       DHCNVC=0.
       DHRAIN=0.
       ELSE
!
!     ##################################################################
!     #                                                                #
!     #  WE MUST CHECK TO SEE IF WE ARE 1 HOUR AFTER ANY OF THE        #
!     #  ACCUMULATION BUCKETS HAVE BEEN EMPTIED.  IF WE ARE AT SUCH A  #
!     #  TIME THEN WE NEED TO SET TO ZERO THE VARIABLES USED TO HOLD   #
!     #  THE PRECEDING HOUR'S VALUES.                                  #
!     #                                                                #     
!     ##################################################################      
!
       TIME=(NTSD-1)*DT
       RESET0=TIME-(NTSD/NPREC)*NPREC*DT
       RESET1=(NPHS-1)*DT+3600.
!
!JLG       IF(MOD(NTSD,NPREC) >= NPHS.AND.RESET0 <= RESET1) THEN
       IF(MOD(NTSD,NPREC) >= TPREC.AND.RESET0 <= RESET1) THEN
       STATPR=0.
       STACPR=0.
       STASNM=0.
       STASNO=0.
       STASRF=0.
       STABRF=0.
       ELSE
       STATPR=ACPREC0*1.E3
       STACPR=CUPREC0*1.E3
       STASNM=ACSNOM0*1.E3
       STASNO=ACSNOW0*1.E3
       STASRF=SSROFF0*1.E3
       STABRF=BGROFF0*1.E3
       ENDIF          
!
       RESET0=TIME-(NTSD/NRDSW)*NRDSW*DT
!
!JLG       IF(MOD(NTSD,NRDSW) >= NPHS .AND. RESET0 <= RESET1) THEN
       IF(MOD(NTSD,NRDSW) >= TRDSW .AND. RESET0 <= RESET1) THEN
       STASWI=0.
       STASWO=0.
       STASWT=0.
       ELSE
       STASWI=ASWIN0
       STASWO=ASWOUT0
       STASWT=ASWTOA0
       ENDIF 
!
       RESET0=TIME-(NTSD/NRDLW)*NRDLW*DT
!
!JLG       IF(MOD(NTSD,NRDLW) >= NPHS .AND. RESET0 <= RESET1) THEN
       IF(MOD(NTSD,NRDLW) >= TRDLW .AND. RESET0 <= RESET1) THEN
       STALWI=0.
       STALWO=0.
       STALWT=0.
       ELSE
       STALWI=ALWIN0
       STALWO=ALWOUT0
       STALWT=-ALWTOA0
       ENDIF
!
       RESET0=TIME-(NTSD/NSRFC)*NSRFC*DT
!
!JLG       IF(MOD(NTSD,NSRFC) >= NPHS .AND. RESET0 <= RESET1) THEN
       IF(MOD(NTSD,NSRFC) >= TSRFC .AND. RESET0 <= RESET1) THEN
       STAEVP=0.
       STAPOT=0.
       STASHX=0.
       STASUB=0.
       STAPCX=0.
       ELSE
       STAEVP=SFCLHX0
       STAPOT=POTFLX0
       STASHX=SFCSHX0
       STASUB=SUBSHX0
       STAPCX=SNOPCX0
       ENDIF
!
       RESET0=TIME-(NTSD/NHEAT)*NHEAT*DT
!
!JLG       IF(MOD(NTSD,NHEAT) >= NCNVC .AND. RESET0 <= RESET1) THEN
       IF(MOD(NTSD,NHEAT) >= THEAT .AND. RESET0 <= RESET1) THEN
       DHCNVC=0.
       DHRAIN=0.
       ELSE
       DO N=1,NUMSTA
       DO L=1,LM
        DHCNVC(L,N)=TCUCN0(N,L)
        DHRAIN(L,N)=TRAIN0(N,L)
       ENDDO
       ENDDO
       ENDIF
!
!     ##################################################################
!
       ENDIF 
!
!     ##################################################################
!     #                                                                #
!     #    FOR ROTATION OF WINDS FROM E-GRID TO GEODETIC ORIENTATION   #
!     #    WE NEED THE TWO QUANTITIES BELOW                            #
!     #                                                                #     
!     ##################################################################      
!
       SINPH0=SIN(TPH0D*DTR)
       COSPH0=COS(TPH0D*DTR)
!
!     ##################################################################
!     #                                                                #
!     #    INITIAL CALCULATIONS/PREPARATIONS. WE LOAD SEVERAL ARRAYS   #
!     #    WITH PROFILE VARIABLES                                      #
!     #                                                                #     
!     ##################################################################      
!
       DO N=1,NUMSTA
        IF(CZMEAN(N) > 0.) THEN
        FACTR(N)=CZEN(N)/CZMEAN(N)
        ELSE
        FACTR(N)=0.
        ENDIF
       ENDDO  
!
!     ##################################################################
!     #                                                                #
!     #    ADJUST SHORTAVE TENDENCIES TO ACCOUNT FOR CHANGE OF SOLAR   #
!     #    POSITION BETWEEN CALLS TO RADIATION                         #
!     #                                                                #     
!     ##################################################################      
!
       DO L=1,LM
       DO N=1,NUMSTA
        RSWTT(N,L)=RSWTT(N,L)*FACTR(N)
       ENDDO
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #                           COMPUTE RTOP                         #
!     #                                                                #     
!     ##################################################################      
!
       DO L=1,LM
       DO N=1,NUMSTA
        APEL=PT+AETA(L)*PDSL1(N)
        RTOP(N,L)=R*T(N,L)*(1.+0.608*Q(N,L))/APEL
       ENDDO
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #                    PDS IS SURFACE PRESSURE                     #
!     #                                                                #     
!     ##################################################################      
!
      PDS=PD+PT
!
!     ##################################################################
!     #                                                                #
!     #               EGRID2 IS THE SURFACE TEMPERATURE                #
!     #                                                                #     
!     ##################################################################      
!
       DO N=1,NUMSTA
        EGRID2(N)=THS(N)*(PDS(IHINDX(N),JHINDX(N))*1.E-5)**CAPA
        IF(ACPREC(N) < 0.) ACPREC(N)=0.
        IF(CUPREC(N) < 0.) CUPREC(N)=0.
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #               SET CYCLE, DATE, AND FORECAST TIME               #
!     #                                                                #     
!     ##################################################################      
!
       IHR=NTSD/NTSPH+0.5
       IYR=IDAT(3)
       IMNTH=IDAT(1)
       IDAY=IDAT(2)
       IFCST=(NTSD-1)*DT
!
       WRITE(LIST,*)' POST PROFILE FOR ',IYR,IMNTH,IDAY,IHR      
!
!     ##################################################################
!     #                                                                #
!     #  SET RTSPH,RTSCU,RTSRA TO 1. OVER THE NUMBER OF TIMES THE      #
!     #  VARIOUS PHYSICS ROUTINES HAVE BEEN CALLED SINCE LAST OUTPUT   #
!     #  OF PROFILER DATA. NECESSARY FOR CORRECT AVERAGING OF          #
!     #  VARIABLES                                                     #
!     #                                                                #     
!     ##################################################################      
!
       IF(APHTIM > 0.) THEN
       RTSPH=1./APHTIM
       ELSE
       RTSPH=1.
       ENDIF
!
       IF(ACUTIM > 0.) THEN
       RTSCU=1./ACUTIM
       ELSE
       RTSCU=1.
       ENDIF
!
       IF(ARATIM > 0.) THEN
       RTSRA=1./ARATIM
       ELSE
       RTSRA=1.
       ENDIF
!
!     ##################################################################
!     #                                                                #
!     #   OUTPUT PROFILE DATA.                                         #
!     #   THE FOLLOWING LOOP IS OVER ALL PROFILE SITES                 #
!     #                                                                #     
!     ##################################################################      
!
       DO 1000 N=1,NUMSTA
!
!     ##################################################################
!     #                                                                #
!     #                       ZERO OUTPUT ARRAY                        #
!     #                                                                #     
!     ##################################################################      
!
       PRODAT=0.
       FPACK=0.
!
!     ##################################################################
!     #                                                                #
!     #  CONSTRUCT HEADER FOR CURRENT PROFILE SITE.  THE HEADER        #
!     #  CONTAINS THE FOLLOWING INFORMATION:  PACKED CYCLE-DATE,       #
!     #  FORECAST TIME, INTEGER STATION ID, STATION LATITUDE, STATION  #
!     #  STATION ELEVATION, NUMBER OF VERTICAL LEVELS IN PROFILE,      #
!     #  NUMBER OF MULTI-LEVEL PARAMETERS, NUMBER OF SINGLE LEVEL      #
!     #  PARAMETERS, TOTAL LENGTH (IN WORDS) OF MULTI- AND SINGLE      #
!     #  LEVEL DATA, PROFILE CLASS FLAG, AND A DUMMY WORD FOR FUTURE   #
!     #  USE                                                           #
!     #                                                                #     
!     ##################################################################      
!
       IH=IHINDX(N)
       JH=JHINDX(N)
       LMHK=LMH(IH,JH)
       NWORD2=2*LMHK
       NWORD3=3*LMHK
       NWORD4=4*LMHK
       NWORD5=5*LMHK
       NWORD6=6*LMHK
       NWORD7=7*LMHK
       NWORD8=8*LMHK
       NWORD9=9*LMHK
       NWORD10=10*LMHK
       NWORD11=11*LMHK
       NWORD12=12*LMHK
       NWORD13=13*LMHK
       NWORD14=14*LMHK
       NWORD15=15*LMHK
       ISTAT=IDSTN(N)
       CISTAT=CIDSTN(N)
!
       FPACK(1)=STNLAT(N)/DTR
       FPACK(2)=-STNLON(N)/DTR
       IF(FPACK(2)<-180.) FPACK(2)=FPACK(2)+360.
       FPACK(3)=FIS(N)*GI
       FPACK(4)=FLOAT(LMHK)
       FPACK(5)=LCL1ML
       FPACK(6)=LCL1SL
       FPACK(7)=9+FPACK(5)*FPACK(4)+FPACK(6)
       FPACK(8)=999.
       FPACK(9)=999.
!
!     ##################################################################
!     #                                                                #
!     #                WIND ROTATION SINES AND COSINES                 #
!     #                                                                #     
!     ##################################################################      
!	   
      DLM=STNLON(N)+TLM0D*DTR
      XX=COSPH0*COS(STNLAT(N))*COS(DLM)                                &
      +SINPH0*SIN(STNLAT(N))
      YY=-COS(STNLAT(N))*SIN(DLM)
      TLON=ATAN(YY/XX)
      ALPHA=ASIN(SINPH0*SIN(TLON)/COS(STNLAT(N)))
      SINALP=SIN(ALPHA)
      COSALP=COS(ALPHA)
!
!     ##################################################################
!     #                                                                #
!     #  EXTRACT PRESSURE AND TEMPERATURE PROFILES.                    #
!     #  EXTRACT/ROTATE U AND V WIND COMPONENT PROFILES.               #
!     #  EXTRACT SPECIFIC HUMIDITY AND TEMPERATURE TENDENCY.           #
!     #  EXTRACT CLOUD WATER, HEATING DUE TO CONVECTION, LARGE         #
!     #  SCALE RAIN, SHORT WAVE RADIATION, LONG WAVE RADIATION,        #
!     #  AND CLOUD FRACTION                                            #
!     #                                                                #     
!     ##################################################################      
!
       DO LV=1,LMHK
        LVL=LMHK-LV+1
        PRODAT(LVL)=PDSL1(N)*AETA(LV)+PT
        PRODAT(LMHK+LVL)=T(N,LV)
!
!     ##################################################################
!     #                                                                #
!     #                          ROTATE WINDS                          #
!     #                                                                #     
!     ##################################################################      
!
        UT=U(N,LV)
        VT=V(N,LV)
        PRODAT(NWORD2+LVL)=UT*COSALP+VT*SINALP
        PRODAT(NWORD3+LVL)=VT*COSALP-UT*SINALP
!
        PRODAT(NWORD4+LVL)=Q(N,LV)
!
        IF(RTOP(N,LV) > 1.E-12)                                        &
        PRODAT(NWORD5+LVL)=OMGALF(N,LV)*CP/(RTOP(N,LV)*DT)
!  
        IF(IW(N,LV) > 0.5) THEN
        PRODAT(NWORD6+LVL)=-CWM(N,LV)
        ELSE
        PRODAT(NWORD6+LVL)=CWM(N,LV)
        ENDIF
!
        PRODAT(NWORD7+LVL)=TCUCN(N,LV)
        PRODAT(NWORD8+LVL)=TRAIN(N,LV)
        PRODAT(NWORD9+LVL)=RSWTT(N,LV)
        PRODAT(NWORD10+LVL)=RLWTT(N,LV)
        PRODAT(NWORD11+LVL)=CCR(N,LV)*100.
!
        IF(LV == 1) THEN
        PRODAT(NWORD12+LVL)=Q2(N,LV)
        ELSE
        PRODAT(NWORD12+LVL)=(Q2(N,LV)+Q2(N,LV-1))*0.5
        ENDIF
!
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #   MODIFY ACCUMLATIONS SO AS TO REPRESENT ACCUMULATED           #
!     #   CHANGE SINCE LAST PROFILE OUTPUT TIME                        #
!     #                                                                #     
!     ##################################################################      
!
       DO LL=1,LMHK
        STADHC(LL)=PRODAT(NWORD7+LL)-DHCNVC(LL,N)
        STADHR(LL)=PRODAT(NWORD8+LL)-DHRAIN(LL,N)
!
        DHCNVC(LL,N)=PRODAT(NWORD7+LL)
        DHRAIN(LL,N)=PRODAT(NWORD8+LL)
!
        IF(MOD(NTSD,NHEAT)<NCNVC) THEN
        DHCNVC(LL,N)=0.
        DHRAIN(LL,N)=0.
        ENDIF
!
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #   EXTRACT SINGLE LEVEL DATA                                    #
!     #   EGRID2 IS SURFACE TEMPERATURE                                #
!     #                                                                #     
!     ##################################################################      
!
       PRODAT(NWORD15+1)=PSLP(N)
       PRODAT(NWORD15+2)=PDS(IH,JH)
       PRODAT(NWORD15+3)=EGRID2(N)
       PRODAT(NWORD15+4)=TLMIN(N)
       PRODAT(NWORD15+5)=TLMAX(N)
       PRODAT(NWORD15+6)=SMSTAV(N)*100.
       PRODAT(NWORD15+7)=ACPREC(N)*1000.
       PRODAT(NWORD15+8)=CUPREC(N)*1000.
       PRODAT(NWORD15+27)=Z0(N)
!
       STAPRX=PRODAT(NWORD15+7)-STATPR(N)
       STACRX=PRODAT(NWORD15+8)-STACPR(N)
!
!     ##################################################################
!     #                                                                #
!     #                          ROTATE WINDS                          #
!     #                                                                #     
!     ##################################################################      
!
       UT=U10(N)
       VT=V10(N)
       PRODAT(NWORD15+28)=UT*COSALP+VT*SINALP
       PRODAT(NWORD15+29)=VT*COSALP-UT*SINALP
!
       PRODAT(NWORD15+30)=TH10(N)
       PRODAT(NWORD15+31)=Q10(N)
       PRODAT(NWORD15+32)=TSHLTR(N)
       PRODAT(NWORD15+33)=QSHLTR(N)
       PRODAT(NWORD15+34)=SFCEXC(N)
       PRODAT(NWORD15+35)=VEGFRC(N)*100.
       PRODAT(NWORD15+36)=CMC(N)*1000.
       PRODAT(NWORD15+37)=SMC(N,1)
       PRODAT(NWORD15+38)=SMC(N,2)
       PRODAT(NWORD15+39)=SMC(N,3)
       PRODAT(NWORD15+40)=SMC(N,4)
       PRODAT(NWORD15+41)=STC(N,1)
       PRODAT(NWORD15+42)=STC(N,2)
       PRODAT(NWORD15+43)=STC(N,3)
       PRODAT(NWORD15+44)=STC(N,4)
       PRODAT(NWORD15+45)=SM(N)+SICE(N)
       PRODAT(NWORD15+46)=CFRACL(N)*100.
       PRODAT(NWORD15+47)=CFRACM(N)*100.
       PRODAT(NWORD15+48)=CFRACH(N)*100.
       PRODAT(NWORD15+49)=SR(N)*100.
       PRODAT(NWORD15+50)=NINT(HBOT(N))
!
       PRODAT(NWORD15+9)=SFCLHX(N)
       PRODAT(NWORD15+10)=POTFLX(N)
       PRODAT(NWORD15+11)=SFCSHX(N)
       PRODAT(NWORD15+12)=SUBSHX(N)
       PRODAT(NWORD15+13)=SNOPCX(N)
       PRODAT(NWORD15+14)=ASWIN (N)
       PRODAT(NWORD15+15)=ASWOUT(N)
       PRODAT(NWORD15+16)=ALWIN (N)
       PRODAT(NWORD15+17)=ALWOUT(N)
       PRODAT(NWORD15+18)=-ALWTOA(N)
       PRODAT(NWORD15+19)=ASWTOA(N)
       PRODAT(NWORD15+20)=ACSNOW(N)*1000.
       PRODAT(NWORD15+21)=SMSTOT(N)*1000.
       PRODAT(NWORD15+22)=SNO(N)*1000.
       PRODAT(NWORD15+23)=ACSNOM(N)*1000.
       PRODAT(NWORD15+24)=SSROFF(N)*1000.
       PRODAT(NWORD15+25)=BGROFF(N)*1000.
       PRODAT(NWORD15+26)=SOILTB(N)
!
!     ##################################################################
!     #                                                                #
!     #       ACCUMULATED CHANGE SINCE LAST PROFILE OUTPUT TIME        #
!     #                                                                #     
!     ##################################################################      
!
       PSFCEVP=PRODAT(NWORD15+9 )-STAEVP(N)
       PPOTEVP=PRODAT(NWORD15+10)-STAPOT(N)
       PSFCSHX=PRODAT(NWORD15+11)-STASHX(N)
       PSFCSUB=PRODAT(NWORD15+12)-STASUB(N)
       PSNOPCX=PRODAT(NWORD15+13)-STAPCX(N)
       PRSWIN=PRODAT(NWORD15+14)-STASWI(N)
       PRSWOUT=PRODAT(NWORD15+15)-STASWO(N)
       PRLWIN=PRODAT(NWORD15+16)-STALWI(N)
       PRLWOUT=PRODAT(NWORD15+17)-STALWO(N)
       PRLWTOA=PRODAT(NWORD15+18)-STALWT(N)
       PRSWTOA=PRODAT(NWORD15+19)-STASWT(N)
       PACSNOW=PRODAT(NWORD15+20)-STASNO(N)
       PACSNOM=PRODAT(NWORD15+23)-STASNM(N)
       PSSROFF=PRODAT(NWORD15+24)-STASRF(N)
       PBGROFF=PRODAT(NWORD15+25)-STABRF(N)
!
!     ##################################################################
!     #                                                                #
!     #     TRANSFER STATION PROFILE DATA TO "PACKED" OUTPUT ARRAY     #
!     #                                                                #     
!     ##################################################################      
!
       NN=0
       NLEN=FPACK(7)
!
       DO NL=10,NLEN
        NN=NL-9
        FPACK(NL)=PRODAT(NN)
       ENDDO
!
!     ##################################################################
!     #                                                                #
!     #   REPLACE ACCUMULATED QUANTITIES WITH ACCUMULATION             #
!     #   SINCE LAST PROFILE OUTPUT TIME                               #
!     #                                                                #     
!     ##################################################################      
!
       DO LL=1,LMHK
        FPACK(9+NWORD7+LL)=STADHC(LL)*RTSCU
        FPACK(9+NWORD8+LL)=STADHR(LL)*RTSRA
       ENDDO
!
       FPACK(9+NWORD15+7)=STAPRX
       FPACK(9+NWORD15+8)=STACRX
       FPACK(9+NWORD15+9)=PSFCEVP*RTSPH
       FPACK(9+NWORD15+10)=PPOTEVP*RTSPH
       FPACK(9+NWORD15+11)=PSFCSHX*RTSPH
       FPACK(9+NWORD15+12)=PSFCSUB*RTSPH
       FPACK(9+NWORD15+13)=PSNOPCX*RTSPH
       FPACK(9+NWORD15+14)=PRSWIN*RTSPH
       FPACK(9+NWORD15+15)=PRSWOUT*RTSPH
       FPACK(9+NWORD15+16)=PRLWIN*RTSPH
       FPACK(9+NWORD15+17)=PRLWOUT*RTSPH
       FPACK(9+NWORD15+18)=PRLWTOA*RTSPH
       FPACK(9+NWORD15+19)=PRSWTOA*RTSPH
       FPACK(9+NWORD15+20)=PACSNOW
       FPACK(9+NWORD15+23)=PACSNOM
       FPACK(9+NWORD15+24)=PSSROFF
       FPACK(9+NWORD15+25)=PBGROFF
!
       IF(RESTRT) THEN
       DO LL=1,LMHK
        FPACK(9+NWORD7+LL)=0.
        FPACK(9+NWORD8+LL)=0.
       ENDDO
       FPACK(9+NWORD15+7)=0.
       FPACK(9+NWORD15+8)=0.
       FPACK(9+NWORD15+9)=0.
       FPACK(9+NWORD15+10)=0.
       FPACK(9+NWORD15+11)=0.
       FPACK(9+NWORD15+12)=0.
       FPACK(9+NWORD15+13)=0.
       FPACK(9+NWORD15+14)=0.
       FPACK(9+NWORD15+15)=0.
       FPACK(9+NWORD15+16)=0.
       FPACK(9+NWORD15+17)=0.
       FPACK(9+NWORD15+18)=0.
       FPACK(9+NWORD15+19)=0.
       FPACK(9+NWORD15+20)=0.
       FPACK(9+NWORD15+23)=0.
       FPACK(9+NWORD15+24)=0.
       FPACK(9+NWORD15+25)=0.
       ENDIF
!
!     ##################################################################
!     #                                                                #
!     #                        WRITE PROFILE DATA                      #
!     #                                                                #     
!     ##################################################################      
!  
       NREC=IFHR*NUMSTA+N
!
       IF(ISTAT == 724930) THEN
       DO NL=1,NLEN
        WRITE(6,*) 'NL,FPACK(NL): ', NL,FPACK(NL)
       ENDDO
       ENDIF
       WRITE(LCLAS1,REC=NREC)IHRST,IDAT,IFCST,ISTAT,CISTAT,            &
       (FPACK(NL),NL=1,NLEN)
!
 1000 CONTINUE
!
!     ##################################################################
!     #                                                                #
!     #          SAVE ARRAY DATA OF THE PRECEDING RESTRT FILE          #
!     #                 READ THE PREVIOUS RESTRT FILE                  #
!     #                                                                #     
!     ##################################################################      
!
      TRAIN0=TRAIN
      TCUCN0=TCUCN
      ACPREC0=ACPREC
      CUPREC0=CUPREC
      ACSNOW0=ACSNOW
      ACSNOM0=ACSNOM
      SSROFF0=SSROFF
      BGROFF0=BGROFF
      SFCSHX0=SFCSHX
      SFCLHX0=SFCLHX
      SUBSHX0=SUBSHX
      SNOPCX0=SNOPCX
      ASWIN0=ASWIN
      ASWOUT0=ASWOUT
      ASWTOA0=ASWTOA
      ALWIN0=ALWIN
      ALWOUT0=ALWOUT
      ALWTOA0=ALWTOA
      POTFLX0=POTFLX
!
!     ##################################################################
!     #                                                                #
!     #                     END OF PROFILE SITE LOOP                   #
!     #                                                                #
!     ##################################################################      
!
      ITAG=ITAG+NINC     
!
  100 CONTINUE
!
!     ##################################################################
!     #                                                                #
!     #                     END PROFILE POSTING CODE                   #
!     #                                                                #     
!     ##################################################################
!
      CLOSE(LCLAS1)
!      STOP
      END PROGRAM POST0
