      SUBROUTINE BCEX(NHB,LRSTRT,IHOUR,PTOP
     1,               IUNWGT,IUDETI,IUBCF)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  BCEX        BC EXTRACTION FOR NESTS
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-05-25
C
C ABSTRACT:  THIS ROUTINE GENERATES THE BOUNDARY CONDITION
C            TENDENCIES FOR THE ETA DOMAINS THAT ARE NESTED
C            WITHIN THE OPERATIONAL DOMAIN AND WRITES THEM 
C            TO SEPARATE FILES FOR EACH OF THE NESTED DOMAINS
C
C PROGRAM HISTORY LOG:
C   99-05-25  T BLACK - ORIGINATOR
C   00-12-05  T BLACK - MODIFIED FOR SLP REDUCTION IN NESTS
C
C USAGE:  CALL BCEX FROM PROGRAM POST0
C
C   INPUT ARGUMENT LIST:
C        NHB - THE UNIT FOR THE NHB FILE
C     LRSTRT - THE UNIT FOR THE RESTRT FILES
C      IHOUR - THE INTEGER FORECAST HOUR
C      IUWGT - THE BASE UNIT NUMBER FOR INTERPOLATION WEIGHTS
C     IUDETI - THE UNIT NUMBER FOR NEST DETA's
C      IUBCF - COUNTER FOR THE OUTPUT UNIT NUMBERS
C
C   OUTPUT ARGUMENT LIST:
C     NONE
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C           HORIZ - INTERPOLATES HORIZONTALLY FROM THE PARENT GRID 
C                   TO THE NEST'S BOUNDARY USING INTERPOLATION
C                   WEIGHTS
C          PTETAE - INTERPOLATES VERTICALLY FROM THE OUTPUT OF
C                   SUBROUTINE HORIZ TO THE INNER GRID'S MODEL
C                   LEVELS
C
C----------------------------------------------------------------------
      INCLUDE "parmeta"
C----------------------------------------------------------------------
C***
C***  HARDWIRE THE TIME TENDENCY INTERVAL TO 1 HOUR
C***
                             P A R A M E T E R
     & (TBOCO=1.)
C----------------------------------------------------------------------
                             L O G I C A L
     & EXBC,RUN
C----------------------------------------------------------------------
C
C***  NEST SPECIFICATION ARRAYS FROM THE BCEXDATA NAMELIST
C
                              R E A L
     & TBCEND(15),TPH0DI(15),TLM0DI(15)
     &,DETAO(LM),RDETA(LM),AETAO(LM)
C----------------------------------------------------------------------
C
C***  THESE ARE THE INTERPOLATION WEIGHT ARRAYS
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:) :: HLATI,HLONI,VLATI,VLONI
C
                            I N T E G E R
     &,ALLOCATABLE,DIMENSION(:,:) :: IOUTHB,JOUTHB,IOUTVB,JOUTVB
     &,                              IONEAR,JONEAR
                            I N T E G E R
     &,ALLOCATABLE,DIMENSION(:,:) :: IOUTHB2,JOUTHB2
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:,:) :: HWGTS,VWGTS,HWGTS2
C----------------------------------------------------------------------
C
C***  ARRAYS FROM THE OUTPUT RESTART FILES
C
                              R E A L
     & PD(IM,JM),RES(IM,JM),USTAR(IM,JM),PSLP(IM,JM)
     &,T(IM,JM),Q(IM,JM),U(IM,JM),V(IM,JM)
     &,Q2(IM,JM),CWM(IM,JM)
C
                              R E A L
     & DUM(IM,JM,4)
C----------------------------------------------------------------------
                              R E A L
     & HTM(IM,JM,LM),VTM(IM,JM,LM)
C----------------------------------------------------------------------
C
C***  THESE ARE THE BOUNDARY VALUES AT THE INNER GRID POINTS AFTER
C***  HORIZONTAL INTERPOLATION FROM THE PARENT GRID
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:) :: PDBO,USTARBO
     &,                              PSLPB,PSLPB2
     &,                              PSLPB_L,PSLPB2_L
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:,:) :: TBO,QBO,UBO,VBO,Q2BO,CWMBO
     &,                                HTMBO
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:,:) :: TSLPB,TSLPB2
     &,                                TSLPB_L,TSLPB2_L
     &,                                PMSLPB,PMSLPB2
     &,                                PMSLPB_L,PMSLPB2_L
C----------------------------------------------------------------------
C
C***  THESE ARE THE BOUNDARY VALUES AT THE INNER GRID POINTS AFTER
C***  VERTICAL INTERPOLATION AT THE CURRENT TIME
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:,:) :: PDBI
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:,:,:) :: TBI,QBI,UBI,VBI,Q2BI,CWMBI
C----------------------------------------------------------------------
C
C***  THESE ARE THE BOUNDARY VALUES AT THE INNER GRID POINTS AFTER
C***  VERTICAL INTERPOLATION AT THE PREVIOUS TIME
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:) :: PDBI0
C
                            R E A L
     &,ALLOCATABLE,DIMENSION(:,:,:) :: TBI0,QBI0,UBI0,VBI0
     &,                                Q2BI0,CWMBI0
C----------------------------------------------------------------------
                              I N T E G E R
     & IMI(15),JMI(15),LMI(15),IDAT(3)
     &,IHE(JM),IVE(JM)
C----------------------------------------------------------------------
                              C H A R A C T E R
     & RSTFIL*50,RESTHR*4,FNAME*80,ENVAR*50,BLANK*4
     &,BCWGTFIL*50,BCNEST*50
C
C----------------------------------------------------------------------
      DATA BLANK/'    '/
C----------------------------------------------------------------------
      NAMELIST /BCEXDATA/
     & EXBC,NBCEX,TBEGP
     &,TPH0DO,TLM0DO
     &,TBCEND
     &,IMI,JMI,LMI,TPH0DI,TLM0DI
C----------------------------------------------------------------------
C**********************************************************************
C-----------------------------------------------------------------------
C
C***  READ THE BCEXDATA NAMELIST
C
      REWIND 17
      READ(17,BCEXDATA)
C
C***  ONLY PROCEED IF EXBC IS ON
C
      IF(.NOT.EXBC)RETURN
C-----------------------------------------------------------------------
C
C***  FIND THE LARGEST BOUNDARY DIMENSIONS OF THE NESTS
C
      KBIMX=0
      LMIMX=0
      DO NB=1,NBCEX
C
        KBI=2*IMI(NB)+JMI(NB)-3
        IF(KBI.GT.KBIMX)THEN
          KBIMX=KBI
        ENDIF
C
        LMI0=LMI(NB)
        IF(LMI0.GT.LMIMX)THEN
          LMIMX=LMI0
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
C
C***  ALLOCATE BOUNDARY ARRAYS FOR THIS LARGEST SIZE
C***  TO HOLD THE HORIZONTAL INTERPOLATIONS
C
      ALLOCATE(PDBO(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(USTARBO(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(PSLPB(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(PSLPB2(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(PSLPB_L(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(PSLPB2_L(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(TSLPB(1:KBIMX,LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(TSLPB2(1:KBIMX,LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(TSLPB_L(1:KBIMX,LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(TSLPB2_L(1:KBIMX,LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(PMSLPB(1:KBIMX,LM,1:NBCEX),STAT=I)
      ALLOCATE(PMSLPB2(1:KBIMX,LM,1:NBCEX),STAT=I)
      ALLOCATE(PMSLPB_L(1:KBIMX,LM,1:NBCEX),STAT=I)
      ALLOCATE(PMSLPB2_L(1:KBIMX,LM,1:NBCEX),STAT=I)
      ALLOCATE(TBO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
      ALLOCATE(QBO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
      ALLOCATE(UBO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
      ALLOCATE(VBO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
      ALLOCATE(Q2BO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
      ALLOCATE(CWMBO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
      ALLOCATE(HTMBO(1:KBIMX,1:LM,1:NBCEX),STAT=I)
C-----------------------------------------------------------------------
C***
C***  ALLOCATE THE INTERPOLATING ARRAYS AND FILL THEM
C***
      ALLOCATE(HLATI(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(HLONI(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(VLATI(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(VLONI(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(IOUTHB(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(JOUTHB(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(IOUTHB2(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(JOUTHB2(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(IOUTVB(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(JOUTVB(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(IONEAR(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(JONEAR(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(HWGTS(1:KBIMX,4,1:NBCEX),STAT=I)
      ALLOCATE(VWGTS(1:KBIMX,4,1:NBCEX),STAT=I)
      ALLOCATE(HWGTS2(1:KBIMX,4,1:NBCEX),STAT=I)
C
      IUNW=IUNWGT
      DO NB=1,NBCEX
cccccc  IUNW=IUNW+1
C
        CLOSE(IUNW)
        WRITE(BCWGTFIL,10)NB
   10   FORMAT('eta_bcwgts.',I2.2)
        OPEN(UNIT=IUNW,FILE=BCWGTFIL,FORM='UNFORMATTED',IOSTAT=IER)
        REWIND IUNW
C
        KBI=2*IMI(NB)+JMI(NB)-3
        KBI2=2*IMI(NB)+JMI(NB)-7
        READ(IUNW)(HLATI(N,NB),N=1,KBI)
     1,           (HLONI(N,NB),N=1,KBI)
     2,           (VLATI(N,NB),N=1,KBI)
     3,           (VLONI(N,NB),N=1,KBI)
        READ(IUNW)(IOUTHB(N,NB),N=1,KBI)
     1,           (JOUTHB(N,NB),N=1,KBI)
     2,           (IOUTVB(N,NB),N=1,KBI)
     3,           (JOUTVB(N,NB),N=1,KBI)
     4,           (IONEAR(N,NB),N=1,KBI)
     5,           (JONEAR(N,NB),N=1,KBI)
        READ(IUNW)((HWGTS(N,N2,NB),N=1,KBI),N2=1,4)
     1,           ((VWGTS(N,N2,NB),N=1,KBI),N2=1,4)
        READ(IUNW)(IOUTHB2(N,NB),N=1,KBI2)
     1,           (JOUTHB2(N,NB),N=1,KBI2)
        READ(IUNW)((HWGTS2(N,N2,NB),N=1,KBI2),N2=1,4)
      ENDDO
C-----------------------------------------------------------------------
C***
C***  READ MASKS NEEDED FROM THE NHB FILE
C***
      REWIND NHB
C
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
C
      DO L=1,LM
        READ(NHB)((HTM(I,J,L),I=1,IM),J=1,JM)
      ENDDO
C
      DO L=1,LM
        READ(NHB)((VTM(I,J,L),I=1,IM),J=1,JM)
      ENDDO
C
      READ(NHB)DY,CPGFV,EN,ENT,R,PT,TDDAMP
     1,        F4D,F4Q,EF4T,DETAO,RDETA,AETAO
C----------------------------------------------------------------------
C
C***  GENERATE THE NAME OF THE CURRENT RESTRT FILE
C
C----------------------------------------------------------------------
      ENVAR=' '
      CALL GETENV("RSTFNL",ENVAR)
      CALL GETENV("tmmark",RESTHR)
      KPATH = INDEX(ENVAR,' ') -1
      IF(KPATH.LE.0) KPATH = LEN(ENVAR)
      print *,'kpath= ',kpath
C
      IF(RESTHR.EQ.'    ')THEN
        WRITE(RSTFIL,20)IHOUR
   20   FORMAT('restrt',I6.6)
      ELSE
        WRITE(RSTFIL,25)IHOUR,RESTHR
   25   FORMAT('restrt',I6.6,'.',a4)
      ENDIF
C
      KRST = INDEX(RSTFIL,' ') -1
      IF(KRST.LE.0) KRST = LEN(RSTFIL)
      print *,'krst= ',krst

C----------------------------------------------------------------------
C***
C***  READ THIS RESTART FILE AND DO THE HORIZONTAL INTERPOLATION
C***  TO THE BOUNDARIES LEVEL BY LEVEL
C***
C----------------------------------------------------------------------
C***
C***  CONVERT SINGLE LEVELS FROM RESTRT FILE TO 
C***  SINGLE LEVELS ON THE INNER BOUNDARY 
C***
      CLOSE(LRSTRT)
      IF(ENVAR(1:4).EQ.BLANK) THEN
        OPEN(UNIT=LRSTRT,FILE=RSTFIL,FORM='UNFORMATTED',IOSTAT=IER)
      ELSE
        FNAME = ENVAR(1:KPATH) // RSTFIL(1:KRST)
        OPEN(UNIT=LRSTRT,FILE=FNAME,FORM='UNFORMATTED',IOSTAT=IER)
      ENDIF
        
      IF(IER.NE.0)WRITE(LIST,*)' LRSTRT1 OPEN UNIT ERROR IER=',IER
C
      REWIND LRSTRT
      READ(LRSTRT)RUN,IDAT,IHRST,NTSD
      READ(LRSTRT)
C
      DO L=1,LM
        READ(LRSTRT)
      ENDDO
C
      READ(LRSTRT)
      READ(LRSTRT)PD,RES
      READ(LRSTRT)
C
      DO L=1,LM
        READ(LRSTRT)T
        READ(LRSTRT)Q
        READ(LRSTRT)U
        READ(LRSTRT)V
        READ(LRSTRT)Q2
        READ(LRSTRT)
        READ(LRSTRT)CWM
        READ(LRSTRT)
        READ(LRSTRT)
C
        CALL HORIZ(KBIMX,NBCEX,L,PTOP
     1,            IMI,JMI,TPH0DO,TLM0DO,TPH0DI,TLM0DI
     2,            IOUTHB,JOUTHB,IOUTVB,JOUTVB
     3,            VLATI,VLONI
     4,            HTM(1,1,L),VTM(1,1,L)
     5,            HWGTS,VWGTS,AETAO
     6,            T,Q,U,V,Q2,CWM
     7,            TBO,QBO,UBO,VBO,Q2BO,CWMBO
     8,            LMIMX,IOUTHB2,JOUTHB2,HWGTS2
     9,            PD,RES,TSLPB_L,TSLPB2_L,PMSLPB_L,PMSLPB2_L)
C
      ENDDO
C
      READ(LRSTRT)
      READ(LRSTRT)
      READ(LRSTRT)DUM,USTAR
      READ(LRSTRT)DUM,PSLP
C
      CLOSE(LRSTRT)      
C***
C***  INTERPOLATE SLP FROM THE OUTER GRID TO THE INNER GRID'S
C***  TWO BOUNDARY ROWS ONLY FOR USE IN GENERATING THE SLP ON
C***  THE INNER GRID IF THIS IS THE FINAL FORECAST HOUR
C***
C***  CALCULATE THE I-INDEX EAST-WEST INCREMENTS
C***
      DO J=1,JM
        IHE(J)=MOD(J+1,2)
      ENDDO
C
      DO 30 NB=1,NBCEX
cccc  IF(NINT(TBCEND(NB)).EQ.IHOUR)THEN
C
        KBI=2*IMI(NB)+JMI(NB)-3
        DO K=1,KBI
          IH=IOUTHB(K,NB)
          JH=JOUTHB(K,NB)
          I1=IH
          J1=JH
          I2=IH+1
          J2=JH
          I3=IH+IHE(JH)
          J3=JH-1
          I4=IH+IHE(JH)
          J4=JH+1
C
          H1=HWGTS(K,1,NB)
          H2=HWGTS(K,2,NB)
          H3=HWGTS(K,3,NB)
          H4=HWGTS(K,4,NB)
C
          PSLPB_L(K,NB)=PSLP(I1,J1)*H1+PSLP(I2,J2)*H2
     1                 +PSLP(I3,J3)*H3+PSLP(I4,J4)*H4
        ENDDO
C
        KBI2=2*IMI(NB)+JMI(NB)-7
        DO K=1,KBI2
          IH=IOUTHB2(K,NB)
          JH=JOUTHB2(K,NB)
          I1=IH
          J1=JH
          I2=IH+1
          J2=JH
          I3=IH+IHE(JH)
          J3=JH-1
          I4=IH+IHE(JH)
          J4=JH+1
C
          H1=HWGTS2(K,1,NB)
          H2=HWGTS2(K,2,NB)
          H3=HWGTS2(K,3,NB)
          H4=HWGTS2(K,4,NB)
C
          PSLPB2_L(K,NB)=PSLP(I1,J1)*H1+PSLP(I2,J2)*H2
     1                  +PSLP(I3,J3)*H3+PSLP(I4,J4)*H4
        ENDDO
C
cccc  ENDIF
C
   30 CONTINUE
C-----------------------------------------------------------------------
C***
C***  ASSIGN TO EACH POINT ON THE NEST BOUNDARIES
C***  THE NEAREST NEIGHBOR VALUES ON THE PARENT GRID
C***  FOR PD, USTAR, AND HTM
C***
      DO NB=1,NBCEX
        KBI=2*IMI(NB)+JMI(NB)-3
C
        DO K=1,KBI
          PDBO(K,NB)=PD(IONEAR(K,NB),JONEAR(K,NB))
          USTARBO(K,NB)=USTAR(IONEAR(K,NB),JONEAR(K,NB))
        ENDDO
C
        DO L=1,LM
          DO K=1,KBI
            HTMBO(K,L,NB)=HTM(IONEAR(K,NB),JONEAR(K,NB),L)
          ENDDO
        ENDDO
C
      ENDDO
C-----------------------------------------------------------------------
C
C***  ALLOCATE TEMPORARY BOUNDARY ARRAYS FOR THE NESTS
C
      ALLOCATE(PDBI0(1:KBIMX,1:NBCEX),STAT=I)
      ALLOCATE(TBI0(1:KBIMX,1:LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(QBI0(1:KBIMX,1:LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(UBI0(1:KBIMX,1:LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(VBI0(1:KBIMX,1:LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(Q2BI0(1:KBIMX,1:LMIMX,1:NBCEX),STAT=I)
      ALLOCATE(CWMBI0(1:KBIMX,1:LMIMX,1:NBCEX),STAT=I)
C----------------------------------------------------------------------
C***
C***  DO THE VERTICAL INTERPOLATION ON THE NEST BOUNDARIES
C***  AT THE CURRENT TIME
C***
      CALL PTETAE(PTOP,KBIMX,LMIMX,NBCEX
     1,           IMI,JMI,LMI
     2,           DETAO,IUDETI
     3,           PDBO,TBO,QBO,UBO,VBO,Q2BO,CWMBO
     4,           HTMBO,USTARBO
     5,           PDBI0,TBI0,QBI0,UBI0,VBI0,Q2BI0,CWMBI0)
C----------------------------------------------------------------------
C***
C***  INTERPOLATE THE TSLPB ARRAYS TO THE NEST VERTICAL DISTRIBUTION.
C***  THESE WILL ONLY BE USED FOR THE NEST SLP REDUCTION.
C***
      CALL VRT_TSLP(TSLPB_L,TSLPB2_L,PMSLPB_L,PMSLPB2_L
     1,             PSLPB_L,PSLPB2_L
     2,             IMI,JMI,LMI,LMIMX,KBIMX,NBCEX,PTOP,IUDETI)
C----------------------------------------------------------------------
C
C***  ALLOCATE FINAL BOUNDARY ARRAYS FOR THE NESTS
C
      ALLOCATE(PDBI(1:KBIMX,2,1:NBCEX),STAT=I)
      ALLOCATE(TBI(1:KBIMX,1:LMIMX,2,1:NBCEX),STAT=I)
      ALLOCATE(QBI(1:KBIMX,1:LMIMX,2,1:NBCEX),STAT=I)
      ALLOCATE(UBI(1:KBIMX,1:LMIMX,2,1:NBCEX),STAT=I)
      ALLOCATE(VBI(1:KBIMX,1:LMIMX,2,1:NBCEX),STAT=I)
      ALLOCATE(Q2BI(1:KBIMX,1:LMIMX,2,1:NBCEX),STAT=I)
      ALLOCATE(CWMBI(1:KBIMX,1:LMIMX,2,1:NBCEX),STAT=I)
C
C***  COPY THE TEMPORARY BOUNDARY FILES INTO SLOT 1
C***  OF THE FINAL ARRAYS
C
      DO NB=1,NBCEX
        KBI=2*IMI(NB)+JMI(NB)-3
        LMI0=LMI(NB)
C
        DO K=1,KBI
          PDBI(K,1,NB)=PDBI0(K,NB)
        ENDDO
C
        DO L=1,LMI0
        DO K=1,KBI
          TBI(K,L,1,NB)=TBI0(K,L,NB)
          QBI(K,L,1,NB)=QBI0(K,L,NB)
          UBI(K,L,1,NB)=UBI0(K,L,NB)
          VBI(K,L,1,NB)=VBI0(K,L,NB)
          Q2BI(K,L,1,NB)=Q2BI0(K,L,NB)
          CWMBI(K,L,1,NB)=CWMBI0(K,L,NB)
        ENDDO
        ENDDO
      ENDDO
C----------------------------------------------------------------------
C
C***  GENERATE THE NAME OF THE PREVIOUS RESTRT FILE
C
      ENVAR=' '
      CALL GETENV("RSTFNL",ENVAR)
      CALL GETENV("tmmark",RESTHR)
      KPATH = INDEX(ENVAR,' ') -1
      IF(KPATH.LE.0) KPATH = LEN(ENVAR)
      print *,'kpath= ',kpath
C
      IHOUR0=IHOUR-1
      IF(RESTHR.EQ.'    ')THEN
        WRITE(RSTFIL,40)IHOUR0
   40   FORMAT('restrt',I6.6)
      ELSE
        WRITE(RSTFIL,45)IHOUR0,RESTHR
   45   FORMAT('restrt',I6.6,'.',a4)
      ENDIF
C
      KRST = INDEX(RSTFIL,' ') -1
      IF(KRST.LE.0) KRST = LEN(RSTFIL)
      print *,'krst= ',krst
C----------------------------------------------------------------------
C***
C***  READ THIS RESTART FILE AND DO THE HORIZONTAL INTERPOLATION
C***  TO THE BOUNDARIES LEVEL BY LEVEL
C***
C----------------------------------------------------------------------
C***
C***  CONVERT SINGLE LEVELS FROM RESTRT FILE TO 
C***  SINGLE LEVELS ON THE INNER BOUNDARY 
C***
      CLOSE(LRSTRT)
      IF(ENVAR(1:4).EQ.BLANK) THEN
        OPEN(UNIT=LRSTRT,FILE=RSTFIL,FORM='UNFORMATTED',IOSTAT=IER)
      ELSE
        FNAME = ENVAR(1:KPATH) // RSTFIL(1:KRST)
        OPEN(UNIT=LRSTRT,FILE=FNAME,FORM='UNFORMATTED',IOSTAT=IER)
      ENDIF
C
      IF(IER.NE.0)WRITE(LIST,*)' LRSTRT0 OPEN UNIT ERROR IER=',IER
C
      REWIND LRSTRT
      READ(LRSTRT)RUN,IDAT,IHRST,NTSD
      READ(LRSTRT)
C
      DO L=1,LM
        READ(LRSTRT)
      ENDDO
C
      READ(LRSTRT)
      READ(LRSTRT)PD,RES
      READ(LRSTRT)
C
      DO L=1,LM
        READ(LRSTRT)T
        READ(LRSTRT)Q
        READ(LRSTRT)U
        READ(LRSTRT)V
        READ(LRSTRT)Q2
        READ(LRSTRT)
        READ(LRSTRT)CWM
        READ(LRSTRT)
        READ(LRSTRT)
C
        CALL HORIZ(KBIMX,NBCEX,L,PTOP
     1,            IMI,JMI,TPH0DO,TLM0DO,TPH0DI,TLM0DI
     2,            IOUTHB,JOUTHB,IOUTVB,JOUTVB
     3,            VLATI,VLONI
     4,            HTM(1,1,L),VTM(1,1,L)
     4,            HWGTS,VWGTS,AETAO
     5,            T,Q,U,V,Q2,CWM
     6,            TBO,QBO,UBO,VBO,Q2BO,CWMBO
     8,            LMIMX,IOUTHB2,JOUTHB2,HWGTS2
     9,            PD,RES,TSLPB,TSLPB2,PMSLPB,PMSLPB2)
C
      ENDDO
C
      READ(LRSTRT)
      READ(LRSTRT)
      READ(LRSTRT)DUM,USTAR
      READ(LRSTRT)DUM,PSLP
C
      CLOSE(LRSTRT)      
C***
C***  INTERPOLATE SLP FROM THE OUTER GRID TO THE INNER GRID'S
C***  TWO BOUNDARY ROWS ONLY FOR USE IN GENERATING THE SLP ON
C***  THE INNER GRID IF THIS IS THE FINAL FORECAST HOUR
C***
C***  CALCULATE THE I-INDEX EAST-WEST INCREMENTS
C***
      DO J=1,JM
        IHE(J)=MOD(J+1,2)
      ENDDO
C
      DO 50 NB=1,NBCEX
cccc  IF(NINT(TBCEND(NB)).EQ.IHOUR)THEN
C
        KBI=2*IMI(NB)+JMI(NB)-3
        DO K=1,KBI
          IH=IOUTHB(K,NB)
          JH=JOUTHB(K,NB)
          I1=IH
          J1=JH
          I2=IH+1
          J2=JH
          I3=IH+IHE(JH)
          J3=JH-1
          I4=IH+IHE(JH)
          J4=JH+1
C
          H1=HWGTS(K,1,NB)
          H2=HWGTS(K,2,NB)
          H3=HWGTS(K,3,NB)
          H4=HWGTS(K,4,NB)
C
          PSLPB(K,NB)=PSLP(I1,J1)*H1+PSLP(I2,J2)*H2
     1               +PSLP(I3,J3)*H3+PSLP(I4,J4)*H4
        ENDDO
C
        KBI2=2*IMI(NB)+JMI(NB)-7
C
        DO K=1,KBI2
          IH=IOUTHB2(K,NB)
          JH=JOUTHB2(K,NB)
          I1=IH
          J1=JH
          I2=IH+1
          J2=JH
          I3=IH+IHE(JH)
          J3=JH-1
          I4=IH+IHE(JH)
          J4=JH+1
C
          H1=HWGTS2(K,1,NB)
          H2=HWGTS2(K,2,NB)
          H3=HWGTS2(K,3,NB)
          H4=HWGTS2(K,4,NB)
C
          PSLPB2(K,NB)=PSLP(I1,J1)*H1+PSLP(I2,J2)*H2
     1                +PSLP(I3,J3)*H3+PSLP(I4,J4)*H4
        ENDDO
C
cccc  ENDIF
C
   50 CONTINUE
C-----------------------------------------------------------------------
C***
C***  ASSIGN TO EACH POINT ON THE NEST BOUNDARIES
C***  THE NEAREST NEIGHBOR VALUES ON THE PARENT GRID
C***  FOR PD AND USTAR
C***
      DO NB=1,NBCEX
        KBI1=2*IMI(NB)+JMI(NB)-3
C
        DO K=1,KBI1
          PDBO(K,NB)=PD(IONEAR(K,NB),JONEAR(K,NB))
          USTARBO(K,NB)=USTAR(IONEAR(K,NB),JONEAR(K,NB))
        ENDDO
C
      ENDDO
C-----------------------------------------------------------------------
C***
C***  DO THE VERTICAL INTERPOLATION ON THE NEST BOUNDARIES
C***  AT THE PREVIOUS TIME
C***
      CALL PTETAE(PTOP,KBIMX,LMIMX,NBCEX
     1,           IMI,JMI,LMI
     2,           DETAO,IUDETI
     3,           PDBO,TBO,QBO,UBO,VBO,Q2BO,CWMBO
     4,           HTMBO,USTARBO
     5,           PDBI0,TBI0,QBI0,UBI0,VBI0,Q2BI0,CWMBI0)
C----------------------------------------------------------------------
C***
C***  INTERPOLATE THE TSLPB ARRAYS TO THE NEST VERTICAL DISTRIBUTION.
C***  THESE WILL ONLY BE USED FOR THE NEST SLP REDUCTION.
C***
      CALL VRT_TSLP(TSLPB,TSLPB2,PMSLPB,PMSLPB2,PSLPB,PSLPB2
     1,             IMI,JMI,LMI,LMIMX,KBIMX,NBCEX,PTOP,IUDETI)
C----------------------------------------------------------------------
C***
C***  COMPUTE TIME TENDENCIES
C***
      RTBOCO=1./(TBOCO*3600.)
      IUBCF1=IUBCF
C
      DO 300 NB=1,NBCEX
C
      KBI=2*IMI(NB)+JMI(NB)-3
      KBI2=2*IMI(NB)+JMI(NB)-7
      DO N=1,KBI
        PDBI(N,2,NB)=(PDBI(N,1,NB)-PDBI0(N,NB))*RTBOCO
        PDBI(N,1,NB)=PDBI0(N,NB)
      ENDDO
C
      LMI1=LMI(NB)
      DO L=1,LMI1
C
        DO N=1,KBI
          TBI(N,L,2,NB)=(TBI(N,L,1,NB)-TBI0(N,L,NB))*RTBOCO
          QBI(N,L,2,NB)=(QBI(N,L,1,NB)-QBI0(N,L,NB))*RTBOCO
          UBI(N,L,2,NB)=(UBI(N,L,1,NB)-UBI0(N,L,NB))*RTBOCO
          VBI(N,L,2,NB)=(VBI(N,L,1,NB)-VBI0(N,L,NB))*RTBOCO
          Q2BI(N,L,2,NB)=(Q2BI(N,L,1,NB)-Q2BI0(N,L,NB))*RTBOCO
          CWMBI(N,L,2,NB)=(CWMBI(N,L,1,NB)-CWMBI0(N,L,NB))*RTBOCO
C
          TBI(N,L,1,NB)=TBI0(N,L,NB)
          QBI(N,L,1,NB)=QBI0(N,L,NB)
          UBI(N,L,1,NB)=UBI0(N,L,NB)
          VBI(N,L,1,NB)=VBI0(N,L,NB)
          Q2BI(N,L,1,NB)=Q2BI0(N,L,NB)
          CWMBI(N,L,1,NB)=CWMBI0(N,L,NB)
        ENDDO
      ENDDO
C----------------------------------------------------------------------
C***
C***  OPEN UNIT TO THE BCNEST FILE
C***
ccccc IUBCF1=IUBCF1+1
      LRECBC=4*(1+(1+6*LMI1)*KBI*2+(KBI+KBI2)*(LMI1+1))
      WRITE(6,*)' BCNEST RECORD LENGTH IS ',LRECBC,' BYTES'
      WRITE(BCNEST,75)NB
   75 FORMAT('bcnest.',I2.2)
      CLOSE(IUBCF1)
      OPEN(UNIT=IUBCF1,FILE=BCNEST,ACCESS='DIRECT',RECL=LRECBC
     1,    IOSTAT=IER)
C
      IF(IHOUR.EQ.1)THEN
        RUN=.TRUE.
        WRITE(IUBCF1,REC=1)RUN,IDAT,IHRST,TBOCO
      ENDIF
C
C***  FOR USE IN IDENTIFYING WHERE WE ARE IN THE NEST'S
C***  BC FILE, INCLUDE THE FORECAST HOUR WHEN WRITING OUT
C***  ALL THE BOUNDARY ARRAYS
C
      FHOUR=IHOUR-TBOCO
      NREC=IHOUR+1
      WRITE(IUBCF1,REC=NREC)FHOUR
     1,                   ((PDBI(K,N,NB),K=1,KBI),N=1,2)
     2,                   (((TBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     3,                   (((QBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     4,                   (((UBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     5,                   (((VBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     6,                  (((Q2BI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     7,                 (((CWMBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     8,                  ((TSLPB(K,L,NB),K=1,KBI),L=1,LMI1)
     9,                 ((TSLPB2(K,L,NB),K=1,KBI2),L=1,LMI1)
     O,                   (PSLPB(K,NB),K=1,KBI)
     1,                  (PSLPB2(K,NB),K=1,KBI2)
C***
C***  SINCE NEW BC CONDITIONS ARE NEEDED FOR THE LAST TIMESTEP
C***  OF THE NEST'S FORECAST, DUPLICATE THE LAST SET
C***  EXCEPT FOR THE APPROPRIATE BOUNDARY SLPs
C***
      IF(NINT(TBCEND(NB)).EQ.IHOUR)THEN
        NREC=NREC+1
        WRITE(IUBCF1,REC=NREC)FHOUR
     1,                     ((PDBI(K,N,NB),K=1,KBI),N=1,2)
     2,                     (((TBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     3,                     (((QBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     4,                     (((UBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     5,                     (((VBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     6,                    (((Q2BI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     7,                   (((CWMBI(K,L,N,NB),K=1,KBI),L=1,LMI1),N=1,2)
     8,                    ((TSLPB_L(K,L,NB),K=1,KBI),L=1,LMI1)
     9,                   ((TSLPB2_L(K,L,NB),K=1,KBI2),L=1,LMI1)
     O,                     (PSLPB_L(K,NB),K=1,KBI)
     1,                    (PSLPB2_L(K,NB),K=1,KBI2)
      ENDIF
C
      CLOSE(IUBCF1)
C
  300 CONTINUE
C----------------------------------------------------------------------
      RETURN
      END
