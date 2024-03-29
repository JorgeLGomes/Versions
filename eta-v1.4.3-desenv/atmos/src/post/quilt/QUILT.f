                  PROGRAM QUILT
C
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C                .      .    .
C MAIN PROGRAM:  QUILT       ONE RESTRT FROM MANY
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-04-29
C
C ABSTRACT:  THIS PROGRAM READS IN THE SMALL RESTRT FILES 
C            WRITTEN BY EACH PROCESSOR AND CREATES THE SINGLE
C            GLOBAL RESTRT FILE APPENDED WITH THE MISCELLANEOUS
C            QUANTITIES NEEDED FOR THE PROFILE JOB
C
C PROGRAM HISTORY LOG:
C   99-04-29  T BLACK - ORIGINATOR
C   99-07-27  E ROGERS - MODIFIED SO THAT SMALL RESTART FILES DO NOT
C                        HAVE TO BE IN THE LOCAL DIRECTORY
C   00-01-20  JIM TUCCILLO - MPI IMPLEMENTATION
C   00-02-20  JIM TUCCILLO - READ A SINGLE DIRECT ACCESS FILE FROM ETAFCST.X
C
C USAGE:  MAIN PROGRAM
C
C   INPUT ARGUMENT LIST:
C     NONE
C
C   OUTPUT ARGUMENT LIST:
C     NONE
C
C   INPUT FILES:  NONE 
C
C   OUTPUT FILES:  NONE
C
C   SUBPROGRAMS CALLED:
C     UNIQUE: 
C            SLP
C
C   EXIT STATES:
C     COND =   0 - NORMAL EXIT
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE : IBM SP
C
C$$$
C
C     THIS CODE ASSUMES THAT NSOIL IS GE TO 4. IF THIS IS NOT TRUE,
C     THE CODE WILL STOP. THE EQUIVALENCING IS THE PROBLEM.
C
C-----------------------------------------------------------------------
      INCLUDE "parmeta"
      INCLUDE "parmsoil"
      INCLUDE "mpif.h"
C-----------------------------------------------------------------------
      INCLUDE "PARA.comm"
      INCLUDE "BUFFER.comm"
C-----------------------------------------------------------------------
                              P A R A M E T E R
     & (NPES=INPES*JNPES,LNIP=IM/INPES,LNJP=JM/JNPES)
                              P A R A M E T E R
     & (LB=2*IM+JM-3)
C-----------------------------------------------------------------------
C
       REAL DUM1(IM,JM),DUM2(IM,JM),DUM3(IM,JM),DUM4(IM,JM)
       REAL DUM5(IM,JM),DUM6(IM,JM),DUM7(IM,JM)
       REAL DUMS(IM,JM,NSOIL)
       EQUIVALENCE ( DUM1(1,1), DUMS(1,1,1) )
       EQUIVALENCE ( DUM2(1,1), DUMS(1,1,2) )
       EQUIVALENCE ( DUM3(1,1), DUMS(1,1,3) )
       EQUIVALENCE ( DUM4(1,1), DUMS(1,1,4) )
C
C-----------------------------------------------------------------------
      REAL, ALLOCATABLE ::
     & PDOMG(:,:),RESOMG(:,:),PD(:,:),RES(:,:),FIS(:,:)
     &,RSWIN(:,:),RSWOUT(:,:),TG(:,:),Z0(:,:),AKMS(:,:)
     &,CZEN(:,:),AKHS(:,:),THS(:,:),QS(:,:),TWBS(:,:)
     &,QWBS(:,:),HBOT(:,:),CFRACL(:,:),THZ0(:,:),QZ0(:,:)
     &,UZ0(:,:),VZ0(:,:),USTAR(:,:),HTOP(:,:),CFRACM(:,:)
     &,SNO(:,:),SI(:,:),CLDEFI(:,:),RF(:,:),PSLP(:,:)
     &,CUPPT(:,:),CFRACH(:,:),SOILTB(:,:),SFCEXC(:,:)
     &,SMSTAV(:,:),SMSTOT(:,:),GRNFLX(:,:),PCTSNO(:,:)
     &,RLWIN(:,:),RADOT(:,:),CZMEAN(:,:),SIGT4(:,:)
     &,U00(:,:),SR(:,:),PREC(:,:),ACPREC(:,:),ACCLIQ(:,:)
     &,CUPREC(:,:),ACFRCV(:,:),ACFRST(:,:),SFCSHX(:,:)
     &,ACSNOW(:,:),ACSNOM(:,:),SSROFF(:,:),BGROFF(:,:)
     &,SFCLHX(:,:),SUBSHX(:,:),SNOPCX(:,:),SFCUVX(:,:)
     &,SFCEVP(:,:),POTEVP(:,:),ASWIN(:,:),ASWOUT(:,:)
     &,ASWTOA(:,:),ALWIN(:,:),ALWOUT(:,:),ALWTOA(:,:)
     &,TH10(:,:),Q10(:,:),U10(:,:),V10(:,:),TSHLTR(:,:)
     &,QSHLTR(:,:),PSHLTR(:,:),CMC(:,:),POTFLX(:,:)
     &,TLMIN(:,:),TLMAX(:,:),ALBEDO(:,:)

      INTEGER, ALLOCATABLE:: HBM2(:,:)
C
      REAL UL(2*LM),DETA(LM)
C
      REAL, ALLOCATABLE ::
     & OMGALF(:,:,:),T(:,:,:),Q(:,:,:),U(:,:,:)
     &,V(:,:,:),Q2(:,:,:),TTND(:,:,:),CWM(:,:,:)
     &,TRAIN(:,:,:),TCUCN(:,:,:)
     &,RSWTT(:,:,:),RLWTT(:,:,:),T0(:,:,:),Q0(:,:,:)
     &,P0(:,:),RSWTOA(:,:),RLWTOA(:,:),SM(:,:)
C
      REAL, ALLOCATABLE ::
     & SMC(:,:,:),STC(:,:,:),SH2O(:,:,:)
                              R E A L
     & PDB(LB,2),TB(LB,LM,2),QB(LB,LM,2),UB(LB,LM,2),VB(LB,LM,2)
Cmp
     &,Q2B(LB,LM,2),CWMB(LB,LM,2)
Cmp
C
C-----------------------------------------------------------------------
      INTEGER IDAT(3)
C
      INTEGER, ALLOCATABLE ::
     & LC(:,:),NCFRCV(:,:),NCFRST(:,:)
C-----------------------------------------------------------------------
                              L O G I C A L
     & RUN,FIRST
C-----------------------------------------------------------------------
                              C H A R A C T E R
     & RSTFIL1*50,RSTFIL2*50,RESTHR*4,LABEL*32
     &,FNAME*80,ENVAR*50,BLANK*4
C
       LOGICAL LME,SIGMA,RESTRT,SINGLRST,SUBPOST,
     &	       NEST,HYDRO,SPLINE
C-----------------------------------------------------------------------
      DATA LRSTRT1/21/,LRSTRT2/61/,NHB/12/,BLANK/'    '/
C-----------------------------------------------------------------------
C
	REAL TSHDE(999999),SPL(LSM)

      NAMELIST /FCSTDATA/
     & TSTART,TEND,TCP,RESTRT,SINGLRST,SUBPOST,NMAP,TSHDE,SPL
     &,NPHS,NCNVC,NRADSH,NRADLH,NTDDMP
     &,TPREC,THEAT,TCLOD,TRDSW,TRDLW,TSRFC
     &,NEST,HYDRO,SPLINE

C-----------------------------------------------------------------------
      real*8 timef
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C      call start()
      CALL MPI_FIRST

C
C	GET SIGMA INFO FROM NHB FILE, SPLINE FROM RESTART later?
C

        open(unit=NHB,file='cnst.file',access='sequential',
     +          form='unformatted')
      REWIND NHB
      READ(NHB)NFCST,NBC,LIST,DT,IDTAD,SIGMA
	close(NHB)

	write(6,*) 'NHB read, SIGMA= ', SIGMA

       REWIND 11
       READ(11,FCSTDATA)
	write(6,*) 'discoverd that NEST,HYDRO,SPLINE= ', NEST,HYDRO,
     +	SPLINE


c     IF(ME.EQ.0)THEN
c       CALL W3TAGB('ETA_QUILT',1999,0267,0082,'NP22')
c     ENDIF
C
      IF(NSOIL.LT.4)THEN
        print *, ' NSOIL IS LESS THAN 4. CHANGE THE EQUIVALENCES'
        print *, ' STOPPING'
        stop
      ENDIF
C
      IF(ME.EQ.0)THEN
        LME=.TRUE.
      ELSE
        LME=.FALSE.
      ENDIF
C
      btim=timef()
C
      ALLOCATE(PDOMG(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RESOMG(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(OMGALF(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(PD(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RES(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(FIS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(T(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(Q(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(U(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(V(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(Q2(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(TTND(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(CWM(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(TRAIN(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(TCUCN(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(RSWIN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RSWOUT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(TG(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(Z0(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(AKMS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CZEN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(AKHS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(THS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(QS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(TWBS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(QWBS(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(HBOT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CFRACL(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(THZ0(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(QZ0(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(UZ0(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(VZ0(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(USTAR(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(HTOP(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SNO(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SI(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CLDEFI(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RF(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(PSLP(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CUPPT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CFRACH(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CFRACM(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SOILTB(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SFCEXC(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SMSTAV(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SMSTOT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(GRNFLX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(PCTSNO(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RLWIN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RADOT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CZMEAN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SIGT4(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(U00(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(LC(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SR(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(PREC(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ACPREC(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ACCLIQ(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(CUPREC(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ACFRCV(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(NCFRCV(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ACFRST(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(NCFRST(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ACSNOW(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ACSNOM(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SSROFF(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(BGROFF(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SFCSHX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SFCLHX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SUBSHX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SNOPCX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SFCUVX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SFCEVP(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(POTEVP(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ASWIN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ASWOUT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ASWTOA(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ALWIN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ALWOUT(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ALWTOA(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(TH10(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(Q10(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(U10(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(V10(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(TSHLTR(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(QSHLTR(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(PSHLTR(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SMC(MY_ISD:MY_IED,MY_JSD:MY_JED,1:NSOIL))
      ALLOCATE(SH2O(MY_ISD:MY_IED,MY_JSD:MY_JED,1:NSOIL))
      ALLOCATE(CMC(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(STC(MY_ISD:MY_IED,MY_JSD:MY_JED,1:NSOIL))
      ALLOCATE(POTFLX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(TLMIN(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(TLMAX(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RSWTT(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(RLWTT(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(RSWTOA(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(RLWTOA(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(P0(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(T0(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(Q0(MY_ISD:MY_IED,MY_JSD:MY_JED,1:LM))
      ALLOCATE(HBM2(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(SM(MY_ISD:MY_IED,MY_JSD:MY_JED))
      ALLOCATE(ALBEDO(MY_ISD:MY_IED,MY_JSD:MY_JED))
C
C-----------------------------------------------------------------------
C***
C***  READ IN THE FIRST FORECAST HOUR, THE TOTAL NUMBER OF OUTPUT TIMES
C***  TO PROCESS, AND THE TIME INCREMENT IN HOURS BETWEEN OUTPUT TIMES
C***
      IF(ME.EQ.0)THEN
        READ(5,*)IHOUR,KOUNT,NINC
      ENDIF

C
      CALL MPI_BCAST(IHOUR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IRTN)
      CALL MPI_BCAST(KOUNT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IRTN)
      CALL MPI_BCAST(NINC ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IRTN)
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IRTN)
C
      IF(ME.EQ.0)WRITE(6,*)IHOUR,KOUNT,NINC
C-----------------------------------------------------------------------
C***
C***  LOOP OVER ALL THE OUTPUT TIMES
C***
C-----------------------------------------------------------------------
      DO 500 NK=1,KOUNT
C-----------------------------------------------------------------------
C***
C***  GENERATE THE NAME OF THE INPUT RESTRT FILE FROM THIS PE
C***
      ENVAR = ' '
      CALL GETENV("tmmark",RESTHR)
      CALL GETENV("RSHOME",ENVAR)
      KPATH = INDEX(ENVAR,' ') -1
      IF(KPATH.LE.0) KPATH = LEN(ENVAR)
      print *,'kpath= ',kpath
C
      IF(RESTHR.EQ.'    ')THEN
c       WRITE(RSTFIL1,50)IHOUR,IPE
c  50   FORMAT('restrt',I3.3,'.',I3.3)
        WRITE(RSTFIL1,50)IHOUR
   50   FORMAT('restrt',I6.6,'.quilt')
      ELSE
c       WRITE(RSTFIL1,55)IHOUR,IPE,RESTHR
c  55   FORMAT('restrt',I3.3,'.',I3.3,'.',a4)
        WRITE(RSTFIL1,55)IHOUR,RESTHR
   55   FORMAT('restrt',I6.6,'.quilt.',a4)
      ENDIF
C
      KRST = INDEX(RSTFIL1,' ') -1
      IF(KRST.LE.0) KRST = LEN(RSTFIL1)
      print *,'krst= ',krst
C***
C***  OPEN UNIT TO THE LOCAL RESTART FILE
C***
      CLOSE(LRSTRT1)
C
      IF(ENVAR(1:4).EQ.BLANK) THEN
         FNAME = RSTFIL1(1:KRST)
      ELSE
         FNAME = ENVAR(1:KPATH) // RSTFIL1(1:KRST)
      END IF
c      OPEN(UNIT=LRSTRT1,FILE=FNAME,FORM='UNFORMATTED',IOSTAT=IER)
       OPEN(UNIT=LRSTRT1,FILE=FNAME,FORM='UNFORMATTED',IOSTAT=IER,
     *   access='direct',recl=4)
       read(LRSTRT1,rec=1) irecl
       close(LRSTRT1)
      if ( irecl/4 .gt. ibufmax ) then
         print *, ' IBUFMAX in parmbuf is too small'
         print *, ' It must be at least ',irecl/4, ' stopping'
         stop
      end if
c
       OPEN(UNIT=LRSTRT1,FILE=FNAME,FORM='UNFORMATTED',IOSTAT=IER,
     *   access='direct',recl=irecl)
      print *, ' IRECL, FNAME = ',irecl, fname
C-----------------------------------------------------------------------
      DO 200 IPE=JSTA(ME), JEND(ME)
C-----------------------------------------------------------------------
C***
      MY_IS_GLB = MY_IS_GLB_A(IPE)
      MY_IE_GLB = MY_IE_GLB_A(IPE)
      MY_JS_GLB = MY_JS_GLB_A(IPE)
      MY_JE_GLB = MY_JE_GLB_A(IPE)
C
      is = MY_IS_GLB
      ie = MY_IE_GLB
      js = MY_JS_GLB
      je = MY_JE_GLB
      len_ch = (ie-is+1) * (je-js+1)
C***
C***  READ FROM THE LOCAL RESTRT FILE INTO THE APPROPRIATE LOCATIONS
C***  OF THE GLOBAL ARRAYS
C***
C
	write(6,*) 'reading from LRSTRT1 ', IPE

      read(LRSTRT1,rec=ipe+1,iostat=ier) (buf(i),i=1,irecl/4)
      if ( ier .ne. 0 ) then
         print *, ' error from direct access read = ',ier
      end if
c
C     EXTRACT RECORD LENGTH
      call decoal(idum,-1)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)RUN,IDAT,IHRST,NTSD,LABEL
      call decoal(run,1)
      call decoal(IDAT,3)
      call decoal(IHRST,1)
      call decoal(NTSD,1)
      call decoal(LABEL,8)
	write(6,*) 'RUN,IDAT,IHRST,NTSD= ', RUN,IDAT,IHRST,NTSD
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((PDOMG(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((RESOMG(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(pdomg(is:ie,js:je),len_ch)
      call decoal(resomg(is:ie,js:je),len_ch)

C-----------------------------------------------------------------------
      DO L=1,LM
c       READ(LRSTRT1)((OMGALF(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                               J=MY_JS_GLB,MY_JE_GLB)
      call decoal(omgalf(is:ie,js:je,l),len_ch)
      ENDDO
C-----------------------------------------------------------------------
c     READ(LRSTRT1)RUN,IDAT,IHRST,NTSD,LABEL
c    1,            FIRST,IOUT,NSHDE
        call decoal(RUN,1)
        call decoal(IDAT,3)
        call decoal(IHRST,1)
        call decoal(NTSD,1)
        call decoal(LABEL,8)
        call decoal(FIRST,1)
        call decoal(IOUT,1)
        call decoal(NSHDE,1)

C-----------------------------------------------------------------------
c     READ(LRSTRT1)((PD(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                       J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((RES(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                        J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((FIS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                        J=MY_JS_GLB,MY_JE_GLB)
      call decoal(pd(is:ie,js:je),len_ch)
      call decoal(res(is:ie,js:je),len_ch)
      call decoal(fis(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
C
C***  DUMMY READ FOR BC ARRAYS
C
      CALL DECOAL(PDB,LB*2)
      CALL DECOAL(TB,LB*LM*2)
      CALL DECOAL(QB,LB*LM*2)
      CALL DECOAL(UB,LB*LM*2)
      CALL DECOAL(VB,LB*LM*2)
      CALL DECOAL(Q2B,LB*LM*2)
      CALL DECOAL(CWMB,LB*LM*2)

c     READ(LRSTRT1)PDB,TB,QB,UB,VB
c     READ(LRSTRT1)
C-----------------------------------------------------------------------
      DO L=1,LM
c       READ(LRSTRT1)((T(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
      call decoal(t(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((Q(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
      call decoal(q(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((U(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
      call decoal(u(is:ie,js:je,l),len_ch)
	if (L .eq. 1) then
	do J=js,je
	do I=is,ie
	write(6,*) 'I,J,U: ', I,J,U(I,J,L)
	enddo
	enddo
	endif
c       READ(LRSTRT1)((V(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
      call decoal(v(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((Q2(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(q2(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((TTND(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                             J=MY_JS_GLB,MY_JE_GLB)
      call decoal(ttnd(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((CWM(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                            J=MY_JS_GLB,MY_JE_GLB)
      call decoal(cwm(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((TRAIN(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                              J=MY_JS_GLB,MY_JE_GLB)
      call decoal(train(is:ie,js:je,l),len_ch)
c       READ(LRSTRT1)((TCUCN(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                              J=MY_JS_GLB,MY_JE_GLB)
      call decoal(tcucn(is:ie,js:je,l),len_ch)
      ENDDO
C-----------------------------------------------------------------------
c     READ(LRSTRT1)RUN,IDAT,IHRST,NTSD,LABEL
c    1,          ((RSWIN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    2                        J=MY_JS_GLB,MY_JE_GLB)
c    3,          ((RSWOUT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    4                         J=MY_JS_GLB,MY_JE_GLB)
c    5,          ((TG(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    6                     J=MY_JS_GLB,MY_JE_GLB)
c    7,          ((Z0(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    8                     J=MY_JS_GLB,MY_JE_GLB)
c    9,          ((AKMS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    O                       J=MY_JS_GLB,MY_JE_GLB)
c    1,          ((CZEN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    2                       J=MY_JS_GLB,MY_JE_GLB)
      call decoal(RUN,1)
      call decoal(IDAT,3)
      call decoal(IHRST,1)
      call decoal(NTSD,1)
      call decoal(LABEL,8)
      call decoal(rswin(is:ie,js:je),len_ch)
      call decoal(rswout(is:ie,js:je),len_ch)
      call decoal(tg(is:ie,js:je),len_ch)
      call decoal(z0(is:ie,js:je),len_ch)
      call decoal(akms(is:ie,js:je),len_ch)
      call decoal(czen(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((AKHS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                         J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((THS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                        J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((QS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                       J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((TWBS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                         J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((QWBS(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                         J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((HBOT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                         J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((CFRACL(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(akhs(is:ie,js:je),len_ch)
      call decoal(ths(is:ie,js:je),len_ch)
      call decoal(qs(is:ie,js:je),len_ch)
      call decoal(twbs(is:ie,js:je),len_ch)
      call decoal(qwbs(is:ie,js:je),len_ch)
      call decoal(hbot(is:ie,js:je),len_ch)
      call decoal(cfracl(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((THZ0(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                         J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((QZ0(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                        J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((UZ0(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                       J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((VZ0(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                         J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((USTAR(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                         J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((HTOP(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                         J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((CFRACM(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(thz0(is:ie,js:je),len_ch)
      call decoal(qz0(is:ie,js:je),len_ch)
      call decoal(uz0(is:ie,js:je),len_ch)
      call decoal(vz0(is:ie,js:je),len_ch)
      call decoal(ustar(is:ie,js:je),len_ch)
      call decoal(htop(is:ie,js:je),len_ch)
      call decoal(cfracm(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((SNO(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                        J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((SI(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                        J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((CLDEFI(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                           J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((RF(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                       J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((PSLP(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                         J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((CUPPT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((CFRACH(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(sno(is:ie,js:je),len_ch)
      call decoal(si(is:ie,js:je),len_ch)
      call decoal(cldefi(is:ie,js:je),len_ch)
      call decoal(rf(is:ie,js:je),len_ch)
      call decoal(pslp(is:ie,js:je),len_ch)
      call decoal(cuppt(is:ie,js:je),len_ch)
      call decoal(cfrach(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((SOILTB(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((SFCEXC(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((SMSTAV(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                           J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((SMSTOT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                           J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((GRNFLX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                           J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((PCTSNO(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(soiltb(is:ie,js:je),len_ch)
      call decoal(sfcexc(is:ie,js:je),len_ch)
      call decoal(smstav(is:ie,js:je),len_ch)
      call decoal(smstot(is:ie,js:je),len_ch)
      call decoal(grnflx(is:ie,js:je),len_ch)
      call decoal(pctsno(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((RLWIN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((RADOT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                          J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((CZMEAN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                           J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((SIGT4(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                          J=MY_JS_GLB,MY_JE_GLB)
      call decoal(rlwin(is:ie,js:je),len_ch)
      call decoal(radot(is:ie,js:je),len_ch)
      call decoal(czmean(is:ie,js:je),len_ch)
      call decoal(sigt4(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((U00(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                        J=MY_JS_GLB,MY_JE_GLB)
c    2,              UL
c    3,            ((LC(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    4                       J=MY_JS_GLB,MY_JE_GLB)
c    5,            ((SR(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    6                       J=MY_JS_GLB,MY_JE_GLB)
      call decoal(u00(is:ie,js:je),len_ch)
      call decoal(ul,2*lm)
      call decoal(lc(is:ie,js:je),len_ch)
      call decoal(sr(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)RUN,IDAT,IHRST,NTSD,LABEL
c    1,            ((PREC(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    2                         J=MY_JS_GLB,MY_JE_GLB)
c    3,            ((ACPREC(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    4                           J=MY_JS_GLB,MY_JE_GLB)
c    5,            ((ACCLIQ(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    6                           J=MY_JS_GLB,MY_JE_GLB)
c    7,            ((CUPREC(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    8                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(RUN,1)
      call decoal(IDAT,3)
      call decoal(IHRST,1)
      call decoal(NTSD,1)
      call decoal(LABEL,8)
      call decoal(prec(is:ie,js:je),len_ch)
      call decoal(acprec(is:ie,js:je),len_ch)
      call decoal(accliq(is:ie,js:je),len_ch)
      call decoal(cuprec(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((ACFRCV(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((NCFRCV(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((ACFRST(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                           J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((NCFRST(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(acfrcv(is:ie,js:je),len_ch)
      call decoal(ncfrcv(is:ie,js:je),len_ch)
      call decoal(acfrst(is:ie,js:je),len_ch)
      call decoal(ncfrst(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((ACSNOW(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((ACSNOM(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((SSROFF(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((BGROFF(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(acsnow(is:ie,js:je),len_ch)
      call decoal(acsnom(is:ie,js:je),len_ch)
      call decoal(ssroff(is:ie,js:je),len_ch)
      call decoal(bgroff(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((SFCSHX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((SFCLHX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((SUBSHX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                           J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((SNOPCX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                           J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((SFCUVX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                           J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((SFCEVP(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((POTEVP(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(SFCSHX(is:ie,js:je),len_ch)
      call decoal(SFCLHX(is:ie,js:je),len_ch)
      call decoal(SUBSHX(is:ie,js:je),len_ch)
      call decoal(SNOPCX(is:ie,js:je),len_ch)
      call decoal(SFCUVX(is:ie,js:je),len_ch)
      call decoal(SFCEVP(is:ie,js:je),len_ch)
      call decoal(POTEVP(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((ASWIN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                          J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((ASWOUT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((ASWTOA(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                           J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((ALWIN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                          J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((ALWOUT(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                           J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((ALWTOA(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(ASWIN(is:ie,js:je),len_ch)
      call decoal(ASWOUT(is:ie,js:je),len_ch)
      call decoal(ASWTOA(is:ie,js:je),len_ch)
      call decoal(ALWIN(is:ie,js:je),len_ch)
      call decoal(ALWOUT(is:ie,js:je),len_ch)
      call decoal(ALWTOA(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)ARDSW,ARDLW,ASRFC,AVRAIN,AVCNVC
      call decoal(ARDSW,1)
      call decoal(ARDLW,1)
      call decoal(ASRFC,1)
      call decoal(AVRAIN,1)
      call decoal(AVCNVC,1)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((TH10(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                         J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((Q10(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                        J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((U10(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                        J=MY_JS_GLB,MY_JE_GLB)
c    6,            ((V10(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    7                        J=MY_JS_GLB,MY_JE_GLB)
c    8,            ((TSHLTR(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    9                           J=MY_JS_GLB,MY_JE_GLB)
c    O,            ((QSHLTR(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((PSHLTR(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                           J=MY_JS_GLB,MY_JE_GLB)
      call decoal(TH10(is:ie,js:je),len_ch)
      call decoal(Q10(is:ie,js:je),len_ch)
      call decoal(U10(is:ie,js:je),len_ch)
      call decoal(V10(is:ie,js:je),len_ch)
      call decoal(TSHLTR(is:ie,js:je),len_ch)
      call decoal(QSHLTR(is:ie,js:je),len_ch)
      call decoal(PSHLTR(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)(((SMC(I,J,N),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB),
c    2                           N=1,NSOIL)
      call decoal(SMC(is:ie,js:je,1:nsoil),len_ch*nsoil)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((CMC(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                        J=MY_JS_GLB,MY_JE_GLB)
      call decoal(CMC(is:ie,js:je),len_ch)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)(((STC(I,J,N),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB),
c    2                           N=1,NSOIL)
      call decoal(STC(is:ie,js:je,1:nsoil),len_ch*nsoil)
      call decoal(SH2O(is:ie,js:je,1:nsoil),len_ch*nsoil)
C-----------------------------------------------------------------------
c     READ(LRSTRT1)((POTFLX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    1                           J=MY_JS_GLB,MY_JE_GLB)
c    2,            ((TLMIN(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    3                          J=MY_JS_GLB,MY_JE_GLB)
c    4,            ((TLMAX(I,J),I=MY_IS_GLB,MY_IE_GLB),
c    5                          J=MY_JS_GLB,MY_JE_GLB)
c    6,            ACUTIM,ARATIM,APHTIM
c    7,            NHEAT,NPHS,NCNVC,NPREC,NRDSW,NRDLW,NSRFC
c    8,            TPH0D,TLM0D,RESTRT
      call decoal(ALBEDO(is:ie,js:je),len_ch)
      call decoal(POTFLX(is:ie,js:je),len_ch)
      call decoal(TLMIN(is:ie,js:je),len_ch)
      call decoal(TLMAX(is:ie,js:je),len_ch)
      call decoal(ACUTIM,1)
      call decoal(ARATIM,1)
      call decoal(APHTIM,1)
      call decoal(NHEAT,1)
      call decoal(NPHS,1)
      call decoal(NCNVC,1)
      call decoal(NPREC,1)
      call decoal(NRDSW,1)
      call decoal(NRDLW,1)
      call decoal(NSRFC,1)
      call decoal(TPH0D,1)
      call decoal(TLM0D,1)
      call decoal(RESTRT,1)
C-----------------------------------------------------------------------
      DO L=1,LM
c       READ(LRSTRT1)((RSWTT(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                              J=MY_JS_GLB,MY_JE_GLB)
c       READ(LRSTRT1)((RLWTT(I,J,L),I=MY_IS_GLB,MY_IE_GLB),
c    1                              J=MY_JS_GLB,MY_JE_GLB)
      call decoal(RSWTT(is:ie,js:je,l),len_ch)
      call decoal(RLWTT(is:ie,js:je,l),len_ch)
      ENDDO

	do l=1,lm
	call decoal(t0(is:ie,js:je,l),len_ch)
	call decoal(q0(is:ie,js:je,l),len_ch)
	enddo

	call decoal(p0(is:ie,js:je),len_ch)
	call decoal(hbot(is:ie,js:je),len_ch)
	call decoal(htop(is:ie,js:je),len_ch)
	call decoal(rswtoa(is:ie,js:je),len_ch)
	call decoal(rlwtoa(is:ie,js:je),len_ch)

C	mp additions in chkout
	
	call decoal(hbm2(is:ie,js:je),len_ch)
	call decoal(sm(is:ie,js:je),len_ch)
	call decoal(spl(1:lsm),lsm)
	call decoal(deta(1:lm),lm)
	call decoal(pt,1)
	call decoal(spline,1)

C	write(6,*) 'deta= ', deta
C	write(6,*) 'pt= ', pt
	write(6,*) 'sigma= ', sigma
	write(6,*) 'decoaled that spline is: ', spline

C-----------------------------------------------------------------------
C***
C***  CLOSE THE LOCAL RESTRT FILE
C***
C-----------------------------------------------------------------------
  200 CONTINUE
      CLOSE(LRSTRT1)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C***  FINISHED ASSEMBLING THE GLOBAL RESTRT FILE
C***  FROM ALL OF THE LOCAL ONES.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C***
C*** BEFORE WRITING OUT THE RESTRT FILE, COMPUTE THE MSLP
C***

	IF(SIGMA)THEN
	if(spline)then
	        print*,'calling slpsigspline'
        call slpsigspline(PD,FIS,T,Q,SPL,LSM
     1,            DETA,PT,PSLP)
       else
        print*,'calling slpsig'
        CALL SLPSIG(PD,FIS,SM,T,Q,CWM,HBM2,U00,SPL,LSM
     1,            UL,DETA,PT,PSLP)
       end if
      ELSE
       print*,'calling regular slp ', NHB
       CALL SLP(NHB,PD,RES,FIS,T,Q,NTSD,PSLP)
      END IF
	write(6,*) 'back from SLP'

C	write(6,*) 'PSLP near N boundary: '
	do J=JM,JM-3,-1
C	write(6,*) (PSLP(I,J),I=15,20)
	enddo
C	write(6,*) 'PSLP near S boundary: '
	do J=3,1,-1
C	write(6,*) (PSLP(I,J),I=15,20)
	enddo

	blrg=-999999999.
	bsml=999999999.

	do J=1,JM
	do I=1,IM
	if (pslp(i,j) .gt. blrg) then
	 blrg=pslp(i,j)
	 ilrg=i
	 jlrg=j
	endif
	if (pslp(i,j) .lt. bsml) then
		bsml=pslp(i,j)
		isml=i
		jsml=j
	endif
	enddo
	enddo

	write(6,*) 'extremes in pslp: ', bsml,blrg

C	write(6,*) 'surrounding large point: ', ilrg,jlrg
	
	do j=jlrg+1,jlrg-1,-1
C	write(6,*) (pslp(i,j),i=ilrg-1,ilrg+1)
	enddo

C	write(6,*) 'surrounding min point: ', isml,jsml
	
	do j=jsml+1,jsml-1,-1
C	write(6,*) (pslp(i,j),i=isml-1,isml+1)
	enddo

C
C-----------------------------------------------------------------------
C***  WRITE OUT THE GLOBAL FILE.
C-----------------------------------------------------------------------
C***
C***  GENERATE THE NAME OF THE GLOBAL OUTPUT RESTRT FILE
C***
      ENVAR=' '
      CALL GETENV("RSTFNL",ENVAR)
      CALL GETENV("tmmark",RESTHR)
      KPATH = INDEX(ENVAR,' ') -1
      IF(KPATH.LE.0) KPATH = LEN(ENVAR)
      print *,'kpath= ',kpath
C
      IF(RESTHR.EQ.'    ')THEN
        WRITE(RSTFIL2,280)IHOUR
  280   FORMAT('restrt',I6.6)
      ELSE
        WRITE(RSTFIL2,285)IHOUR,RESTHR
  285   FORMAT('restrt',I6.6,'.',a4)
      ENDIF
C
      KRST = INDEX(RSTFIL2,' ') -1
      IF(KRST.LE.0) KRST = LEN(RSTFIL2)
      print *,'krst= ',krst
C***
C***  OPEN UNIT TO THE GLOBAL RESTART FILE
C***
      CLOSE(LRSTRT2)
C
      IF(ENVAR(1:4).EQ.BLANK) THEN
       OPEN(UNIT=LRSTRT2,FILE=RSTFIL2,FORM='UNFORMATTED',IOSTAT=IER)
      ELSE
       FNAME = ENVAR(1:KPATH) // RSTFIL2(1:KRST)
       OPEN(UNIT=LRSTRT2,FILE=FNAME,FORM='UNFORMATTED',IOSTAT=IER)
      ENDIF
C-----------------------------------------------------------------------
      IF ( LME ) WRITE(LRSTRT2)RUN,IDAT,IHRST,NTSD,LABEL,ihour
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)PDOMG,RESOMG
      CALL COLLECT(PDOMG,DUM1)
      CALL COLLECT(RESOMG,DUM2)
      IF ( LME ) WRITE(LRSTRT2) DUM1, DUM2
C-----------------------------------------------------------------------
      DO L=1,LM
cwas    WRITE(LRSTRT2)((OMGALF(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(OMGALF(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
      ENDDO
C-----------------------------------------------------------------------
      IF ( LME ) WRITE(LRSTRT2)RUN,IDAT,IHRST,NTSD,LABEL,ihour,
     1              FIRST,IOUT,NSHDE
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)PD,RES,FIS
      CALL COLLECT(PD,DUM1)
      CALL COLLECT(RES,DUM2)
      CALL COLLECT(FIS,DUM3)
      IF ( LME ) WRITE(LRSTRT2) DUM1, DUM2, DUM3
C-----------------------------------------------------------------------
Cwrong      IF ( LME ) WRITE(LRSTRT2)PDB,TB,QB,UB,VB
      IF ( LME ) WRITE(LRSTRT2)PDB,TB,QB,UB,VB,Q2B,CWMB
C-----------------------------------------------------------------------
      DO L=1,LM
cwas    WRITE(LRSTRT2)((T(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(T(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((Q(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(Q(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((U(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(U(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((V(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(V(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((Q2(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(Q2(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((TTND(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(TTND(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((CWM(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(CWM(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((TRAIN(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(TRAIN(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((TCUCN(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(TCUCN(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
      ENDDO
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)RUN,IDAT,IHRST,NTSD,LABEL
cwas 1,            RSWIN,RSWOUT,TG,Z0,AKMS,CZEN
      CALL COLLECT(RSWIN,DUM1)
      CALL COLLECT(RSWOUT,DUM2)
      CALL COLLECT(TG,DUM3)
      CALL COLLECT(Z0,DUM4)
      CALL COLLECT(AKMS,DUM5)
      CALL COLLECT(CZEN,DUM6)
      IF ( LME ) WRITE(LRSTRT2)RUN,IDAT,IHRST,NTSD,LABEL,ihour
     1,           DUM1,DUM2,DUM3,DUM4,DUM5,DUM6
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)AKHS,THS,QS,TWBS,QWBS,HBOT,CFRACL
      CALL COLLECT(AKHS,DUM1)
      CALL COLLECT(THS,DUM2)
      CALL COLLECT(QS,DUM3)
      CALL COLLECT(TWBS,DUM4)
      CALL COLLECT(QWBS,DUM5)
      CALL COLLECT(HBOT,DUM6)
      CALL COLLECT(CFRACL,DUM7)
      IF ( LME ) WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)THZ0,QZ0,UZ0,VZ0,USTAR,HTOP,CFRACM
      CALL COLLECT(THZ0,DUM1)
      CALL COLLECT(QZ0,DUM2)
      CALL COLLECT(UZ0,DUM3)
      CALL COLLECT(VZ0,DUM4)
      CALL COLLECT(USTAR,DUM5)
      CALL COLLECT(HTOP,DUM6)
      CALL COLLECT(CFRACM,DUM7)
      IF ( LME ) WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)SNO,SI,CLDEFI,RF,PSLP,CUPPT,CFRACH
      CALL COLLECT(SNO,DUM1)
      CALL COLLECT(SI,DUM2)
      CALL COLLECT(CLDEFI,DUM3)
      CALL COLLECT(RF,DUM4)
      CALL COLLECT(PSLP,DUM5)
      CALL COLLECT(CUPPT,DUM6)
      CALL COLLECT(CFRACH,DUM7)
      IF ( LME ) WRITE(LRSTRT2) DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)SOILTB,SFCEXC,SMSTAV,SMSTOT,GRNFLX,PCTSNO
      CALL COLLECT(SOILTB,DUM1)
      CALL COLLECT(SFCEXC,DUM2)
      CALL COLLECT(SMSTAV,DUM3)
      CALL COLLECT(SMSTOT,DUM4)
      CALL COLLECT(GRNFLX,DUM5)
      CALL COLLECT(PCTSNO,DUM6)
      IF ( LME )WRITE(LRSTRT2) DUM1,DUM2,DUM3,DUM4,DUM5,DUM6
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)RLWIN,RADOT,CZMEAN,SIGT4
      CALL COLLECT(RLWIN,DUM1)
      CALL COLLECT(RADOT,DUM2)
      CALL COLLECT(CZMEAN,DUM3)
      CALL COLLECT(SIGT4,DUM4)
      IF ( LME ) WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)U00,UL,LC,SR
      CALL COLLECT(U00,DUM1)
      CALL COLLECT(LC,DUM2)
      CALL COLLECT(SR,DUM3)
      IF ( LME ) WRITE(LRSTRT2)DUM1,UL,DUM2,DUM3
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)RUN,IDAT,IHRST,NTSD,LABEL
cwas 1,             PREC,ACPREC,ACCLIQ,CUPREC
      CALL COLLECT(PREC,DUM1)
      CALL COLLECT(ACPREC,DUM2)
      CALL COLLECT(ACCLIQ,DUM3)
      CALL COLLECT(CUPREC,DUM4)
      IF ( LME ) WRITE(LRSTRT2)RUN,IDAT,IHRST,NTSD,LABEL,ihour
     1,             DUM1,DUM2,DUM3,DUM4
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)ACFRCV,NCFRCV,ACFRST,NCFRST
      CALL COLLECT(ACFRCV,DUM1)
      CALL COLLECT(NCFRCV,DUM2)
      CALL COLLECT(ACFRST,DUM3)
      CALL COLLECT(NCFRST,DUM4)
      IF ( LME ) WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)ACSNOW,ACSNOM,SSROFF,BGROFF
      CALL COLLECT(ACSNOW,DUM1)
      CALL COLLECT(ACSNOM,DUM2)
      CALL COLLECT(SSROFF,DUM3)
      CALL COLLECT(BGROFF,DUM4)
      IF ( LME ) WRITE(LRSTRT2) DUM1,DUM2,DUM3,DUM4
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)SFCSHX,SFCLHX,SUBSHX,SNOPCX,SFCUVX,SFCEVP,POTEVP
      CALL COLLECT(SFCSHX,DUM1)
      CALL COLLECT(SFCLHX,DUM2)
      CALL COLLECT(SUBSHX,DUM3)
      CALL COLLECT(SNOPCX,DUM4)
      CALL COLLECT(SFCUVX,DUM5)
      CALL COLLECT(SFCEVP,DUM6)
      CALL COLLECT(POTEVP,DUM7)
      IF ( LME ) WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)ASWIN,ASWOUT,ASWTOA,ALWIN,ALWOUT,ALWTOA
      CALL COLLECT(ASWIN,DUM1)
      CALL COLLECT(ASWOUT,DUM2)
      CALL COLLECT(ASWTOA,DUM3)
      CALL COLLECT(ALWIN,DUM4)
      CALL COLLECT(ALWOUT,DUM5)
      CALL COLLECT(ALWTOA,DUM6)
      IF ( LME )WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4,DUM5,DUM6
C-----------------------------------------------------------------------
      IF ( LME ) WRITE(LRSTRT2)ARDSW,ARDLW,ASRFC,AVRAIN,AVCNVC
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)TH10,Q10,U10,V10,TSHLTR,QSHLTR,PSHLTR
      CALL COLLECT(TH10,DUM1)
      CALL COLLECT(Q10,DUM2)
      CALL COLLECT(U10,DUM3)
      CALL COLLECT(V10,DUM4)
      CALL COLLECT(TSHLTR,DUM5)
      CALL COLLECT(QSHLTR,DUM6)
      CALL COLLECT(PSHLTR,DUM7)
      IF ( LME )WRITE(LRSTRT2)DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)SMC
      DO L = 1, NSOIL
         CALL COLLECT(SMC(:,:,L), DUMS(:,:,L))
      END DO
      IF ( LME ) WRITE(LRSTRT2) DUMS
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)CMC
      CALL COLLECT(CMC,DUM1)
      IF ( LME ) WRITE(LRSTRT2) DUM1
C-----------------------------------------------------------------------
cwas  WRITE(LRSTRT2)STC
      DO L = 1, NSOIL
         CALL COLLECT(STC(:,:,L), DUMS(:,:,L))
      END DO
      IF ( LME ) WRITE(LRSTRT2) DUMS
      DO L=1,NSOIL
        CALL COLLECT(SH2O(:,:,L), DUMS(:,:,L))
      ENDDO
      IF(LME)WRITE(LRSTRT2) DUMS
C
      CALL COLLECT(ALBEDO,DUM1)
      IF(LME)WRITE(LRSTRT2) DUM1

C-----------------------------------------------------------------------
C***
C***  AUXILIARY QUANTITIES NEEDED FOR PROFILE JOB BUT NOT BY THE POST
C***
cwas  WRITE(LRSTRT2)POTFLX,TLMIN,TLMAX
cwas 1,             ACUTIM,ARATIM,APHTIM
cwas 2,             NHEAT,NPHS,NCNVC,NPREC,NRDSW,NRDLW,NSRFC
cwas 3,             TPH0D,TLM0D,RESTRT
      CALL COLLECT(POTFLX,DUM1)
      CALL COLLECT(TLMIN,DUM2)
      CALL COLLECT(TLMAX,DUM3)
      IF ( LME ) WRITE(LRSTRT2) DUM1, DUM2, DUM3
     1,             ACUTIM,ARATIM,APHTIM
     2,             NHEAT,NPHS,NCNVC,NPREC,NRDSW,NRDLW,NSRFC
     3,             TPH0D,TLM0D,RESTRT
C-----------------------------------------------------------------------
      DO L=1,LM
cwas    WRITE(LRSTRT2)((RSWTT(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(RSWTT(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
cwas    WRITE(LRSTRT2)((RLWTT(I,J,L),I=1,IM),J=1,JM)
        CALL COLLECT(RLWTT(:,:,L),DUM1)
        IF ( LME ) WRITE(LRSTRT2) DUM1
      ENDDO


Cmp
	DO L=1,LM
	CALL COLLECT(T0(:,:,L),DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1
	CALL COLLECT(Q0(:,:,L),DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1
	ENDDO
Cmp

	CALL COLLECT(P0,DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1
	CALL COLLECT(HBOT,DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1
	CALL COLLECT(HTOP,DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1
	CALL COLLECT(RSWTOA,DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1
	CALL COLLECT(RLWTOA,DUM1)
	IF ( LME ) WRITE(LRSTRT2) DUM1

Cmp


C-----------------------------------------------------------------------
C***
C***  CLOSE THE GLOBAL RESTRT FILE
C***
      CLOSE(LRSTRT2)
C-----------------------------------------------------------------------
      IHOUR=IHOUR+NINC
  500 CONTINUE
C-----------------------------------------------------------------------
      tot_tim=timef()-btim
      if(me.eq.0)then
        write(6,*)' tot_tim=',tot_tim*1.e-3
C        call summary()
      endif
C
c     IF(ME.EQ.0)THEN
c       CALL W3TAGE('ETA_QUILT')
c     ENDIF
C
      CALL MPI_LAST
C
      STOP0
      END
C***********************************************************************
      subroutine decoal(a,len_ch)
      real a(*)
      include "BUFFER.comm"
      if ( len_ch .lt. 0 ) then
         ip = 0
      end if
      do i = 1, abs(len_ch)
         ip = ip + 1
         a(i) = buf(ip)
      end do
      end

