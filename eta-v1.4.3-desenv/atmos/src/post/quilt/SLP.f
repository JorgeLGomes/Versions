      SUBROUTINE SLP(NHB,PDx,RESx,FISx,Tx,Qx,NTSD,PSLPx)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  SLP         SLP REDUCTION
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-04-22
C
C ABSTRACT:  THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE
C            REDUCTION USING EITHER THE MESINGER RELAXATION
C            METHOD OR THE STANDARD NCEP REDUCTION
C
C PROGRAM HISTORY LOG:
C   99-04-22  T BLACK - ORIGINATOR
C   00-01-10  JIM TUCCILLO - MPI VERSION
C
C USAGE:  CALL SLP FROM PROGRAM POST0
C
C   INPUT ARGUMENT LIST:
C     NHB  - UNIT NUMBER FOR READING THE NHB FILE
C     PD   - SFC PRESSURE MINUS PTOP
C     RES  - RECIPROCAL OF ETA AT THE GROUND
C     FIS  - SURFACE GEOPOTENTIAL
C     T    - TEMPERATURE 
C     Q    - SPECIFIC HUMIDITY
C     NTSD - THE TIMESTEP
C
C   OUTPUT ARGUMENT LIST:
C     PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C             NONE
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C-----------------------------------------------------------------------
      INCLUDE "parmeta"
      INCLUDE "PARA.comm"
      INCLUDE 'mpif.h'
C-----------------------------------------------------------------------
                            P A R A M E T E R
     & (LP1=LM+1,IMJM=IM*JM-JM/2,JM2=JM-2)
                            P A R A M E T E R
     & (NRLX1=250,NRLX2=100,KSLPD=1
     &, OVERRC=1.75,AD05=OVERRC*0.05,CFT0=OVERRC-1.
     &, ROG=287.04/9.8)
C-----------------------------------------------------------------------
                            R E A L
     & PD(IM,JM),RES(IM,JM),FIS(IM,JM)
C 
       REAL, ALLOCATABLE :: HTM(:,:,:)
                            R E A L
     & PSLP(IM,JM),TTV(IM,JM)
     &,PDSL1(IM,JM),PBI(IM,JM),SLPX(IM,JM)
                            R E A L
     & DETA(LM),RDETA(LM),AETA(LM),F4Q2(LM),ETA(LP1),DFL(LP1)
C-----------------------------------------------------------------------
                            I N T E G E R
     & KMNTM,LMH(IM,JM),LMV(IM,JM)
                            I N T E G E R
     & IMNT(IMJM),JMNT(IMJM)
                            I N T E G E R
     & IHE(JM),IHW(JM),IVE(JM),IVW(JM)
C-----------------------------------------------------------------------
                            L O G I C A L
     & SIGMA,STDRD,NO_REDUCE
C
       REAL PDx(MY_ISD:MY_IED,MY_JSD:MY_JED)
       REAL RESx(MY_ISD:MY_IED,MY_JSD:MY_JED)
       REAL FISx(MY_ISD:MY_IED,MY_JSD:MY_JED)
       REAL PSLPx(MY_ISD:MY_IED,MY_JSD:MY_JED)
       REAL Tx(MY_ISD:MY_IED,MY_JSD:MY_JED,LM)
       REAL Qx(MY_ISD:MY_IED,MY_JSD:MY_JED,LM)
C
       REAL, ALLOCATABLE :: T(:,:,:), Q (:,:,:)
C
       REAL DUM ( IM, JM )
       EQUIVALENCE ( DUM, SLPX )
C
C-----------------------------------------------------------------------

C
      CALL COLLECT(PDx,PD)
      CALL COLLECT(RESx,RES)
      CALL COLLECT(FISx,FIS)
C
      IF ( ME .EQ. 0 ) THEN
C***
C***  READ IN THE ARRAYS AND CONSTANTS THAT ARE NEEDED
C***
	open(unit=NHB,file='cnst.file',access='sequential',
     +		form='unformatted')
      REWIND NHB
C
	write(6,*) 'to NHB file '
      READ(NHB)NFCST,NBC,LIST,DT,IDTAD,SIGMA
	write(6,*) 'from NHB file: ', NFCST, NBC, DT, IDTAD,SIGMA
      READ(NHB)LMH
      READ(NHB)LMV
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
      READ(NHB)
C
      STDRD=.FALSE.
      NO_REDUCE=.FALSE.
	write(6,*) 'LMH,LMV sample: ', LMH(10,10),LMV(10,10)
      IF(SIGMA) STDRD=.TRUE.
C
      DO L=1,LM
cwas    READ(NHB)((HTM(I,J,L),I=1,IM),J=1,JM)
        READ(NHB)((DUM(I,J),I=1,IM),J=1,JM)
        DO J = 1, JM
           DO I = 1, IM
              IF ( DUM(I,J) .LT. 0.5 ) GOTO 666
           END DO
        END DO
C       IF WE GET TO HERE, WE HAVE ALL ATMOSPHERE
        GOTO 667
666     CONTINUE
C       IF WE GET HERE, WE HAVE FOUND A NON_ATM POINT
        LHMNT = L
        GOTO 668
667     CONTINUE
      ENDDO
      LHMNT = LM+1
      GO TO 669
668   CONTINUE
      BACKSPACE NHB
      ALLOCATE(HTM(IM,JM,LHMNT:LM))
      DO L = LHMNT, LM
         READ(NHB)((HTM(I,J,L),I=1,IM),J=1,JM)
      END DO
669   CONTINUE
C
      DO L=1,LM
        READ(NHB)
      ENDDO
C
      READ(NHB)DY,CPGFV,EN,ENT,R,PT,TDDAMP  
     1,        F4D,F4Q,EF4T,DETA,RDETA,AETA,F4Q2,ETA,DFL
C-----------------------------------------------------------------------
	write(6,*) 'past NHB reads'
C***
C***  CALCULATE THE I-INDEX EAST-WEST INCREMENTS
C***
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
C!$OMP parallel do 
      DO J=1,JM
      DO I=1,IM
        PSLP(I,J)=0.
        TTV(I,J)=0.
      ENDDO
      ENDDO
C
C!$OMP parallel do 
      doout110: DO J=1,JM
      doin110: DO I=1,IM
      PDSL1(I,J)=RES(I,J)*PD(I,J)
      PSLP(I,J)=PD(I,J)+PT
      PBI (I,J)=PSLP(I,J)
      END DO doin110
      END DO doout110
C
C***  CALCULATE SEA LEVEL PRESSURE FOR PROFILES (AND POSSIBLY
C***  FOR POSTING BY POST PROCESSOR).
C
C***  "STDRD" REFERS TO THE "STANDARD" SLP REDUCTION SCHEME.
C***  THIS IS THE ONLY SCHEME AT PRESENT AVAILABLE FOR A SIGMA=.TRUE.
C***  ETA MODEL RUN.
C
cwas  IF(SIGMA)THEN
cwas    STDRD=.TRUE.
cwas    GO TO 400
cwas  ENDIF
cwas  IF(STDRD)GO TO 400
C--------------------------------------------------------------------
C***
C***  WE REACH THIS LINE IF WE WANT THE MESINGER ETA SLP REDUCTION
C***  BASED ON RELAXATION TEMPERATURES.  THE FIRST STEP IS TO
C***  FIND THE HIGHEST LAYER CONTAINING MOUNTAINS.
C***
c     DO 210 L=LM,1,-1
C
c     DO J=1,JM
c     DO I=1,IM
c       IF(HTM(I,J,L).LT.0.5)GO TO 210
c     ENDDO
c     ENDDO
c
c     LHMNT=L+1
c     GO TO 220
c 210 CONTINUE
C
  220 IF(LHMNT.EQ.LP1)THEN
cwas    GO TO 430
        NO_REDUCE = .TRUE.
      ENDIF
C***
C
      END IF !   END OF IF TEST OM ME
C
	write(6,*) 'STDRD pre MPI_BCAST ', STDRD
      CALL MPI_BCAST(NO_REDUCE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(STDRD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LHMNT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
	write(6,*) 'STDRD post MPI_BCAST ', STDRD
C
C     COLLECT T AND Q
C
      IF ( .NOT. NO_REDUCE) THEN
C
      ALLOCATE(T(IM,JM,LHMNT-1:LM))
      IF ( .NOT. STDRD ) ALLOCATE(Q(IM,JM,LHMNT-1:LM))
      DO K = LHMNT-1, LM
         CALL COLLECT(Tx(:,:,K),T(:,:,K))
         IF ( .NOT. STDRD ) THEN
            CALL COLLECT(Qx(:,:,K),Q(:,:,K))
         END IF
      END DO
C
      END IF
C
	write(6,*) 'to this if test on ME, ME= ', ME
      IF ( ME .EQ. 0 ) THEN
C
	write(6,*) 'right here, NO_REDUCE= ', NO_REDUCE
      IF ( NO_REDUCE ) GOTO 430
	write(6,*) ' STDRD check...goto 400 ', STDRD
      IF ( STDRD ) GOTO 400
C
C***  NOW GATHER THE ADDRESSES OF ALL THE UNDERGROUND POINTS.
C***
C!$OMP parallel do private(kmn,kount)
c     DO 250 L=LHMNT,LM
c     KMN=0
c     KMNTM(L)=0
c     KOUNT=0
c     DO 240 J=3,JM-2
c     DO 240 I=2,IM-1
c     KOUNT=KOUNT+1
cwas  IMNT(KOUNT,L)=0
cwas  JMNT(KOUNT,L)=0
c     IF(HTM(I,J,L).GT.0.5)GO TO 240
c     KMN=KMN+1
c     IMNT(KMN,L)=I
c     JMNT(KMN,L)=J
c 240 CONTINUE
c     KMNTM(L)=KMN
c 250 CONTINUE
C
C***  AS THE FIRST GUESS, SET THE UNDERGROUND TEMPERATURES EQUAL
C***  TO THE VALUE GIVING PD/ETAS+PT AS SEA LEVEL PRESSURE.
C
c     IF(NTSD.EQ.1)THEN
c       KMM=KMNTM(LM)
c!$OMP parallel do private(i,j,lmap1,tgss),shared(t)
c       DO 260 KM=1,KMM
c       I=IMNT(KM,LM)
c       J=JMNT(KM,LM)
c       TGSS=FIS(I,J)/(R*ALOG((PDSL1(I,J)+PT)/(PD(I,J)+PT)))
c       LMAP1=LMH(I,J)+1
c       DO 260 L=LMAP1,LM
c       T(I,J,L)=TGSS
c 260   CONTINUE
c     ENDIF
c
c     IF(NTSD.EQ.1)THEN
        LL = LM
        DO J=3,JM-2
        DO I=2,IM-1
        IF(HTM(I,J,LL).LE.0.5) THEN
        TGSS=FIS(I,J)/(R*ALOG((PDSL1(I,J)+PT)/(PD(I,J)+PT)))
        LMAP1=LMH(I,J)+1
        DO 260 L=LMAP1,LM
        T(I,J,L)=TGSS
  260   CONTINUE
        END IF
       END DO 
       END DO
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
C!$OMP parallel do private(dposp,i,j,kmm,lmst,pbin,phbi,phti,prbin,prtin,
C!$OMP*                    ptin,tinit,trtv,ttv)
      DO 300 L=LHMNT,LM
C
      KMN=0
      KMNTM=0
      doout240: DO J=3,JM-2
      doin240: DO  I=2,IM-1
      IF(HTM(I,J,L).GT.0.5) CYCLE doin240
      KMN=KMN+1
      IMNT(KMN)=I
      JMNT(KMN)=J
      END DO doin240
      END DO doout240
      KMNTM=KMN

C
      TTV(1:IM,1:JM)=T(1:IM,1:JM,L)
C
C     FOR GRID BOXES NEXT TO MOUNTAINS REPLACE ttv BY AN "EQUIVALENT"
C     TV, ONE WHICH CORRESPONDS TO THE CHANGE IN P BETWEEN REFERENCE
C     INTERFACE GEOPOTENTIALS, INSTEAD OF BETWEEN LAYER INTERFACES
C
      DO J=3,JM-2
      DO I=2,IM-1
        IF(HTM(I,J,L).GT.0.5.AND.
     2     HTM(I+IHW(J),J-1,L)*HTM(I+IHE(J),J-1,L)
     3    *HTM(I+IHW(J),J+1,L)*HTM(I+IHE(J),J+1,L)
     4    *HTM(I-1     ,J  ,L)*HTM(I+1     ,J  ,L)
     5    *HTM(I       ,J-2,L)*HTM(I       ,J+2,L).LT.0.5)THEN
          LMST=LMH(I,J)
C***
C***  FIND P AT THE REFERENCE INTERFACE GEOPOTENTIAL AT THE BOTTOM
C***
          PBIN=PT+PD(I,J)
          PHBI=DFL(LMST+1)
          DO LI=LMST,1,-1
            PTIN=PBIN-DETA(LI)*PD(I,J)*RES(I,J)
            TRTV=2.*R*T(I,J,LI)*(1.+0.608*Q(I,J,LI))
            PHTI=PHBI+TRTV*(PBIN-PTIN)/(PBIN+PTIN)
            IF(PHTI.GE.DFL(L+1))GO TO 273
            PBIN=PTIN
            PHBI=PHTI
          ENDDO
  273     DPOSP=(PHTI-DFL(L+1))/TRTV
          PRBIN=(1.+DPOSP)/(1.-DPOSP)*PTIN
C***
C***  FIND P AT THE REFERENCE INTERFACE GEOPOTENTIAL AT THE TOP
C***
          PBIN=PT+PD(I,J)
          PHBI=DFL(LMST+1)
C
          DO LI=LMST,1,-1
            PTIN=PBIN-DETA(LI)*PD(I,J)*RES(I,J)
            TRTV=2.*R*T(I,J,LI)*(1.+0.608*Q(I,J,LI))
            PHTI=PHBI+TRTV*(PBIN-PTIN)/(PBIN+PTIN)
            IF(PHTI.GE.DFL(L))GO TO 275
            PBIN=PTIN
            PHBI=PHTI
          ENDDO
C
  275     DPOSP=(PHTI-DFL(L))/TRTV
          PRTIN=(1.+DPOSP)/(1.-DPOSP)*PTIN
C
          TTV(I,J)=(DFL(L)-DFL(L+1))/(2.*R)*(PRBIN+PRTIN)/(PRBIN-PRTIN)
        ENDIF
      ENDDO
      ENDDO
C
      KMM=KMNTM
C
      DO 285 N=1,NRLX
c     IF(L.EQ.LM)DLTMX=0.
      DO 280 KM=1,KMM
      I=IMNT(KM)
      J=JMNT(KM)
      TINIT=TTV(I,J)
      TTV(I,J)=AD05*(4.*(TTV(I+IHW(J),J-1)+TTV(I+IHE(J),J-1)
     1                  +TTV(I+IHW(J),J+1)+TTV(I+IHE(J),J+1))
     2                  +TTV(I-1,J)       +TTV(I+1,J)
     3                  +TTV(I,J-2)       +TTV(I,J+2))
     4                  -CFT0*TTV(I,J)
C
c     IF(L.EQ.LM)THEN
c       DLTT=ABS(TTV(I,J)-TINIT)
c       IF(DLTT.GT.DLTMX)DLTMX=DLTT
c     ENDIF
  280 CONTINUE
C
c     IF(L.EQ.LM)WRITE(51,2802)N,DLTMX
c2802 FORMAT(' ',I5,E10.3)
C
  285 CONTINUE
C
      DO 290 KM=1,KMM
      I=IMNT(KM)
      J=JMNT(KM)
      T(I,J,L)=TTV(I,J)
  290 CONTINUE
  300 CONTINUE
C----------------------------------------------------------------
C***
C***  CALCULATE THE SEA LEVEL PRESSURE AS PER THE NEW SCHEME.
C***
C     VALUES FOR IMNT AND JMNT ARE FOR LAYER LM - THIS IS WHAT WE WANT
      KMM=KMNTM
C!$OMP parallel do private(dposp,i,j,lmap1,pbin,ptin),shared(pslp)
      DO 320 KM=1,KMM
      I=IMNT(KM)
      J=JMNT(KM)
      LMAP1=LMH(I,J)+1
      PBIN=PT+PD(I,J)
C
      DO L=LMAP1,LM
        PTIN=PBIN
        DPOSP=(DFL(L)-DFL(L+1))/(2.*R*T(I,J,L))
        PBIN=(1.+DPOSP)/(1.-DPOSP)*PTIN
        if (J .eq. 198 .and. ( I .ge. 67 .and. I .le. 67) ) then
	write(6,263) L,T(I,J,L),PTIN,PBIN,DPOSP
	write(6,262) DFL(L),DFL(L+1),(1.+DPOSP)/(1.-DPOSP)
	endif
      ENDDO
C

  262   format('DFL(L),DFL(L+1),factor ++  ', f9.2,1x,f9.2,1x,f9.5)
  263   format('L, T, PTIN, PBIN, DPOSP :: ', I3,1x,
     +                                          f6.2,1x,2(f7.0,1x),f8.5)
C

      PSLP(I,J)=PBIN
  320 CONTINUE
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
      write(6,*) 'doing standard reduction!!!'
      doout410: DO J=1,JM
      doin410: DO I=1,IM
      IF(FIS(I,J).GE.1.)THEN
        LMA=LMH(I,J)
        ALPP1=ALOG(PDSL1(I,J)*ETA(LMA+1)+PT)
        SLOP=0.0065*ROG*T(I,J,LMA)
        IF(SLOP.LT.0.50)THEN
          SLPP=ALPP1+FIS(I,J)/(R*T(I,J,LMA))
        ELSE
          TTT=-(ALOG(PDSL1(I,J)*ETA(LMA)+PT)+ALPP1)
     1         *SLOP*0.50+T(I,J,LMA)
          SLPP=(-TTT+SQRT(TTT*TTT+2.*SLOP*
     1          (FIS(I,J)/R+
     2          (TTT+0.50*SLOP*ALPP1)*ALPP1)))/SLOP
        ENDIF
        PSLP(I,J)=EXP(SLPP)
      ENDIF
      END DO doin410
      END DO doout410
C
C****************************************************************
C     AT THIS POINT WE HAVE A SEA LEVEL PRESSURE FIELD BY
C     EITHER METHOD.  5-POINT AVERAGE THE FIELD ON THE E-GRID.
C****************************************************************
C
  430 CONTINUE
C
C!$OMP parallel do 
      SLPX(1:IM,1:JM)=PSLP(1:IM,1:JM)
C
      DO 480 KS=1,KSLPD
C
C!$OMP parallel do private(ihh2)
      doout460: DO J=3,JM2
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
C
C!$OMP parallel do
      DO J=1,JM
      DO I=1,IM
        PSLP(I,J)=SLPX(I,J)
      ENDDO
      ENDDO
C
  480 CONTINUE
C
      END IF    ! END IF OF ME
C
      CALL DIST(PSLP,PSLPx)
      IF ( .NOT. STDRD .AND. .NOT. NO_REDUCE ) THEN
         DO K = LHMNT, LM
            CALL DIST(T(:,:,K),Tx(:,:,K))
         END DO
      END IF
      IF (ALLOCATED(T)) DEALLOCATE(T)
      IF (ALLOCATED(Q)) DEALLOCATE(Q)
      IF (ALLOCATED(HTM)) DEALLOCATE(HTM)
C----------------------------------------------------------------
	write(6,*) 'leaving 4'
      RETURN
      END
