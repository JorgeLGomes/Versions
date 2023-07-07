C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)

      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*8   NEMO,NEMS,NEMT,NEMF
      CHARACTER*1   TAB
      DIMENSION     NEMS(250),IRPS(250),KNTS(250)
      LOGICAL       REP

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IF(ITAB.LE.0 .OR. ITAB.GT.NTBD(LUN)) GOTO 900

      REP  = .FALSE.

C  CLEAR THE RETURN VALUES
C  -----------------------

      NSEQ = 0

      DO I=1,250
      NEMS(I) = ' '
      IRPS(I) = 0
      KNTS(I) = 0
      ENDDO

C  PARSE THE TABLE D ENTRY
C  -----------------------

      NEMO = TABD(ITAB,LUN)(7:14)
      IDSC = IDND(ITAB,LUN)
      CALL UPTDD(ITAB,LUN,0,NDSC)

      IF(IDSC.LT.IFXY('300000')) GOTO 901
      IF(IDSC.GT.IFXY('363255')) GOTO 901
C     IF(NDSC.LE.0             ) GOTO 902

      DO J=1,NDSC
      IF(NSEQ+1.GT.250) GOTO 903
      CALL UPTDD(ITAB,LUN,J,IDSC)
      CALL NUMTAB(LUN,IDSC,NEMT,TAB,IRET)
      IF(TAB.EQ.'R') THEN
         IF(REP) GOTO 904
         REP = .TRUE.
         IF(IRET.LT.0) THEN
            IRPS(NSEQ+1) = 1
            KNTS(NSEQ+1) = ABS(IRET)
         ELSEIF(IRET.GT.0) THEN
            IRPS(NSEQ+1) = IRET
         ENDIF
      ELSEIF(TAB.EQ.'F') THEN
         IF(.NOT.REP) GOTO 904
         IRPS(NSEQ+1) = IRET
         REP = .FALSE.
      ELSEIF(TAB.EQ.'D'.OR.TAB.EQ.'C') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         NEMS(NSEQ) = NEMT
      ELSEIF(TAB.EQ.'B') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         IF(NEMT(1:1).EQ.'.') THEN
            CALL UPTDD(ITAB,LUN,J+1,IDSC)
            CALL NUMTAB(LUN,IDSC,NEMF,TAB,IRET)
            CALL RSVFVM(NEMT,NEMF)
            IF(TAB.NE.'B') GOTO 906
         ENDIF
         NEMS(NSEQ) = NEMT
      ELSE
         GOTO 905
      ENDIF
      ENDDO

      RETURN
900   CALL BORT('NEMTBD - ITAB NOT IN TABLE D   '                )
901   CALL BORT('NEMTBD - BAD DESCRIPTOR VALUE: '          //NEMO)
902   CALL BORT('NEMTBD - ZERO LENGTH SEQUENCE: '          //NEMO)
903   CALL BORT('NEMTBD - TOO MANY DESCRIPTORS IN SEQ: '   //NEMO)
904   CALL BORT('NEMTBD - REPLICATOR OUT OF ORDER IN SEQ: '//NEMO)
905   CALL BORT('NEMTBD - BAD DESCRIPTOR IN SEQUENCE: '    //NEMO)
906   CALL BORT('NEMTBD - FOLLOWING VALUE NOT FROM TABLEB:'//NEMF)
      END
