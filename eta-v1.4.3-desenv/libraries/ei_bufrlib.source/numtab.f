C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NUMTAB(LUN,IDN,NEMO,TAB,IRET)

      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*6   ADN30,CID
      CHARACTER*3   TYPS
      CHARACTER*1   REPS,TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      NEMO = ' '
      IRET = 0
      TAB = ' '

C  LOOK FOR A REPLICATOR OR A REPLICATOR FACTOR
C  --------------------------------------------

      IF(IDN.GE.IDNR(1,1) .AND. IDN.LE.IDNR(1,2)) THEN
         TAB  = 'R'
         IRET = -MOD(IDN,256)
         RETURN
      ENDIF

      DO I=2,5
      IF(IDN.EQ.IDNR(I,1)) THEN
         TAB  = 'R'
         IRET = I
         RETURN
      ELSEIF(IDN.EQ.IDNR(I,2)) THEN
         TAB  = 'F'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE D
C  -----------------------

      DO I=1,NTBD(LUN)
      IF(IDN.EQ.IDND(I,LUN)) THEN
         NEMO = TABD(I,LUN)(7:14)
         TAB  = 'D'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE B
C  -----------------------

      DO I=1,NTBB(LUN)
      IF(IDN.EQ.IDNB(I,LUN)) THEN
         NEMO = TABB(I,LUN)(7:14)
         TAB  = 'B'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE C
C  -----------------------

      CID = ADN30(IDN,6)
      CID = CID(1:3)
      IF(CID.EQ.'201' .OR. CID.EQ.'202' .OR. CID.EQ.'206') THEN
         NEMO = ADN30(IDN,6)
         TAB  = 'C'
         IRET = MOD(ID,10)
         RETURN
      ENDIF

      RETURN
      END
