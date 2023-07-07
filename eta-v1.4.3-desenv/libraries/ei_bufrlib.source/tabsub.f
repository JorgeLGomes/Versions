C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TABSUB(LUN,NEMO)

      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /TABCCC/ ICDW,ICSC,ICRV

      CHARACTER*10 TAG
      CHARACTER*8  NEMO,NEMS,NEM
      CHARACTER*3  TYP
      CHARACTER*1  TAB
      DIMENSION    NEM(250,10),IRP(250,10),KRP(250,10)
      DIMENSION    DROP(10),JMP0(10),NODL(10),NTAG(10,2)
      LOGICAL      DROP

      DATA MAXLIM /10/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE MNEMONIC
C  ------------------

      CALL NEMTAB(LUN,NEMO,IDN,TAB,ITAB)
      IF(TAB.NE.'D') GOTO 900

C  STORE A SUBSET NODE AND JUMP/LINK THE TREE
C  ------------------------------------------

      CALL INCTAB(NEMO,'SUB',NODE)
      JUMP(NODE) = NODE+1
      JMPB(NODE) = 0
      LINK(NODE) = 0
      IBT (NODE) = 0
      IRF (NODE) = 0
      ISC (NODE) = 0

      CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,1),IRP(1,1),KRP(1,1))
      NTAG(1,1) = 1
      NTAG(1,2) = NSEQ
      JMP0(1)   = NODE
      LIMB      = 1

      ICDW = 0
      ICSC = 0
      ICRV = 0

C  THIS LOOP RESOLVES ENTITIES IN A SUBSET BY EMULATING RECURSION
C  --------------------------------------------------------------

1     DO N=NTAG(LIMB,1),NTAG(LIMB,2)

      NTAG(LIMB,1) = N+1
      NODL(LIMB)   = NTAB+1
      DROP(LIMB)   = N.EQ.NTAG(LIMB,2)

      CALL NEMTAB(LUN,NEM(N,LIMB),IDN,TAB,ITAB)
      NEMS = NEM(N,LIMB)

C  SPECIAL TREATMENT FOR CERTAIN OPERATOR DESCRIPTORS (TAB=C)
C  ----------------------------------------------------------

      IF(TAB.EQ.'C') THEN
         NODL(LIMB) = NTAB
         READ(NEMS,'(3X,I3)') IYYY
         IF(ITAB.EQ.1) THEN
            ICDW = IYYY-128
            IF(IYYY.EQ.0) ICDW = 0
         ELSEIF(ITAB.EQ.2) THEN
            ICSC = IYYY-128
            IF(IYYY.EQ.0) ICSC = 0
         ENDIF
      ELSE
         IREP = IRP(N,LIMB)
         IKNT = KRP(N,LIMB)
         JUM0 = JMP0(LIMB)
         CALL TABENT(LUN,NEMS,TAB,ITAB,IREP,IKNT,JUM0)
      ENDIF

      IF(TAB.EQ.'D') THEN
         LIMB = LIMB+1
         IF(LIMB.GT.MAXLIM) GOTO 901
         CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,LIMB),IRP(1,LIMB),KRP(1,LIMB))
         NTAG(LIMB,1) = 1
         NTAG(LIMB,2) = NSEQ
         JMP0(LIMB)   = NTAB
         GOTO 1
      ELSEIF(DROP(LIMB)) THEN
2        LINK(NODL(LIMB)) = 0
         LIMB = LIMB-1
         IF(LIMB.EQ.0 ) THEN
            IF(ICDW.NE.0) GOTO 902
            IF(ICSC.NE.0) GOTO 903
            RETURN
         ENDIF
         IF(DROP(LIMB)) GOTO 2
         LINK(NODL(LIMB)) = NTAB+1
         GOTO 1
      ELSEIF(TAB.NE.'C') THEN
         LINK(NODL(LIMB)) = NTAB+1
      ENDIF

      ENDDO

      CALL BORT('TABSUB - SHOULD NOT GET HERE               ')
900   CALL BORT('TABSUB - SUBSET NODE NOT IN TABLE D: '//NEMO)
901   CALL BORT('TABSUB - TOO MANY LIMBS                    ')
902   CALL BORT('TABSUB - CHANGE DATA WIDTH OPERATOR NOT CANCELED')
903   CALL BORT('TABSUB - CHANGE DATA SCALE OPERATOR NOT CANCELED')
      END
