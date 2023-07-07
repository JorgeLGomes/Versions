      SUBROUTINE SNO8GET(SNODEP,IYR,IMN,IDY,  INSNOAF,LIST,LUSAF)
C
C USAGE: Program to read in USAF gribbed snow/ice data on 1/8th bedient
C polar stereographic grid.  The gribbing program has already reversed
C the coordinates in the N-S direction so now the (1,1) point is at
C the 'lower left' corner.
C
      INTEGER KPDS(25), KGDS(22), JPDS(25), JGDS(22)
      REAL SNODEP(512,512), ICE(512,512)  
      LOGICAL BITMAP(512,512), LUSAF
C
C     PROCESS USAF N.H. SNOW/ICE ANAL (FROM AFGWC VIA SPN), GRIBBED
C     USAF SNOW/ICE ANAL IS DAILY AND HIGH RES (45 KM)             
C                                                                  
C  COMBINE THE TWO SEPARATE ICE/SNOW FIELDS INTO ONE FIELD.  
C  SNOW IS IN METERS.  ICE FLAG OF '1' IS CONVERTED TO A VALUE
C  OF '11' IN SNODEP FIELD.
C                                                                  
      JPDS = -1
	 CALL BAOPEN(INSNOAF,'snowdepth.grb',IRETBA)
Cmp      CALL GETGBOPL(INSNOAF,0,512*512,0,JPDS,JGDS,KF,KNUM,KPDS,KGDS,
      CALL GETGB(INSNOAF,0,512*512,0,JPDS,JGDS,KF,KNUM,KPDS,KGDS,
     &     BITMAP,ICE,IRET1)
      WRITE(6,*) 'AFTER GETGB FOR AF ICE, IRET=', IRET1
C
      JPDS = -1
Cmp      CALL GETGBOPL(INSNOAF,0,512*512,1,JPDS,JGDS,KF,KNUM,KPDS,KGDS,
      CALL GETGB(INSNOAF,0,512*512,1,JPDS,JGDS,KF,KNUM,KPDS,KGDS,
     &     BITMAP,SNODEP,IRET2)
      WRITE(6,*) 'AFTER GETGB FOR AF SNOW, IRET=', IRET2
C
      IYR2D = KPDS(8)
      IMN = KPDS(9)
      IDY = KPDS(10)
      ICENT = KPDS(21)
      IF(IYR2D.LT.100) THEN
       IYR = (ICENT - 1) * 100 + IYR2D
      ELSE
       IYR = ICENT * 100
      ENDIF
C
      LUSAF = IRET1.EQ.0 .AND. IRET2.EQ.0
      WRITE(6,*) 'LUSAF=', LUSAF
      IF (.NOT.LUSAF) RETURN
C
      DO 20 I = 1, 512
        DO 10 J = 1, 512
          IF (ABS(ICE(I,J)-1.) .LT. 0.0001) SNODEP(I,J) = 11.
 10     CONTINUE
 20   CONTINUE
      RETURN
      END
