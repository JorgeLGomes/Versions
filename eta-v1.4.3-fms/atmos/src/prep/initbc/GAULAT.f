      SUBROUTINE GAULAT(GAUL,K)
C
      DIMENSION A(500)
C      DIMENSION GAUL(1)
      DIMENSION GAUL(K)

	real(8) XZ,AVSP,SP,PK,PKMRK,PKM1,PKM2
C
C      ESP=1.E-14
      ESP=1.E-12
      C=(1.E0-(2.E0/3.14159265358979E0)**2)*0.25E0
      FK=K
      KK=K/2
      CALL BSSLZ1(A,KK)
      do30: DO IS=1,KK
      XZ=COS(A(IS)/SQRT((FK+0.5E0)**2+C))
      ITER=0
   10 PKM2=1.E0
      PKM1=XZ
      ITER=ITER+1
      IF(ITER.GT.10) GO TO 70
      do20: DO N=2,K
      FN=N
      PK=((2.E0*FN-1.E0)*XZ*PKM1-(FN-1.E0)*PKM2)/FN
      PKM2=PKM1
      PKM1=PK
      END DO do20
      PKM1=PKM2
      PKMRK=(FK*(PKM1-XZ*PK))/(1.E0-XZ**2)
      SP=PK/PKMRK
      XZ=XZ-SP
      AVSP=ABS(SP)
      IF(AVSP.GT.ESP) GO TO 10
      A(IS)=XZ
      END DO do30
      IF(K.EQ.KK*2) GO TO 50
      A(KK+1)=0.E0
      PK=2.E0/FK**2
      do40: DO N=2,K,2
      FN=N
      PK=PK*FN**2/(FN-1.E0)**2
      END DO do40
   50 CONTINUE
      DO 60 N=1,KK
      L=K+1-N
      A(L)=-A(N)
   60 CONTINUE
C
      RADI=180./(4.*ATAN(1.))
      DO 211 N=1,K
      GAUL(N)=90.-ACOS(A(N))*RADI
C	if (mod(N,3) .eq. 0) write(6,*) 'N, GAULAT: ', N, GAUL(N)
  211 CONTINUE
C
C     PRINT *,'GAUSSIAN LAT (DEG) FOR JMAX=',K
C     PRINT *,(GAUL(N),N=1,K)
C
      RETURN
   70 WRITE(6,6000)
 6000 FORMAT(//5X,14HERROR IN GAUAW//)
      STOP
      END
