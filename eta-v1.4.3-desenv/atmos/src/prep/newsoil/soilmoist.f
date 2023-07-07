      PROGRAM KMVEGT
      include 'parmetasoil'
c      PARAMETER  (IM=149,JM=311,IMJM=IM*JM-JM/2)
c      PARAMETER  (IMA=5000,JMA=5000,IMJMA=IMA*JMA)
      PARAMETER  (IMA=5025,JMA=5400,IMJMA=IMA*JMA)
      PARAMETER  (NTYP = 18)
      PARAMETER
     1 (WBD=-(IM-1)*DLMD, SBD=-((JM-1)*0.5)*DPHD)
      PARAMETER
     1 (IM1=IM-1,NINC=2*IM-1,DTR=1.745329E-2)
      PARAMETER
     1 (WB=WBD*DTR,SB=SBD*DTR,DLM=DLMD*DTR,DPH=DPHD*DTR)
c      DIMENSION NSIBV(NTYP),NEGRD(IMJM,NTYP)
c      DIMENSION NSIBG(6,6) ,NEGRDG(IMJM,   6)
      REAL      ALAT,ALON
      DIMENSION HGT(IM,JM),SM(IM,JM),SMK(IMJM)
      real vegtyp(IM,JM)
      INTEGER   NEGRD(IMJM),INDICE(IM,JM)
      INTEGER IK,IN,INIJ,FINIJ,INII,FINII,IJ,II
      REAL NTY,SUMNTY(IMJM),nmax(imjm)

      DATA MISS/0/,NZERO/0/
C
C *** ALAT AND ALON ARE THE GEODETIC LATITUE 
C *** AND LONGITUDE OF THE 1KM LAT/LON GRID. 
C *** FOR EACH PAIR OF ALAT/ALON, FIND A K-INDEX ON ETA
C *** E GRID H POINT. THEN AT EVERY H POINT, 
C *** WE HAVE THE NUMBERS OF OB FOR THE 
C *** 20 VEG TYPES AND THE NUMBERS OF OB FOR THE 6 VEG 
C *** TYPE GROUPS. AT LAST, DETERMINE THE DOMINANT VEG 
C *** TYPE IN THE DOMAINT VEG GROUP ON EVERY H GRID POINT
C
       DO K = 1,IMJM
          NEGRD(K) = NZERO
       ENDDO
C***
C     READ SEA-LAND MASK FOR 40KM ETA-GRID (1d)
      read(10) hgt,sm
C       open(76,file='mask_dan1.asc',status='unknown')

      IMT=2*IM-1
      DO K=1,IMJM
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX.LE.IM)THEN
          I=IX
          J=2*JX-1
        ELSE
          I=IX-IM
          J=2*JX
        ENDIF
        SMK(K)=SM(I,J)
C	write(76,*)k,smk(k),i,j
      ENDDO
cJLG        write(67)SMK
C	close(76)

C
C     READ MASK in 2-D WITH 1KM RESOLUTION 
C     READ SOIL TYPE MASK in 1-D WITH 1KM RESOLUTION 
C     (INTEGER VARIABLE)
C
	NNN=1
	NNOB=1
c
       dxy=0.01   !resolucao do mapa de umidade. (1km)
c
c       open(76,file='newmoist_dan1.asc',status='unknown')
c            write(76,*)'inicio',WBD,IM,DLMD,SBD,JM,DPHD
c       open(77,file='newmoist_dan2.asc',status='unknown')
        n=0
        do 210 nj=1,jma
            alat = -40.0 + (nj-1)*dxy
          do 210 ni=1,ima
              alon = -82.0 + (ni-1)*dxy
              n=n+1
C
            read(33,*)NTY
C            write(77,*)NTY
 6122       format(f10.4)
c
	IF (NTY .LT. 0.0 .OR. NTY .GT. 1.0 )then 
           NNN=NNN+1
C            write(77,*)'saindo1',NTY
	   GO TO 210
	ENDIF
C
C    BE WARE THE WEST IS POSITIVE.
C
C    The comment above is assuming the given longitudes
C     were positive. But since we have given the western
C     boundary as negative, we comment the line below 
C     in which the longitude was set to negative
C
         FLAT = ALAT
         FLON = -ALON
CJLG	 write(64,*) FLAT, FLON, NTY
C
C        FIND KINDEX OF OBSERVATION
C
        ALN =   FLON
        APH =   FLAT
C
        CALL TOM( APH, ALN, KINDEX,WB,SB,DLM,DPH,TLM,TPH)
C         write(76,*)K,NTY
        IF(KINDEX .LE. 0 .OR. KINDEX .GT. IMJM) THEN
C            WRITE(6,6002) KINDEX, APH, ALN, ALAT
            NNOB=NNOB+1
            GOTO 210
	ELSE
c            write(76,*)'pasou1',K,NTY,dtr
	    tlmd=tlm/dtr
	    tphd=tph/dtr
c            write(76,*)'pasou1',K,NTY,tlm,tlm/dtr,wbd
	    if (abs(tlmd).gt.abs(wbd)) go to 210 ! Chou: out of domain    
c            write(76,*)'pasou2',K,NTY,tlm,dtr,tlm/dtr,wbd
        ENDIF
6002    FORMAT(1x,'REJECTED:    KINDEX, APH, ALN, ALAT',I5,1x,3F11.6) 
        K= KINDEX
        IF(SMK(K) .EQ. 1.0) GO TO 210
        NEGRD(K)=NEGRD(K) + 1
        SUMNTY(K)=SUMNTY(K)+ NTY
            write(76,*)'pasou 2',K,NTY,SMK(K),NEGRD(K),SUMNTY(K)
210	CONTINUE
c        close(76)
c
c  fim do primeiro loop   
c
       open(75,file='newsoilmoist_prep_dan.asc',status='unknown')
       write(75,*)'imjm=',imjm
       DO 220 K=1,IMJM
        IF (SMK(K) .EQ. 1.0) GO TO 220
	IF (NEGRD(K) .LE. 0) THEN
		WRITE(6,*)'PROBLEMA',K,NEGRD(K)
		NMAX(K)=-99.00
c		STOP
C		GO TO 220
	ELSE
		NMAX(K)= SUMNTY(K)/NEGRD(K)
	ENDIF
		

        write(75,*)k, Nmax(k)
 220   CONTINUE

c      WRITE(32) NMAX


c
c Convert veg_map from k-index into 2-dim index im,jm
c
       IMM1=IM-1
         K=0
         DO J=1,JM
         DO I=1,IMM1+MOD(J,2)
            K=K+1
            VEGTYP(I,J)=NMAX(K)
	    INDICE(I,J)=K
        write(75,*)I,J,K,VEGTYP(I,J),NMAX(K)
         ENDDO
         ENDDO
         DO J=2,JM-1,2
            VEGTYP(IM,J)=VEGTYP(IMM1,J)
         ENDDO
        write(75,*)
c        write(75,*)'etatopo j=1,j=jm'
c        j=1
c	write(75,429) (SMK(indice(I,J)),I=1,IM,1)
c	write(75,429) (VEGTYP(I,J),I=1,IM,1)
c        j=2
c	write(75,429) (SMK(indice(I,J)),I=1,IM,1)
c	write(75,429) (VEGTYP(I,J),I=1,IM,1)
c        j=jm-1
c        write(75,429) (SMK(indice(I,J)),I=1,IM,1)
c	write(75,429) (VEGTYP(I,J),I=1,IM,1)
c        j=jm
c        write(75,429) (SMK(indice(I,J)),I=1,IM,1)
c	write(75,429) (VEGTYP(I,J),I=1,IM,1)
c        write(75,*)	
c        write(75,*)'umidade 2d'

cJLG       write(65,*) 'sample veg type values '
        DO J=JM,1,-1
        write(75,429) (VEGTYP(I,J),I=1,IM,1)
        ENDDO
429   format(300(f10.4,x))
      close(75)
c
c CORRIGE DIFERENCAS ENTRE A MASCARA DE SOLOS E A UMIDADE
c

        IMM1=IM-1
c        K=0
	DO J=1,JM
	DO I=1,IMM1+MOD(J,2)
c        	K=K+1
		IF(VEGTYP(I,J).LE.0.001)THEN
		IF(SMK(INDICE(I,J)) .NE. 1.0) THEN			
		WRITE(6,*)'PROBLEMA',i,j,INDICE(I,J),VEGTYP(I,J),
     &		                   NMAX(INDICE(I,J)),SMK(INDICE(I,J))
			IN = -1
			IK = 0
			DOWHILE(IN.LT.0)
				IK=IK+1
				INIJ= MAX(J-IK,1)
				FINIJ= MIN(J+IK,JM)
				INII=MAX(I-IK,1)
				FINII=MIN(I+IK,IMM1+MOD(J,2))
				
				IJ=INIJ
				DOWHILE(IJ.LE.FINIJ.AND.IN.LT.0)
				II = INII
				DOWHILE(II.LE.FINII.AND.IN.LT.0)
				IF(SMK(INDICE(II,IJ)) .NE. 1.0) THEN			
				IF(VEGTYP(II,IJ).GT.0.0)THEN
					VEGTYP(I,J)= VEGTYP(II,IJ)
					NMAX(INDICE(I,J))= VEGTYP(II,IJ)
					IN=1
		WRITE(6,*)'PROBLEMA SOL',INDICE(I,J),IK,IN,VEGTYP(II,IJ)
				ENDIF
				ENDIF
c		WRITE(6,*)'PROBLEMA NOSOL',INDICE(I,J),IK,IN,VEGTYP(II,IJ)
				II = II+1
				ENDDO
				IJ=IJ+1
				ENDDO
			ENDDO
		ELSE
			VEGTYP(I,J)=0.0
			NMAX(INDICE(I,J))=0.0
		
		ENDIF
		ENDIF
	ENDDO
	ENDDO
	DO J=1,JM
	DO I=1,IMM1+MOD(J,2)
		IF(SMK(INDICE(I,J)).EQ.1.0)THEN
			VEGTYP(I,J)=1.0
			NMAX(INDICE(I,J))=1.0
		ENDIF
	ENDDO
	ENDDO

        DO J=2,JM-1,2
           VEGTYP(IM,J)=VEGTYP(IMM1,J)
        ENDDO
	
	DO J=1,2
	DO I=1,IM
		VEGTYP(I,J)=0.0
		NMAX(INDICE(I,J))=0.0
	ENDDO
	ENDDO
	DO J=JM-1,JM
	DO I=1,IM
		VEGTYP(I,J)=0.0
		NMAX(INDICE(I,J))=0.0
	ENDDO
	ENDDO
	


       open(75,file='newsoilmoist_corregido_prep_dan.asc',
     * status='unknown')
       DO 720 K=1,IMJM
c        IF (SMK(K) .EQ. 1.0) GO TO 720
        write(75,*)k, Nmax(k),SMK(K)
c
 720   CONTINUE

cJLG       write(65,*) 'sample veg type values '
        DO J=JM,1,-1
        write(75,431) (VEGTYP(I,J),I=1,IM,1)
        ENDDO
431   format(300(f10.4,x))
      close(75)

       open(75,file='newsoilmoist_erros_prep_dan.asc',
     * status='unknown')
       	DO K=1,IMJM
        	IF (SMK(K) .NE. 1.0.and.nmax(k).le.0.0)then
        	write(75,*)'terra',k, Nmax(k),SMK(K)
		endif
c        	IF(SMK(K) .EQ. 1.0.and.nmax(k).gt.0.0)then
c        	write(75,*)'mar',k, Nmax(k),SMK(K)
c		endif
	ENDDO
	close(75)

      WRITE(32) NMAX

      WRITE(31) vegtyp

 1062 FORMAT(1x,I5,6(1x,I4), 6x,I4)
 1063 FORMAT(1x,I5,20I4)
C 
 9999 WRITE(6,*)  'Data length; ',N-1
      WRITE(6,*)  'Missing data number: ', NNN-1
      WRITE(6,*)  'Rejected data number: ',NNOB-1
C
c      close(75)
      END
C
C     THE ROUTINES BELOW ARE FROM THE TLL.F FILES
C
      SUBROUTINE ETRAN(ELAT,ELON,RLAT,RLON,PHI0,AM0)
      PI=3.14159265
C     ASSUMES DEGREES SENT AND RECEIVED
      DEGRAD=PI/180.
      A=SIN(ELAT*DEGRAD)*COS(PHI0)+COS(ELAT*DEGRAD)*SIN(PHI0)*
     &COS(ELON*DEGRAD)
      RLAT=ASIN(A)
      B=COS(ELAT*DEGRAD)*COS(ELON*DEGRAD)/(COS(RLAT)*COS(PHI0))
      F=B-TAN(RLAT)*TAN(PHI0)
      IF (F.GE.1.0) F=1.0
      C=ACOS(F)
      IF (ELON.LT.0.0) THEN
      RLON=(AM0+C)/DEGRAD
      ELSE
      RLON=(AM0-C)/DEGRAD
      END IF
      RLAT=RLAT/DEGRAD
      RETURN
      END
C
      SUBROUTINE TOM(APHD,ALMD,KMIN,WB,SB,DLM,DPH,TLM,TPH)
C     ******************************************************************
C     *                                                                *
C     *     ROUTINE FOR CALCULATING THE K INDEX OF THE NEAREST         *
C     *       HEIGHT POINT ON THE E GRID TO THAT OF A POINT            *
C     *                GIVEN IN GEODETIC COORDINATES                   *
C     *  LAT POS NORTH    LON POS WEST                                 *
C     ******************************************************************
C                             P A R A M E T E R
C     1 (IM=149,DLMD=.21634615387500000000
C     2  ,DPHD=.20192307712500000000)
CJLG     3  ,WBD=-5.913462,SBD=-4.173077)
      include 'parmetasoil'
      INTEGER IM1,NINC
                             P A R A M E T E R
     1 (DTR=1.745329E-2,NMAX=4)
                             D I M E N S I O N
     1  AH(NMAX),KH(NMAX)
C----------------------------------------------------------------------
      IM1=IM-1
      NINC=2*IM-1
      AIM1=REAL(IM1)
      RDLM=1./DLM
      RDPH=1./DPH
C
C***  CONVERT LOCATION FROM GEODETIC TO TRANSFORMED LAT/LON
C
      CALL TLLRAD(APHD,ALMD,TLM,TPH)
C
C***  X1,Y1 ARE THE COORDINATES WITH RESPECT TO THE X,Y TRANSFORMED
C***  SYSTEM
C
      X1=(TLM-WB)*RDLM
      Y1=(TPH-SB)*RDPH
C
C***  X2,Y2 ARE THE COORDINATES WITH RESPECT TO THE X',Y' TRANSFORMED
C***  SYSTEM (TRANSLATED BY AIM1 ALONG Y')'
C
      X2=.5*( X1+Y1)
      Y2=.5*(-X1+Y1)+AIM1
C
C***  I2,J2 ARE THE COORDINATES OF THE SOUTHERN HEIGHT POINT IN
C***  A CLUSTER OF FOUR
C
      I2=INT(X2)
      J2=INT(Y2)
C
      P=X2-I2
      Q=Y2-J2
      PQ=P*Q
C
      JR=J2-IM1
      I3=I2-JR
      J3=I2+JR
C***
C***  INDEX K CORRESPONDING TO POINT (I2,J2) (SOUTHERNMOST OF THE
C***  FOUR SURROUNDING HEIGHT POINTS)
C***
      K=J3*IM-J3/2+(I3+2)/2
C***
C***  FIND WHICH OF THE FOUR H POINTS IS NEAREST
C***
      AH(1)=P*Q
      AH(2)=(1.-P)*Q
      AH(3)=P*(1.-Q)
      AH(4)=(1.-P)*(1.-Q)
      KH(1)=K
      KH(2)=K+IM
      KH(3)=K+IM1
      KH(4)=K+NINC
C
      KMIN=KH(1)
      NM=1
  100 NL=NMAX-NM
      DO 200 NN=1,NL
      NX=NM+NN
      IF(AH(NM).GT.AH(NX))THEN
        KMIN=KH(NX)
        NM=NX
        IF(NX.LT.NMAX)GO TO 100
      ENDIF
  200 CONTINUE
C----------------------------------------------------------------------
      RETURN
      END
C
C----------------------------------------------------------------------
      SUBROUTINE TLLRAD(APHD,ALMD,TLM,TPH)
C     ******************************************************************
C     *                                                                *
C     *  ROUTINE FOR CONVERTING FROM GEODETIC COORDINATES (RADIANS)    *
C     *            TO TRANSFORMED COORDINATES (RADIANS)                *
C     *                                                                *
C     ******************************************************************
       include 'parmetasoil'
                             P A R A M E T E R
C
     1 (DTR=1.745329E-2)

C Alterado o sinal
      TLM0=-TLM0D*DTR ! Chou: tlmo is positive if long is west
      TPH0=TPH0D*DTR
C---------------------------------------------------------------------

      STPH0=SIN(TPH0)
      CTPH0=COS(TPH0)
      ALM=ALMD*DTR
      APH=APHD*DTR
      RELM=ALM-TLM0
      SRLM=SIN(RELM)
      CRLM=COS(RELM)
C
      SPH=SIN(APH)
      CPH=COS(APH)
      CC=CPH*CRLM
C
      ANUM=CPH*SRLM
      DENOM=CTPH0*CC+STPH0*SPH
C----------------------------------------------------------------------
CJLG        TLM=ATAN2(ANUM,DENOM)
      TLM=-ATAN2(ANUM,DENOM)
      TPH=ASIN(CTPH0*SPH-STPH0*CC)
C----------------------------------------------------------------------
      RETURN
      END

