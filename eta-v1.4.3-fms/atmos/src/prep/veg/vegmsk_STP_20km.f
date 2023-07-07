      PROGRAM KMVEGT
      include 'parmeta'
      CHARACTER arqtopo*16 
      PARAMETER  (IMA=9600,JMA=9600,IMJMA=IMA*JMA)
      PARAMETER  (NTYP=12)
      PARAMETER
     1 (WBD=-(IM-1)*DLMD, SBD=-((JM-1)*0.5)*DPHD)
      PARAMETER
     1 (IM1=IM-1,NINC=2*IM-1,DTR=.01745329251994329576)
      PARAMETER
     1 (WB=WBD*DTR,SB=SBD*DTR,DLM=DLMD*DTR,DPH=DPHD*DTR)
      DIMENSION NSIBV(NTYP),NEGRD(IMJM,NTYP)
      DIMENSION NSIBG(6,6) ,NEGRDG(IMJM,6)
      REAL      ALAT,ALON
      DIMENSION HGT(IM,JM),SM(IM,JM),SMK(IMJM)
      INTEGER*4 vegtyp(IM,JM),NMAX(IMJM)
      INTEGER   NTY
C
      DATA NSIBV/1,2,3,4,5,6,7,8,9,10,11,12/
      DATA (NSIBG(1,I),I=1,6)/1,2,3,4,5,6/
      DATA (NSIBG(2,I),I=1,6)/7,12,0,0,0,0/
      DATA (NSIBG(3,I),I=1,6)/8,9,10,0,0,0/
      DATA (NSIBG(4,I),I=1,6)/11,0,0,0,0,0/
      DATA (NSIBG(5,I),I=1,6)/0,0,0,0,0,0/
      DATA (NSIBG(6,I),I=1,6)/0,0,0,0,0,0/
C
c      DATA MISS/-99/,NZERO/0/
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
C       DO K=1,IMJM
C          NMAX(K)=0
C        DO N=1,NTYP
C          NEGRD(K,N)=0
C          NEGRDG(K,N)=0
C        ENDDO
C       ENDDO 
          NMAX=0    
	  NEGRD=0
	  NEGRDG=0

C***
C     READ SEA-LAND MASK FOR ETA-GRID (1d)
      open(10,form='unformatted',access='sequential')
      read(10) hgt,sm

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
      ENDDO
cJLG        write(67)SMK


C
C     READ MASK in 2-D WITH 1KM RESOLUTION 
C     READ VEGETATION TYPE MASK in 1-D WITH 1KM RESOLUTION 
C     (INTEGER VARIABLE)
C
       NNN=1
       NNOB=1
c
       dxy=0.004166666666667   !resolucao do mapa de veg. (0.5km)
c
        n=0
        do 210 nj=1,jma
            alat = -19.9979166666667 + (nj-1)*dxy
          do 210 ni=1,ima
              alon = -19.9979166666667 + (ni-1)*dxy
              n=n+1
C
Chou change unit
            open(20,form='formatted')
	    read(20,6122)NTY
 6122       format(i2)
c
         IF (NTY .LE. 0 .OR. NTY .GT. NTYP )then 
           NNN=NNN+1
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
        IF(KINDEX .LE. 0 .OR. KINDEX .GT. IMJM) THEN
C            WRITE(6,6002) KINDEX, APH, ALN, ALAT
            NNOB=NNOB+1
            GOTO 210
         ELSE
          tlmd=tlm/dtr
          tphd=tph/dtr
          if (abs(tlmd).gt.abs(wbd)) go to 210 ! Chou: out of domain    
        ENDIF
 6002  FORMAT(1x,'REJECTED:    KINDEX, APH, ALN, ALAT',I5,1x,3F11.6) 
       K= KINDEX
       IF(SMK(K) .GT. 0.95) GO TO 210
C
C     COUNT THE NUMBERS OF EVERY VEG TYPE IN ONE GRID BOX
C
        DO 500 NV = 1, NTYP
          IF(NTY .eq. NSIBV(NV) ) THEN
            NEGRD(K,NV)=NEGRD(K,NV) + 1
        ENDIF
 500    CONTINUE
C
C    GROUP THE 20 VEG TYPES INTO 6 AND COUNT THE NUMBER OF EACH GROUP
C
      DO 620 NG  = 1, 6
      DO 600 NGG = 1, 6
         IF (NTY .EQ. NSIBG(NG,NGG)) THEN
            NEGRDG(K,NG)= NEGRDG(K,NG)+1
         ENDIF
 600   CONTINUE
 620   CONTINUE
 210   CONTINUE
c
c  fim do primeiro loop   
c
       DO 220 K=1,IMJM
        IF(SMK(K) .ge. 0.95) GO TO 220
          NTYY=0
C     DETERMINE THE DOMINANT GROUP IN THE GRID BOX(exclude water)
C
      MMAX = -99999
      DO 650 NG=1,6     
        IF (NEGRDG(K,NG) .LE.MMAX .OR. NEGRDG(K,NG) .EQ. NZERO) THEN
            GO TO 650
        ENDIF
      MMAX=NEGRDG(K,NG)
        NTYY=NG
 650  CONTINUE
	
Carol	
c     	IF (NEGRDG(K,5).gt.0) ntyy=5
						
      IF(NTYY .EQ. 0) THEN
C	  write(47,*) 'No available data at the grid point ',K
          GO TO 220
       ENDIF
C
C      NOW DETERMINE THE DOMINANT VEG TYPE IN THE DOMINANT GROUP
C
      MMAX = -99999
      NGG=1
       DO 302 NV = 1, NTYP
        IF (NSIBV(NV) .NE. NSIBG(NTYY,NGG))  GO TO 302
        IF (NGG .GT. 6) THEN
        GO TO 220
        ELSEIF (NEGRD(K,NV) .LT. MMAX .OR. NEGRD(K,NV) .EQ. NZERO) THEN
            NGG=NGG+1
         GO TO 302
        ELSE
         NGGG=NGG
            NVV=NV
         NMAX(K)=NSIBG(NTYY,NGGG)
         NGG=NGG+1
         MMAX=NEGRD(K,NV)
         GO TO 302
        ENDIF
  302  CONTINUE
c	write(48,*) 'K=',K,' NTYY= ',NTYY,' NVV= ',NVV,' NGG= ',NGGG
 220   CONTINUE
C
      WRITE(32) real(NMAX)
CJLG      WRITE(61,1061) (NMAX(K),K=1,IMJM)
CJLG 1061 FORMAT(1x,149(1x,I4))

           write(65,*) 'sample veg type values '
         DO J=1,imjm
         if (nmax(j).eq.0. and. smk(j).lt.0)  write(65,*) j,nmax(k) 
         ENDDO


c
c Convert veg_map from k-index into 2-dim index im,jm
c
       IMM1=IM-1
         K=0
         DO J=1,JM
         DO I=1,IMM1+MOD(J,2)
            K=K+1
            VEGTYP(I,J)=NMAX(K)
         ENDDO
         ENDDO
         DO J=2,JM-1,2
            VEGTYP(IM,J)=VEGTYP(IMM1,J)
         ENDDO

cJLG       write(65,*) 'sample veg type values '
cJLG        DO J=JM,1,-1
cJLG        write(65,429) (VEGTYP(I,J),I=1,IM,1)
cJLG        ENDDO
cJLG  429   format(149(I2,x))


      WRITE(31) vegtyp
C
      DO 720 NN = 1,IMJM
        if (NMAX(NN) .eq. MISS) GO TO 720
c        WRITE(63,1062) NN, (NEGRDG(NN,N1),N1=1,6), NMAX(NN)
c        WRITE(62,1063) NN, (NEGRD(NN,MM),MM=1,20)
 720  CONTINUE
 1062 FORMAT(1x,I5,6(1x,I4), 6x,I4)
 1063 FORMAT(1x,I5,20I4)
C 
 9999 WRITE(6,*)  'Data length; ',N-1
      WRITE(6,*)  'Missing data number: ', NNN-1
      WRITE(6,*)  'Rejected data number: ',NNOB-1
C
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
      include 'parmeta'
      INTEGER IM1,NINC
                             P A R A M E T E R
     1 (DTR=.01745329251994329576,NMAX=4)
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
       include 'parmeta'
                             P A R A M E T E R
C
     1 (DTR=.01745329251994329576)

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
