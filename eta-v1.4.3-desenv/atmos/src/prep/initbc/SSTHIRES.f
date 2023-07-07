      SUBROUTINE SSTHIRES (SST,SM,GLAT,GLON,IDAT,LIST,DTR)
!
      IMPLICIT REAL (A-H, O-Z)
!
      INCLUDE "parmeta"
!
!hires
	REAL INCR_LAT,INCR_LON,ILAT1,ILON1,ILAT2,ILON2
!	PARAMETER(INCR=0.25) !mel
        PARAMETER(INCR_LAT=0.0833333)
        PARAMETER(INCR_LON=0.0833333)
!hires

      PARAMETER  (H90=90.0,H360=360.0,D5=5.E-1,D00=0.0,H1=1.0)
!
      INTEGER IDATE(4),IDAT(3),MONTH(12)
      DIMENSION SSTLL(4321,2160),SALTLK(12),SALTLA(2),SALTLO(2) !mel
!
      DIMENSION  SST(IM,JM), SM(IM,JM), GLAT(IM,JM), GLON(IM,JM)
      DIMENSION  SST2(IM,JM)
!
      DATA   INSST/39/
      DATA   INDXST/0/

      DATA MONTH/31,28,31,30,31,30,31,31,30,31,30,31/
!
      DATA SALTLK/273.38,274.27,278.50,283.01,287.33,293.41
     1,           297.13,297.73,294.97,289.58,282.31,275.67/
!
!     CORNERS OF SALT LAKE LATITUDE/LONGITUDE BOX
!     in degrees---> 40.0     42.0            111.0    114.0
      DATA SALTLA/0.698132,0.733038/,SALTLO/1.937315,1.989675/
!
      IOUTUPRT = LIST

!JLG 20081104 - Removed reading of sst grib. Now read a binary sst converted from sst grib2 with wgrib2
 
      CALL GRIDST(INSST,SSTLL)


!----  INTERPOLATE 0.25-DEG GLOBAL SST TO ETA GRID  -------
!
!-CP NOTE:  THIS SUBROUTINE AND INTERPOLATION ALGORITHM ASSUME
!-CP A 0.25-DEG GLOBAL SST FIELD IN THE FOLLOWING FORMAT:  
!
!	I=1 at 0.125 E, I=1440 at 359.875E, I=1441 at 0.125 E
!	J=1 at 89.875 S, J=720 at 89.875 N 
!
!	Old 1 degree data
!	I=1 AT 0.5 E,  I=2 AT 1.5 E, ... , I=360 at 0.5W
!	J=1 AT 89.5S, J=2 AT 88.5 S, ..., J=180 at 89.5N
!
!	Old 0.5 degree data
!	I=1 at 0.25 E, I=720 at 359.75E, I=721 at 0.25 E
!	J=1 at 89.75 S, J=360 at 89.75 N 
!
!-CP  
!-CP In the interpolation algorithm below, glon is positive westward,
!-CP from 0 to 360, with 0 at the greenwich meridian.  Elon is positive 
!-CP eastward, thus the need to subtract glon from 360 to get the index
!-CP of the correct oisst point.  If your input 1 deg SST field is in
!-CP a different indexing scheme, you will need to change the algorithm
!-CP below - see "grdeta.oldoi"

      DO J=1,JM
      DO I=1,IM
	
      ELAT=H90+GLAT(I,J)/DTR
      ELON=H360-GLON(I,J)/DTR

      IF(ELON.GT.H360)ELON=ELON-H360

	DIF=ELON-INT(ELON)
CCCC	write(6,*) 'DIF_LON: ',DIF
CCCC    0.0833333358/2 = 0.0416666679
CCCC    0.0833333358*1.5=0.12500000370

        IF (DIF .lt. INCR_LON) ILON1=INT(ELON)+(0.5*INCR_LON)
	IF (DIF .ge. INCR_LON .and. DIF .lt. (2*INCR_LON)) THEN 
	ILON1=INT(ELON)+(1.5*INCR_LON)
	ENDIF 
	IF (DIF .ge. (2*INCR_LON) .and. DIF .lt. (3*INCR_LON)) THEN 
	ILON1=INT(ELON)+(2.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (3*INCR_LON) .and. DIF .lt. (4*INCR_LON)) THEN 
	ILON1=INT(ELON)+(3.5*INCR_LON)
	ENDIF 
	IF (DIF .ge. (4*INCR_LON) .and. DIF .lt. (5*INCR_LON)) THEN 
	ILON1=INT(ELON)+(4.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (5*INCR_LON) .and. DIF .lt. (6*INCR_LON)) THEN 
	ILON1=INT(ELON)+(5.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (6*INCR_LON) .and. DIF .lt. (7*INCR_LON)) THEN 
	ILON1=INT(ELON)+(6.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (7*INCR_LON) .and. DIF .lt. (8*INCR_LON)) THEN 
	ILON1=INT(ELON)+(7.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (8*INCR_LON) .and. DIF .lt. (9*INCR_LON)) THEN 
	ILON1=INT(ELON)+(8.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (9*INCR_LON) .and. DIF .lt. (10*INCR_LON)) THEN 
	ILON1=INT(ELON)+(9.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (10*INCR_LON) .and. DIF .lt. (11*INCR_LON)) THEN 
	ILON1=INT(ELON)+(10.5*INCR_LON)
	ENDIF
	IF (DIF .ge. (11*INCR_LON)) THEN 
	ILON1=INT(ELON)+(11.5*INCR_LON)
	ENDIF

!old      IF(ILON1.EQ.D00)ILON1=360.
      IF(ILON1.LE.D00)ILON1=360.
      ILON2=ILON1+INCR_LON


!-MP	New approach sets ILAT1, ILON1 to point on SST grid that is
!-MP	SW of the Eta Grid point.

        DIF=ELAT-INT(ELAT)
CCCC	write(6,*) 'DIF_LAT: ',DIF
CCCC    0.0833719298/2=0.0416859649
CCCC    0.0833719298*1.5=0.12505789470
	IF (DIF .lt. INCR_LAT) ILAT1=INT(ELAT)+(0.5*INCR_LAT)
	IF (DIF .ge. INCR_LAT .and. DIF .lt. (2*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(1.5*INCR_LAT)
	ENDIF 
	IF (DIF .ge. (2*INCR_LAT) .and. DIF .lt. (3*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(2.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (3*INCR_LAT) .and. DIF .lt. (4*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(3.5*INCR_LAT)
	ENDIF 
	IF (DIF .ge. (4*INCR_LAT) .and. DIF .lt. (5*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(4.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (5*INCR_LAT) .and. DIF .lt. (6*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(5.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (6*INCR_LAT) .and. DIF .lt. (7*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(6.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (7*INCR_LAT) .and. DIF .lt. (8*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(7.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (8*INCR_LAT) .and. DIF .lt. (9*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(8.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (9*INCR_LAT) .and. DIF .lt. (10*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(9.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (10*INCR_LAT) .and. DIF .lt. (11*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(10.5*INCR_LAT)
	ENDIF
	IF (DIF .ge. (11*INCR_LAT)) THEN 
	ILAT1=INT(ELAT)+(11.5*INCR_LAT)
	ENDIF

      ILAT2=ILAT1+INCR_LAT

!hires,notsure      W1=ELON-ILON1+D5
      W1=ELON-ILON1+INCR_LON/2.
      IF(W1.LT.D00)W1=W1+H360
!hires,notsure      W2=ELAT-ILAT1+D5
      W2=ELAT-ILAT1+INCR_LAT/2.
      AR1=W1*W2
      AR2=W1*(H1-W2)
      AR3=(H1-W1)*(H1-W2)
      AR4=(H1-W1)*W2
	LON1INDX=12*(ILON1) 
	LON2INDX=12*(ILON2)
	LAT1INDX=12*(ILAT1) 
	LAT2INDX=12*(ILAT2) 

	if (mod (I,20) .eq. 0 .and. mod(J,20) .eq. 0) then
	write(6,*) 'weights: ',AR1,AR2,AR3,AR4
	write(6,*) 'ILAT1,ILON1,ELAT,ELON: ', ILAT1,ILON1,ELAT,ELON
	write(6,*) '------------------------------------------'
	write(6,*) 'corresponding indices: ', LAT1INDX,LON1INDX,
     +					      LAT2INDX,LON2INDX
	endif

	if (LON1INDX .lt. 1 .or. LON1INDX .gt. 4321) then 
	write(6,*) 'out of bounds on index!!', LON1INDX
	endif
	SST2(I,J)=AR1*SSTLL(LON2INDX,LAT2INDX)+
     +	 	 AR2*SSTLL(LON2INDX,LAT1INDX)+
     +	 	 AR3*SSTLL(LON1INDX,LAT1INDX)+
     +	 	 AR4*SSTLL(LON1INDX,LAT2INDX)
	SST(I,J)=SSTLL(LON2INDX,LAT2INDX)
      ENDDO
      ENDDO

CGSM			 write(15) SST			
CGSM			 write(18) SST2

!***
!***  INSERT TEMPERATURES FOR THE GREAT SALT LAKE
!***
      ID1=IDAT(1)
      ID2=IDAT(2)
      MARG0=ID1-1
      IF(MARG0.LT.1)MARG0=12
      MNTH0=MONTH(MARG0)
      MNTH1=MONTH(ID1)
      IF(ID2.LT.15)THEN
        NUMER=ID2+MNTH0-15
        DENOM=MNTH0
        IARG1=MARG0
        IARG2=ID1
      ELSE
        NUMER=ID2-15
        DENOM=MNTH1
        IARG1=ID1
        IARG2=ID1+1
        IF(IARG2.GT.12)IARG2=1
      ENDIF
      FRAC=NUMER/DENOM
      DO J=1,JM
      DO I=1,IM
        IF(GLAT(I,J).GT.SALTLA(1).AND.GLAT(I,J).LT.SALTLA(2))THEN
          IF(GLON(I,J).GT.SALTLO(1).AND.GLON(I,J).LT.SALTLO(2))THEN
            IF(SM(I,J).GT.0.5)
     1        SST(I,J)=SALTLK(IARG1)+
     2                (SALTLK(IARG2)-SALTLK(IARG1))*FRAC
          ENDIF
        ENDIF
      ENDDO
      ENDDO
C
      RETURN
C
 4500 CONTINUE
C              ERROR OCCURRED WHEN INPUTING SST FROM GRIB.
      WRITE (IOUTUPRT, 4550) INSST
 4550 FORMAT ('0', 'ERROR OCCURRED WHEN READING IN SST        ',
     1             'ON UNIT', I3, ' GRIB ' /
     2        ' ', 'EXECUTION TERMINATING.')
C
      STOP 222
C
      END





