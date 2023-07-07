      SUBROUTINE SNOHIRES (SI,   SM,GLAT,GLON,LIST,LUSAF,DTR)
C
C  INPUTS:  SM (ETA LAND/SEA MASK)
C            GLAT (LATITUDE ARRAY OF ETA GRID)
C            GLON (LONGITUDE ARRAY OF ETA GRID)
C            LIST (UNIT NUMBER OF PRINTOUT)
C            LUSAF (LOGICAL - TRUE:USE IMS SNOW AND USAF SNOW, 
C                            FALSE: USE IMS SNOW ONLY)
C            DTR   (DEGREES TO RADIANS CONVERSION FACTOR)
C
C  OUTPUTS: SI (SNOWDEPTH AND SEA-ICE ON ETA GRID, SEE CALLING
C                SUBROUTINE CNSTS FOR MORE DETAILS ON SI OUTPUT)
C
      IMPLICIT REAL (A-H, O-Z)
C
      INCLUDE "parmeta"
C
      LOGICAL  LUSAF
C
      REAL         AFSNO   (512,512)
      INTEGER      MSKAF   (512,512)
C
      REAL         SCVH    (1024,1024)
      INTEGER      MSKSCVH (1024,1024)
C
      REAL         SI(IM,JM), SM(IM,JM), GLAT(IM,JM), GLON(IM,JM)
C
C  INPUT UNITS (45-NESDIS DAILY SNOW/ICE, 41-USAF DAILY SNOW/ICE,
C                43-NESDIS 1/16-BEDIENT LAND/SEA MASK,
C                42-USAF   1/8-BEDIENT LAND/SEA MASK)
C 
      DATA   INSNOIMS/45/, INSNOAF/41/, INMSKAF/42/, INMSKIMS/43/
C
C  SET LOWER THRESHOLD FOR SNODEPTH IN TENTHS OF INCHES (MININUM
C  SNOW DEPTH IS 1.5" WHERE THERE IS SNOW COVERAGE).
C
      DATA     IDEPTH /15/
C
C**************************  BEGIN EXECUTION ***********************
C
C     SPECIFY THE UNIT NUMBER OF THE PRINTER.
      IOUTUPRT = LIST
C
C     SPECIFY PARAMETERS OF THE N.H. 1024X1024 IMS 1/16-MESH GRID
C
C     LOCATION OF THE POLE:
      XPNMC8 = 513.
      YPNMC8 = 512.
C     N -- THE NUMBER OF GRID INTERVALS FROM THE POLE TO THE EQUATOR:
      ENNMC8 = 16.0E0 * 31.2043316E0
C     THE LONGITUDINAL ROTATION ANGLE OF THE GRID
      ALNMC8 = 10.E0
C     THE ORIENTATION WEST LONGITUDE OF THE GRID:
      ORIENT8 = 80.E0
C
C  SPECIFY PARAMETERS OF THE N.H. 512 X 512 GRID TYPE 88 (USAF GRID).
C
C  LOCATION OF POLE:
      XPNMCAF = 257.E0
      YPNMCAF = 256.E0
C  GRID MESH LENGTH AT 60N IN KM (FYI, BUT NOT NEEDED)
C.... XMESHL  = 47.625E0
C  NUMBER OF GRID INTERVALS FROM POLE TO EQUATOR
      ENNMCAF = 8.0 * 31.2043316E0
C  THE LONGITUDE ROTATION ANGLE OF THE GRID
      ALNMCAF = 10.E0
C  THE ORIENTATION WEST LONGITUDE OF THE GRID
      ORIENTAF  = 80.E0
C
C  CALL TO READ NESDIS/IMS 1/16-BEDIENT DAILY N.H. SNOW/ICE VIA SNO16GET
C  
      write(6,*) 'calling SNO16GET'
      CALL SNO16GET(SCVH,IYEAR,IMONTH,IDAY,  INSNOIMS,LIST)
      write(6,*) 'back from SNO16GET'
Ctest
Ctest      write(60,7701) IYEAR, IMONTH, IDAY
Ctest  7701 format('IMS DATA : ', 3i2.2)
Ctest
C
C At present time the IMS file has a 2-digit year.  Add century:
C Now IYEAR from SNO16GET is 4-digits so we don't need to do this
c     IF (IYEAR .GT. 90) THEN
c        IYEAR = IYEAR + 1900
c     ELSE
c        IYEAR = IYEAR + 2000
c     ENDIF
C
C Calculate Julian date:
      NDAY = IW3JDN(IYEAR,IMONTH,IDAY)
      WRITE (IOUTUPRT, 1) INSNOIMS, IYEAR, IMONTH, IDAY, NDAY 
    1 FORMAT(1H ,' DAILY NESDIS SNOW READ IN VIA UNIT NO.=',I3/
     2       1H ,' FOR VALID DATE OF Y-M-D =',I4,'-',I2,'-',I2/
     3       1H ,' WHICH IS JULIAN DAY =', I10) 
C
C   READ THE NESDIS-IMS LAND-WATER MASK
C                                          (SEA=0,LAND=1)
      DO 10 J = 1, 1024
         READ(INMSKIMS,'(80I1)') (MSKSCVH(I,J),I=1,1024)
 10   CONTINUE
C
      CALL PRINTIMS (SCVH, MSKSCVH)
C
C  PRINT SPECIFIED USAF SNOW PROCESSING FLAG
      WRITE(IOUTUPRT,2211)   LUSAF
 2211 FORMAT(//1H ,5X,'SUBROUTINE SNOHIRES                    '/
     1     1H ,1X,'  WILL TRY TO PROCESS USAF SNOWCOVER: LUSAF = ',L2)
C
      IF (LUSAF) THEN
C
C  READ AIR FORCE 1/8-BEDIENT DAILY N.H. SNOW/ICE VIA SNO8GET
C
	write(6,*) 'calling SNO8GET'
      CALL SNO8GET(AFSNO,IYEAR,IMONTH,IDAY,   INSNOAF,LIST,LUSAF)
Ctest
Ctest      write(60,7702) IYEAR, IMONTH, IDAY
Ctest  7702 format('USAF DATA: ', 3i2.2)
Ctest
C
C  IF I/O ERROR ENCOUNTERED IN SNO8GET, THEN LUSAF IS RETURNED FALSE
C  
        IF (.NOT. LUSAF) THEN
          WRITE(IOUTUPRT,56)
   56     FORMAT(1H ,'WARNING: FILE ERR IN USAF SNOW ANAL')
        ELSE 
C  The USAF snow header has a 2-digit year.  Add century:
C  Now IYEAR from SNO8GET is 4-digits so we don't need to do this
c         IF (IYEAR .GT. 90) THEN
c           IYEAR = IYEAR + 1900
c         ELSE
c           IYEAR = IYEAR + 2000
c         ENDIF
C
C Calculate Julian date:
          NAFDAY = IW3JDN(IYEAR,IMONTH,IDAY)
C
          WRITE(IOUTUPRT,2) INSNOAF, IYEAR, IMONTH, IDAY, NAFDAY 
    2     FORMAT(1H ,' DAILY USAF SNOW READ IN VIA UNIT NO.=',I3/
     1           1H ,' FOR VALID DATE OF Y-M-D ='I4,'-',I2,'-',I2/
     2           1H ,' WHICH IS JULIAN DAY =', I10) 
C
           IF (NAFDAY .LT. NDAY-4) THEN
             LUSAF = .FALSE.
             WRITE(IOUTUPRT,2644) NDAY,NAFDAY
 2644        FORMAT(/1H ,'******  WARNING   ****** WARNING ******'/
     1            1H ,'DATE OF USAF SNOW ANAL MORE THAN 4 DAYS OLD'/
     1            1H ,' -  WILL FALL-BACK TO NESDIS IMS SNOW COVER'/
     2            1H ,43X,'   NESDIS SNOW ANAL JULDAY=',I7/
     3            1H ,43X,'     USAF SNOW ANAL JULDAY=',I7)
C 
           ENDIF
        ENDIF
      ENDIF
C
C    NOTE:  UPON RETURN FROM CALL SNO8GET ABOVE, USAF SNOW/ICE FIELD
C            HAS FOLLOWING PHYSICAL RANGES:
C         - VALUES OVER SEA POINTS ARE 0.0 OR 11.0 (SEA-ICE FLAG)
C         - VALUES OVER LAND/COAST ARE 0.0 OR POS DEPTH IN METERS
C         - SNOWDEPTH OVER LAND IS ACTUAL, NOT WATER EQUIVALENT YET
C
C   READ THE USAF AFGWC LAND/COAST/SEA MASK
C                           (SEA=1,LAND=2,COASTAL-LAND=4,OFFWORLD=9)
      READ(INMSKAF)  MSKAF
C
C---------- I/O OF PRIMARY INPUT FIELDS IS COMPLETE -------------
C                        INIT RADIANS TO DEGREES
      RTD = 1./DTR
C
C IDEPTH IS CRITERION OF SNOWDEPTH THRESHOLD IN TENTHS OF
C INCHES BELOW WHICH USAF SNOWDEPTH WILL BE ASSUMED ZERO.
C HERE CONVERT IDEPTH TO METERS
C
       DEPTH = FLOAT(IDEPTH) * 2.54E-3   
C
       IF ( LUSAF) THEN
         WRITE(IOUTUPRT,2321) LUSAF,IDEPTH
 2321    FORMAT(1H //' USAF SNODEP ANAL WILL BE USED, LUSAF=',L2/
     1    1H , 35X,'SNODEPTH THRESHOLD (TENTHS OF INCHES) =',I3)
         CALL PRINTAF (AFSNO,MSKAF)
      ELSE
         WRITE(IOUTUPRT,2322) LUSAF
 2322    FORMAT(1H //' USAF SNODEP ANALYSIS WILL BE IGNORED'/
     1          1H , 35X,'LOGICAL FLAG LUSAF=',L2)
      ENDIF
C
C----------INITIALIZE SNOW/ICE ARRAYS TO ZERO ON ETA GRID------------
C
      SI = 0.0
C--------------------------------------------------------------------
C
C ****** NOW BEGIN MAJOR LOOP OVER ALL ETA GRIDS AND POINTS *******
C
         DO J=1,JM
         DO I=1,IM
C
C--------------- DETERMINE LAT/LON OF ETA GRID POINT -------------
C                    (HERE LONG WILL BE EAST LONG)
C
      YYLAT = GLAT(I,J)*RTD
      XLONG = 360. - GLON(I,J)*RTD
C
C  WHERE ETA DOMAIN SOUTH OF 22 N LAT (INCLUDING ANY S.H.),
C  WE KEEP DEFAULT ZERO SNOW/ICE 
C
      IF (YYLAT.LT.22.0E0) GO TO 4300
C----------------------------------------------------------------
C
C  DETERMINE LOCATION OF ETA POINT ON THE 1024 X 1024 NESDIS/IMS GRID
C
         RM= ENNMC8*COS(YYLAT * DTR) / (1.E0 + SIN(YYLAT * DTR))
         RAD = (XLONG - ALNMC8) * DTR
         X = XPNMC8 + RM * COS(RAD)
         Y = YPNMC8 + RM * SIN(RAD)
C
      IS  = INT(X)
      IP1 = IS + 1
      JS  = INT(Y)
      JP1 = JS + 1
C
C--------------------------------------------------------------------
C.......FIRST INTERPOLATE NESDIS SNOW/ICE AS PRIMARY DEFAULT......
C
C     THE VALUE OF SNOW COVER OR SEA ICE ON THE 1024 X 1024
C     NESDIS/IMS GRID IS 1 FOR SNOW/ICE POINTS, 0 FOR SNOW/ICE-FREE.
C     WE UTILIZE AN 1024 X 1024 LAND MASK
C     FROM IMS TO DISTINGUISH SEA ICE FROM SNOW POINTS
C     (ACTUALLY, THE DATA WE GOT FROM IMS HAVE DIFFERENT VALUES FOR
C     SNOW AND FOR ICE.  BUT SINCE WE WANT TO USE AS MUCH OF THE ORIGINAL
C     PROGRAM [WRITTEN FOR THE 1/2-MESH SAB SNOW/ICE] AS POSSIBLE, WE GIVE
C     SNOW COVER AND SEA ICE THE SAME VALUE (=1) AND LET THE MASK DO THE
C     JOB.
C
C  NOW USE ETA AND IMS LAND-SEA MASKS TO ENSURE ONLY LAND
C  POINTS ARE INTERPOLATED TO LAND POINTS (TO DETERMINE SNOW)
C  AND ONLY SEA POINTS ARE INTERPOLATED TO SEA POINTS (FOR ICE)
C  (NESDIS/IMS LAND MASK: SEA=0,LAND=1, WHILE THE ETA MASK IS THE
C  OTHER WAY ROUND).
C
      ILAND = 1
      IF( SM(I,J) .GT. 0.9 ) ILAND=0
C
      IPOINT = NINT(X)
      JPOINT = NINT(Y)
      IF ( MSKSCVH(IPOINT,JPOINT) .EQ. ILAND ) THEN
        SI(I,J) = SCVH(IPOINT,JPOINT)
        GO TO 3351
      ENDIF
C
C  NEAREST NEIGHBOR NOT SAME SFC TYPE, SO USE ALL 4 SURROUNDING POINTS
C
      KOUNT = 0
C
      XRATIO = X - REAL(IS)
      YRATIO = Y - REAL(JS)
C
      AREA11 = (1.0E0 - XRATIO) * (1.0E0 - YRATIO)
      AREA21 = XRATIO * (1.0E0 - YRATIO)
      AREA12 = (1.0E0 - XRATIO) * YRATIO
      AREA22 = XRATIO * YRATIO
C
      IF( MSKSCVH(IS, JS) .EQ. ILAND) THEN
         KOUNT  = KOUNT + 1
         AREA   = AREA11
         IPOINT = IS
         JPOINT = JS
      END IF
C
      IF( MSKSCVH(IS, JP1) .EQ. ILAND ) THEN
         KOUNT = KOUNT +1
         IF (KOUNT .EQ. 1) THEN
            IPOINT = IS
            JPOINT = JP1
         ELSEIF (AREA12 .GT. AREA) THEN
            AREA   = AREA12
            IPOINT = IS
            JPOINT = JP1
         END IF
      END IF
C
      IF( MSKSCVH(IP1, JS) .EQ. ILAND ) THEN
         KOUNT = KOUNT + 1
         IF (KOUNT .EQ. 1) THEN
            AREA   = AREA21
            IPOINT = IP1
            JPOINT = JS
         ELSEIF (AREA21 .GT. AREA) THEN
            AREA   = AREA21
            IPOINT = IP1
            JPOINT = JS
         END IF
      END IF
C
      IF( MSKSCVH(IP1, JP1) .EQ. ILAND ) THEN
         KOUNT = KOUNT + 1
         IF (KOUNT .EQ. 1) THEN
            AREA   = AREA22
            IPOINT = IP1
            JPOINT = JP1
         ELSEIF (AREA22 .GT. AREA) THEN
            AREA   = AREA22
            IPOINT = IP1
            JPOINT = JP1
         END IF
      END IF
C
C     DETERMINE SNO/ICE USING NEAREST NEIGHBOR WITH SAME SFC TYPE 
C
      IF(KOUNT .GT. 0) THEN
         SI(I,J) = SCVH(IPOINT,JPOINT)
      ELSE
C
C         NO IMMEDIATELY SURROUNDING POINTS IN THE 1024 X 1024 FIELD OF
C         SNOW/ICE HAVE THE SAME LAND-SEA TYPE AS THE ETA POINT.  THE
C         ETA POINT MAY BE SMALL ISLAND OR LAKE OR SMALL BAY OR PENNIN.
C         (INVARIABLY A SMALL LAKE IN ETA GRID)
C         SO EXPAND SEARCH RADIUS AND TAKE FIRST SFC TYPE MATCH
C
          IPOINT = NINT(X)
          JPOINT = NINT(Y)
C
          do3346: DO LL=1,4
           JPE = MIN (1024, JPOINT+LL)
           JPB = MAX (1 , JPOINT-LL)
           IPE = MIN (1024, IPOINT+LL)
           IPB = MAX (1 , IPOINT-LL)
C
             doout2346: DO MK=JPB,JPE
             doin2346: DO  NK=IPB,IPE
               IF (MSKSCVH(NK,MK) .EQ. ILAND) THEN
               SI(I,J) = SCVH(NK,MK)
               GO TO 3351
               ENDIF

            END DO doin2346
            END DO doout2346
          END DO do3346
C
C  NO LAND/SEA MASK MATCHES FOUND, SO 
C     A) NORTH OF 55N, WE ASSIGN SNOW/ICE IRRESPECTIVE OF SFC TYPE
C     B) SOUTH OF 55N, WE KEEP A PRIORI ZERO DEFAULT
C   (THE "B" OPTION BEST FOR WARMER LATS OF U.S., WHERE THIS CONDITION 
C   IS VIRTUALLY ALWAYS A SMALL ETA LAKE WITH NO COUNTERPART WATER 
C   NEARBY IN THE NESDIS/IMS GRID, E.G., SALT LAKE, WHERE WE MUST  
C   AVOID GETTING SEA-ICE OWING TO SURROUNDING SNOW COVER)
C
       WRITE (IOUTUPRT, 3347) I,J,YYLAT,XLONG,IPOINT,JPOINT,ILAND
 3347   FORMAT(1H ,'**NO IMS MSK MATCH,',
     1         'ETA-I,J,ELAT,ELON,IMS-I,IMS-J,LND:',2I6,2F7.2,2I3,I2)
C
             IF (YYLAT .GE. 55.0 ) THEN
               SI(I,J) = SCVH(IPOINT,JPOINT) 
             ELSE
               SI(I,J) = 0.0
             ENDIF 
      ENDIF
C
 3351 CONTINUE
C 
      IF (.NOT. LUSAF) GO TO 4300
C
C
C-------------- BEGIN USAF SNOW/ICE INTERPOLATION------------------
c  If current USAF snow/ice anal was successfully read (lusaf=true), 
c  add the USAF information as follows:
C
c  The quality of NESDIS/IMS sea-ice cover and snow coverage (especially
c  in areas with small amounts of snow) are better than the usaf data.
c  We do not use usaf sea-ice cover data at all.  Over land, the
c  presence/absence of snow is determined by the IMS snow coverage, i.e., 
c  at a given location, if the IMS data indicate no snow, then we assume 
c  there is no snow, no matter what the USAF data say.  if the IMS data
c  show snow but the USAF data have no snow or less than 1.5" of snow at
c  the location, we assume there is a 1.5" of snow.  if both the IMS data
c  and USAF data indicate snow, and the USAF snow depth is more than 1.5",
c  then we use the USAF snow depth.
c
c     On the 512 x 512 USAF grid, the data values are 0.0 for
c     no snow or ice, 11.0 for ice points, and positive but
c     less than 11.0 for snow points.
C
C--------------------------------------IF ETA SEA POINT, SKIP USAF ANL
      IF ( SM(I,J) .GT. 0.9 ) GO TO 4300
C
C If the IMS data indicate no snow on this point, also skip USAF ANL:
      IF (SI(I,J).LT.0.001) GO TO 4300
C
C Otherwise, set snow depth to be the lower threshold, 1.5":
      SI(I,J) = DEPTH
C
C-------------------------------------------------------------------
C  THIS IS AN ETA LAND POINT, SO APPLY USAF SNOWDEPTH ANAL
C
C  DETERMINE LOCATION OF ETA POINT ON THE 512 X 512 USAF GRID
C
         RM= ENNMCAF*COS(YYLAT * DTR) / (1.E0 + SIN(YYLAT * DTR))
         RAD = (XLONG - ALNMCAF) * DTR
         X = XPNMCAF + RM * COS(RAD)
         Y = YPNMCAF + RM * SIN(RAD)
C
      IS  = INT(X)
      IP1 = IS + 1
      JS  = INT(Y)
      JP1 = JS + 1
C
C-----IF OUTSIDE OF USAF GRID DOMAIN (I.E. S.H.) WE KEEP ZERO DEFAULT--
C
      IF ((IS .LT. 1) .OR. (IS .GT. 511) .OR. (JS .LT. 1)
     1         .OR. (JS .GT. 511))  THEN
        GO TO 4300
      ENDIF
C--------------------------------------------------------------------
C
C  NOW USE ETA AND USAF LAND-SEA MASK TO ENSURE ONLY LAND POINTS ARE 
C  INTERPOLATED TO LAND POINTS (TO DETERMINE SNOW)
C   (USAF LAND MASK: SEA=1,LAND=2,COASTAL-LAND=4,OFFWORLD=9)
C   NOTE: IN REACHING THIS STAGE, WE HAVE ALREADY INSURED WE ARE ON 
C   ETA LAND POINT AND IN N.H., I.E. NOT MSKAF=9 (I.E. NOT OFFWORLD)
C  
      ILAND = 2
C
      IPOINT = NINT(X)
      JPOINT = NINT(Y)
      IF ( MSKAF(IPOINT,JPOINT) .GE. ILAND) THEN
        SI(I,J) = AFSNO(IPOINT,JPOINT)
        SI(I,J) = AMAX1(DEPTH,SI(I,J))
        GO TO 4300
      ENDIF
C
      KOUNT = 0
C
      XRATIO = X - REAL(IS)
      YRATIO = Y - REAL(JS)
C
      AREA11 = (1.0E0 - XRATIO) * (1.0E0 - YRATIO)
      AREA21 = XRATIO * (1.0E0 - YRATIO)
      AREA12 = (1.0E0 - XRATIO) * YRATIO
      AREA22 = XRATIO * YRATIO
C
      IF( MSKAF(IS, JS) .GE. ILAND) THEN
         KOUNT  = KOUNT + 1
         AREA   = AREA11
         IPOINT = IS
         JPOINT = JS
      END IF
C
      IF( MSKAF(IS, JP1) .GE. ILAND ) THEN
         KOUNT = KOUNT +1
         IF (KOUNT .EQ. 1) THEN
            IPOINT = IS
            JPOINT = JP1
         ELSEIF (AREA12 .GT. AREA) THEN
            AREA   = AREA12
            IPOINT = IS
            JPOINT = JP1
         END IF
      END IF
C
      IF( MSKAF(IP1, JS) .GE. ILAND ) THEN
         KOUNT = KOUNT + 1
         IF (KOUNT .EQ. 1) THEN
            AREA   = AREA21
            IPOINT = IP1
            JPOINT = JS
         ELSEIF (AREA21 .GT. AREA) THEN
            AREA   = AREA21
            IPOINT = IP1
            JPOINT = JS
         END IF
      END IF
C
      IF( MSKAF(IP1, JP1) .GE. ILAND ) THEN
         KOUNT = KOUNT + 1
         IF (KOUNT .EQ. 1) THEN
            AREA   = AREA22
            IPOINT = IP1
            JPOINT = JP1
         ELSEIF (AREA22 .GT. AREA) THEN
            AREA   = AREA22
            IPOINT = IP1
            JPOINT = JP1
         END IF
      END IF
C
C     DETERMINE SNO/ICE CONSIDERING THE NUMBER OF POINTS SURROUNDING
C     ETA GRID POINT WITH THE SAME LAND-SEA MASK FLAG
C
      IF (KOUNT .GT. 0) THEN
          SI(I,J) = AFSNO(IPOINT,JPOINT)
          SI(I,J) = AMAX1(DEPTH,SI(I,J))
C
      ELSE
C
C         NO IMMEDIATELY SURROUNDING POINTS IN THE 512 X 512 FIELD OF
C         SNO/ICE HAVE THE SAME LAND-SEA MASK AS THE ETA POINT.
C         THE ETA POINT MAY BE AN ISLAND OR LAKE OR BAY OR PENNINSULA.
C         SO EXPAND SEARCH RADIUS AND TAKE FIRST MASK FLAG MATCH
C
          IPOINT = NINT(X)
          JPOINT = NINT(Y)
C
          do7346: DO LL=1,7
           JPE = MIN (512, JPOINT+LL)
           JPB = MAX (1 ,  JPOINT-LL)
           IPE = MIN (512, IPOINT+LL)
           IPB = MAX (1 ,  IPOINT-LL)
C
             doout6346: DO MK=JPB,JPE
             doin6346:  DO NK=IPB,IPE
               IF (MSKAF(NK,MK) .GE. ILAND) THEN
               SI(I,J) = AFSNO(NK,MK)
               SI(I,J) = AMAX1(DEPTH,SI(I,J))
               GO TO 4300
               ENDIF
             END DO doin6346
             END DO doout6346
         END DO do7346
C
C   NO LAND MASK MATCH FOUND, SO WE PRINT WARNING AND STAY
C   WITH EARLIER VALUE DETERMINED FROM NESDIS/IMS ANAL, WHICH 
C   WE CONVERT FROM COVER FLAG (0,1) TO DEFAULT DEPTH OF 1.5".
C
          SI(I,J) = DEPTH
          WRITE (IOUTUPRT, 7347) I,J,YYLAT,XLONG,IS,JS,SI(I,J)
 7347      FORMAT(1H ,'*** WARNING ***..NO USAF LAND MSK MATCH ',
     1 ' AT ETA-I,J,ELAT,ELON,USAF-I,J:',2I6,2F7.2,2I3/
     2 1H ,' DEFAULT TO NESDIS VALUE OF ',F5.1,' TIMES .10 M')
C
      ENDIF
C     
C
C******************** END MAJOR ETA GRID POINT LOOP ****************
C
 4300 CONTINUE
         ENDDO
         ENDDO
C
C  PRINT SAMPLE OF SNOW/ICE ON ETA GRID
C
      CALL PRINTETA (SI,SM)
Ctest CALL PRINTYL(SI,SM)

      RETURN

      END
