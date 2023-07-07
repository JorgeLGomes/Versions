	SUBROUTINE VFRAC(im,jm,glatr,glonr,TLM0D,TPH0D,DLMD,DPHD,SM,
     +  SICE,VEGFRC)

	real*4 glat(IM,JM), glon(IM,JM)
	real*4 glatr(IM,JM), glonr(IM,JM)
	real*4 TLM0D,TPH0D,DLMD,DPHD
	real*4 sm(im,jm),sice(im,jm),vegfrc(im,jm)
	real*4 dum2(im,jm),dum2b(im,jm)
	integer*4 JULD,MON2,MON1,day1,day2,JULM(13)
        real*4 wght1,wght2,rday,dtr
        INTEGER IDAT(3)
Chou 20080202
        LOGICAL LVEG
	REAL*4 VEGFRM(im,jm)     ! monthly veg greeness	

      data JULM/0,31,59,90,120,151,181,212,243,273,304,334,365/

       INTEGER KPDS(200),KGDS(200),JPDS(200),JGDS(200),KF,KNUM
C      REAL*4 FPAR1(2500,1250),FPAR2(2500,1250),FSVE

       REAL FPAR1(15680,10600),FPAR2(15680,10600),FSVE
       INTEGER RECORDLEN

C       LOGICAL BITMAP(2500,1250)

        common /mytime/idat

C    FOURTH, MONTHLY GREEN VEG FRACTION, INTERPOLATE TO DAY OF YEAR.
C    THE MONTHLY FIELDS ASSUMED VALID AT 15TH OF MONTH.
C    VALUES RANGE FROM 0.001 TO 1.0 OVER LAND, 0.0 OVER WATER.
C    INTERPOLATE IN SPACE.
C
C    READ FROM A FILE THAT HAS 13 RECORDS, ONE PER MONTH WITH
C    JANUARY OCCURRING BOTH AS RECORD 1 AND AGAIN AS RECORD 13,
C    THE LATTER TO SIMPLIFY TIME INTERPOLATION FOR DAYS
C    BETWEEN DEC 16 AND JAN 15. WE TREAT JAN 1 TO JAN 15
C    AS JULIAN DAYS 366 TO 380 BELOW, I.E WRAP AROUND YEAR.
C *** THE FOLLOWING PART WAS REVISED BY F. CHEN 7/96 TO REFLECT
C        A NEW NESDIS VEGETATION FRACTION PRODUCT (FIVE-YEAR
C        CLIMATOLOGY WITH 0.144 DEGREE RESOLUTION
C        FROM 89.928S, 180W TO 89.928N, 180E)
C
       REWIND 38
C
C ****  DO TIME INTERPOLATION ****
       JULD=JULM(IDAT(1))+IDAT(2)
       IF(JULD.LE.15) JULD=JULD+365
       MON2=IDAT(1)
       IF(IDAT(2).GT.15) MON2=MON2+1
       IF(MON2.EQ.1) MON2=13
       MON1=MON2-1
C **** ASSUME DATA VALID AT 15TH OF MONTH
       DAY2=JULM(MON2)+15
       DAY1=JULM(MON1)+15
       RDAY=JULD
       WGHT1=(DAY2-RDAY)/(DAY2-DAY1)
       WGHT2=(RDAY-DAY1)/(DAY2-DAY1)
C
       RECORDLEN = 4*15680
       OPEN(38,FILE='FCOVER.bin',
     &               ACCESS='direct',FORM='UNFORMATTED',RECL=RECORDLEN)
       DO J=1,10600
          READ(38,REC=((10600*(MON1-1))+J)) FPAR1(1:15680,J)
       ENDDO
      
       DO J=1,10600
          READ(38,REC=((10600*(MON2-1))+J)) FPAR2(1:15680,J)
       ENDDO
       CLOSE(38)
C
       DO JJ=1,10600
        DO I=1,15680
         FSVE=WGHT1*FPAR1(I,JJ)+WGHT2*FPAR2(I,JJ)
         FPAR1(I,JJ)=FSVE
        END DO
       END DO
c
C ** SPACE INTERPOLATION OF FPAR1 TO E GRID
        dtr=acos(-1.)/180.
        glat=glatr/dtr
        glon=glonr/dtr
       WRITE(6,*) 'lat and lon for PUTVEG ', glat(im/2,jm/2),
     &             glon(im/2,jm/2)
       CALL PUTVEG(TLM0D,TPH0D,DLMD,DPHD,GLAT,GLON,
     &             FPAR1,SM,SICE,vegfrc)
       WRITE(6,*) 'veg fracs ', vegfrc(im/2,jm/2),vegfrc(im/2+5,jm/2+5)

Chou  CREATE FILE of MONTHLY VEG GREENESS INTERPOLATEd TO Egrid in IM,JM
CJLG   VGREEN_12mo.dat done in vgreen module
       LVEG=.FALSE.  

       IF(LVEG) then

       OPEN(83,file='VGREEN_12MO.dat',form='unformatted',
     &        status='unknown')


       OPEN(38,FILE='FCOVER.bin',
     &               ACCESS='direct',FORM='UNFORMATTED',RECL=RECORDLEN)

       JMS=0
C       JME=10600

       DO imo=1,12
          
       DO J= 1 ,10600
          READ(38,REC=J+JMS) FPAR1(1:15680,J)
       ENDDO
       
C interpolate to Egrid.
       CALL PUTVEG(TLM0D,TPH0D,DLMD,DPHD,GLAT,GLON,
     &             FPAR1,SM,SICE,vegfrm)
     
       CALL PUTVEG(TLM0D,TPH0D,DLMD,DPHD,GLAT,GLON,
     &             FPAR1,SM,SICE,vegfrm)

       write(83) vegfrm

       JMS=JMS+10600
C       JME=JME+10600      
       ENDDO  ! END imo Loop
 
       CLOSE(38)
      
       ENDIF
Chou

       RETURN
       END
