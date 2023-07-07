PROGRAM VFRAC

    USE PARMCONF
    USE F77KINDS

    REAL   (KIND=R4KIND)   , DIMENSION(15680,10600)    :: FPAR1,FPAR2
    INTEGER(KIND=R4KIND)   , PARAMETER                 :: RECORDLEN = 4*15680

    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: GLAT
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: GLON
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: GLATR
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: GLONR
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: HGT
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: SM
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: SICE
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: VEGFRC    
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: VEGFRM  ! monthly veg greeness	

        integer*4 JULD,MON2,MON1,day1,day2,JULM(13)
        real*4 wght1,wght2,rday,dtr

!       A NEW VEGETATION FRACTION PRODUCT (1999-2014
!       CLIMATOLOGY WITH 0.00892857 DEGREE RESOLUTION
!       Fuster et al (2018)  
!

!   MAKE SICE ZERO
    SICE = 0

!   READ SEA-LAND MASK ETA-GRID
    READ(10) HGT,SM
!    
    DO J=JM,1,-1
!        WRITE(1200,'(921I2)')(INT(SM(I,J)),I=1,IM)
    ENDDO
!    CLOSE(1200)

    CALL GEOGRAPHIC(GLATR, GLONR,IM,JM,DLMD,DPHD,TLM0D,TPH0D)

    REWIND 38

!Chou/Andre  CREATE FILE of MONTHLY VEG GREENESS INTERPOLATEd TO Egrid in IM,JM

       OPEN(83,FILE='VGREEN_12MO.dat',form='unformatted',status='unknown')

       OPEN(38,FILE='FCOVER.bin',ACCESS='direct',FORM='UNFORMATTED',RECL=RECORDLEN)

       JMS=0

        DTR=ACOS(-1.)/180.
        GLAT=GLATR/DTR
        GLON=GLONR/DTR
       
       DO J= 1 ,10600
          READ(38,REC=J+JMS) FPAR1(1:15680,J)
       ENDDO
       
!       CALL PUTVEG(TLM0D,TPH0D,DLMD,DPHD,GLAT,GLON,FPAR1,SM,SICE,vegfrm)

       DO imo=1,12
          
       DO J= 1 ,10600
          READ(38,REC=J+JMS) FPAR1(1:15680,J)
       ENDDO
!	
! SPACE INTERPOLATION OF FPAR1 TO E GRID

       WRITE(6,*) 'lat and lon for PUTVEG ', GLAT(IM/2,JM/2), GLON(IM/2,JM/2)

! INTERPOLATE to Egrid.
       CALL PUTVEG(GLAT,GLON,FPAR1,SM,SICE,vegfrm)

       WRITE(6,*) 'veg fracs ', VEGFRC(IM/2,JM/2),VEGFRC(IM/2+5,JM/2+5)
     
!       CALL PUTVEG(TLM0D,TPH0D,DLMD,DPHD,GLAT,GLON,FPAR1,SM,SICE,vegfrm)

       WRITE(83) VEGFRM

       JMS=JMS+10600

       ENDDO  ! END imo Loop
 
       CLOSE(38)
      
!Chou/Andre

END
