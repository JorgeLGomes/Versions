    SUBROUTINE CED_IJ(GLAT,GLON,XINT_I4,YINT_I4                                         &
               &      ,XYDEC_R8,XDEC_R8,YDEC_R8                                         &
               &      ,NX,NY,LAT_P1,LON_P1,LAT_PNY,LON_PNX,INCR_LON,INCR_LAT)
    USE CONSTANTS
    USE DIAGNOSTIC
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND)                          , INTENT(INOUT)           ::  GLAT
    REAL   (KIND=R4KIND)                          , INTENT(INOUT)           ::  GLON
    INTEGER(KIND=I4KIND)                          , INTENT(INOUT)           ::  XINT_I4
    INTEGER(KIND=I4KIND)                          , INTENT(INOUT)           ::  YINT_I4
    REAL   (KIND=R4KIND)                          , INTENT(INOUT)           ::  XYDEC_R8
    REAL   (KIND=R4KIND)                          , INTENT(INOUT)           ::  XDEC_R8
    REAL   (KIND=R4KIND)                          , INTENT(INOUT)           ::  YDEC_R8
    INTEGER(KIND=I8KIND)                          , INTENT(IN)              ::  NX
    INTEGER(KIND=I8KIND)                          , INTENT(IN)              ::  NY
    REAL   (KIND=R4KIND)                          , INTENT(IN)              ::  LAT_P1
    REAL   (KIND=R4KIND)                          , INTENT(IN)              ::  LON_P1
    REAL   (KIND=R4KIND)                          , INTENT(IN)              ::  LON_PNX
    REAL   (KIND=R4KIND)                          , INTENT(IN)              ::  LAT_PNY
    REAL   (KIND=R4KIND)                          , INTENT(IN)              ::  INCR_LON
    REAL   (KIND=R4KIND)                          , INTENT(IN)              ::  INCR_LAT
!
    REAL   (KIND=R4KIND)                                                    ::  XPTS_R8
    REAL   (KIND=R4KIND)                                                    ::  YPTS_R8

    REAL   (KIND=R4KIND)                                                    ::  TMP1
    REAL   (KIND=R4KIND)                                                    ::  XMAX
    REAL   (KIND=R4KIND)                                                    ::  XMIN
    REAL   (KIND=R4KIND)                                                    ::  YMAX
    REAL   (KIND=R4KIND)                                                    ::  YMIN
    REAL   (KIND=R4KIND)                                                    ::  FILL
    REAL   (KIND=R4KIND)                                                    ::  HI
    REAL   (KIND=R4KIND)                                                    ::  HJ
    INTEGER(KIND=I4KIND)                                                    ::  NRET
!
      fill=-999
!
!      HI=(-1.)**ISCAN      ! Compilation problem
!      HJ=(-1.)**(1-JSCAN)  ! Compilation problem
! 
      HI= 1
      HJ= 1
      XMIN=0
      XMAX=NX+1
      IF(NX.EQ.NINT(360/ABS(INCR_LON))) XMAX=NX+2
! 
      YMIN=0
      YMAX=NY+1
! 
      NRET=0
!
      IF(ABS(GLON).LE.360.AND.ABS(GLAT).LE.90) THEN
!
         TMP1=MOD(HI*(GLON-LON_P1)+3600,360.)
         XPTS_R8=(1+HI*TMP1* INCR_LON**(-1))
!              
         TMP1=(GLAT-LAT_P1)
         YPTS_R8=1+TMP1*INCR_LAT**(-1)
!
         IF(XPTS_R8.GE.XMIN.AND.XPTS_R8.LE.XMAX.AND. &
      &     YPTS_R8.GE.YMIN.AND.YPTS_R8.LE.YMAX) THEN
            NRET=NRET+1
         ELSE
            XPTS_R8=FILL
            YPTS_R8=FILL
         ENDIF
      ELSE
         XPTS_R8=FILL
         YPTS_R8=FILL
      ENDIF
!
      if (XPTS_R8 .lt. 1. .or. XPTS_R8 .gt. NX .or. &
      &   YPTS_R8 .lt. 1. .or. YPTS_R8 .gt. NY) then
          if (XPTS_R8 .lt. 1) XPTS_R8=XPTS_R8+NX
          if (XPTS_R8 .gt. NX) XPTS_R8=XPTS_R8-NX
      endif        
      XDEC_R8=XPTS_R8-FLOOR(XPTS_R8)
      XINT_I4=INT(XPTS_R8) 
      YDEC_R8=YPTS_R8-FLOOR(YPTS_R8)
      YINT_I4=INT(YPTS_R8)  
      XYDEC_R8=XDEC_R8*YDEC_R8
    END SUBROUTINE CED_IJ
       

