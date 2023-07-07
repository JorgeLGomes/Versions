    MODULE SCRTCH
!>--------------------------------------------------------------------------------------------------
!> MODULE SCRTCH
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!>
!> DRIVER     : TABLE
!>--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! MODULE TABCOM CONTAINS QUANTITIES PRECOMPUTED IN SUBROUTINE TABLE FOR USE IN THE LONGWAVE     
! RADIATION PROGRAM:
! EM1     -  E1 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! EM1WDE  -  E1 FUNCTION, EVALUATED OVER THE 160-560 CM-1 INTERVAL
! TABLE1  -  E2 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! TABLE2  -  TEMPERATURE DERIVATIVE OF TABLE1
! TABLE3  -  MASS DERIVATIVE OF TABLE1
! EM3     -  E3 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! SOURCE  -  PLANCK FUNCTION, EVALUATED AT SPECIFIED TEMPS. FOR BANDS USED IN CTS CALCULATIONS
! DSRCE   -  TEMPERATURE DERIVATIVE OF SOURCE
! IND     -  INDEX, WITH VALUE IND(I) = I. USED IN FST88
! INDX2   -  INDEX VALUES USED IN OBTAINING "LOWER TRIANGLE" ELEMENTS OF AVEPHI,ETC., IN FST88
! KMAXV   -  INDEX VALUES USED IN OBTAINING "UPPER TRIANGLE" ELEMENTS OF AVEPHI,ETC., IN FST88
! KMAXVM  -  KMAXV(L), USED FOR DO LOOP INDICES
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : NBLX, NBLW
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(28,180)                                                     ::&
    & SUM     , PERTSM  , SUM3    , SUMWDE
!
    REAL   (KIND=R4KIND), DIMENSION(28,NBLX)                                                    ::&
    & SRCWD
!
    REAL   (KIND=R4KIND), DIMENSION(28,NBLW)                                                    ::&
    & SRC1NB  , DBDTNB
!
    REAL   (KIND=R4KIND), DIMENSION(181)                                                        ::&
    & ZMASS   , ZROOT
!
    REAL   (KIND=R4KIND), DIMENSION(28)                                                         ::&
    & SC      , DSC     , XTEMV   , TFOUR   , FORTCU  , X       , X1      ,                       &
    & SRCS    , SUM4    , SUM6    , SUM7    , SUM8    , SUM4WD  , R1      ,                       &
    & R2      , S2      , T3      , R1WD
!
    REAL   (KIND=R4KIND), DIMENSION(180)                                                        ::&
    & X2      , EXPO    , FAC
!
    REAL   (KIND=R4KIND), DIMENSION(30)                                                         ::&
    & CNUSB   , DNUSB
!
    REAL   (KIND=R4KIND), DIMENSION(NBLW)                                                       ::&
    & ALFANB  , AROTNB
!
    REAL   (KIND=R4KIND), DIMENSION(NBLW)                                                       ::&
    & ANB     , BNB     , CENTNB  , DELNB   , BETANB
!
    END MODULE SCRTCH
