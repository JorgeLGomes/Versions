    MODULE TABCOM
!>--------------------------------------------------------------------------------------------------
!> MODULE TABCOM
!>
!> ABSTRACT:
!> MODULE TABCOM CONTAINS QUANTITIES PRECOMPUTED IN SUBROUTINE TABLE FOR USE IN THE LONGWAVE 
!> RADIATION PROGRAM:
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!> 
!> DRIVER     : E1E290
!>              E290
!>              E2SPEC
!>              E3V88
!>              FST88
!>              RADFS
!>              TABLE 
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : IMAX, LP1, LP1V, NBLY
!
    IMPLICIT NONE
!
    SAVE
!--------------------------------------------------------------------------------------------
! EM1    - E1 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! EM1WDE - E1 FUNCTION, EVALUATED OVER THE 160-560 CM-1 INTERVAL
! TABLE1 - E2 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! TABLE2 - TEMPERATURE DERIVATIVE OF TABLE1
! TABLE3 - MASS DERIVATIVE OF TABLE1
! EM3    - E3 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS
! SOURCE - PLANCK FUNCTION, EVALUATED AT SPECIFIED TEMPS. FOR BANDS USED IN CTS CALCULATIONS
! DSRCE  - TEMPERATURE DERIVATIVE OF SOURCE
! IND    - INDEX, WITH VALUE IND(I)=I. USED IN FST88
! INDX2  - INDEX VALUES USED IN OBTAINING "LOWER TRIANGLE" ELEMENTS OF AVEPHI, ETC., IN FST88
! KMAXV  - INDEX VALUES USED IN OBTAINING "UPPER TRIANGLE" ELEMENTS OF AVEPHI, ETC., IN FST88
! KMAXVM - KMAXV(L),USED FOR DO LOOP INDICES
!--------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), DIMENSION(IMAX)                                                       ::&
    & IND
!
    INTEGER(KIND=I4KIND), DIMENSION(LP1V)                                                       ::&
    & INDX2
!
    INTEGER(KIND=I4KIND), DIMENSION(LP1)                                                        ::&
    & KMAXV
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & KMAXVM  
!
    REAL   (KIND=R4KIND), DIMENSION(28, 180)                                                    ::&
    & EM1     , EM1WDE  , TABLE1  , TABLE2  , TABLE3  , EM3
!
    REAL   (KIND=R4KIND), DIMENSION(28, NBLY)                                                   ::&
    & SOURCE  , DSRCE
!-----------------------------------  
! VARIABLES EQUIVALENCED FROM E1E290
!-----------------------------------  
    REAL   (KIND=R4KIND), DIMENSION(5040)                                                       ::&
    & T1      , T2      , T4      , EM1V    , EM1VW 
!
    EQUIVALENCE (EM1V(1),    EM1(1,1)), (EM1VW(1), EM1WDE(1,1))
    EQUIVALENCE (  T1(1), TABLE1(1,1)), (   T2(1), TABLE2(1,1))
    EQUIVALENCE (  T4(1), TABLE3(1,1)) 
!----------------------------------  
! VARIABLES EQUIVALENCED FROM E3V88
!----------------------------------
    REAL   (KIND=I4KIND), DIMENSION(5040)                                                       ::&
    & EM3V
!
    EQUIVALENCE (EM3V(1) , EM3(1,1))  
!
    END MODULE TABCOM
