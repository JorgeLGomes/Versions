    MODULE INPUT
!>--------------------------------------------------------------------------------------------------
!> MODULE INPUT
!>
!> ABSTRACT:
!> THIS MODULE IS FOR THE INPUT DATA WHEN IT IS DERIVED INTERNALLY. 
!> IN RADMN, THESE DATA ARE TRANSFERRED INTO THE EAL INPUT MODULE, RADISW.
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!>
!> DRIVER     : GFDLRD
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : L, LP1    
!
    IMPLICIT NONE
!
    SAVE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ICH     , ICM     , ICT     , ICB
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CH      , CM      , CL      , EMCH    , EMCM    , EMCL
!
    REAL   (KIND=R4KIND), DIMENSION(L)                                                          ::&
    & RR      , QQO3
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & DTEMP   , PPRESS 
!------------
! FROM GFDLRD
!------------
    DATA CH  , CM  , CL  , ICH, ICM, ICT, ICB / .159,  .07,  .269, 5, 11, 12, 14/
    DATA EMCH, EMCM, EMCL                     /1.   , 1.  , 1.                  /
!
    END MODULE INPUT


