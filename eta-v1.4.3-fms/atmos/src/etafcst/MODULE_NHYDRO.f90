    MODULE NHYDRO
!>--------------------------------------------------------------------------------------------------
!> MODULE NHYDRO
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : BOCOH
!>              BOCOHF
!>              CHKOUT
!>              DIVHOA
!>              DIVHOAST
!>              EBU
!>              EPS
!>              GOSSIP
!>              HADZ
!>              INIT
!>              INITS
!>              PDTEDT
!>              PGCOR
!>              VADZ
!>              VTADV
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : LM, IDIM1, IDIM2, JDIM1, JDIM2
!
    IMPLICIT NONE
!
    SAVE
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & HYDRO   , SPLINE
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                               ::&
    & DWDT    , PDWDT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM+1)                             ::&
    & PINT    ,                                                                                   &
    & W       ,                                                                                   &
    & Z     
!
    END MODULE NHYDRO                             
