    MODULE IMASS_TABLES
!>--------------------------------------------------------------------------------------------------
!> MODULE IMASS_TABLES
!>
!> USE MODULES: F77KINDS
!>              IACCR_TABLES
!>
!> DRIVER     : GSMCOLUMN
!>              GSMCONST
!>--------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------
! IMASS_TABLES - LOOKUP TABLES FOR MASS CONTENT OF ICE
!-----------------------------------------------------
!
!---------------------------------------------------- 
! MASSI  - INTEGRATED QUANTITY ASSOCIATED W/ ICE MASS
!---------------------------------------------------- 
    USE F77KINDS
    USE IACCR_TABLES, ONLY : MDIMIN, MDIMAX
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(MDIMIN:MDIMAX)                                              ::&
    & MASSI
!
    END MODULE IMASS_TABLES
