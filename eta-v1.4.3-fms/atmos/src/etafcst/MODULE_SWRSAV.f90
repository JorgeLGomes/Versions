    MODULE SWRSAV
!>--------------------------------------------------------------------------------------------------
!> MODULE SWRSAV
!>
!> USE MODULES: F77KINDS
!>              RDPARM
!>
!> DRIVER     : GFDLRD
!>              RADTN
!>              GFDLRD
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE RDPARM  , ONLY : NB
!
    IMPLICIT NONE
!
    SAVE
!
    REAL   (KIND=R4KIND), DIMENSION(NB)                                                         ::&
    & ABCFF   , PWTS
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CFCO2   , CFO3    , REFLO3  , RRAYAV
!-------------------------------------------------------------------------------------------------- 
! SPECIFICATION OF DATA STATEMENTS:
! ABCFF - ABSORPTION COEFFICIENTS FOR BANDS IN K-DISTRIBUTION.
! ORIGINALLY GIVEN BY LACIS AND HANSEN, REVISED BY RAMASWAMY
! PWTS - CORRESPONDING WEIGHTS ASSIGNED TO BANDS IN THE K-DISTRIBUTION
! REFLO3, RRAYAV - REFLECTION COEFFICIENTS GIVEN BY LACIS AND HANSEN TO ACCOUNT FOR EFFECTS OF 
!                  RAYLEIGH SCATTERING IN THE VISIBLE FREQUENCIES (BAND 1)
! CFCO2, CFO3 = CONVERSION FACTORS FROM GM / CM ** 2 TO CM-ATM(STP)
!
! THE FOLLOWING ARE THE COEFFICIENTS FOR THE 12-BAND SHORTWAVE RADIATION CODE, SPECIFIED BY
! RAMASWAMY.
!--------------------------------------------------------------------------------------------------
    DATA ABCFF /                                                                                  &
    & 2*4.0E-5, .002   , .035 , .377 , 1.95 , 9.40 , 44.6 , 190. , 989. , 2706.  , 39011./
!
    DATA PWTS /                                                                                   &
    & .5000  , .121416, .0698, .1558, .0631, .0362, .0243, .0158, .0087, .001467, .002342, .001075/
!--------------------------------------------------------------------------------------------------
! THE ORIGINAL 9-BAND LACIS-HANSEN COEFFICIENTS ARE GIVEN HERE; IT THE USER INSISTS ON USING THESE
! VALUES, SHE MUST ALSO CHANGE THE PARAMETER NB FROM 12 TO 9.
! THIS PARAMETER IS DEFINED IN RDPARM. H. NO OTHER CHANGES ARE REQUIRED
!--------------------------------------------------------------------------------------------------
    DATA CFCO2  / 508.96   /
    DATA CFO3   / 466.64   /
    DATA REFLO3 /   1.9    /
    DATA RRAYAV /   0.144  /
!
    END MODULE SWRSAV

