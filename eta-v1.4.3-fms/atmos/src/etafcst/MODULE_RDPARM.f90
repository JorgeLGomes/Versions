    MODULE RDPARM
!>--------------------------------------------------------------------------------------------------
!> MODULE RDPARM
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : GFDLRD
!>              GOSSIP
!>              GRADFS
!>              MODULE_BANDTA
!>              MODULE_BDCOMB
!>              MODULE_CO2BD2
!>              MODULE_CO2BD3
!>              MODULE_CO2BD4
!>              MODULE_CO2BD5
!>              MODULE_CTLBLK
!>              MODULE_DYNAM
!>              MODULE_INPUT
!>              MODULE_MAPOT
!>              MODULE_PHYS
!>              MODULE_SCRTCH
!>              MODULE_SWRSAV
!>              MODULE_TABCOM
!>              MODULE_TBLTMP
!>              RADFS
!>              SWR93
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
    SAVE
!------------------------------------------------------------------------------------------
! PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
! IMAX   =  NO. POINTS ALONG THE LAT. CIRCLE USED IN CALCS.
! L      =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
! NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE BANDTA FOR DEFINITION
! NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS
! NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE BDCOMB FOR DEFINITION
! INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
! NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS.
!         
! NB, KO2 ARE SHORTWAVE PARAMETERS; OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS.
!------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: L      =  LM
    INTEGER(KIND=I4KIND), PARAMETER :: IMAX   =  IM
    INTEGER(KIND=I4KIND), PARAMETER :: NBLW   = 163          
    INTEGER(KIND=I4KIND), PARAMETER :: NBLX   =  47
    INTEGER(KIND=I4KIND), PARAMETER :: NBLY   =  15
!
    INTEGER(KIND=I4KIND), PARAMETER :: NBLM   = NBLY - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1    = L    + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP2    = L    + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LP3    = L    + 3
!
    INTEGER(KIND=I4KIND), PARAMETER :: LM1    = L    - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM2    = L    - 2
    INTEGER(KIND=I4KIND), PARAMETER :: LM3    = L    - 3
!
    INTEGER(KIND=I4KIND), PARAMETER :: LL     =   2  * L
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1   = LL   + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LLP2   = LL   + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LLP3   = LL   + 3
!
    INTEGER(KIND=I4KIND), PARAMETER :: LLM1   = LL   - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LLM2   = LL   - 2
    INTEGER(KIND=I4KIND), PARAMETER :: LLM3   = LL   - 3
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M   = LP1  * LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M1  = LP1M - 1
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1V   = LP1  * (1 + 2 * L / 2)
    INTEGER(KIND=I4KIND), PARAMETER :: LP121  = LP1  * NBLY
!
    INTEGER(KIND=I4KIND), PARAMETER :: LL3P   =   3  *  L + 2
!
    INTEGER(KIND=I4KIND), PARAMETER :: NB     =  12
    INTEGER(KIND=I4KIND), PARAMETER :: NB1    = NB    - 1
!
    INTEGER(KIND=I4KIND), PARAMETER :: INLTE  =   3
    INTEGER(KIND=I4KIND), PARAMETER :: INLTEP = INLTE + 1 
!
    INTEGER(KIND=I4KIND), PARAMETER :: NNLTE  =  56
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1I   = IMAX  * LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1I  = IMAX  * LLP1
    INTEGER(KIND=I4KIND), PARAMETER :: LL3PI  = IMAX  * LL3P
!
    INTEGER(KIND=I4KIND), PARAMETER :: KO2    =  12
    INTEGER(KIND=I4KIND), PARAMETER :: KO21   = KO2   + 1 
    INTEGER(KIND=I4KIND), PARAMETER :: KO2M   = KO2   - 1
!
    END MODULE RDPARM
