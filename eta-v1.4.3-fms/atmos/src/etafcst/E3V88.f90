    SUBROUTINE E3V88(EMV, TV, AV)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE E3V88
!>
!> SUBPROGRAM: E3V88 - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> SUBROUTINE E3V88 COMPUTES NEARBY LAYER TRANSMISSIVITIES FOR H2O USING A TABLE LOOKUP OF THE PRE-
!> COMPUTED E3 FUNCTION (DESCRIBED IN REF. (4)).
!> 
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> TV  - 
!> AV  -
!>
!> OUTPUT ARGUMENT LIST:
!> EM3 -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              HCON
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              PHYCON
!>              TABCOM
!>              TEMPCOM
!>              TOPO
!>  
!> DRIVER     : FST88
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE HCON
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE PHYCON
    USE TABCOM
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
!--------------------------------------------------------------------------------------------------
! PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
! IMAX   =  NO. POINTS ALONG THE LAT. CIRCLE USED IN CALCS.
! L      =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS 
! NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE BANDTA FOR DEFINITION              
! NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS
! NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE BDCOMB FOR DEFINITION  
! INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
! NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS. NB,KO2 ARE SHORTWAVE PARAMETERS;
!           OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS.
!--------------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: L      = LM
    INTEGER(KIND=I4KIND), PARAMETER :: NCOL   = IMAX
    INTEGER(KIND=I4KIND), PARAMETER :: NBLW   = 163
    INTEGER(KIND=I4KIND), PARAMETER :: NBLX   = 47
    INTEGER(KIND=I4KIND), PARAMETER :: NBLM   = NBLY  - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP2    = L     + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LP3    = L     + 3
    INTEGER(KIND=I4KIND), PARAMETER :: LM1    = L     - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM2    = L     - 2
    INTEGER(KIND=I4KIND), PARAMETER :: LM3    = L     - 3
    INTEGER(KIND=I4KIND), PARAMETER :: LL     = 2     * L
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1   = LL    + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LLP2   = LL    + 2
    INTEGER(KIND=I4KIND), PARAMETER :: LLP3   = LL    + 3
    INTEGER(KIND=I4KIND), PARAMETER :: LLM1   = LL    - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LLM2   = LL    - 2
    INTEGER(KIND=I4KIND), PARAMETER :: LLM3   = LL    - 3
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M   = LP1   * LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M1  = LP1M  - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP121  = LP1   * NBLY
    INTEGER(KIND=I4KIND), PARAMETER :: LL3P   = 3     * L    + 2
    INTEGER(KIND=I4KIND), PARAMETER :: NB     = 12
    INTEGER(KIND=I4KIND), PARAMETER :: INLTE  = 3
    INTEGER(KIND=I4KIND), PARAMETER :: INLTEP = INLTE + 1
    INTEGER(KIND=I4KIND), PARAMETER :: NNLTE  = 56
    INTEGER(KIND=I4KIND), PARAMETER :: LP1I   = IMAX  * LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1I  = IMAX  * LLP1
    INTEGER(KIND=I4KIND), PARAMETER :: LL3PI  = IMAX  * LL3P
    INTEGER(KIND=I4KIND), PARAMETER :: NB1    = NB    - 1
    INTEGER(KIND=I4KIND), PARAMETER :: KO2    = 12
    INTEGER(KIND=I4KIND), PARAMETER :: KO21   = KO2   + 1
    INTEGER(KIND=I4KIND), PARAMETER :: KO2M   = KO2   - 1
!--------------------------------------------------------------------------------------------------
! PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
! IMAX   =  NO. POINTS SENT TO RADFS
! L      =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL 
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
! NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE BANDTA FOR DEFINITION                    
! NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS
! NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE BDCOMB FOR DEFINITION               
! INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
! NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS. NB,KO2 ARE SHORTWAVE PARAMETERS;
!           OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS.
!
! COMMON BLOCK TABCOM CONTAINS QUANTITIES PRECOMPUTED IN SUBROUTINE
! TABLE FOR USE IN THE LONGWAVE RADIATION PROGRAM:
! EM1     =  E1 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS                    
! EM1WDE  =  E1 FUNCTION, EVALUATED OVER THE 160-560 CM-1 INTERVAL               
! TABLE1  =  E2 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS            
! TABLE2  =  TEMPERATURE DERIVATIVE OF TABLE1
! TABLE3  =  MASS DERIVATIVE OF TABLE1 
! EM3     =  E3 FUNCTION, EVALUATED OVER THE 0-560 AND 1200-2200 CM-1 INTERVALS               
! SOURCE  =  PLANCK FUNCTION, EVALUATED AT SPECIFIED TEMPS. FOR BANDS USED IN CTS CALCULATIONS               
! DSRCE   =  TEMPERATURE DERIVATIVE OF SOURCE
! IND     =  INDEX, WITH VALUE IND(I)=I. USED IN FST88
! INDX2   =  INDEX VALUES USED IN OBTAINING "LOWER TRIANGLE" ELEMENTS OF AVEPHI,ETC., IN FST88             
! KMAXV   =  INDEX VALUES USED IN OBTAINING "UPPER TRIANGLE" ELEMENTS OF AVEPHI,ETC., IN FST88                  
! KMAXVM  =  KMAXV(L),USED FOR DO LOOP INDICES
!--------------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                                          ::&
    & IT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                                          ::&  
    & WW1     , DT      ,                                                                         &
    & WW2     , DU
!------------------------------------------------------ 
! THE FOLLOWING ARRAYS ARE EQUIVALENCED TO VTEMP ARRAYS
!------------------------------------------------------ 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                                          ::&
    & FXO     , FYO     , TMP3
!-------------------------------------- 
! DIMENSIONS OF ARRAYS IN ARGUMENT LIST
!--------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                    , INTENT(INOUT)       ::& 
    & EMV
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP1)                    , INTENT(IN)          ::& 
    & TV      , AV
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K 
!
    EQUIVALENCE (FXO,WW1), (FYO,WW2), (IT,TMP3)
!-------------------------------------------------------------------------
! THE FOLLOWING LOOP REPLACES A DOUBLE LOOP OVER I (1-IMAX) AND K (1-LLP1)
!-------------------------------------------------------------------------
    DO 203 K=1,LLP1
        DO 203 I=MYIS,MYIE
            FXO (I,K) = AINT(TV(I,K)    * HP1)
            TMP3(I,K) = LOG10(AV(I,K))  + H16E1
            DT  (I,K) = TV(I,K)   - TEN * FXO(I,K)
            FYO (I,K) = AINT(TMP3(I,K)  * TEN)
            DU  (I,K) = TMP3(I,K) - HP1 * FYO(I,K)
!--------------------------------------------------------------------------------------------------
! OBTAIN INDEX FOR TABLE LOOKUP; THIS VALUE WILL HAVE TO BE DECREMENTED BY 9 TO ACCOUNT FOR TABLE 
! TEMPS STARTING AT 100K.
!--------------------------------------------------------------------------------------------------
            IT (I,K) = FXO(I,K) + FYO(I,K) * H28E1
            WW1(I,K) = TEN - DT(I,K)
            WW2(I,K) = HP1 - DU(I,K)
            EMV(I,K) = WW1(I,K) * WW2(I,K) * EM3V(IT(I,K) -  9)                                   &
            &        + WW2(I,K) *  DT(I,K) * EM3V(IT(I,K) -  8)                                   &
            &        + WW1(I,K) *  DU(I,K) * EM3V(IT(I,K) + 19)                                   &
            &        +  DT(I,K) *  DU(I,K) * EM3V(IT(I,K) + 20)
203 END DO
!
    RETURN
!
    END SUBROUTINE E3V88
