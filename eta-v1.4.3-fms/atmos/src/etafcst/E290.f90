    SUBROUTINE E290(EMISSB, EMISS, AVEPHI, KLEN, FXOE2, DTE2)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE E290
!>
!> SUBPROGRAM: E290 - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> SUBROUTINE E290 COMPUTES THE EXCHANGE TERMS IN THE FLUX EQUATION FOR LONGWAVE RADIATION FOR ALL 
!> TERMS EXCEPT THE EXCHANGE WITH THE TOP OF THE ATMOSPHERE. THE METHOD IS A TABLE LOOKUP ON A PRE-
!> COMPUTED E2 FUNCTION (DEFINED IN REF. (4)).
!> CALCULATIONS ARE DONE IN THE FREQUENCY RANGE:
!> 1) 0-560, 1200-2200 CM-1   FOR Q(APPROX)
!> MOTIVATION FOR THESE CALCULATIONS IS IN REFERENCES (1) AND (4).
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
!> AVEPHI -
!> KLEN   - 
!> FXOE2  - 
!> DTE2   -  
!>
!> OUTPUT ARGUMENT LIST:
!> EMISSB - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> EMISS  -
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
!--------------------------------------------------------------------------------------------------
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
! COMMON BLOCK TABCOM CONTAINS QUANTITIES PRECOMPUTED IN SUBROUTINE
! TABLE FOR USE IN THE LONGWAVE RADIATION PROGRAM:
! EM1     =  E1 FUNCTION, EVALUATED OVER THE 0-560 AND  1200-2200 CM-1 INTERVALS     
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
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(INOUT)       ::&
    & EMISSB     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(INOUT)       ::&
    & EMISS
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & AVEPHI
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & IVAL    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & DT      , FYO     , DU
!---------------------------------------- 
! TMP3 MAY BE EQUIVALENCED TO DT IN VTEMP
!---------------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & TMP3
!------------------------------- 
! VARIABLES IN THE ARGUMENT LIST
!------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & FXOE2   , DTE2
!
    EQUIVALENCE (TMP3,DT)
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & KLEN 
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K     		 
!--------------------------------------------------------------------------------------------------
! FIRST WE OBTAIN THE EMISSIVITIES AS A FUNCTION OF TEMPERATURE (INDEX FXO) AND WATER AMOUNT 
! (INDEX FYO). THIS PART OF THE CODE THUS GENERATES THE E2 FUNCTION.
!
! CALCULATIONS FOR VARYING KP (FROM KP=K+1 TO LP1, INCLUDING SPECIAL CASE: RESULTS ARE IN EMISS
!--------------------------------------------------------------------------------------------------
    DO 132 K=1,LP2-KLEN
        DO 132 I=MYIS,MYIE
            TMP3(I,K) = LOG10(AVEPHI(I,KLEN+K-1)) + H16E1
            FYO (I,K) = AINT(TMP3(I,K)  * TEN)
            DU  (I,K) = TMP3(I,K) - HP1 * FYO(I,K)
            FYO (I,K) = H28E1           * FYO(I,K)
            IVAL(I,K) =  FYO(I,K) +       FXOE2(I,KLEN+K-1)
            EMISS(I,KLEN+K-1) = T1(IVAL(I,K)) + DU(I,K)                                           &
    &                         * T2(IVAL(I,K)) + DTE2(I,KLEN+K-1)                                  &
    &                         * T4(IVAL(I,K))
132 END DO
!---------------------------------------------------------------------------------------------- 
! THE SPECIAL CASE EMISS(I,L) (LAYER KP) IS OBTAINED NOW BY AVERAGING THE VALUES FOR L AND LP1:
!---------------------------------------------------------------------------------------------- 
    DO 134 I=MYIS,MYIE
        EMISS(I,L) = HAF * (EMISS(I,L) + EMISS(I,LP1))
134 END DO
!------------------------------------------------------- 
! NOTE THAT EMISS(I,LP1) IS NOT USEFUL AFTER THIS POINT.
!------------------------------------------------------- 
!
!--------------------------------------------------------------------------------------------------
! CALCULATIONS FOR KP=KLEN AND VARYING K; RESULTS ARE IN EMISSB.
! IN THIS CASE, THE TEMPERATURE INDEX IS UNCHANGED, ALWAYS BEING FXO(I,KLEN-1); THE WATER INDEX 
! CHANGES, BUT IS SYMMETRICAL WITH THAT FOR THE VARYING KP CASE.
! NOTE THAT THE SPECIAL CASE IS NOT INVOLVED HERE. (FIXED LEVEL) K VARIES FROM (KLEN+1) TO LP1; 
! RESULTS ARE IN EMISSB(I,(KLEN) TO L)
!--------------------------------------------------------------------------------------------------
    DO 142 K=1,LP1-KLEN
        DO 142 I=MYIS,MYIE
              DT(I,K) = DTE2(I,KLEN-1)
            IVAL(I,K) = FYO(I,K) + FXOE2(I,KLEN-1)
142 END DO
!
    DO 234 K=1,LP1-KLEN
        DO 234 I=MYIS,MYIE
            EMISSB(I,KLEN+K-1) = T1(IVAL(I,K)) + DU(I,K)                                          &
    &                          * T2(IVAL(I,K)) + DT(I,K)                                          &
    &                          * T4(IVAL(I,K))
234 END DO
!
    RETURN
!
    END SUBROUTINE E290
