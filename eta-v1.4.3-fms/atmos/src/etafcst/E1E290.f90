    SUBROUTINE E1E290(G1, G2, G3, G4, G5, EMISS, FXOE1, DTE1, FXOE2, DTE2, AVEPHI, TEMP, T)
!--------------------------------------------------------------------------------------------------
!> SUBROUTINE E1E290
!>
!> SUBPROGRAM: E1E290 - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> SUBROUTINE E1E290 COMPUTES THE EXCHANGE TERMS IN THE FLUX EQUATION FOR LONGWAVE RADIATION FOR ALL
!> TERMS EXCEPT THE EXCHANGE WITH THE TOP OF THE ATMOSPHERE.
!> THE METHOD IS A TABLE LOOKUP ON A PRE-COMPUTED E2 FUNCTION (DEFINED IN REF. (4)).
!> THE E1 FUNCTION  CALCULATIONS (FORMERLY DONE IN SUBROUTINE E1V88 COMPUTE THE FLUX RESULTING FROM
!> THE EXCHANGE OF PHOTONS BETWEEN A LAYER AND THE TOP OF THE ATMOSPHERE.
!> THE METHOD IS A TABLE LOOKUP ON A PRE-COMPUTED E1 FUNCTION.
!> CALCULATIONS ARE DONE IN TWO FREQUENCY RANGES:
!> 1)   0-560, 1200-2200 CM-1   FOR Q(APPROX)
!> 2) 160-560            CM-1   FOR Q(APPROX,CTS).
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
!> FXOE1  - 
!> FXOE2  - 
!> DTE1   - 
!> DTE2   - 
!> AVEPHI - 
!> TEMP   - 
!> T      - 
!>
!> OUTPUT ARGUMENT LIST:
!> G1     -
!> G2     - 
!> G3     - 
!> G4     - 
!> G5     - 
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
!
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
!
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
! 
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
!
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
! KMAXV   =  INDEX VALUES USED IN OBTAINING "UPPER TRIANGLE" ELEMENTS OF AVEPHI,ETC.,IN FST88
! KMAXVM  =  KMAXV(L),USED FOR DO LOOP INDICES
!--------------------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & TEMP    , T       , AVEPHI  
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LL3P)                                          ::&
    & IT1 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & FYO     ,                                                                                   &
    & DU      , WW1     ,  WW2    ,                                                               &
    & TMP3    
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(INOUT)       ::&
    & EMISS
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & TMP5    , TMP9
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & IVAL
!------------------------------- 
! VARIABLES IN THE ARGUMENT LIST
!------------------------------- 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & FXOE1   , DTE1    ,                                                                         &
    & FXOE2   , DTE2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(INOUT)       ::&
    & G1      , G3      , G4
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                       , INTENT(OUT)         ::&
    & G2      , G5
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K       , KP 
!--------------------------------------------------------------------------------------------------
! FIRST WE OBTAIN THE EMISSIVITIES AS A FUNCTION OF TEMPERATURE (INDEX FXO) AND WATER AMOUNT 
! (INDEX FYO). 
! THIS PART OF THE CODE THUS GENERATES THE E2 FUNCTION. THE FXO INDICES HAVE BEEN OBTAINED IN FST88
! FOR CONVENIENCE.
!
! THIS SUBROUTINE EVALUATES THE K=1 CASE ONLY
!
! THIS LOOP REPLACES LOOPS GOING FROM I=1,IMAX AND KP=2,LP1 
! PLUS THE SPECIAL CASE FOR THE LP1TH LAYER
!--------------------------------------------------------------------------------------------------
    DO 132 K=1,LP1
        DO 132 I=MYIS,MYIE
            TMP3 (I,K) = LOG10(AVEPHI(I,K)) + H16E1
            FYO  (I,K) = AINT(TMP3(I,K) * TEN)
            DU   (I,K) = TMP3(I,K) - HP1 * FYO(I,K)
            FYO  (I,K) = H28E1 * FYO(I,K)
            IVAL (I,K) = FYO(I,K) + FXOE2(I,K)
            EMISS(I,K) = T1(IVAL(I,K)) + DU(I,K) * T2(IVAL(I,K)) + DTE2(I,K) * T4(IVAL(I,K))
132 END DO
!--------------------------------------------------------------------------------------------------
! THE SPECIAL CASE EMISS(I,L) (LAYER KP) IS OBTAINED NOW BY AVERAGING THE VALUES FOR L AND LP1:
!--------------------------------------------------------------------------------------------------
    DO 134 I=MYIS,MYIE
        EMISS(I,L) = HAF * (EMISS(I,L) + EMISS(I,LP1))
134 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATIONS FOR THE KP=1 LAYER ARE NOT PERFORMED, AS THE RADIATION CODE ASSUMES THAT THE TOP 
! FLUX LAYER (ABOVE THE TOP DATA LEVEL) IS ISOTHERMAL, AND HENCE CONTRIBUTES NOTHING TO THE FLUXES 
! AT OTHER LEVELS.
!
! THE FOLLOWING IS THE CALCULATION FOR THE E1 FUNCTION, FORMERLY DONE IN SUBROUTINE E1V88. 
! THE MOVE TO E1E288 IS DUE TO THE SAVINGS IN OBTAINING INDEX VALUES (THE TEMP. INDICES HAVE BEEN 
! OBTAINED IN FST88, WHILE THE U-INDICES ARE OBTAINED IN THE E2 CALCS., WITH K=1).
!
! FOR TERMS INVOLVING TOP LAYER, DU IS NOT KNOWN; IN FACT, WE USE INDEX 2 TO REPERSENT INDEX 1 IN 
! PREV. CODE. THIS MEANS THAT THE IT1 INDEX 1 AND LLP1 HAS TO BE CALCULATED SEPARATELY. 
! THE INDEX LLP2 GIVES THE SAME VALUE AS 1; IT CAN BE OMITTED.
!--------------------------------------------------------------------------------------------------
    DO 208 I=MYIS,MYIE
        IT1(I,1) = FXOE1(I,1)
        WW1(I,1) = TEN - DTE1(I,1)
        WW2(I,1) = HP1
208 END DO
!
    DO 209 K=1,L
        DO 209 I=MYIS,MYIE
            IT1(I,K+1)     = FYO(I,K) + FXOE1(I,K+1)
            IT1(I,LP2+K-1) = FYO(I,K) + FXOE1(I,K)
            WW1(I,K+1)     = TEN - DTE1(I,K+1)
            WW2(I,K+1)     = HP1 - DU(I,K)
209 END DO
!
    DO 211 KP=1,L
        DO 211 I=MYIS,MYIE
            IT1(I,KP+LLP1) = FYO(I,KP) + FXOE1(I,1)
211 END DO
!------------------------------------------------------ 
! G3(I,1) HAS THE SAME VALUES AS G1 (AND DID ALL ALONG)
!------------------------------------------------------ 
    DO 230 I=MYIS,MYIE
        G1(I,1) = WW1(I,1) * WW2(I,1) * EM1V(IT1(I,1)) + WW2(I,1) * DTE1(I,1) * EM1V(IT1(I,1)+1)
        G3(I,1) =  G1(I,1)
230 END DO
!
    DO 240 K=1,L
        DO 240 I=MYIS,MYIE
            G1(I,K+1)  =  WW1(I,K+1) *  WW2(I,K+1)  *  EM1V(IT1(I,K+1))                           &
    &                  +  WW2(I,K+1) * DTE1(I,K+1)  *  EM1V(IT1(I,K+1) +  1)                      &
    &                  +  WW1(I,K+1) *   DU(I,K)    *  EM1V(IT1(I,K+1) + 28)                      &
    &                  + DTE1(I,K+1) *   DU(I,K)    *  EM1V(IT1(I,K+1) + 29)
!
            G2(I,K)   =   WW1(I,K)   *  WW2(I,K+1)  *  EM1V(IT1(I,K+LP2-1))                       &
    &                 +   WW2(I,K+1) * DTE1(I,K)    *  EM1V(IT1(I,K+LP2-1) +  1)                  &
    &                 +   WW1(I,K)   *   DU(I,K)    *  EM1V(IT1(I,K+LP2-1) + 28)                  &
    &                 +  DTE1(I,K)   *   DU(I,K)    *  EM1V(IT1(I,K+LP2-1) + 29)
240 END DO
!
    DO 241 KP=2,LP1
        DO 241 I=MYIS,MYIE
            G3(I,KP)  =   WW1(I,1)   *  WW2(I,KP)   *  EM1V(IT1(I,LL+KP))                         &
    &                 +   WW2(I,KP)  * DTE1(I,1)    *  EM1V(IT1(I,LL+KP) +  1)                    &
    &                 +   WW1(I,1)   *   DU(I,KP-1) *  EM1V(IT1(I,LL+KP) + 28)                    &
    &                 +  DTE1(I,1)   *   DU(I,KP-1) *  EM1V(IT1(I,LL+KP) + 29)
241 END DO
!
    DO 244 I=MYIS,MYIE
        G4(I,1) = WW1(I,1) * WW2(I,1) * EM1VW(IT1(I,1)) + WW2(I,1) * DTE1(I,1) * EM1VW(IT1(I,1)+1)
244 END DO
!
    DO 242 K=1,L
        DO 242 I=MYIS,MYIE
            G4(I,K+1) =   WW1(I,K+1) *  WW2(I,K+1)  *  EM1VW(IT1(I,K+1))                          &
    &                 +   WW2(I,K+1) * DTE1(I,K+1)  *  EM1VW(IT1(I,K+1) +  1)                     &
    &                 +   WW1(I,K+1) *   DU(I,K)    *  EM1VW(IT1(I,K+1) + 28)                     &
    &                 +  DTE1(I,K+1) *   DU(I,K)    *  EM1VW(IT1(I,K+1) + 29)
!
            G5(I,K)   =   WW1(I,K)   *  WW2(I,K+1)  *  EM1VW(IT1(I,K+LP2-1))                      &
    &                 +   WW2(I,K+1) * DTE1(I,K)    *  EM1VW(IT1(I,K+LP2-1) +  1)                 &
    &                 +   WW1(I,K)   *   DU(I,K)    *  EM1VW(IT1(I,K+LP2-1) + 28)                 &
    &                 +  DTE1(I,K)   *   DU(I,K)    *  EM1VW(IT1(I,K+LP2-1) + 29)
242 END DO
!
    RETURN
!
    END SUBROUTINE E1E290
