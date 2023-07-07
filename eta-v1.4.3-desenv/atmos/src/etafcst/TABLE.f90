    SUBROUTINE TABLE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TABLE
!>
!> SUBPROGRAM: TABLE - COMPUTES TABLE ENTRIES USED IN THE LONGWAVE RADIA PROGRAM
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT:
!> SUBROUTINE TABLE COMPUTES TABLE ENTRIES USED IN THE LONGWAVE RADIA PROGRAM. ALSO CALCULATED ARE
!> INDICES USED IN STRIP-MINING AND FOR SOME PRE-COMPUTABLE FUNCTIONS.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????       - ORIGINATOR
!> 18-01-15  LUCCI       - MODERNIZATION OF THE CODE, INCLUDING:
!>                         * F77 TO F90/F95
!>                         * INDENTATION & UNIFORMIZATION CODE
!>                         * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                         * DOCUMENTATION WITH DOXYGEN
!>                         * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              HCON
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              PHYCON
!>              RNDDTA
!>              SCRTCH
!>              TABCOM
!>              TBLTMP
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : CONRAD
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
    USE RNDDTA
    USE SCRTCH
    USE TABCOM
    USE TBLTMP
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
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
! NB, KO2 ARE SHORTWAVE PARAMETERS; OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS.
!------------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: L      = LM
    INTEGER(KIND=I4KIND), PARAMETER :: NCOL   = IMAX
    INTEGER(KIND=I4KIND), PARAMETER :: NBLM   = NBLY  -  1
    INTEGER(KIND=I4KIND), PARAMETER :: LP2    = L     +  2 
    INTEGER(KIND=I4KIND), PARAMETER :: LP3    = L     +  3
    INTEGER(KIND=I4KIND), PARAMETER :: LM1    = L     -  1
    INTEGER(KIND=I4KIND), PARAMETER :: LM2    = L     -  2
    INTEGER(KIND=I4KIND), PARAMETER :: LM3    = L     -  3
    INTEGER(KIND=I4KIND), PARAMETER :: LL     =  2    *  L
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1   = LL    +  1 
    INTEGER(KIND=I4KIND), PARAMETER :: LLP2   = LL    +  2
    INTEGER(KIND=I4KIND), PARAMETER :: LLP3   = LL    +  3
    INTEGER(KIND=I4KIND), PARAMETER :: LLM1   = LL    -  1
    INTEGER(KIND=I4KIND), PARAMETER :: LLM2   = LL    -  2
    INTEGER(KIND=I4KIND), PARAMETER :: LLM3   = LL    -  3
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M   = LP1   *  LP1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1M1  = LP1M  -  1 
    INTEGER(KIND=I4KIND), PARAMETER :: LP121  = LP1   * NBLY
    INTEGER(KIND=I4KIND), PARAMETER :: LL3P   =  3    *  L + 2
    INTEGER(KIND=I4KIND), PARAMETER :: NB     = 12
    INTEGER(KIND=I4KIND), PARAMETER :: INLTE  =  3 
    INTEGER(KIND=I4KIND), PARAMETER :: INLTEP = INLTE + 1 
    INTEGER(KIND=I4KIND), PARAMETER :: NNLTE  = 56
    INTEGER(KIND=I4KIND), PARAMETER :: LP1I   = IMAX  * LP1 
    INTEGER(KIND=I4KIND), PARAMETER :: LLP1I  = IMAX  * LLP1
    INTEGER(KIND=I4KIND), PARAMETER :: LL3PI  = IMAX  * LL3P
    INTEGER(KIND=I4KIND), PARAMETER :: NB1    = NB    - 1
    INTEGER(KIND=I4KIND), PARAMETER :: KO2    = 12
    INTEGER(KIND=I4KIND), PARAMETER :: KO21   = KO2   + 1 
    INTEGER(KIND=I4KIND), PARAMETER :: KO2M   = KO2   - 1
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , N       , ICNT    , I1      , I2      , I2E     , J       , JP      , IA      ,   &
    & NSB     , NSUBDS         
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & CENT    , DEL     , BDLO    , BDHI    , ANU     , C1
!--------------------------------------------------------------
! COMPUTE LOCAL QUANTITIES AND AO3,BO3,AB15 FOR NARROW-BANDS...
!--------------------------------------------------------------
    DO 101 N=1,NBLW
           ANB(N) = ARNDM(N)
           BNB(N) = BRNDM(N)
        CENTNB(N) = HAF       * (BANDLO(N) + BANDHI(N))
         DELNB(N) = BANDHI(N) -  BANDLO(N)
        BETANB(N) = BETAD(N)
101 END DO
!
    AB15(1) = ANB(57) * BNB(57)
    AB15(2) = ANB(58) * BNB(58)
!------------------  
! FOR WIDE BANDS...
!------------------ 
    AB15WD = AWIDE * BWIDE
!--------------------------------- 
! COMPUTE INDICES: IND,INDX2,KMAXV
!--------------------------------- 
    DO 111 I=1,IMAX
        IND(I) = I
111 END DO
!
    ICNT = 0
!
    DO 113 I1=1,L
        I2E = LP1 - I1
        DO 115 I2=1,I2E
            ICNT        = ICNT + 1
            INDX2(ICNT) = LP1 * (I2 - 1) + LP2 * I1
    115 END DO
113 END DO
!
    KMAXV(1) = 1
!
    DO 117 I=2,L
        KMAXV(I) = KMAXV(I-1) + (LP2-I)
117 END DO
!
    KMAXVM = KMAXV(L)
!------------------------------- 
! COMPUTE RATIOS OF CONT. COEFFS
!--------------------------------- 
    SKC1R = BETAWD    / BETINW
    SKO3R = BETAD(61) / BETINW
    SKO2D = ONE       / BETINW
!-------------------------------------------------------------------------------------------------- 
! BEGIN TABLE COMPUTATIONS HERE 
! COMPUTE TEMPS, MASSES FOR TABLE ENTRIES
! NOTE: THE DIMENSIONING AND INITIALIZATION OF XTEMV AND OTHER ARRAYS WITH DIMENSION OF 28 IMPLY
! A RESTRICTION OF MODEL TEMPERATURES FROM 100K TO 370K.
! THE DIMENSIONING OF ZMASS,ZROOT AND OTHER ARRAYS WITH DIMENSION OF 180 IMPLY A RESTRICTION OF
! MODEL H2O AMOUNTS SUCH THAT OPTICAL PATHS ARE BETWEEN 10**-16 AND 10**2, IN CGS UNITS.
!--------------------------------------------------------------------------------------------------
    ZMASS(1) = H1M16
!
    DO 201 J=1,180
        JP = J+1
        ZROOT(J) = SQRT(ZMASS(J))
        ZMASS(JP) = ZMASS(J) * H1P25892
201 END DO
!
    DO 203 I=1,28
         XTEMV(I) = HNINETY + TEN      * I
         TFOUR(I) =           XTEMV(I) * XTEMV(I) * XTEMV(I) * XTEMV(I)
        FORTCU(I) = FOUR    * XTEMV(I) * XTEMV(I) * XTEMV(I)
203 END DO
!------------------------------------------------------------------------------------------------------------------------------- 
! THE COMPUTATION OF SOURCE,DSRCE IS NEEDED ONLY FOR THE COMBINED WIDE-BAND CASE.TO OBTAIN THEM, THE SOURCE MUST BE COMPUTED FOR 
! EACH OF THE (NBLX) WIDE BANDS(=SRCWD) THEN COMBINED (USING IBAND) INTO SOURCE.
!-------------------------------------------------------------------------------------------------------------------------------
    DO 205 N=1,NBLY
        DO 205 I=1,28
            SOURCE(I,N) = ZERO
205 END DO
!
    DO 207 N=1,NBLX
        DO 207 I=1,28
             SRCWD(I,N) = ZERO
207 END DO
!------------------------ 
! BEGIN FREQ. LOOP (ON N)
!------------------------
    DO 211 N=1,NBLX
        IF (N <= 46) THEN
!------------------------
! THE 160-1200 BAND CASES
!------------------------
            CENT = CENTNB(N+16)
            DEL  =  DELNB(N+16)
            BDLO = BANDLO(N+16)
            BDHI = BANDHI(N+16)
        END IF
!
        IF (N == NBLX) THEN
!------------------------
! THE 2270-2380 BAND CASE
!------------------------
            CENT = CENTNB(NBLW)
            DEL  =  DELNB(NBLW)
            BDLO = BANDLO(NBLW)
            BDHI = BANDHI(NBLW)
        END IF
!---------------------------------------------------------------------------------------------------------------------------------
! FOR PURPOSES OF ACCURACY, ALL EVALUATIONS OF PLANCK FCTNS ARE MADE ON 10 CM-1 INTERVALS, THEN SUMMED INTO THE (NBLX) WIDE BANDS.
!---------------------------------------------------------------------------------------------------------------------------------
        NSUBDS = (DEL - H1M3) / 10 + 1
!
        DO 213 NSB=1,NSUBDS
            IF (NSB /= NSUBDS) THEN
                CNUSB(NSB) = TEN * (NSB-1) + BDLO + FIVE
                DNUSB(NSB) = TEN
            ELSE
                CNUSB(NSB) = HAF * (TEN*(NSB-1) + BDLO + BDHI)
                DNUSB(NSB) = BDHI - (TEN * (NSB-1) + BDLO)
            END IF
!
            C1 = (H37412M5) * CNUSB(NSB) ** 3
!------------------------ 
! BEGIN TEMP. LOOP (ON I)
!------------------------
            DO 215 I=1,28
                    X(I) = H1P4387 * CNUSB(NSB) / XTEMV(I)
                   X1(I) = EXP(X(I))
                 SRCS(I) = C1 / (X1(I) - ONE)
                SRCWD(I,N) = SRCWD(I,N) + SRCS(I) * DNUSB(NSB)
        215 END DO
    213 END DO
!
211 END DO
!------------------------------------------------------------------------------
! THE FOLLOWING LOOPS CREATE THE COMBINED WIDE BAND QUANTITIES SOURCE AND DSRCE
!------------------------------------------------------------------------------
    DO 221 N=1,40
        DO 221 I=1,28
            SOURCE(I,IBAND(N)) = SOURCE(I,IBAND(N)) + SRCWD(I,N)
221 END DO
!
    DO 223 N=9,NBLY
        DO 223 I=1,28
            SOURCE(I,N) = SRCWD(I,N+32)
223 END DO
!
    DO 225 N=1,NBLY
        DO 225 I=1,27
            DSRCE(I,N) = (SOURCE(I+1,N) - SOURCE(I,N)) * HP1
225 END DO
!
    DO 231 N=1,NBLW
        ALFANB(N) = BNB(N) * ANB(N)
        AROTNB(N) =  SQRT(ALFANB(N))
231 END DO
!--------------------------------------------------------------------------------------------------
! FIRST COMPUTE PLANCK FCTNS (SRC1NB) AND DERIVATIVES (DBDTNB) FOR USE IN TABLE EVALUATIONS. THESE
! ARE DIFFERENT FROM SOURCE,DSRCE BECAUSE DIFFERENT FREQUENCY PTS ARE USED IN EVALUATION, 
! THE FREQ. RANGES ARE DIFFERENT, AND THE DERIVATIVE ALGORITHM IS DIFFERENT.
!--------------------------------------------------------------------------------------------------
    DO 301 N=1,NBLW
        CENT = CENTNB(N)
         DEL =  DELNB(N)
!--------------------------------------------------------------------------------------------------
! NOTE: AT PRESENT, THE IA LOOP IS ONLY USED FOR IA=2. THE LOOP STRUCT IS KEPT SO THAT IN THE 
! FUTURE, WE MAY USE A QUADRATURE SCHEME FOR THE PLANCK FCTN EVALUATION, RATHER THAN USE THE 
! MID-BAND FREQUENCY.
!--------------------------------------------------------------------------------------------------
        DO 303 IA=1,3
            ANU = CENT + HAF * (IA-2) * DEL
             C1 = (H37412M5) * ANU * ANU * ANU + H1M20
!----------------- 
! TEMPERATURE LOOP
!-----------------
            DO 305 I=1,28
                  X(I) = H1P4387 * ANU / XTEMV(I)
                 X1(I) = EXP(X(I))
                 SC(I) = C1 / ((X1(I) - ONE) + H1M20)
                DSC(I) = SC(I) * SC(I) * X(I) * X1(I) / (XTEMV(I) *C1)
        305 END DO
!
            IF (IA == 2) THEN
                DO 307 I=1,28
                    SRC1NB(I,N) = DEL * SC(I)
                    DBDTNB(I,N) = DEL * DSC(I)
            307 END DO
            END IF
!
    303 END DO
!
301 END DO
!--------------------------------------------------------------------------------------------------
! NEXT COMPUTE R1,R2,S2,AND T3-COEFFICIENTS USED FOR E3 FUNCTION WHEN THE OPTICAL PATH IS LESS
! THAN 10-4. IN THIS CASE, WE ASSUME A DIFFERENT DEPENDENCE ON (ZMASS). ALSO OBTAIN R1WD, 
! WHICH IS R1 SUMMED OVER THE 160-560 CM-1 RANGE
!--------------------------------------------------------------------------------------------------
    DO 311 I=1,28
          SUM4(I) = ZERO
          SUM6(I) = ZERO
          SUM7(I) = ZERO
          SUM8(I) = ZERO
        SUM4WD(I) = ZERO
311 END DO
!
    DO 313 N=1,NBLW
        CENT = CENTNB(N)
!------------------------------------------------------------------------------------- 
! PERFORM SUMMATIONS FOR FREQ. RANGES OF 0-560,1200-2200 CM-1 FOR SUM4, SUM6,SUM7,SUM8
!------------------------------------------------------------------------------------- 
        IF (CENT < 560. .OR. CENT > 1200. .AND. CENT <= 2200.) THEN
            DO 315 I=1,28
                SUM4(I) = SUM4(I) + SRC1NB(I,N)
                SUM6(I) = SUM6(I) + DBDTNB(I,N)
                SUM7(I) = SUM7(I) + DBDTNB(I,N) * AROTNB(N)
                SUM8(I) = SUM8(I) + DBDTNB(I,N) * ALFANB(N)
        315 END DO
        END IF
!----------------------------------------------------------------------  
! PERFORM SUMMATIONS OVER 160-560 CM-1 FREQ RANGE FOR E1 CALCS (SUM4WD)
!----------------------------------------------------------------------
        IF (CENT > 160. .AND. CENT < 560.) THEN
            DO 316 I=1,28
                SUM4WD(I) = SUM4WD(I) + SRC1NB(I,N)
        316 END DO
        END IF
!
313 END DO
!
    DO 317 I=1,28
          R1(I) =   SUM4(I) /  TFOUR(I)
          R2(I) =   SUM6(I) / FORTCU(I)
          S2(I) =   SUM7(I) / FORTCU(I)
          T3(I) =   SUM8(I) / FORTCU(I)
        R1WD(I) = SUM4WD(I) /  TFOUR(I)
317 END DO
!
    DO 401 J=1,180
        DO 401 I=1,28
               SUM(I,J) = ZERO
            PERTSM(I,J) = ZERO
              SUM3(I,J) = ZERO
            SUMWDE(I,J) = ZERO
401 END DO
!---------------------- 
! FREQUENCY LOOP BEGINS 
!----------------------
    DO 411 N=1,NBLW
        CENT = CENTNB(N)
!--------------------------------------------------------------
! PERFORM CALCULATIONS FOR FREQ. RANGES OF 0-560,1200-2200 CM-1
!--------------------------------------------------------------
        IF (CENT < 560. .OR. CENT > 1200. .AND. CENT <= 2200.) THEN
            DO 413 J=1,180
                  X2(J) = AROTNB(N) * ZROOT(J)
                EXPO(J) = EXP(-X2(J))
        413 END DO
!
            DO 415 J=1,180
                IF (X2(J) >= HUNDRED) THEN
                    EXPO(J) = ZERO
                END IF
        415 END DO
!
            DO 417 J=121,180
                FAC(J) = ZMASS(J) * (ONE - (ONE + X2(J)) * EXPO(J)) / (X2(J) * X2(J))
        417 END DO
!
            DO 419 J=1,180
                DO 419 I=1,28
                       SUM(I,J) =    SUM(I,J) + SRC1NB(I,N) * EXPO(J)
                    PERTSM(I,J) = PERTSM(I,J) + DBDTNB(I,N) * EXPO(J)
        419 END DO
!
            DO 421 J=121,180
                DO 421 I=1,28
                    SUM3(I,J) = SUM3(I,J) + DBDTNB(I,N) * FAC(J)
        421 END DO
        END IF
!-----------------------------------------------------------------
! COMPUTE SUM OVER 160-560 CM-1 RANGE FOR USE IN E1 CALCS (SUMWDE)
!-----------------------------------------------------------------
        IF (CENT > 160. .AND. CENT < 560.) THEN
            DO 420 J=1,180
                DO 420 I=1,28
                    SUMWDE(I,J) = SUMWDE(I,J) + SRC1NB(I,N) * EXPO(J)
        420 END DO
        END IF
!
411 END DO
!
    DO 431 J=1,180
        DO 431 I=1,28
               EM1(I,J) =    SUM(I,J) /  TFOUR(I)
            TABLE1(I,J) = PERTSM(I,J) / FORTCU(I)
431 END DO
!
    DO 433 J=121,180
        DO 433 I=1,28
            EM3(I,J) = SUM3(I,J) / FORTCU(I)
433 END DO
!
    DO 441 J=1,179
        DO 441 I=1,28
            TABLE2(I,J) = (TABLE1(I,J+1) - TABLE1(I,J)) * TEN
441 END DO
!
    DO 443 J=1,180
        DO 443 I=1,27
            TABLE3(I,J) = (TABLE1(I+1,J) - TABLE1(I,J)) * HP1
443 END DO
!
    DO 445 I=1,28
        TABLE2(I,180) = ZERO
445 END DO
!
    DO 447 J=1,180
        TABLE3(28,J)  = ZERO
447 END DO
!
    DO 449 J=1,2
        DO 449 I=1,28
            EM1(I,J)  = R1(I)
449 END DO
!
    DO 451 J=1,120
        DO 451 I=1,28
            EM3(I,J)  = R2(I) / TWO-S2(I) * SQRT(ZMASS(J)) / THREE + T3(I) * ZMASS(J) / EIGHT
451 END DO
!
    DO 453 J=121,180
        DO 453 I=1,28
            EM3(I,J)  = EM3(I,J) / ZMASS(J)
453 END DO
!------------------------------------------------------------------------------------------
! NOW COMPUTE E1 TABLES FOR 160-560 CM-1 BANDS ONLY. WE USE R1WD AND SUMWDE OBTAINED ABOVE.
!------------------------------------------------------------------------------------------
    DO 501 J=1,180
        DO 501 I=1,28
            EM1WDE(I,J) = SUMWDE(I,J) / TFOUR(I)
501 END DO
!
    DO 503 J=1,2
        DO 503 I=1,28
            EM1WDE(I,J) = R1WD(I)
503 END DO
!
    RETURN
!
    END SUBROUTINE TABLE
