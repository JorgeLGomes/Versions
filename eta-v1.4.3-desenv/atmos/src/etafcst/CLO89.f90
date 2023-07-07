    SUBROUTINE CLO89(CLDFAC, CAMT, NCLDS, KBTM, KTOP)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE CLO89
!>
!> SUBPROGRAM: CLO89 - ?????
!> PROGRAMMER: Q. ZHAO
!> ORG: W/NP2
!> DATE: 95-03-22
!>
!> ABSTRACT: 
!> COMPUTES CLOUD TRANSMISSION FUNCTIONS FOR THE LONGWAVE CODE,USING CODE WRITTEN BY BERT
!> KATZ (301-763-8161) AND MODIFIED BY DAN SCHWARZKOPF IN DECEMBER, 1988.
!>
!> PROGRAM HISTORY LOG:
!> 95-03-22  Q. ZHAO  - ORIGINATOR            
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT  ARGUMENT LIST:
!> CAMT   -
!> KTOP   -
!> KBTM   -
!> NCLDS  -
!>
!> OUTPUT ARGUMENT LIST:
!> CLDFAC -
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
!>              TEMPCOM
!>              TOPO
!> 
!> DRIVER     : RADFS
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
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
!----------------------------------------------------------------------------------------
! PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
! IMAX   =  NO. POINTS ALONG THE LAT. CIRCLE USED IN CALCS
! L      =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL
! NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
!
! NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE BANDTA FOR DEFINITION
! NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS
! NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE BDCOMB FOR DEFINITION
! INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
! NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS.
!
! NB,KO2 ARE SHORTWAVE PARAMETERS; OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS
!----------------------------------------------------------------------------------------
    INTEGER(KIND=I4KIND), PARAMETER :: L      = LM
    INTEGER(KIND=I4KIND), PARAMETER :: IMAX   = IM
    INTEGER(KIND=I4KIND), PARAMETER :: NCOL   = IMAX
    INTEGER(KIND=I4KIND), PARAMETER :: NBLW   = 163
    INTEGER(KIND=I4KIND), PARAMETER :: NBLX   = 47
    INTEGER(KIND=I4KIND), PARAMETER :: NBLY   = 15
    INTEGER(KIND=I4KIND), PARAMETER :: NBLM   = NBLY  - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LP1    = L     + 1
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
    INTEGER(KIND=I4KIND), PARAMETER :: LP1V   = LP1   * (1 + 2 * L / 2)
    INTEGER(KIND=I4KIND), PARAMETER :: LP121  = LP1   * NBLY
    INTEGER(KIND=I4KIND), PARAMETER :: LL3P   = 3     * L + 2
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
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IQ      , ITOP    , JTOP    , IP      , IR      , J       , I       , K1      , K2      ,   &
    & KB      , K       , KP      , KT      , NC
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & XCLD
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(IN)          ::&
    & NCLDS
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & KTOP    , KBTM
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & CAMT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1, LP1)                , INTENT(INOUT)       ::&
    & CLDFAC
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                                        ::&
    & CLDROW
!
    REAL   (KIND=R4KIND), DIMENSION(LP1, LP1, 64)                                               ::&
    & CLDIPT
!
    DO 1 IQ=MYIS,MYIE,64
        ITOP = IQ + 63
        IF (ITOP > MYIE) ITOP = MYIE
        JTOP = ITOP - IQ + 1
        DO 11 IP=1,JTOP
            IR = IQ + IP - 1
!
            IF (NCLDS(IR) == 0) THEN
                DO 25 J=1,LP1
                    DO 25 I=1,LP1
                        CLDIPT(I,J,IP) = 1.
                 25 END DO
            END IF
!
            IF (NCLDS(IR) >= 1) THEN
                XCLD = 1. - CAMT(IR,2)
                K1   = KTOP(IR,2) + 1
                K2   = KBTM(IR,2)
!
                DO 27 J=1,LP1
                    CLDROW(J) = 1.
             27 END DO
!
                DO 29 J=1,K2
                    CLDROW(J) = XCLD
             29 END DO
!
                KB = MAX(K1,K2+1)
!
                DO 33 K=KB,LP1
                    DO 33 KP=1,LP1
                        CLDIPT(KP,K,IP) = CLDROW(KP)
             33 END DO
!
                DO 37 J=1,LP1
                    CLDROW(J) = 1.
             37 END DO
!
                DO 39 J=K1,LP1
                    CLDROW(J) = XCLD
             39 END DO
!
                KT = MIN(K1-1,K2)
!
                DO 43 K=1,KT
                    DO 43 KP=1,LP1
                        CLDIPT(KP,K,IP) = CLDROW(KP)
            43  END DO
!    
                IF (K2+1 <= K1-1) THEN
!
                    DO 31 J=K2+1,K1-1
                        DO 31 I=1,LP1
                            CLDIPT(I,J,IP) = 1.
                 31 END DO
!
                ELSE IF (K1 <= K2) THEN
!
                    DO 32 J=K1,K2
                        DO 32 I=1,LP1
                            CLDIPT(I,J,IP) = XCLD
                 32 END DO
!
                END IF
!    
            END IF
!
            IF (NCLDS(IR) >= 2) THEN
                DO 21 NC=2,NCLDS(IR)
                    XCLD = 1. - CAMT(IR,NC+1)
                    K1 = KTOP(IR,NC+1) + 1
                    K2 = KBTM(IR,NC+1)
!
                    DO 47 J=1,LP1
                        CLDROW(J) = 1.
                 47 END DO
!
                    DO 49 J=1,K2
                        CLDROW(J) = XCLD
                 49 END DO
!
                    KB = MAX(K1,K2+1)
!
                    DO 53 K=KB,LP1
                        DO 53 KP=1,LP1
                            CLDIPT(KP,K,IP) = CLDIPT(KP,K,IP) * CLDROW(KP)
                 53 END DO
!
                    DO 57 J=1,LP1
                        CLDROW(J) = 1.
                 57 END DO
!
                    DO 59 J=K1,LP1
                        CLDROW(J) = XCLD
                 59 END DO
!
                    KT = MIN(K1-1,K2)
                    DO 63 K=1,KT
                        DO 63 KP=1,LP1
                            CLDIPT(KP,K,IP) = CLDIPT(KP,K,IP) * CLDROW(KP)
                 63 END DO
!
                    IF (K1 <= K2) THEN
                        DO 52 J=K1,K2
                            DO 52 I=1,LP1
                                CLDIPT(I,J,IP) = CLDIPT(I,J,IP) * XCLD
                     52 END DO
                    END IF

             21 END DO
            END IF
!
     11 END DO
!
        DO 71 J=1,LP1
            DO 71 I=1,LP1
                DO 71 IP=1,JTOP
                    IR = IQ + IP - 1
                    CLDFAC(IR,I,J) = CLDIPT(I,J,IP)
     71 END DO
!
  1 END DO
!
    RETURN
!
    END SUBROUTINE CLO89
