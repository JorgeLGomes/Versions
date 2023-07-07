    SUBROUTINE VDIFQ(LMHK, KTM, DTQ2, Q2, EL, Z)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VDIFQ
!>
!> SUBPROGRAM: VDIFQ - VERTICAL DIFFUSION
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT: 
!> ?????
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
!> LMHK  - MASS POINT MODEL SURFACE (LOWEST MODEL HEIGHT)
!> KTM   -
!> DTQ2  - PHYSICS TINE STEP
!> EL    - 
!> Z     - INTERFACE HEIGHT
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> Q2    - TURBULENT KINETIC ENERGY
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : TURBL
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA
!
    IMPLICIT NONE
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1  = LM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM1  = LM - 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM2  = LM - 2
!
    REAL   (KIND=R4KIND), PARAMETER :: ESQ  = 0.20
    REAL   (KIND=R4KIND), PARAMETER :: ELZ0 = 0.
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(INOUT)       ::&
    & Q2
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(IN)          ::&
    & EL
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(IN)          ::&
    & Z
!
    REAL   (KIND=R4KIND), DIMENSION(LM2)                                                        ::&
    & CM      , CR      , RSQ2    , AKQ     , DTOZ    
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMHK    , KTM  
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DTQ2 
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMHM    , LMHP    , L       , KT      , IVI     , LMH2
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DTDIF   , CF      , DTOZS   , AKQS    
!
    DTDIF = DTQ2 / FLOAT(KTM)
    LMHM  = LMHK - 1
    LMH2  = LMHK - 2
    LMHP  = LMHK + 1
!
    DO 300 KT=1,KTM
!
        DO 100 L=1,LMH2
            DTOZ(L) = (DTDIF + DTDIF) / (Z(L) - Z(L+2))
             AKQ(L) = SQRT((Q2(L) + Q2(L+1))*0.5) * (EL(L)+EL(L+1)) * 0.5 * ESQ / (Z(L+1)-Z(L+2))
              CR(L) = -DTOZ(L) * AKQ(L)
    100 END DO
!    
          CM(1) = DTOZ(1) * AKQ(1) + 1.
        RSQ2(1) =   Q2(1)
!
        DO 110 L=2,LMH2
              CF      = -DTOZ(L)   * AKQ(L-1) / CM(L-1)
              CM(L) = -  CR(L-1) * CF + (AKQ(L-1) + AKQ(L)) * DTOZ(L) + 1.
            RSQ2(L) = -RSQ2(L-1) * CF + Q2(L)
    110 END DO
!
        DTOZS = (DTDIF + DTDIF) / (Z(LMHM) - Z(LMHP))
        AKQS  = SQRT((Q2(LMHM) + Q2(LMHK))*0.5) * (EL(LMHM)+ELZ0) * 0.5 * ESQ / (Z(LMHK)-Z(LMHP))
!    
        CF    = -DTOZS * AKQ(LMH2) / CM(LMH2)
!
        Q2(LMHM) = (DTOZS * AKQS * Q2(LMHK) - RSQ2(LMH2) * CF+Q2(LMHM)) / ((AKQ(LMH2) + AKQS)     & 
    &            *  DTOZS - CR(LMH2) * CF+1.)
!    
        DO 120 IVI=1,LMH2
            L     = LMHM - IVI
            Q2(L) = (-CR(L) * Q2(L+1) + RSQ2(L)) / CM(L)
    120 END DO
!
300 END DO
!
    RETURN
!
    END SUBROUTINE VDIFQ

