    SUBROUTINE VDIFH(LMHK, KTM, DTQ2, THZ0, QZ0, AKHS, CT, CKLQ, T, Q, AKH, APE, Z, EVPR)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VDIFH
!>
!> SUBPROGRAM: VDIFH - VERTICAL DIFFUSION OF MASS VARIABLES
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
!> THZ0  - POTENTIAL TEMPERATURE AT TOP OF VISCOUS LAYER
!> QZ0   - SPECIF HUMIDITY AT TOP OF VISCOUS LAYER
!> AKHS  - SURFACE COEF. FOR T AND Q DIVIDED BY DELTA Z
!> CT    -
!> CKLQ  - MASK VALUE
!> AKH   - SURFACE EXCHANGE COEF. FOR T AND Q (TRANSPOSED)
!> APE   - EXNER FUNCTION (TRANSPOSED)
!> Z     - INTERFACE HEIGHT (TRANSPOSED)
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> T     - (TRANSPOSED)
!> Q     - (TRANSPOSED)
!> EVPR  - 
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : TURBL
!>
!> CALLS      : ----- 
!>--------------------------------------------------------------------------------------------------	
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LM1 = LM - 1
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(INOUT)       ::&
    & T       , Q
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::&
    & APE
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(IN)          ::&
    & Z
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(IN)          ::&
    & AKH      
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                                        ::&
    & CM      , CR      , RST     , RSQ     , DTOZ    , AKCT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMHK    , KTM      
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DTQ2    , THZ0    , QZ0     , AKHS    , CT      , CKLQ
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & EVPR
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMHM    , LMHP    , L       , KT   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DTDIF   , DTOZL   , CF      , DTOZS   , AKHH    , AKQS    , CMB     , CMTB    , RSTB    ,   &
    & RSQB    , RCML    , CMQB
!
    DTDIF = DTQ2 / FLOAT(KTM)
    LMHM  = LMHK-1
    LMHP  = LMHK+1
!
    EVPR = 0.E0 !DULE
!
    DO 100 L=1,LMHM
        DTOZ(L) =  DTDIF    /  (Z(L) - Z(L+1))
          CR(L) = -DTOZ(L)  * AKH(L)
        AKCT(L) =  AKH(L)   *  (Z(L) - Z(L+2)) * 0.5 * CT
100 END DO
!
    CM(1) = DTOZ(1) * AKH(1) + 1.
!
    DO 300 KT=1,KTM
!
        RST(1) = -AKCT(1) * DTOZ(1) + T(1) * APE(1)
        RSQ(1) =     Q(1)
!
        DO 110 L=2,LMHM
            DTOZL  =  DTOZ(L)
            CF     = -DTOZL     * AKH(L-1) /    CM(L-1)
            CM(L)  = -  CR(L-1) * CF       +  (AKH(L-1) +  AKH(L)) * DTOZL + 1.
            RST(L) = - RST(L-1) * CF       + (AKCT(L-1) - AKCT(L)) * DTOZL + T(L) * APE(L)
            RSQ(L) = - RSQ(L-1) * CF       +     Q(L)
    110 END DO
!
        DTOZS   =  DTDIF     / (Z(LMHK) - Z(LMHP))
        AKHH    =  AKH(LMHM)
!    
        CF      = - DTOZS    * AKHH     / CM(LMHM)
        AKQS    =   AKHS     * CKLQ
!    
        CMB     =   CR(LMHM) * CF
        CMTB    = - CMB      + (AKHH    + AKHS)       * DTOZS + 1.
        CMQB    = - CMB      + (AKHH    + AKQS)       * DTOZS + 1.
!    
        RSTB    = -RST(LMHM) * CF       + (AKCT(LMHM) - AKHS  * CT) * DTOZS + T(LMHK) * APE(LMHK)
        RSQB    = -RSQ(LMHM) * CF       + Q(LMHK)
!
        T(LMHK) = ( DTOZS    * AKHS     * THZ0 + RSTB) / (APE(LMHK)  * CMTB)
        EVPR    =   EVPR     - Q(LMHK)  * DETA(LMHK)                           !DULE
        Q(LMHK) = ( DTOZS    * AKQS     * QZ0          + RSQB)       / CMQB
        EVPR    =   EVPR     + Q(LMHK)  * DETA(LMHK)                           !DULE
!
        DO 120 L=LMHM,1,-1
            RCML = 1.        / CM(L) 
            T(L) = (-CR(L)   *  T(L+1) *  APE(L+1) + RST(L)) * RCML / APE(L)
            EVPR =  EVPR     -  Q(L)   * DETA(L)                              !DULE
            Q(L) = (-CR(L)   *  Q(L+1) +  RSQ(L))  * RCML
            EVPR =  EVPR     +  Q(L)   * DETA(L)                              !DULE
    120 END DO
!
300 END DO
!
    RETURN
!
    END SUBROUTINE VDIFH
