    SUBROUTINE DIFCOF(LMHK, GM, GH, EL, Q2, Z, AKM, AKH)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE DIFCOF
!> 
!> SUBPROGRAM: DIFCOF - ?????
!> PROGRAMMER: ?????
!> ORG: ??????
!> DATE: ??-??-??       
!>     
!> ABSTRACT:
!> ?????
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
!> LMHK  - MASS POINT MODEL SURFACE (LOWEST MODEL HEIGHT)
!> GM    - WIND SHEAR
!> GH    - VERTICAL GRADIENT OS POTENTIAL TEMPERATURE
!> EL    - 
!> Q2    - TURBULENT KINETIC ENERGY (TRANSPOSED)
!> Z     - INTERFACE HEIGHT (TRANSPOSED) 
!>
!> OUTPUT ARGUMENT LIST:
!> AKM   - SURFACE EXCHANGE COEFFICIENTS FOR U AND V
!> AKH   - SURFACE EXCHANGE COEFFICIENTS FOR T AND Q  
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
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
!
!--------------------------------- 
! LEVEL 2.5 DIFFUSION COEFFICIENTS
!---------------------------------
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
    INTEGER(KIND=I4KIND), PARAMETER :: LP1 = LM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM1 = LM - 1
!
    REAL   (KIND=R4KIND), PARAMETER :: G    =  9.8
    REAL   (KIND=R4KIND), PARAMETER :: BETA =  1.  / 270. 
    REAL   (KIND=R4KIND), PARAMETER :: BTG  = BETA * G
    REAL   (KIND=R4KIND), PARAMETER :: PRT  =  1.0
    REAL   (KIND=R4KIND), PARAMETER :: GAM1 =  0.2222222222222222222 
!---------------------------------------------------------------------
!DCR: Nakanishi and Niino, 2009
    REAL   (KIND=R4KIND), PARAMETER :: A1   =  1.18
    REAL   (KIND=R4KIND), PARAMETER :: A2   =  0.665
    REAL   (KIND=R4KIND), PARAMETER :: B1   =  24.0
    REAL   (KIND=R4KIND), PARAMETER :: B2   =  15.0
    REAL   (KIND=R4KIND), PARAMETER :: C1   =  0.137
!-------------------------------------------------------------------------
!J2001    REAL   (KIND=R4KIND), PARAMETER :: A1   =  0.659888514560862645
!J2001    REAL   (KIND=R4KIND), PARAMETER :: A2   =  0.6574209922667784586
!J2001    REAL   (KIND=R4KIND), PARAMETER :: B1   = 11.87799326209552761
!J2001    REAL   (KIND=R4KIND), PARAMETER :: B2   =  7.226971804046074028
!J2001    REAL   (KIND=R4KIND), PARAMETER :: C1   =  0.000830955950095854396
!----------------------------------------------------------
! N     (G=9.8,BETA=1./270.,BTG=BETA*G
! N     ,PRT=1.0,GAM1=0.2222222222222222222
! N     ,A1=0.3310949523016403346,A2=0.8273378704055731278
! N     ,B1=5.959709141429526024,B2=3.626088092074591135
! N     ,C1=-0.3330651924968952113)
!----------------------------------------------------------
! Y     (G=9.8,BETA=1./270.,BTG=BETA*G
! Y     ,PRT=0.8,GAM1=0.2222222222222222222
! Y     ,A1=0.9222222350809054114,A2=0.7350190142719400952
! Y     ,B1=16.60000023145629741,B2=10.10000014082581951
! Y     ,C1=0.0805318118080613468)
!----------------------------------------------------------
!
!--------------------------------------------
! COEFFICIENTS FOR THE SM AND SH DETERMINANTS
!--------------------------------------------
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & BSMH = -3. * A1 * A2 * (3. * A2 + 3. * B2 * C1 + 12. * A1 * C1 - B2) * BTG
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & CESM =       A1 * (1. - 3. * C1)                                     
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & BSHM = 18. * A1 * A1 * A2 * C1                                          
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & BSHH =  9. * A1 * A2 * A2 * BTG                                                             
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & CESH = A2              
!---------------------------------------------
! COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!---------------------------------------------
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & ADNM = 18. * A1 * A1 * A2 * (B2 - 3. * A2) * BTG
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::& 
    & ADNH =  9. * A1 * A2 * A2 * (12. *  A1 + 3. * B2) * BTG * BTG
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & BDNM =  6. * A1 * A1
!
    REAL   (KIND=R4KIND), PARAMETER                                                             ::&
    & BDNH =  3. * A2 *  (7. * A1 + B2) * BTG
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::&
    & Q2
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(IN)          ::&
    & GM      , GH      , EL      
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(INOUT)       ::&
    & AKM     , AKH
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(IN)          ::&
    & Z
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMHK    
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMHM    , LMHP    , L       
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & ELL     , ELOQ2   , ELOQ4   , GML     , GHL     , ADEN    , BDEN    , CDEN    , BESM    ,   &
    & BESH    , RDEN    , ESM     , ESH     , RDZ     , Q1L     , ELQDZ 
!
    LMHM = LMHK - 1
    LMHP = LMHK + 1
!
    DO 100 L=1,LMHM
        ELL = EL(L)
!
        ELOQ2 = ELL   * ELL   / Q2(L)
        ELOQ4 = ELOQ2 * ELOQ2
!
        GML = GM(L)
        GHL = GH(L)
!---------------------------------------------
! COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!---------------------------------------------
        ADEN = (ADNM * GML + ADNH * GHL) * GHL
        BDEN =  BDNM * GML + BDNH * GHL
        CDEN = 1.
!------------------------------------ 
! COEFFICIENTS FOR THE SM DETERMINANT
!------------------------------------
        BESM = BSMH * GHL
!------------------------------------ 
! COEFFICIENTS FOR THE SH DETERMINANT
!------------------------------------
        BESH = BSHM * GML + BSHH * GHL
!--------------- 
! 1./DENOMINATOR
!---------------
        RDEN = 1. / (ADEN * ELOQ4 + BDEN * ELOQ2 + CDEN)
!---------- 
! SM AND SH
!----------
        ESM = (BESM * ELOQ2 + CESM) * RDEN
        ESH = (BESH * ELOQ2 + CESH) * RDEN
!----------------------- 
! DIFFUSION COEFFICIENTS
!----------------------- 
        RDZ    = 2. / (Z(L) - Z(L+2))
        Q1L    = SQRT(Q2(L))
        ELQDZ  = ELL * Q1L * RDZ
        AKM(L) = ELQDZ * ESM
        AKH(L) = ELQDZ * ESH
!
100 END DO
!
    RETURN
!
    END SUBROUTINE DIFCOF
