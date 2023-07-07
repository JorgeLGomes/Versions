    SUBROUTINE SFCDIF(LMHK , SM   , THS  , QS   , UZ0  , VZ0  , THZ0 , QZ0  , USTAR, WSTAR, Z0   ,&
    &                 ZEFF , AKMS , AKHS , HPBL , CT   , U10  , V10  , TH02 , TH10 , Q02  , Q10  ,&
    &                 TH100, Q100 , U100 , V100 , ULM  , VLM  , T    , Q    , APE  , Z    , PD   ,&
    &                 PT   , TLM  , VEGTYP)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SFCDIF
!>
!> SUBROUTINE: SFCDIF - SURFACE LAYER
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> AMMENDED TO USE THE "EFFECTIVE ROUGHNESS" OF MASON (1986, SEE GEORGELIN ET AL., MWR JULY 1994)
!>
!> PROGRAM HISTORY LOG:
!> 96-03-28  BLACK      - ORIGINATOR
!> 97-06-??  MEYS       - MODIFIED FOR DISTRIBUTED MEMORY
!> 99-07-06  BLACK      - FULL ARRAY RATHER THAN JUST EDGES
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LMHK  - 
!> SM    - 
!> THS   - 
!> QS    -
!> ZEFF  -
!> HPBL  - 
!> VLM   -
!> T     -
!> Q     - 
!> Z     -
!> PD    - 
!> PT    - 
!> TLM   - 
!>
!> OUTPUT ARGUMENT LIST:
!> CT    - 
!> U10   - 
!> V10   -
!> TH02  -
!> TH10  - 
!> Q02   -
!> Q10   -
!> TH100 -
!> U100  -
!> V100  - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> UZ0   - 
!> VZ0   - 
!> THZ0  - 
!> QZ0   -
!> USTAR - 
!> Z0    -
!> AKMS  - 
!> AKHS  - 
!> ULM   - 
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MOMENTO
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>     
!> DRIVER     : TURBL
!>
!> CALLS      : -----
!--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MOMENTO
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
    USE Z0DATA
!
    IMPLICIT NONE
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: LP1   = LM + 1
!
    REAL   (KIND=R4KIND), PARAMETER :: WWST   = 1.2
    REAL   (KIND=R4KIND), PARAMETER :: WWST2  = WWST * WWST
    REAL   (KIND=R4KIND), PARAMETER :: G      = 9.8
    REAL   (KIND=R4KIND), PARAMETER :: USTFC  = 0.018 / G
    REAL   (KIND=R4KIND), PARAMETER :: VKRM   = 0.40 
    REAL   (KIND=R4KIND), PARAMETER :: RIC    = 0.183
    REAL   (KIND=R4KIND), PARAMETER :: RFC    = 0.191
    REAL   (KIND=R4KIND), PARAMETER :: FHNEU  = 0.8
    REAL   (KIND=R4KIND), PARAMETER :: RRIC   = 1.0   / RIC
    REAL   (KIND=R4KIND), PARAMETER :: RFAC   = RIC   / (FHNEU * RFC * RFC)
    REAL   (KIND=R4KIND), PARAMETER :: EXCM   = 0.001 
    REAL   (KIND=R4KIND), PARAMETER :: BETA   = 1.    / 270.
    REAL   (KIND=R4KIND), PARAMETER :: BTG    = BETA  * G 
    REAL   (KIND=R4KIND), PARAMETER :: ELFC   = VKRM  * BTG
    REAL   (KIND=R4KIND), PARAMETER :: CNV    = 0.608 * G / BTG
    REAL   (KIND=R4KIND), PARAMETER :: WOLD   =  .15
    REAL   (KIND=R4KIND), PARAMETER :: WNEW   = 1.    - WOLD 
!
    INTEGER(KIND=I4KIND), PARAMETER :: ITRMX  = 5  
!        
    REAL   (KIND=R4KIND), PARAMETER :: PIHF   =   3.14159265 / 2. 
    REAL   (KIND=R4KIND), PARAMETER :: PIFR   =   3.14159265 / 4. 
    REAL   (KIND=R4KIND), PARAMETER :: EPSU2  =   1.E-4
    REAL   (KIND=R4KIND), PARAMETER :: EPSUST =   0.07 
    REAL   (KIND=R4KIND), PARAMETER :: EPSIT  =   1.E-4
    REAL   (KIND=R4KIND), PARAMETER :: EPSA   =   1.E-8
    REAL   (KIND=R4KIND), PARAMETER :: ZTMIN  =  -5.
    REAL   (KIND=R4KIND), PARAMETER :: ZTMAX  =   1.
!BHS    REAL   (KIND=R4KIND), PARAMETER :: ZTMAX  =   2.
    REAL   (KIND=R4KIND), PARAMETER :: GLKBS  =  30.0
    REAL   (KIND=R4KIND), PARAMETER :: GLKBR  =  10.0
    REAL   (KIND=R4KIND), PARAMETER :: GRRS   =  GLKBR / GLKBS
    REAL   (KIND=R4KIND), PARAMETER :: VISC   =   1.5E-5
    REAL   (KIND=R4KIND), PARAMETER :: TVISC  =   2.1E-5
    REAL   (KIND=R4KIND), PARAMETER :: QVISC  =   2.1E-5
    REAL   (KIND=R4KIND), PARAMETER :: RVISC  =   1. / VISC 
    REAL   (KIND=R4KIND), PARAMETER :: RTVISC =   1. / TVISC
    REAL   (KIND=R4KIND), PARAMETER :: RQVISC =   1. / QVISC
    REAL   (KIND=R4KIND), PARAMETER :: SQPR   =   0.84
    REAL   (KIND=R4KIND), PARAMETER :: SQSC   =   0.84
    REAL   (KIND=R4KIND), PARAMETER :: ZQRZT  =  SQSC / SQPR
    REAL   (KIND=R4KIND), PARAMETER :: USTR   =   0.225  
    REAL   (KIND=R4KIND), PARAMETER :: USTC   =   0.7
    REAL   (KIND=R4KIND), PARAMETER :: FZT1   =  RVISC  * TVISC * SQPR
    REAL   (KIND=R4KIND), PARAMETER :: FZQ1   =  RTVISC * QVISC * ZQRZT
    REAL   (KIND=R4KIND), PARAMETER :: FZQ2   =  RTVISC * QVISC * ZQRZT
    REAL   (KIND=R4KIND), PARAMETER :: M      =  30.0
    REAL   (KIND=R4KIND), PARAMETER :: ZTFC   =   1.0
!    REAL   (KIND=R4KIND), PARAMETER :: CZIL   =    .2000
    REAL   (KIND=R4KIND), PARAMETER :: SQVISC = 258.2 
!    REAL   (KIND=R4KIND), PARAMETER :: ZILFC  = -CZIL * VKRM * SQVISC
    REAL   (KIND=R4KIND), PARAMETER :: PQ0    = 379.90516
    REAL   (KIND=R4KIND), PARAMETER :: A2     =  17.2693882
    REAL   (KIND=R4KIND), PARAMETER :: A3     = 273.16
    REAL   (KIND=R4KIND), PARAMETER :: A4     =  35.86
    REAL   (KIND=R4KIND), PARAMETER :: CAPA   =   0.28589641E0
    REAL   (KIND=R4KIND), PARAMETER :: H1M5   =   1.E-5
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(IN)          ::&
    & T       , Q       , APE
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(IN)          ::&
    & Z
!
    REAL   (KIND=R4KIND), DIMENSION(4)                                    , INTENT(IN)          ::&
    & ZEFF
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMHK    , VEGTYP
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & SM      , THS     , QS      , HPBL    , VLM     , PD      , PT      , TLM   
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & WSTAR   , CT      , U10     , V10     , TH02    , TH10    , Q02     , Q10     , Q100    ,   &
    & U100    , V100    , TH100
!
    REAL   (KIND=R4KIND)                                                  , INTENT(INOUT)       ::&
    & UZ0     , VZ0     , THZ0    , QZ0     , USTAR   , Z0      , AKMS    , AKHS    , ULM  
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMHP    , ML      , MH      , ITR
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & PSLMU   , ZZ      , PSLMS   , PSLHU   , PSPMU   , XX      , PSPMS   , YY      , PSPHU   ,   &
    & PSPHS   , THLM    , QLM     , RR      , SMALL   , CZIV    , FZUI1   , FZT2    , PSLHS   ,   &
    & FZU1    , ZU      , WGHT    , RWGH    , ZT      , ZQ      , ALPHA   , X       , WLOW    ,   &
    & ZSL     , RDZ     , CXCH    , DTHV    , DU2     , RIB     , BTBH    , WSTAR2  , ZSLU    ,   &
    & ZSLT    , RLOGU   , RLOGT   , BTGH    , RLMO    , ZETALT  , ZETALU  , ZETAU   , ZETAT   ,   &
    & PSMZ    , SIMM    , PSHZ    , SIMH    , USTARK  , RLMN    , RLMP    , RLMA    , XLU4    ,   &
    & XLT4    , XU4     , XT4     , XLU     , XLT     , XU      , XT      , HSFLX   , HLFLX   ,   &
    & AKMS100 , AKHS02  , AKHS100 , ZU10    , ZU100   , ZT02    , ZT10    , ZT100   , RLNU10  ,   &
    & RLNU100 , RLNT02  , RLNT10  , RLNT100 , ZTAU10  , ZTAU100 , ZTAT02  , ZTAT100 , SIMM10  ,   &
    & SIMM100 , XLU104  , AKHS10  , ZTAT10  , SIMH02  , SIMH10  , SIMH100 , XLU1004 , XLT024  ,   &
    & XLT104  , XLT1004 , XLU10   , XLU100  , XLT02   , XLT10   , XLT100  , PDS     , TERM1   ,   &
    & PSHLTR  , T02     , QSAT2   , T10     , QSAT10  , P100    , T100    , QSAT100 , U10E    ,   &
    & V10E    , U100E   , V100E   , ZUUZ    , ZTAU    , XU104   , XU10    , EKMS10  , EKMS100 ,   &
    & ZILFC
!
    PSLMU(ZZ) = -0.96 * ALOG(1.0 - 4.5 * ZZ)
    PSLMS(ZZ) =  ZZ * RRIC - 2.076 * (1. - 1. / (ZZ + 1.))
    PSLHU(ZZ) = -0.96 * ALOG(1.0 - 4.5 * ZZ)
    PSLHS(ZZ) =  ZZ * RFAC - 2.076 * (1. - 1. / (ZZ + 1.))
    PSPMU(XX) = -2. * ALOG((XX + 1.) * 0.5) - ALOG((XX * XX + 1.) * 0.5) + 2. * ATAN(XX) - PIHF
    PSPHU(XX) = -2. * ALOG((XX * XX + 1.) * 0.5)

    PSPMS(YY) =  5. * YY
    PSPHS(YY) =  5. * YY
! BHS
!    PSPMS(YY) = YY + 0.667 * (YY - 5. / 0.35) * EXP(-0.35 * YY) + 9.5238
!    PSPHS(YY) = (1. + YY * 2. / 3.) ** 1.5 + 0.667 * (YY - 5. / 0.35) * EXP(-0.35 * YY) + 9.5238 -1.
!
!DCR
!Czil Dinamico
!DCR    CZIL = 10**(-0.4*(Z0_DATA(VEGTYP)/0.07))
    ZILFC  = -CZIL_DATA(VEGTYP) * VKRM * SQVISC
!    ZILFC  = -(0.2) * VKRM * SQVISC
!DCR 
   LMHP = LMHK + 1
!
    THLM = T(LMHK) * APE(LMHK)
    QLM  = Q(LMHK)
!
    Z0 = (1. - SM) * Z0 + SM * AMAX1(USTFC * USTAR * USTAR, 1.59E-5)
!-----------------  
! VISCOUS SUBLAYER 
!-----------------
    IF (SM > 0.5 .AND. USTAR < USTC) THEN
!----------------------------
! MALAGUTTI CLASS EXPERIMENTS
!----------------------------
        RR    = USTAR * Z0 * 1 / VISC
        SMALL = 11 / (M * (RR ** 0.25))
        CZIV  = SMALL * GLKBS
        FZU1  = CZIV * VISC
        FZT2  = CZIV * GRRS * TVISC * SQPR
!
        IF (USTAR < USTR) THEN
!        
            ZU   = FZU1 * SQRT(SQRT(Z0 * USTAR * RVISC)) / USTAR
            WGHT = AKMS * ZU * RVISC
            RWGH = WGHT / (WGHT + 1.)
            UZ0  = (ULM * RWGH + UZ0) * 0.5
            VZ0  = (VLM * RWGH + VZ0) * 0.5
!        
            ZT   = FZT1 * ZU
            WGHT = AKHS * ZT * RTVISC
            THZ0 = ((WGHT * THLM + THS) / (WGHT + 1.) + THZ0) * 0.5
!        
            ZQ   = FZQ1 * ZT
            WGHT = AKHS * ZQ * RQVISC
            QZ0  = ((WGHT * QLM + QS)   / (WGHT + 1.) + QZ0)  * 0.5
!        
        END IF
!
        IF (USTAR >= USTR .AND. USTAR < USTC) THEN
!        
            ZU  = Z0
            UZ0 = 0.
            VZ0 = 0.
!        
            ZT   = FZT2 * SQRT(SQRT(Z0 * USTAR * RVISC)) / USTAR
            WGHT = AKHS * ZT * RTVISC
            THZ0 = ((WGHT * THLM + THS) / (WGHT + 1.) + THZ0) * 0.5
!        
            ZQ   = FZQ2 * ZT
            WGHT = AKHS * ZQ * RQVISC
            QZ0  = ((WGHT * QLM + QS)   / (WGHT + 1.) + QZ0)  * 0.5
        END IF
!
    ELSE
!
        ZU = Z0
!
        IF (SM <= 0.5) THEN
            IF (ULM == 0.) ULM = EPSU2
!
            ALPHA = ABS(ATAN(VLM / ULM) + PIHF - EPSA)
            X     = ALPHA / PIFR
            ML    = 1 + X
            ML    = MIN(4,ML)
            MH    = 1 + MOD(ML,4)
            WLOW  = X - ML + 1
            ZU    = WLOW * ZEFF(ML) + (1. - WLOW) * ZEFF(MH)
        END IF
!
        UZ0 = 0.
        VZ0 = 0.
!    
        ZT   = Z0
        THZ0 = THS
!    
        ZQ  = Z0
        QZ0 = QS
!
    END IF
!
    ZSL = (Z(LMHK) - Z(LMHP)) * 0.5
!
    ZU = AMIN1(ZU,0.5*ZSL)
!
    RDZ = 1. / ZSL
    CXCH = EXCM * RDZ
!
    IF (SM > 0.5) THEN
        DTHV = (0.608 * QLM + 1.) * THLM - (0.608 * QZ0 + 1.) * THZ0
    ELSE
        DTHV = (QLM - QZ0) * CNV + THLM - THZ0
        ZT = Z0 * ZTFC
    END IF
!
    DU2 = AMAX1((ULM - UZ0) ** 2 + (VLM-VZ0) ** 2, EPSU2)
!
    RIB = BTG * DTHV * ZSL / DU2
!----------------------------    
! BELJARS CORRECTION OF USTAR 
!---------------------------- 
    BTGH   = BTG * HPBL
    WSTAR2 = WWST2 * ABS(BTGH * AKHS * DTHV) ** (2. / 3.)
    USTAR  = AMAX1(SQRT(AKMS * SQRT(DU2 + WSTAR2)), EPSUST)
!---------------------------- 
! MALAGUTTI CLASS EXPERIMENTS
!---------------------------- 
    RR    = USTAR * Z0 * 1 / VISC
    SMALL = 11 / (M * (RR ** 0.25))
    CZIV  = SMALL * GLKBS
    FZU1  = CZIV * VISC
    FZT2  = CZIV * GRRS * TVISC * SQPR
!--------------------------  
! ZILITINKEVITCH FIX FOR ZT 
!--------------------------
    IF (SM < 0.5) ZT = EXP(ZILFC * SQRT(USTAR * Z0)) * Z0
!
    IF (SM > 0.5 .AND. RIB >= RIC) THEN
!
        AKMS = AMAX1( VISC * RDZ, CXCH)
        AKHS = AMAX1(TVISC * RDZ, CXCH)
!
    ELSE
!
        ZSLU = ZSL + ZU
        ZSLT = ZSL + ZT
!    
        RLOGU = ALOG(ZSLU / ZU)
        RLOGT = ALOG(ZSLT / ZT)
!    
        RLMO = ELFC * AKHS * DTHV / USTAR ** 3
!----------------- 
! SEA POINTS FIRST
!-----------------
        IF (SM > 0.5) THEN
            DO 100 ITR=1,ITRMX
!-------------------------------
! 1./MONIN-OBUKKHOV LENGTH-SCALE 
!-------------------------------
                ZETALT = AMAX1(ZSLT*RLMO,ZTMIN)
                RLMO   = ZETALT / ZSLT
                ZETALU = ZSLU   * RLMO            
                ZETAU  = ZU     * RLMO
                ZETAT  = ZT     * RLMO
!---------------------- 
! LL FUNCTIONS OVER SEA 
!----------------------
                IF (RLMO < 0.) THEN
                    PSMZ =          PSLMU(ZETAU)
                    SIMM =          PSLMU(ZETALU) - PSMZ + RLOGU
                    PSHZ =          PSLHU(ZETAT)
                    SIMH = FHNEU * (PSLHU(ZETALT) - PSHZ + RLOGT)
                ELSE
                    PSMZ =         PSLMS(ZETAU)
                    SIMM =         PSLMS(ZETALU) - PSMZ  + RLOGU
                    PSHZ =          PSLHS(ZETAT)
                    SIMH = FHNEU * (PSLHS(ZETALT) - PSHZ + RLOGT)
                END IF
!------------------------------ 
! BELJAARS CORRECTION FOR USTAR
!------------------------------
                USTAR = AMAX1(SQRT(AKMS * SQRT(DU2 + WSTAR2)), EPSUST)
!---------------------------- 
! MALAGUTTI CLASS EXPERIMENTS
!----------------------------
                RR    = USTAR * Z0 * 1 / VISC
                SMALL = 11 / (M * (Rr ** 0.25))
                CZIV  = SMALL * GLKBS
                FZU1  = CZIV * VISC
                FZT2  = CZIV * GRRS * TVISC * SQPR
!
                USTARK = USTAR * VKRM
                AKMS   = AMAX1(USTARK / SIMM, CXCH)
                AKHS   = AMAX1(USTARK / SIMH, CXCH)
!
                WSTAR2 = WWST2 * ABS(BTGH * AKHS * DTHV) ** (2. / 3.)
                RLMN   = ELFC * AKHS * DTHV / USTAR ** 3
!
                RLMP = RLMO
                RLMA = RLMO * WOLD + RLMN * WNEW
!
                RLMO = RLMA
!
        100 END DO
!
        110 CONTINUE
!----------------------------
! END OF SEA POINT PROCESSING 
!----------------------------
        ELSE
!---------------- 
! NOW LAND POINTS 
!---------------- 
            DO 200 ITR=1,ITRMX
!-------------------------------
! 1./MONIN-OBUKKHOV LENGTH-SCALE
!-------------------------------
                ZETALT = AMAX1(ZSLT * RLMO, ZTMIN)
                RLMO   = ZETALT/ZSLT
                ZETALU = ZSLU * RLMO   
                ZETAU  = ZU   * RLMO
                ZETAT  = ZT   * RLMO
!-----------------------------------------------
! PAULSON 1970 FUNCTIONS OVER LAND W RAD. SKIN T
!-----------------------------------------------
                IF (RLMO < 0.) THEN
                    XLU4 = 1. -16. * ZETALU
                    XLT4 = 1. -16. * ZETALT
                    XU4  = 1. -16. * ZETAU
                    XT4  = 1. -16. * ZETAT
!                
                    XLU = SQRT(SQRT(XLU4))
                    XLT = SQRT(SQRT(XLT4))
                    XU  = SQRT(SQRT(XU4 ))
                    XT  = SQRT(SQRT(XT4 ))
!                
                    PSMZ = PSPMU(XU )
                    SIMM = PSPMU(XLU) - PSMZ + RLOGU
                    PSHZ = PSPHU(XT )
                    SIMH = PSPHU(XLT) - PSHZ + RLOGT
                ELSE
                    ZETAU  = AMIN1(ZETAU , ZTMAX)
                    ZETAT  = AMIN1(ZETAT , ZTMAX)
                    ZETALU = AMIN1(ZETALU, ZTMAX)
                    ZETALT = AMIN1(ZETALT, ZTMAX)
                    PSMZ   = PSPMS(ZETAU )
                    SIMM   = PSPMS(ZETALU) - PSMZ + RLOGU
                    PSHZ   = PSPHS(ZETAT )
                    SIMH   = PSPHS(ZETALT) - PSHZ + RLOGT
                END IF
!------------------------------ 
! BELJAARS CORRECTION FOR USTAR 
!------------------------------ 
                USTAR = AMAX1(SQRT(AKMS * SQRT(DU2 + WSTAR2)), EPSUST)
!-------------------------- 
! ZILITINKEVITCH FIX FOR ZT
!--------------------------  
                ZT    = EXP(ZILFC * SQRT(USTAR * Z0)) * Z0
                ZSLT  = ZSL + ZT
                RLOGT = ALOG(ZSLT / ZT)
!
                USTARK = USTAR * VKRM
                AKMS   = AMAX1(USTARK / SIMM, CXCH)
                AKHS   = AMAX1(USTARK / SIMH, CXCH)
!
                WSTAR2 = WWST2 * ABS(BTGH * AKHS * DTHV) ** (2./3.)
                RLMN   = ELFC * AKHS * DTHV / USTAR ** 3
!
                RLMP = RLMO
                RLMA = RLMO * WOLD + RLMN * WNEW
                RLMO = RLMA
!
        200 END DO
!
        210 CONTINUE
!---------------------------------------------------- 
! END OF LAND POINT PROCESSING AND SEA-LAND BRANCHING
!----------------------------------------------------
        END IF
!------------------------------------------ 
! END OF TURBULENCE-NO TURBULENCE BRANCHING
!------------------------------------------ 
    END IF
!-------------------- 
! COUNTERGRADIENT FIX 
!-------------------- 
    CT = 0.
!-----------------  
! DIAGNOSTIC BLOCK
!-----------------
    WSTAR = SQRT(WSTAR2) / WWST
!
    UMFLX = AKMS * (ULM  - UZ0 )
    VMFLX = AKMS * (VLM  - VZ0 )
    HSFLX = AKHS * (THLM - THZ0)
    HLFLX = AKHS * (QLM  - QZ0 )
!
    IF (SM > 0.5 .AND. RIB >= RIC) THEN
!
        AKMS10 = AMAX1( VISC / 10., 4. * CXCH)
!---------   
! SM v100M
!---------
        AKMS100 = AMAX1( VISC / 100., 4. * CXCH)
        AKHS02  = AMAX1(TVISC /  02., 4. * CXCH)    
        AKHS10  = AMAX1(TVISC /  10., 4. * CXCH)
        AKHS100 = AMAX1(TVISC / 100., 4. * CXCH)
!
    ELSE
!
        ZU10 = ZU + 10.
!---------   
! SM v100M
!---------
        ZU100   = ZU + 100.
        ZT02    = ZT + 02.
        ZT10    = ZT + 10.
        ZT100   = ZT + 100.
!
        RLNU10  = ALOG(ZU10  / ZU)
        RLNU100 = ALOG(ZU100 / ZU)
        RLNT02  = ALOG(ZT02  / ZT)
        RLNT10  = ALOG(ZT10  / ZT)
        RLNT100 = ALOG(ZT100 / ZT)
!
        ZTAU10  = ZU10  * RLMP
        ZTAU100 = ZU100 * RLMP
        ZTAT02  = ZT02  * RLMP
        ZTAT10  = ZT10  * RLMP
        ZTAT100 = ZT100 * RLMP
!---------------------- 
! LL FUNCTIONS OVER SEA 
!---------------------- 
        IF (SM > 0.5) THEN
!
            IF (RLMP < 0.) THEN
                SIMM10  = PSLMU(ZTAU10)  - PSMZ + RLNU10
!---------   
! SM v100M
!---------
                SIMM100 = PSLMU(ZTAU100) - PSMZ + RLNU100
                SIMH02  = FHNEU * (PSLHU(ZTAT02)  - PSHZ + RLNT02)
                SIMH10  = FHNEU * (PSLHU(ZTAT10)  - PSHZ + RLNT10)
                SIMH100 = FHNEU * (PSLHU(ZTAT100) - PSHZ + RLNT100)
!
            ELSE
!
                SIMM10 = PSLMS(ZTAU10)   - PSMZ + RLNU10
!---------   
! SM v100M
!---------
                SIMM100 = PSLMS(ZTAU100) - PSMZ + RLNU100
                SIMH02  = FHNEU * (PSLHS(ZTAT02)  - PSHZ + RLNT02)
                SIMH10  = FHNEU * (PSLHS(ZTAT10)  - PSHZ + RLNT10)
                SIMH100 = FHNEU * (PSLHS(ZTAT100) - PSHZ + RLNT100)
            END IF
!----------------------------------------------- 
! PAULSON 1970 FUNCTIONS OVER LAND W RAD. SKIN T 
!----------------------------------------------- 
        ELSE
!
            IF (RLMP < 0.) THEN
                XLU104 = 1. -16. * ZTAU10
!---------   
! SM v100M
!---------
                XLU1004 = 1. -16. * ZTAU100
                XLT024  = 1. -16. * ZTAT02
                XLT104  = 1. -16. * ZTAT10
                XLT1004 = 1. -16. * ZTAT100
!
                XLU10   = SQRT(SQRT(XLU104 ))
                XLU100  = SQRT(SQRT(XLU1004))
                XLT02   = SQRT(SQRT(XLT024 ))
                XLT10   = SQRT(SQRT(XLT104 ))
                XLT100  = SQRT(SQRT(XLT1004))
!
                SIMM10  = PSPMU(XLU10 ) - PSMZ + RLNU10
                SIMM100 = PSPMU(XLU100) - PSMZ + RLNU100
                SIMH02  = PSPHU(XLT02 ) - PSHZ + RLNT02
                SIMH10  = PSPHU(XLT10 ) - PSHZ + RLNT10
                SIMH100 = PSPHU(XLT100) - PSHZ + RLNT100
            ELSE
                ZTAU10 = AMIN1(ZTAU10,ZTMAX)
!---------   
! SM v100M
!---------
                ZTAU100 = AMIN1(ZTAU100, ZTMAX)
                ZTAT02  = AMIN1(ZTAT02 , ZTMAX)
                ZTAT10  = AMIN1(ZTAT10 , ZTMAX)
                ZTAT100 = AMIN1(ZTAT100, ZTMAX)     
!
                SIMM10  = PSPMS(ZTAU10 ) - PSMZ + RLNU10
                SIMM100 = PSPMS(ZTAU100) - PSMZ + RLNU100
                SIMH02  = PSPHS(ZTAT02 ) - PSHZ + RLNT02
                SIMH10  = PSPHS(ZTAT10 ) - PSHZ + RLNT10
                SIMH100 = PSPHS(ZTAT100) - PSHZ + RLNT100
            END IF
!
        END IF
!
        AKMS10 = AMAX1(USTARK / SIMM10, CXCH)
!---------   
! SM v100M
!---------
        AKMS100 = AMAX1(USTARK / SIMM100, CXCH)
        AKHS02  = AMAX1(USTARK / SIMH02 , CXCH)
        AKHS10  = AMAX1(USTARK / SIMH10 , CXCH)
        AKHS100 = AMAX1(USTARK / SIMH100, CXCH)
    END IF
!
    U10 = UMFLX / AKMS10 + UZ0
!---------   
! SM v100M
!---------
    U100  = UMFLX / AKMS100 +  UZ0
    V100  = VMFLX / AKMS100 +  VZ0
    V10   = VMFLX / AKMS10  +  VZ0
    TH02  = HSFLX / AKHS02  + THZ0
    TH10  = HSFLX / AKHS10  + THZ0
    TH100 = HSFLX / AKHS100 + THZ0
!--------------------------------------------------------------------------------------------------  
! GSM - CHANGED THIS SECTION IN RESPONSE TO PROBLEM WITH 2-M DEW POINT OCCASIONALLY BEING GREATER
!       THAN 2-M TEMPERATURE AND SIMILAR PROBLEM AT 10-M. NOW, A SATURATION Q IS CALCULATED AT 
!       EACH LEVEL, AND THE Q IS CONSTRAINED TO BE NO HIGHER THAN THE SATURATION VALUE.
!--------------------------------------------------------------------------------------------------
    PDS    = PD + PT
    TERM1  = -0.068283 / TLM
    PSHLTR = PDS * EXP(TERM1)
    T02    = TH02 * (PSHLTR * H1M5) ** CAPA
    QSAT2  = PQ0 / PSHLTR * EXP(A2 * (T02 - A3) / (T02 - A4))
    Q02    = HLFLX / AKHS02 + QZ0
!
    IF (Q02 < 0.) THEN
        IF (QLM > 0.) THEN
            Q02 = QLM
        ELSE
            Q02 = 0.0001
        END IF
    END IF
!
    IF (Q02 > QSAT2) THEN
        Q02 = QSAT2
    END IF
!
    T10    = TH10 * (PSHLTR * H1M5) ** CAPA
    QSAT10 = PQ0 / PSHLTR * EXP(A2 * (T10 - A3) / (T10 - A4))
    Q10    = HLFLX / AKHS10 + QZ0
!
    IF (Q10 < 0.) THEN
        IF (QLM > 0.) THEN
            Q10 = QLM
        ELSE
            Q10 = 0.0001
        END IF
    END IF
!
    IF (Q10 > QSAT10) THEN
        Q10 = QSAT10
    END IF
!---------   
! SM v100M
!---------
    P100    = PDS * EXP(-100.0 * G / (287.04 * TLM))
    T100    = TH100 * (P100 * H1M5) ** CAPA
    QSAT100 = PQ0 / P100 * EXP(A2 * (T100 - A3) / (T100 - A4))
    Q100    = HLFLX / AKHS100 + QZ0
!
    IF (Q100 < 0.) THEN
        IF (QLM > 0.) THEN
            Q100 = QLM
        ELSE
            Q100 = 0.0001
        END IF
    END IF
!
    IF (Q100 > QSAT100) THEN
        Q100 = QSAT100
    END IF
!---------   
! SM v100M
!---------
!
!------------------------------
! NEW CLACULATION OF 10-M WINDS
!------------------------------ 
    U10E = U10
    V10E = V10
!---------   
! SM v100M
!---------
    U100E = U100
    V100E = V100
!
    IF (SM < 0.5)  THEN
        ZUUZ = AMIN1(ZU * 0.50, 0.10)
        ZU   = AMAX1(ZU * 0.10, ZUUZ)
        ZU10 = ZU + 10.
!---------   
! SM v100M
!---------
        ZU100   = ZU + 100.
        RLNU10  = ALOG(ZU10 / ZU)
        RLNU100 = ALOG(ZU100 / ZU)
        ZTAU    = ZU    * RLMP
        ZTAU10  = ZU10  * RLMP
        ZTAU100 = ZU100 * RLMP
!
        IF (RLMP < 0) THEN
            XLU104 = 1. -16. * ZTAU10
!---------   
! SM v100M
!---------
            XLU1004 = 1. -16. * ZTAU100
            XU104   = 1. -16. * ZTAU
            XLU10   = SQRT(SQRT(XLU104))
            XLU100  = SQRT(SQRT(XLU1004))
            XU10    = SQRT(SQRT(XU104))
            SIMM10  = PSPMU(XLU10)  - PSPMU(XU10) + RLNU10
            SIMM100 = PSPMU(XLU100) - PSPMU(XU10) + RLNU100
        ELSE
            ZTAU10 = AMIN1(ZTAU10, ZTMAX)
!---------   
! SM v100M
!---------
            ZTAU100 = AMIN1(ZTAU100, ZTMAX)
            SIMM10  = PSPMS(ZTAU10)  - PSPMS(ZTAU) + RLNU10
            SIMM100 = PSPMS(ZTAU100) - PSPMS(ZTAU) + RLNU100
        END IF
!
        EKMS10 = AMAX1(USTARK / SIMM10, CXCH)
!---------   
! SM v100M
!---------
        EKMS100 = AMAX1(USTARK / SIMM100, CXCH)
        U10E    = UMFLX / EKMS10  + UZ0
        V10E    = VMFLX / EKMS10  + VZ0
        U100E   = UMFLX / EKMS100 + UZ0
        V100E   = VMFLX / EKMS100 + VZ0
    END IF
!
      U10=U10E
      V10=V10E
!CGSM v100m
      U100=U100E
      V100=V100E
!CGSM v100m
!
    RETURN
!
    END SUBROUTINE SFCDIF
