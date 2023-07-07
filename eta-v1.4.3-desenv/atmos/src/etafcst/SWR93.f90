    SUBROUTINE SWR93(FSWC ,      HSWC,      UFSWC,     DFSWC,     FSWL,      HSWL,      UFSWL,     &
    &                DFSWL,     PRESS,     COSZRO,    TAUDAR,     RH2O,     RRCO2,     SSOLAR,     &
    &                  QO3,     NCLDS,     KTOPSW,    KBTMSW,     CAMT,       CRR,        CTT,     &
    &                 ALVB,      ALNB,       ALVD,      ALND,    GDFVB,     GDFNB,      GDFVD,     &
    &                GDFND)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SWR93
!>
!> SUBPROGRAM: SWR93 - CALCULATE SURFACE CONDITIONS
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT:
!> ?????
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????       - ORIGINATOR
!> 92-03-06  K. CAMPANA  - INCLUDE HPCON, PARMC CHANGED TO HCON, RDPARM
!> 93-02-??  Y. HOU      - INPUTS 12 BANDS CLD REFLECTTANCE AND TRANSMITTANCE CRR, CTT TO REPLACE
!>                         CIRAB, CIRRF, CUVRF
!> 18-01-15  LUCCI       - MODERNIZATION OF THE CODE, INCLUDING:
!>                         * F77 TO F90/F95
!>                         * INDENTATION & UNIFORMIZATION CODE
!>                         * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                         * DOCUMENTATION WITH DOXYGEN
!>                         * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> PRESS  -
!> COSZRO -
!> TAUDAR -
!> RH2O   -
!> RRCO2  -
!> SSOLAR -
!> QO3    -
!> NCLDS  -
!> KTOPSW -
!> KBTMSW -
!> CRR    -
!> CTT    - 
!> CAMT   -
!> ALVB   -
!> ALVD   -
!> ALNB   -
!> ALND   -
!>
!> OUTPUT ARGUMENT LIST:
!> FSWC   -
!> HSWC   -
!> UFSWC  -
!> DFSWC  -
!> FSWL   -
!> HSWL   -
!> UFSWL  -
!> DFSWL  -
!> GDFVB  -
!> GDFVD  -
!> GDFNB  -
!> GDFND  -
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
!>              RDPARM
!>              SWRSAV
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
    USE RDPARM
    USE SWRSAV
    USE TEMPCOM
    USE TOPO
! 
    IMPLICIT NONE
!
#include "sp.h"
!-------------------------------------------------------------------------------------
! PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
! L =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL
! NB IS A SHORTWAVE PARAMETER; OTHER QUANTITIES ARE DERIVED FROM THE ABOVE PARAMETERS.
! VARIABLES AS IN ARGUMENT LIST
!-------------------------------------------------------------------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(INOUT)       ::&
    & FSWC    , HSWC    , FSWL    , HSWL    , UFSWC   , DFSWC   , UFSWL   , DFSWL
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(INOUT)       ::&
    & GDFVB   , GDFNB   , GDFVD   , GDFND
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & PRESS   , CAMT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, L)                       , INTENT(IN)          ::&
    & RH2O    , QO3      
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, NB, LP1)                 , INTENT(IN)          ::&
    & CRR     , CTT 
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, LP1)                     , INTENT(IN)          ::&
    & KBTMSW  , KTOPSW
!  
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(IN)          ::&
    & COSZRO  , TAUDAR  , ALVB    , ALNB    , ALVD    , ALND
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(IN)          ::&
    & NCLDS    
!---------------- 
! LOCAL VARIABLES
!----------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & PP      , DP      , PR2     , DU      , DUCO2   , DUO3    , FF      , FFCO2   , FFO3    ,   &
    & XAMT
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & RRAY    , SECZ    , REFL    , TMP1    , REFL2   , CCMAX
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, NB)                                            ::&
    & DFNTOP
! 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & UD      , UR      , UDCO2   , URCO2   , UDO3    , URO3    , TDCO2   , TUCO2   , TDO3    ,   &
    & TUO3
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LLP2)                                          ::&
    & UCO2    , UO3     , TCO2    , TO3
!    
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & DFN     , UFN     , CR      , TTD     , TTU     , CT      , PPTOP   , DPCLD
!-----------------------------
! EQUIVALENCED LOCAL VARIABLES
!-----------------------------
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LP1)                                           ::&
    & TTUB1   , TUCL1   , TTDB1   , TDCL1    , TDCL2  , UFNTRN  , UFNCLU  , TCLU    , DFNTRN  ,   &
    & DFNCLU  , TCLD    , ALFA    , ALFAU 
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I      , N        , K       , IP      , KCLDS   , NNCLDS  , JTOP    , J1      , J2      ,   &
    & J3     , KK 
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & RRCO2   , SSOLAR  , HTEMP   , TEMPF   , TEMPG   , DENOM
!    
    EQUIVALENCE (UDO3 , UO3 (IDIM1,1   ), DFNCLU)
    EQUIVALENCE (URO3 , UO3 (IDIM1,LP2 ), UFNCLU)
    EQUIVALENCE (UDCO2, UCO2(IDIM1,1   ), TCLD  )
    EQUIVALENCE (URCO2, UCO2(IDIM1,LP2 ), TCLU  )
    EQUIVALENCE (TDO3 , TO3 (IDIM1,1   ), DFNTRN)
    EQUIVALENCE (TUO3 , TO3 (IDIM1,LP2 ), UFNTRN)
    EQUIVALENCE (TDCO2, TCO2(IDIM1,1   ))
    EQUIVALENCE (TUCO2, TCO2(IDIM1,LP2 ))
    EQUIVALENCE (FF   , ALFA )
    EQUIVALENCE (FFCO2, ALFAU)
    EQUIVALENCE (FFO3 , TTDB1)
    EQUIVALENCE (DU   , TTUB1)
    EQUIVALENCE (DUCO2, TUCL1)
    EQUIVALENCE (DUO3 , TDCL1)
    EQUIVALENCE (PR2  , TDCL2)
!--------------------------------------------------------------------------------------------------
! CALCULATE SECANT OF ZENITH ANGLE (SECZ), FLUX PRESSURES(PP), LAYER WIDTH (DP) AND PRESSURE 
! SCALING FACTOR (PR2).
!--------------------------------------------------------------------------------------------------
    DO 100 I=MYIS,MYIE
        SECZ(I)     = H35E1 / SQRT(H1224E3 * COSZRO(I) * COSZRO(I) + ONE)
          PP(I,1  ) = ZERO
          PP(I,LP1) =       PRESS(I,LP1)
        TMP1(I)     = ONE / PRESS(I,LP1)
100 END DO
!
    DO 110 K=1,LM1
        DO 110 I=MYIS,MYIE
            PP(I,K+1) = HAF * (PRESS(I,K+1) + PRESS(I,K))
110 END DO
!
    DO 120 K=1,L
        DO 120 I=MYIS,MYIE
            DP (I,K) = PP(I,K+1) -  PP(I,K)
            PR2(I,K) = HAF       * (PP(I,K) + PP(I,K+1))
120 END DO
!
    DO 130 K=1,L
        DO 130 I=MYIS,MYIE
            PR2(I,K) = PR2(I,K) * TMP1(I)
130 END DO
!--------------------------------------------------------------- 
! CALCULATE ENTERING FLUX AT THE TOP FOR EACH BAND(IN CGS UNITS)
!---------------------------------------------------------------
    DO 140 N=1,NB
        DO 140 IP=MYIS,MYIE
            DFNTOP(IP,N) = SSOLAR * H69766E5 * COSZRO(IP) * TAUDAR(IP) * PWTS(N)
140 END DO
!----------------------------------------------------------------------------
! EXECUTE THE LACIS-HANSEN REFLECTIVITY PARAMETERIZATION FOR THE VISIBLE BAND
!----------------------------------------------------------------------------
    DO 150 I=MYIS,MYIE
        RRAY(I) = HP219   / (ONE + HP816   * COSZRO(I))
        REFL(I) = RRAY(I) + (ONE - RRAY(I))* (ONE - RRAYAV) * ALVB(I) / (ONE - ALVD(I) * RRAYAV)
150 END DO
!
    DO 155 I=MYIS,MYIE
        RRAY (I) = 0.104   / (ONE + 4.8     * COSZRO(I))
        REFL2(I) = RRAY(I) + (ONE - RRAY(I))* (ONE - 0.093) * ALVB(I) / (ONE - ALVD(I) * 0.093)
155 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATE PRESSURE-WEIGHTED OPTICAL PATHS FOR EACH LAYER IN UNITS OF CM-ATM. PRESSURE WEIGHTING 
! IS USING PR2.
! DU= VALUE FOR H2O;DUCO2 FOR CO2;DUO3 FOR O3.
!--------------------------------------------------------------------------------------------------
    DO 160 K=1,L
        DO 160 I=MYIS,MYIE
            DU   (I,K) =  GINV  * RH2O (I,K) * DP (I,K) * PR2(I,K)
            DUCO2(I,K) = (RRCO2 * GINV       * CFCO2)   * DP (I,K) * PR2(I,K)
            DUO3 (I,K) = (GINV  * CFO3)      * QO3(I,K) * DP (I,K)
160 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATE CLEAR SKY SW FLUX
!
! OBTAIN THE OPTICAL PATH FROM THE TOP OF THE ATMOSPHERE TO THE FLUX PRESSURE. ANGULAR FACTORS ARE 
! NOW INCLUDED. 
! UD=DOWNWARD PATH FOR H2O,WIGTH UR THE UPWARD PATH FOR H2O. CORRESPONDING QUANTITIES FOR CO2,O3 
! ARE UDCO2/URCO2 AND UDO3/URO3.
!--------------------------------------------------------------------------------------------------
    DO 200 IP=MYIS,MYIE
        UD   (IP,1) = ZERO
        UDCO2(IP,1) = ZERO
        UDO3 (IP,1) = ZERO
200 END DO
!
    DO 210 K=2,LP1
        DO 210 I=MYIS,MYIE
            UD   (I,K) = UD   (I,K-1) + DU   (I,K-1) * SECZ(I)
            UDCO2(I,K) = UDCO2(I,K-1) + DUCO2(I,K-1) * SECZ(I)
            UDO3 (I,K) = UDO3 (I,K-1) + DUO3 (I,K-1) * SECZ(I)
210 END DO
!
    DO 220 IP=MYIS,MYIE
        UR   (IP,LP1) = UD   (IP,LP1)
        URCO2(IP,LP1) = UDCO2(IP,LP1)
        URO3 (IP,LP1) = UDO3 (IP,LP1)
220 END DO
!
    DO 230 K=L,1,-1
        DO 230 IP=MYIS,MYIE
            UR   (IP,K) = UR   (IP,K+1) + DU   (IP,K) *  DIFFCTR
            URCO2(IP,K) = URCO2(IP,K+1) + DUCO2(IP,K) *  DIFFCTR
            URO3 (IP,K) = URO3 (IP,K+1) + DUO3 (IP,K) * O3DIFCTR
230 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATE CO2 ABSORPTIONS. THEY WILL BE USED IN NEAR INFRARED BANDS.SINCE THE ABSORPTION AMOUNT 
! IS GIVEN (IN THE FORMULA USED BELOW, DERIVED FROM SASAMORI) IN TERMS OF THE TOTAL SOLAR FLUX, AND
! THE ABSORPTION IS ONLY INCLUDED IN THE NEAR IR (50 PERCENT OF THE SOLAR SPECTRUM), THE 
! ABSORPTIONS ARE MULTIPLIED BY 2. SINCE CODE ACTUALLY REQUIRES TRANSMISSIONS, THESE ARE THE VALUES
! ACTUALLY STORED IN TCO2.
!--------------------------------------------------------------------------------------------------
    DO 240 K=1,LL
        DO 240 I=MYIS,MYIE
            TCO2(I,K+1) = ONE - TWO * (H235M3 * EXP(HP26 * LOG(UCO2(I,K+1) + H129M2)) - H75826M4)
240 END DO
!--------------------------------------------------------------------------------------------------
! NOW CALCULATE OZONE ABSORPTIONS. THESE WILL BE USED IN THE VISIBLE BAND.JUST AS IN THE CO2 CASE, 
! SINCE THIS BAND IS 50 PERCENT OF THE SOLAR SPECTRUM,THE ABSORPTIONS ARE MULTIPLIED BY 2. THE 
! TRANSMISSIONS ARE STORED IN TO3.
!--------------------------------------------------------------------------------------------------
    HTEMP = H1036E2 * H1036E2 * H1036E2
    DO 250 K=1,LL
        DO 250 I=MYIS,MYIE
            TO3(I,K+1) = ONE - TWO    * UO3(I,K+1)  * (H1P082 * EXP(HMP805 * LOG(ONE    + H1386E2 &
    &                  * UO3(I,K+1))) + H658M2      / (ONE    + HTEMP      * UO3(I,K+1)           &
    &                  * UO3(I,K+1)   * UO3(I,K+1)) + H2118M2 / (ONE       + UO3(I,K+1)           &
    &                  * (H42M2       + H323M4                             * UO3(I,K+1))))
!
250 END DO
!--------------------------------- 
! START FREQUENCY LOOP (ON N) HERE
!--------------------------------- 
!
!------------------------------------------------ 
! BAND 1 (VISIBLE) INCLUDES O3 AND H2O ABSORPTION
!------------------------------------------------ 
    DO 260 K=1,L
        DO 260 I=MYIS,MYIE
            TTD(I,K+1) = EXP(HM1EZ  * MIN(FIFTY,ABCFF(1) * UD(I,K+1)))
            TTU(I,K  ) = EXP(HM1EZ  * MIN(FIFTY,ABCFF(1) * UR(I,K  )))
            DFN(I,K+1) = TTD(I,K+1) * TDO3(I,K+1)
            UFN(I,K  ) = TTU(I,K  ) * TUO3(I,K  )
260 END DO
!
    DO 270 I=MYIS,MYIE
        DFN(I,1  ) = ONE
        UFN(I,LP1) = DFN(I,LP1)
270 END DO
!----------------------------------------------------------------------------------- 
! SCALE VISIBLE BAND FLUXES BY SOLAR FLUX AT THE TOP OF THE ATMOSPHERE (DFNTOP(I,1)) 
! DFSW/UFSW WILL BE THE FLUXES, SUMMED OVER ALL BANDS
!-----------------------------------------------------------------------------------
    DO 280  K=1,LP1
        DO 280  I=MYIS,MYIE
            DFSWL(I,K) = DFN (I,K) * DFNTOP(I,1)
            UFSWL(I,K) = REFL(I)   * UFN   (I,K) * DFNTOP(I,1)
280 END DO
!
    DO 285 I=MYIS,MYIE
        GDFVB(I) = DFSWL(I,LP1)      * EXP(-0.15746 *      SECZ(I))
        GDFVD(I) = ((ONE - REFL2(I)) * DFSWL(I,LP1) - (ONE-ALVB(I)) * GDFVB(I)) / (ONE - ALVD(I))
        GDFNB(I) = ZERO
        GDFND(I) = ZERO
285 END DO
!--------------------------------------------------------------------------------------------------
! NOW OBTAIN FLUXES FOR THE NEAR IR BANDS. THE METHODS ARE THE SAME AS FOR THE VISIBLE BAND, EXCEPT
! THAT THE REFLECTION AND TRANSMISSION COEFFICIENTS (OBTAINED BELOW) ARE DIFFERENT, AS RAYLEIGH 
! SCATTERING NEED NOT BE CONSIDERED.
!--------------------------------------------------------------------------------------------------
    DO 350 N=2,NB
        IF (N == 2) THEN
!--------------------------------------------------------------------------------------------------
! THE WATER VAPOR TRANSMISSION FUNCTION FOR BAND 2 IS EQUAL TO THAT OF BAND 1 (SAVED AS TTD,TTU) 
! BAND 2-9 (NEAR-IR) INCLUDES O3, CO2 AND H2O ABSORPTION
!--------------------------------------------------------------------------------------------------
            DO 290 K=1,L
                DO 290 I=MYIS,MYIE
                    DFN(I,K+1) = TTD(I,K+1) * TDCO2(I,K+1)
                    UFN(I,K  ) = TTU(I,K  ) * TUCO2(I,K  )
        290 END DO
!
        ELSE
!---------------------------------------------------------------------------------- 
! CALCULATE WATER VAPOR TRANSMISSION FUNCTIONS FOR NEAR INFRARED BANDS. 
! INCLUDE CO2 TRANSMISSION (TDCO2/TUCO2), WHICH IS THE SAME FOR ALL INFRARED BANDS.
!----------------------------------------------------------------------------------
            DO 300 K=1,L
                DO 300 I=MYIS,MYIE
                    DFN(I,K+1) = EXP(HM1EZ * MIN(FIFTY, ABCFF(N) * UD(I,K+1))) * TDCO2(I,K+1)
                    UFN(I,K  ) = EXP(HM1EZ * MIN(FIFTY, ABCFF(N) * UR(I,K  ))) * TUCO2(I,K  )
        300 END DO
!
        END IF
!--------------------------------------------------------------------------------------------------
! AT THIS POINT,INCLUDE DFN(1),UFN(LP1), NOTING THAT DFN(1)=1 FOR ALL BANDS, AND THAT 
! UFN(LP1)=DFN(LP1) FOR ALL BANDS.
!--------------------------------------------------------------------------------------------------
        DO 310 I=MYIS,MYIE
            DFN(I,1)   = ONE
            UFN(I,LP1) = DFN(I,LP1)
    310 END DO
!------------------------------------------------------------------------------- 
! SCALE THE PREVIOUSLY COMPUTED FLUXES BY THE FLUX AT THE TOP AND SUM OVER BANDS
!-------------------------------------------------------------------------------
        DO 320 K=1,LP1
            DO 320 I=MYIS,MYIE
                DFSWL(I,K) = DFSWL(I,K)            + DFN(I,K) * DFNTOP(I,N)
                UFSWL(I,K) = UFSWL(I,K) + ALNB(I)  * UFN(I,K) * DFNTOP(I,N)
    320 END DO
!
        DO 330 I=MYIS,MYIE
            GDFNB(I) = GDFNB(I) + DFN(I,LP1) * DFNTOP(I,N)
    330 END DO
!
350 END DO
!
    DO 360 K=1,LP1
        DO 360 I=MYIS,MYIE
            FSWL(I,K) = UFSWL(I,K) - DFSWL(I,K)
360 END DO
!
    DO 370 K=1,L
        DO 370 I=MYIS,MYIE
            HSWL(I,K) = RADCON * (FSWL(I,K+1) - FSWL(I,K  )) / DP(I,K)
370 END DO
!------------------------------- 
! END OF FREQUENCY LOOP (OVER N)
!------------------------------- 
!
!----------------------------- 
! CALCULATE CLOUDY SKY SW FLUX
!----------------------------- 
    KCLDS=NCLDS(MYIS)
    DO 400 I=MYIS1,MYIE
        KCLDS = MAX(NCLDS(I), KCLDS)
400 END DO
!
    DO 410 K=1,LP1
        DO 410 I=MYIS,MYIE
            DFSWC(I,K) = DFSWL(I,K)
            UFSWC(I,K) = UFSWL(I,K)
            FSWC (I,K) = FSWL (I,K)
410 END DO
!
    DO 420 K=1,L
        DO 420 I=MYIS,MYIE
            HSWC(I,K) = HSWL(I,K)
420 END DO
!
    IF (KCLDS == 0)  RETURN
!
    DO 430 K=1,LP1
        DO 430 I=MYIS,MYIE
            XAMT(I,K) = CAMT(I,K)
430 END DO
!
    DO 470 I=MYIS,MYIE
        NNCLDS   = NCLDS(I)
        CCMAX(I) = ZERO
!
        IF (NNCLDS <= 0) GOTO 470
!
        CCMAX(I) = ONE
!
        DO 450 K=1,NNCLDS
            CCMAX(I) = CCMAX(I) * (ONE - CAMT(I,K+1))
    450 END DO
!
        CCMAX(I) = ONE - CCMAX(I)
!
        IF (CCMAX(I) > ZERO) THEN
            DO 460 K=1,NNCLDS
                XAMT(I,K+1) = CAMT(I,K+1) / CCMAX(I)
        460 END DO
        END IF
!
470 END DO
!
    DO 480 K=1,LP1
        DO 480 I=MYIS,MYIE
            FF   (I,K) =  DIFFCTR
            FFCO2(I,K) =  DIFFCTR
            FFO3 (I,K) = O3DIFCTR
480 END DO
!
    DO 490 IP=MYIS,MYIE
        JTOP = KTOPSW(IP, NCLDS(IP) + 1)
        DO 490 K=1,JTOP
            FF   (IP,K) = SECZ(IP)
            FFCO2(IP,K) = SECZ(IP)
            FFO3 (IP,K) = SECZ(IP)
490 END DO
!
    DO 500 I=MYIS,MYIE
        RRAY(I) = HP219   / (ONE + HP816 * COSZRO(I))
        REFL(I) = RRAY(I) + (ONE - RRAY(I)) * (ONE - RRAYAV) * ALVD(I)                  &
    &           /           (ONE - ALVD(I)  *        RRAYAV)
!
500 END DO
!
    DO 510 IP=MYIS,MYIE
        UD   (IP,1) = ZERO
        UDCO2(IP,1) = ZERO
        UDO3 (IP,1) = ZERO
510 END DO
!
    DO 520 K=2,LP1
        DO 520 I=MYIS,MYIE
            UD   (I,K) = UD   (I,K-1) + DU   (I,K-1) * FF   (I,K)
            UDCO2(I,K) = UDCO2(I,K-1) + DUCO2(I,K-1) * FFCO2(I,K)
            UDO3 (I,K) = UDO3 (I,K-1) + DUO3 (I,K-1) * FFO3 (I,K)
520 END DO
!
    DO 530 IP=MYIS,MYIE
        UR   (IP,LP1) = UD   (IP,LP1)
        URCO2(IP,LP1) = UDCO2(IP,LP1)
        URO3 (IP,LP1) = UDO3 (IP,LP1)
530 END DO
!
    DO 540 K=L,1,-1
        DO 540 IP=MYIS,MYIE
            UR   (IP,K) = UR   (IP,K+1) + DU   (IP,K) *  DIFFCTR
            URCO2(IP,K) = URCO2(IP,K+1) + DUCO2(IP,K) *  DIFFCTR
            URO3 (IP,K) = URO3 (IP,K+1) + DUO3 (IP,K) * O3DIFCTR
540 END DO
!
    DO 550 K=1,LL
        DO 550 I=MYIS,MYIE
            TCO2(I,K+1) = ONE - TWO * (H235M3 * EXP(HP26 * LOG(UCO2(I,K+1) + H129M2)) - H75826M4)
550 END DO
!
    DO 560 K=1,LL
        DO 560 I=MYIS,MYIE
            TO3(I,K+1) = ONE - TWO * UO3(I,K+1) * (H1P082 * EXP(HMP805 * LOG(ONE + H1386E2        &
    &                  * UO3(I,K+1))) + H658M2  / (ONE + HTEMP * UO3(I,K+1) * UO3(I,K+1)          &
    &                  * UO3(I,K+1))  + H2118M2 / (ONE         + UO3(I,K+1) * (H42M2 + H323M4     &
    &                  * UO3(I,K+1))))
560 END DO
!--------------------------------------------------------------------------------------------------
! THE FIRST CLOUD IS THE GROUND; ITS PROPERTIES ARE GIVEN BY REFL (THE TRANSMISSION (0) IS  
! IRRELEVANT FOR NOW!).
!--------------------------------------------------------------------------------------------------
    DO 570 I=MYIS,MYIE
        CR(I,1) = REFL(I)
570 END DO
!--------------------------------------------------------------------------------------------------
! OBTAIN CLOUD REFLECTION AND TRANSMISSION COEFFICIENTS FOR REMAINING CLOUDS (IF ANY) IN THE 
! VISIBLE BAND THE MAXIMUM NO OF CLOUDS IN THE ROW (KCLDS) IS USED. THIS CREATES EXTRA WORK 
! (MAY BE REMOVED IN A SUBSEQUENT UPDATE).
!--------------------------------------------------------------------------------------------------
    DO 581 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 581
!
        DO 580 KK=2,KCLDS+1
            CR(I,KK) =             CRR(I,1,KK)  * XAMT(I,KK)
            CT(I,KK)= ONE - (ONE - CTT(I,1,KK)) * XAMT(I,KK)
    580 END DO
581 END DO
!--------------------------------------------------------------------------------------------------
! OBTAIN THE PRESSURE AT THE TOP,BOTTOM AND THE THICKNESS OF "THICK" CLOUDS (THOSE AT LEAST 2 
! LAYERS THICK). THIS IS USED LATER IS OBTAINING FLUXES INSIDE THE THICK CLOUDS, FOR ALL FREQUENCY
! BANDS.
!--------------------------------------------------------------------------------------------------
    DO 591 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 591
!
        DO 590 KK=1,KCLDS
            IF ((KBTMSW(I,KK+1) - 1) > KTOPSW(I,KK+1)) THEN
                PPTOP(I,KK) =                      PP(I,KTOPSW(I,KK+1))
                DPCLD(I,KK) = ONE / (PPTOP(I,KK) - PP(I,KBTMSW(I,KK+1)))
            END IF
    590 END DO
591 END DO
!
    DO 600 K=1,L
        DO 600 I=MYIS,MYIE
            TTDB1(I,K+1) = EXP(HM1EZ    * MIN(FIFTY, ABCFF(1) * UD(I,K+1)))
            TTUB1(I,K  ) = EXP(HM1EZ    * MIN(FIFTY, ABCFF(1) * UR(I,K  )))
            TTD  (I,K+1) = TTDB1(I,K+1) * TDO3(I,K+1)
            TTU  (I,K  ) = TTUB1(I,K  ) * TUO3(I,K  )
    600 END DO
!
    DO 610 I=MYIS,MYIE
        TTD(I,1)=ONE
        TTU(I,LP1)=TTD(I,LP1)
610 END DO
!--------------------------------------------------------------------------------------------------
! FOR EXECUTION OF THE CLOUD LOOP, IT IS NECESSARY TO SEPARATE OUT TRANSMISSION FCTNS AT THE TOP 
! AND BOTTOM OF THE CLOUDS, FOR EACH BAND N. THE REQUIRED QUANTITIES ARE:
! TTD(I,KTOPSW(I,K),N)  K RUNS FROM 1 TO NCLDS(I)+1:
! TTU(I,KTOPSW(I,K),N)  K RUNS FROM 1 TO NCLDS(I)+1:
! TTD(I,KBTMSW(I,K),N)  K RUNS FROM 1 TO NCLDS(I)+1:
! AND INVERSES OF THE FIRST TWO. THE ABOVE QUANTITIES ARE STORED IN TDCL1, TUCL1, TDCL2, AND 
! DFNTRN, UFNTRN, RESPECTIVELY, AS THEY  HAVE MULTIPLE USE IN THE PGM. FOR FIRST CLOUD LAYER 
! (GROUND) TDCL1,TUCL1 ARE KNOWN:
!--------------------------------------------------------------------------------------------------
    DO 620 I=MYIS,MYIE
        TDCL1 (I,1) =       TTD   (I,LP1)
        TUCL1 (I,1) =       TTU   (I,LP1)
        TDCL2 (I,1) =       TDCL1 (I,1  )
        DFNTRN(I,1) = ONE / TDCL1 (I,1  )
        UFNTRN(I,1) =       DFNTRN(I,1  )
620 END DO
!
    DO 631 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 631
!
        DO 630 KK=2,KCLDS+1
            TDCL1(I,KK) = TTD(I,KTOPSW(I,KK))
            TUCL1(I,KK) = TTU(I,KTOPSW(I,KK))
            TDCL2(I,KK) = TTD(I,KBTMSW(I,KK))
    630 END DO
631 END DO
!----------------- 
! COMPUTE INVERSES
!-----------------
    DO 641 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 641
!
        DO 640 KK=2,KCLDS
            DFNTRN(I,KK) = ONE / TDCL1(I,KK)
            UFNTRN(I,KK) = ONE / TUCL1(I,KK)
    640 END DO
641 END DO
!-------------------------------------------------------------------------------------------------- 
! COMPUTE THE TRANSMISSIVITY FROM THE TOP OF CLOUD (K+1) TO THE TOP OF CLOUD (K). THE CLOUD 
! TRANSMISSION (CT) IS INCLUDED. THIS QUANTITY IS CALLED TCLU (INDEX K). ALSO, OBTAIN THE 
! TRANSMISSIVITY FROM THE BOTTOM OF CLOUD (K+1) TO THE TOP OF CLOUD (K)(A PATH ENTIRELY OUTSIDE 
! CLOUDS). THIS QUANTITY IS CALLED TCLD (INDEX K).
!--------------------------------------------------------------------------------------------------
    DO 651 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 651
!
        DO 650 KK=1,KCLDS
            TCLU(I,KK) = TDCL1(I,KK) * DFNTRN(I,KK+1) * CT(I,KK+1)
            TCLD(I,KK) = TDCL1(I,KK) /  TDCL2(I,KK+1)
    650 END DO
651 END DO
!--------------------------------------------------------------------------------------------------
! THE FOLLOWING IS THE RECURSION RELATION FOR ALFA: THE REFLECTION COEFFICIENT FOR A SYSTEM 
! INCLUDING THE CLOUD IN QUESTION AND THE FLUX COMING OUT OF THE CLOUD SYSTEM INCLUDING ALL CLOUDS 
! BELOW THE CLOUD IN QUESTION. ALFAU IS ALFA WITHOUT THE REFLECTION OF THE CLOUD IN QUESTION
!--------------------------------------------------------------------------------------------------
    DO 660 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 660
!
        ALFA (I,1) = CR(I,1)
        ALFAU(I,1) = ZERO
660 END DO
!--------------------------------------------------- 
! AGAIN,EXCESSIVE CALCULATIONS-MAY BE CHANGED LATER!
!--------------------------------------------------- 
    DO 671 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 671
!
        DO 670 KK=2,KCLDS+1
            ALFAU(I,KK) =  TCLU(I,KK-1) * TCLU(I,KK-1) * ALFA(I,KK-1) / (ONE - TCLD(I,KK-1)       &
    &                   *  TCLD(I,KK-1) * ALFA(I,KK-1) *   CR(I,KK))
!
            ALFA (I,KK) = ALFAU(I,KK)   +   CR(I,KK)
    670 END DO
671 END DO
!--------------------------------------------------------------------------------------------------
! CALCULATE UFN AT CLOUD TOPS AND DFN AT CLOUD BOTTOMS NOTE THAT UFNCLU(I,KCLDS+1) GIVES THE UPWARD
! FLUX AT THE TOP OF THE HIGHEST REAL CLOUD (IF NCLDS(I)=KCLDS). IT GIVES THE FLUX AT THE TOP OF 
! THE ATMOSPHERE IF NCLDS(I) < KCLDS. IN THE FIRST CASE, TDCL1 EQUALS THE TRANSMISSION FCTN TO THE 
! TOP OF THE HIGHEST CLOUD, AS WE WANT. IN THE SECOND CASE, TDCL1=1, SO UFNCLU EQUALS ALFA. 
! THIS IS ALSO CORRECT.
!--------------------------------------------------------------------------------------------------
    DO 680 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 680
!
        UFNCLU(I,KCLDS+1) =  ALFA(I,KCLDS+1) * TDCL1(I,KCLDS+1)
        DFNCLU(I,KCLDS+1) = TDCL1(I,KCLDS+1)
680 END DO
!--------------------------------------------------------------------- 
! THIS CALCULATION IS THE REVERSE OF THE RECURSION RELATION USED ABOVE
!--------------------------------------------------------------------- 
    DO 691 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 691
!
        DO 690 KK=KCLDS,1,-1
            UFNCLU(I,KK) = UFNCLU(I,KK+1) * ALFAU(I,KK+1) / (ALFA(I,KK+1) * TCLU(I,KK))
            DFNCLU(I,KK) = UFNCLU(I,KK)   /  ALFA(I,KK)
    690 END DO
691 END DO
!
    DO 701 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 701
!
        DO 700 KK=1,KCLDS+1
            UFNTRN(I,KK) = UFNCLU(I,KK) * UFNTRN(I,KK)
            DFNTRN(I,KK) = DFNCLU(I,KK) * DFNTRN(I,KK)
    700 END DO
701 END DO
!-----------------------------------------------------------------  
! CASE OF KK=1( FROM THE GROUND TO THE BOTTOM OF THE LOWEST CLOUD)
!----------------------------------------------------------------- 
    DO 720 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 720
!
        J2 = KBTMSW(I,2)
        DO 710 K=J2,LP1
            UFN(I,K) = UFNTRN(I,1) * TTU(I,K)
            DFN(I,K) = DFNTRN(I,1) * TTD(I,K)
    710 END DO
720 END DO
!--------------------------- 
! REMAINING LEVELS (IF ANY!)
!--------------------------- 
    DO 760 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 760
!
        DO 755 KK=2,KCLDS+1
            J1 = KTOPSW(I,KK)
            J2 = KBTMSW(I,KK+1)
!
            IF (J1 == 1) GOTO 755
!
            DO 730 K=J2,J1
                UFN(I,K) = UFNTRN(I,KK) * TTU(I,K)
                DFN(I,K) = DFNTRN(I,KK) * TTD(I,K)
        730 END DO
!--------------------------------------------------------------------------------------------------  
! FOR THE THICK CLOUDS, THE FLUX DIVERGENCE THROUGH THE CLOUD LAYER IS ASSUMED TO BE CONSTANT. THE 
! FLUX DERIVATIVE IS GIVEN BY TEMPF (FOR THE UPWARD FLUX) AND TEMPG (FOR THE DOWNWARD FLUX).
!--------------------------------------------------------------------------------------------------
            J3 = KBTMSW(I,KK)
!
            IF ((J3-J1) > 1) THEN
                TEMPF = (UFNCLU(I,KK) - UFN(I,J3)) * DPCLD(I,KK-1)
                TEMPG = (DFNCLU(I,KK) - DFN(I,J3)) * DPCLD(I,KK-1)
                DO 740 K=J1+1,J3-1
                    UFN(I,K) = UFNCLU(I,KK) + TEMPF * (PP(I,K) - PPTOP(I,KK-1))
                    DFN(I,K) = DFNCLU(I,KK) + TEMPG * (PP(I,K) - PPTOP(I,KK-1))
            740 END DO
            END IF
!
    755 END DO
!
760 END DO
!
    DO 770 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 770
!
        DO 771 K=1,LP1
            DFSWC(I,K) = DFN(I,K) * DFNTOP(I,1)
            UFSWC(I,K) = UFN(I,K) * DFNTOP(I,1)
    771 END DO
770 END DO
!
    DO 780 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 780
!
        TMP1 (I) = ONE     - CCMAX(I)
        GDFVB(I) = TMP1(I) * GDFVB(I)
        GDFNB(I) = TMP1(I) * GDFNB(I)
        GDFVD(I) = TMP1(I) * GDFVD(I) + CCMAX(I) * DFSWC(I,LP1)
780 END DO
!-------------------------------------------------------------------------------------------------- 
! NOW OBTAIN FLUXES FOR THE NEAR IR BANDS. THE METHODS ARE THE SAME AS FOR THE VISIBLE BAND, EXCEPT
! THAT THE REFLECTION AND TRANSMISSION COEFFICIENTS ARE DIFFERENT, AS RAYLEIGH SCATTERING NEED NOT 
! BE CONSIDERED.
!-------------------------------------------------------------------------------------------------- 
    DO 991 N=2,NB
!
        DO 791 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 791
!
            DO 790 K=1,KCLDS+1
                CR(I,K) =              CRR(I,N,K)  * XAMT(I,K)
                CT(I,K) = ONE - (ONE - CTT(I,N,K)) * XAMT(I,K)
        790 END DO
    791 END DO
!
        IF (N == 2) THEN
!-------------------------------------------------------------------------------------------------- 
! THE WATER VAPOR TRANSMISSION FUNCTION FOR BAND 2 IS EQUAL TO THAT OF BAND 1
! (SAVED AS TTDB1,TTUB1)
!--------------------------------------------------------------------------------------------------
            DO 800 I=MYIS,MYIE
                KCLDS = NCLDS(I)
!
                IF (KCLDS == 0) GOTO 800
!
                DO 801 KK=2,LP1
                    TTD(I,KK) = TTDB1(I,KK) * TDCO2(I,KK)
            801 END DO
!
                DO 802 KK=1,L
                    TTU(I,KK) = TTUB1(I,KK) * TUCO2(I,KK)
            802 END DO
!
        800 END DO
!
        ELSE
!
            DO 810 I=MYIS,MYIE
                KCLDS = NCLDS(I)
!
                IF (KCLDS == 0) GOTO 810
!
                DO 811 KK=2,LP1
                    TTD(I,KK) = EXP(HM1EZ * MIN(FIFTY, ABCFF(N) * UD(I,KK))) * TDCO2(I,KK)
            811 END DO
!
                DO 812 KK=1,L
                    TTU(I,KK) = EXP(HM1EZ * MIN(FIFTY, ABCFF(N) * UR(I,KK))) * TUCO2(I,KK)
            812 END DO
!
        810 END DO
!
        END IF
!--------------------------------------------------------------------------------------------------
! AT THIS POINT,INCLUDE TTD(1),TTU(LP1), NOTING THAT TTD(1)=1 FOR ALL BANDS, AND THAT
! TTU(LP1)=TTD(LP1) FOR ALL BANDS.
!--------------------------------------------------------------------------------------------------
        DO 820 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 820
!
            TTU(I,LP1) = TTD(I,LP1)
            TTD(I,1)   = ONE
    820 END DO
!--------------------------------------------------------------------------------------------------
! FOR EXECUTION OF THE CLOUD LOOP, IT IS NECESSARY TO SEPARATE OUT TRANSMISSION FCTNS AT THE TOP 
! AND BOTTOM OF THE CLOUDS, FOR EACH BAND N. THE REQUIRED QUANTITIES ARE:
! TTD(I,KTOPSW(I,K),N)  K RUNS FROM 1 TO NCLDS(I)+1:
! TTD(I,KBTMSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+1:
! TTU(I,KTOPSW(I,K),N)  K RUNS FROM 1 TO NCLDS(I)+1:
! AND INVERSES OF THE ABOVE. THE ABOVE QUANTITIES ARE STORED IN TDCL1, TDCL2, TUCL1, AND DFNTRN,
! UFNTRN, RESPECTIVELY, AS THEY HAVE MULTIPLE USE IN THE PGM. FOR FIRST CLOUD LAYER (GROUND) 
! TDCL1, TUCL1 ARE KNOWN:
!--------------------------------------------------------------------------------------------------
        DO 830 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 830
!
            TDCL1 (I,1) =       TTD   (I,LP1)
            TUCL1 (I,1) =       TTU   (I,LP1)
            TDCL2 (I,1) =       TDCL1 (I,1)
            DFNTRN(I,1) = ONE / TDCL1 (I,1)
            UFNTRN(I,1) =       DFNTRN(I,1)
    830 END DO
!
        DO 841 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 841
!
            DO 840 KK=2,KCLDS+1
                TDCL1(I,KK) = TTD(I,KTOPSW(I,KK))
                TUCL1(I,KK) = TTU(I,KTOPSW(I,KK))
                TDCL2(I,KK) = TTD(I,KBTMSW(I,KK))
        840 END DO
!
    841 END DO
!
        DO 851 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 851
!
            DO 850 KK=2,KCLDS+1
                DFNTRN(I,KK) = ONE / TDCL1(I,KK)
                UFNTRN(I,KK) = ONE / TUCL1(I,KK)
        850 END DO
!
    851 END DO
!
        DO 861 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 861
!
            DO 860 KK=1,KCLDS
                TCLU(I,KK) = TDCL1(I,KK) * DFNTRN(I,KK+1) * CT(I,KK+1)
                TCLD(I,KK) = TDCL1(I,KK) /  TDCL2(I,KK+1)
        860 END DO
    861 END DO
!--------------------------------------------------------------------------------------------------
! THE FOLLOWING IS THE RECURSION RELATION FOR ALFA: THE REFLECTION COEFFICIENT FOR A SYSTEM 
! INCLUDING THE CLOUD IN QUESTION AND THE FLUX COMING OUT OF THE CLOUD SYSTEM INCLUDING ALL CLOUDS 
! BELOW THE CLOUD IN QUESTION.
!--------------------------------------------------------------------------------------------------
        DO 870 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 870
!
            ALFA (I,1) = CR(I,1)
            ALFAU(I,1) = ZERO
    870 END DO
!--------------------------------------------------- 
! AGAIN, EXCESSIVE CALCULATIONS-MAY BE CHANGED LATER
!---------------------------------------------------
        DO 881 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 881
!
            DO 880 KK=2,KCLDS+1
                ALFAU(I,KK) = TCLU (I,KK-1) * TCLU(I,KK-1) * ALFA(I,KK-1) / (ONE - TCLD(I,KK-1)   &
    &                       * TCLD (I,KK-1) * ALFA(I,KK-1) *   CR(I,KK))
!                
                ALFA (I,KK) = ALFAU(I,KK)   +   CR(I,KK)
        880 END DO
!
    881 END DO
!-------------------------------------------------------------------------------------------------- 
! CALCULATE UFN AT CLOUD TOPS AND DFN AT CLOUD BOTTOMS NOTE THAT UFNCLU(I,KCLDS+1) GIVES THE UPWARD
! FLUX AT THE TOP OF THE HIGHEST REAL CLOUD (IF NCLDS(I)=KCLDS). IT GIVES THE FLUX AT THE TOP OF 
! THE ATMOSPHERE IF NCLDS(I) < KCLDS. IT THE FIRST CASE, TDCL1  EQUALS THE TRANSMISSION FCTN TO THE
! TOP OF THE HIGHEST CLOUD, AS WE WANT. IN THE SECOND CASE, TDCL1=1, SO UFNCLU EQUALS ALFA. 
! THIS IS ALSO CORRECT.
!-------------------------------------------------------------------------------------------------- 
        DO 890 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF(KCLDS == 0) GOTO 890
!
            UFNCLU(I,KCLDS+1) = ALFA (I,KCLDS+1) * TDCL1(I,KCLDS+1)
            DFNCLU(I,KCLDS+1) = TDCL1(I,KCLDS+1)
    890 END DO
!
        DO 901 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 901
!
            DO 900 KK=KCLDS,1,-1
!--------------------------------------------------------------------------------------------------            
! FERRIER, 6/17/02:  EMERGENCY CHANGE TO ELIMINATE PROBLEMATIC FEATURES OF UNREALISTICALLY SMALL 
! CLOUD AMOUNTS
!--------------------------------------------------------------------------------------------------   
                DENOM = ALFA(I,KK+1) * TCLU(I,KK)
                IF (DENOM > 1.E-15) THEN
                    UFNCLU(I,KK) = UFNCLU(I,KK+1) * ALFAU(I,KK+1) / DENOM
                ELSE
                    UFNCLU(I,KK) = 0.
                END IF
!
                DFNCLU(I,KK) = UFNCLU(I,KK)   / ALFA (I,KK)
        900 END DO
    901 END DO
!----------------------------------------------------- 
! NOW OBTAIN DFN AND UFN FOR LEVELS BETWEEN THE CLOUDS
!----------------------------------------------------- 
        DO 911 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 911
!
            DO 910 KK=1,KCLDS+1
                UFNTRN(I,KK) = UFNCLU(I,KK) * UFNTRN(I,KK)
                DFNTRN(I,KK) = DFNCLU(I,KK) * DFNTRN(I,KK)
        910 END DO
    911 END DO
!
        DO 930 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 930
!
            J2 = KBTMSW(I,2)
            DO 920 K=J2,LP1
                UFN(I,K) = UFNTRN(I,1) * TTU(I,K)
                DFN(I,K) = DFNTRN(I,1) * TTD(I,K)
        920 END DO
    930 END DO
!
        DO 970  I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 970
!
            DO 965  KK=2,KCLDS+1
                J1 = KTOPSW(I,KK)
                J2 = KBTMSW(I,KK+1)
!
                IF (J1 == 1) GOTO 965
!
                DO 940 K=J2,J1
                    UFN(I,K) =     UFNTRN(I,KK) * TTU(I,K)
                    DFN(I,K) =     DFNTRN(I,KK) * TTD(I,K)
            940 END DO
!
                J3 = KBTMSW(I,KK)
                IF ((J3-J1) > 1) THEN
                    TEMPF   =     (UFNCLU(I,KK) - UFN(I,J3)) * DPCLD(I,KK-1)
                    TEMPG   =     (DFNCLU(I,KK) - DFN(I,J3)) * DPCLD(I,KK-1)
                    DO 950 K=J1+1,J3-1
                        UFN(I,K) = UFNCLU(I,KK) + TEMPF      * (PP(I,K) - PPTOP(I,KK-1))
                        DFN(I,K) = DFNCLU(I,KK) + TEMPG      * (PP(I,K) - PPTOP(I,KK-1))
                950 END DO
                END IF
!
        965 END DO
!
    970 END DO
!
        DO 980 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 980
!
            DO 981 K=1,LP1
                DFSWC(I,K) = DFSWC(I,K) + DFN(I,K) * DFNTOP(I,N)
                UFSWC(I,K) = UFSWC(I,K) + UFN(I,K) * DFNTOP(I,N)
        981 END DO
!
    980 END DO
!
        DO 990 I=MYIS,MYIE
            KCLDS = NCLDS(I)
!
            IF (KCLDS == 0) GOTO 990
!
            GDFND(I) = GDFND(I) + CCMAX(I) * DFN(I,LP1) * DFNTOP(I,N)
     990 END DO
!
991 END DO
!
    DO 992 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 992
!
        DO 993 K=1,LP1
            DFSWC(I,K) = TMP1(I) * DFSWL(I,K) + CCMAX(I) * DFSWC(I,K)
            UFSWC(I,K) = TMP1(I) * UFSWL(I,K) + CCMAX(I) * UFSWC(I,K)
    993 END DO
992 END DO
!
    DO 994 I=MYIS,MYIE
        KCLDS = NCLDS(I)
!
        IF (KCLDS == 0) GOTO 994
!
        DO 995 KK=1,LP1
            FSWC(I,KK) = UFSWC(I,KK) - DFSWC(I,KK)
    995 END DO
994 END DO
!
    DO 996 I=MYIS,MYIE
        KCLDS=NCLDS(I)
!
        IF (KCLDS == 0) GOTO 996
!
        DO 997 KK=1, L
            HSWC(I,KK) = RADCON * (FSWC(I,KK+1) - FSWC(I,KK)) / DP(I,KK)
    997 END DO
996 END DO
!
    RETURN
!
    END SUBROUTINE SWR93
