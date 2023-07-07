    SUBROUTINE GOSSIP
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE GOSSIP
!>
!> SUBPROGRAM: GOSSIP - EXCHANGE OF FIELDS BETWEEN PROCESSORS
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 97-08-30
!>
!> ABSTRACT:
!> GOSSIP EXCHANGES MANY FIELDS BETWEEN PROCESSORS IN ORDER TO FILL THE HALOES
!>
!> PROGRAM HISTORY LOG:
!> 97-05-??  MEYS       - ORIGINATOR
!> 98-10-23  BLACK      - MODIFIED FOR CURRENT VERSION OF MODEL
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: ACMCLD
!>              ACMCLH
!>              ACMPRE
!>              ACMRDL
!>              ACMRDS 
!>              ACMSFC
!>              CLDWTR
!>              CNVCLD
!>              CONTIN
!>              DYNAM
!>              F77KINDS 
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              PHYS
!>              PVRBLS
!>              RDPARM
!>              SOIL
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>              Z0EFFT
!>
!> DRIVER     : EBU
!>
!> CALLS      : EXCH
!>--------------------------------------------------------------------------------------------------
    USE ACMCLD
    USE ACMCLH
    USE ACMPRE
    USE ACMRDL
    USE ACMRDS
    USE ACMSFC
    USE CLDWTR
    USE CNVCLD
    USE CONTIN
    USE DYNAM
    USE EXCHM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE NHYDRO
    USE PARMETA
    USE PARMSOIL
    USE PARM_TBL
    USE PHYS
    USE PVRBLS
    USE RDPARM   , ONLY : LP1
    USE SOIL
    USE TEMPCOM
    USE TOPO
    USE VRBLS
    USE Z0EFFT
!
    IMPLICIT NONE   
!
    INCLUDE "EXCHM.h"
    INCLUDE "mpif.h"
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & J
!
#include "sp.h"   
!--------------- 
! THE NHB ARRAYS
!---------------    
    CALL EXCH(LMH    , 1    , 5 , 5)
    CALL EXCH(LMV    , 1    , 5 , 5)
    CALL EXCH(HBM2   , 1    , 5 , 5)
    CALL EXCH(HBM3   , 1    , 5 , 5)
    CALL EXCH(VBM2   , 1    , 5 , 5)
    CALL EXCH(VBM3   , 1    , 5 , 5)
    CALL EXCH(SM     , 1    , 5 , 5)
    CALL EXCH(SICE   , 1    , 5 , 5)
    CALL EXCH(HTM    , LM   , 5 , 5)
    CALL EXCH(VTM    , LM   , 5 , 5)
    CALL EXCH(DX     , 1    , 5 , 5)
    CALL EXCH(WPDAR  , 1    , 5 , 5)
    CALL EXCH(CPGFU  , 1    , 5 , 5)
    CALL EXCH(CURV   , 1    , 5 , 5)
    CALL EXCH(FCP    , 1    , 5 , 5)
    CALL EXCH(FDIV   , 1    , 5 , 5)
    CALL EXCH(FAD    , 1    , 5 , 5)
    CALL EXCH(F      , 1    , 5 , 5)
    CALL EXCH(DDMPU  , 1    , 5 , 5)
    CALL EXCH(DDMPV  , 1    , 5 , 5)
    CALL EXCH(GLAT   , 1    , 5 , 5)
    CALL EXCH(GLON   , 1    , 5 , 5)
    CALL EXCH(MXSNAL , 1    , 5 , 5)
    CALL EXCH(EPSR   , 1    , 5 , 5)
    CALL EXCH(TG     , 1    , 5 , 5)
    CALL EXCH(GFFC   , 1    , 5 , 5)
    CALL EXCH(SST    , 1    , 5 , 5)
    CALL EXCH(ALBASE , 1    , 5 , 5)
    CALL EXCH(HDAC   , 1    , 5 , 5)
    CALL EXCH(HDACV  , 1    , 5 , 5)
    CALL EXCH(IVGTYP , 1    , 5 , 5)
    CALL EXCH(ISLTYP , 1    , 5 , 5)
    CALL EXCH(ISLOPE , 1    , 5 , 5)
    CALL EXCH(VEGFRC , 1    , 5 , 5)
!-----------------------  
! THE RESTRT FILE ARRAYS
!----------------------- 
    CALL EXCH (OMGALF, LM   , 5 , 5)
    CALL EXCH (PD    , 1    , 5 , 5)
    CALL EXCH (RES   , 1    , 5 , 5)
    CALL EXCH (FIS   , 1    , 5 , 5)
    CALL EXCH (T     , LM   , 5 , 5)
    CALL EXCH (U     , LM   , 5 , 5)
    CALL EXCH (V     , LM   , 5 , 5)
    CALL EXCH (Q     , LM   , 5 , 5)
    CALL EXCH (Q2    , LM   , 5 , 5)
    CALL EXCH (CWM   , LM   , 5 , 5)
    CALL EXCH (TRAIN , LM   , 5 , 5)
    CALL EXCH (TCUCN , LM   , 5 , 5)
    CALL EXCH (RSWIN , 1    , 5 , 5)
    CALL EXCH (RSWOUT, 1    , 5 , 5)
    CALL EXCH (TG    , 1    , 5 , 5)
    CALL EXCH (Z0    , 1    , 5 , 5)
    CALL EXCH (AKMS  , 1    , 5 , 5)
    CALL EXCH (CZEN  , 1    , 5 , 5)
    CALL EXCH (AKHS  , 1    , 5 , 5)
    CALL EXCH (THS   , 1    , 5 , 5)
    CALL EXCH (QS    , 1    , 5 , 5)
    CALL EXCH (TWBS  , 1    , 5 , 5)
    CALL EXCH (QWBS  , 1    , 5 , 5)
    CALL EXCH (HBOT  , 1    , 5 , 5)
    CALL EXCH (CFRACL, 1    , 5 , 5)
    CALL EXCH (THZ0  , 1    , 5 , 5)
    CALL EXCH (QZ0   , 1    , 5 , 5)
    CALL EXCH (UZ0   , 1    , 5 , 5)
    CALL EXCH (VZ0   , 1    , 5 , 5)
    CALL EXCH (USTAR , 1    , 5 , 5)
    CALL EXCH (HTOP  , 1    , 5 , 5)
    CALL EXCH (CFRACM, 1    , 5 , 5)
    CALL EXCH (SNO   , 1    , 5 , 5)
    CALL EXCH (SI    , 1    , 5 , 5)
    CALL EXCH (CLDEFI, 1    , 5 , 5)
    CALL EXCH (RF    , 1    , 5 , 5)
    CALL EXCH (CUPPT , 1    , 5 , 5)
    CALL EXCH (CFRACH, 1    , 5 , 5)
    CALL EXCH (SOILTB, 1    , 5 , 5)
    CALL EXCH (SFCEXC, 1    , 5 , 5)
    CALL EXCH (SMSTAV, 1    , 5 , 5)
    CALL EXCH (SMSTOT, 1    , 5 , 5)
    CALL EXCH (GRNFLX, 1    , 5 , 5)
    CALL EXCH (PCTSNO, 1    , 5 , 5)
    CALL EXCH (RLWIN , 1    , 5 , 5)
    CALL EXCH (RADOT , 1    , 5 , 5)
    CALL EXCH (CZMEAN, 1    , 5 , 5)
    CALL EXCH (SIGT4 , 1    , 5 , 5)
    CALL EXCH (U00   , 1    , 5 , 5)
    CALL EXCH (LC    , 1    , 5 , 5)
    CALL EXCH (SR    , 1    , 5 , 5)
    CALL EXCH (PREC  , 1    , 5 , 5)
    CALL EXCH (ACPREC, 1    , 5 , 5)
    CALL EXCH (ACCLIQ, 1    , 5 , 5)
    CALL EXCH (CUPREC, 1    , 5 , 5)
    CALL EXCH (ACFRCV, 1    , 5 , 5)
    CALL EXCH (NCFRCV, 1    , 5 , 5)
    CALL EXCH (ACFRST, 1    , 5 , 5)
    CALL EXCH (NCFRST, 1    , 5 , 5)
    CALL EXCH (ACSNOW, 1    , 5 , 5)
    CALL EXCH (ACSNOM, 1    , 5 , 5)
    CALL EXCH (SSROFF, 1    , 5 , 5)
    CALL EXCH (BGROFF, 1    , 5 , 5)
    CALL EXCH (SFCSHX, 1    , 5 , 5)
    CALL EXCH (SFCLHX, 1    , 5 , 5)
    CALL EXCH (SUBSHX, 1    , 5 , 5)
    CALL EXCH (SNOPCX, 1    , 5 , 5)
    CALL EXCH (SFCUVX, 1    , 5 , 5)
    CALL EXCH (SFCEVP, 1    , 5 , 5)
    CALL EXCH (POTEVP, 1    , 5 , 5)
    CALL EXCH (ASWIN , 1    , 5 , 5)
    CALL EXCH (ASWOUT, 1    , 5 , 5)
    CALL EXCH (ASWTOA, 1    , 5 , 5)
    CALL EXCH (ALWIN , 1    , 5 , 5)
    CALL EXCH (ALWOUT, 1    , 5 , 5)
    CALL EXCH (ALWTOA, 1    , 5 , 5)
    CALL EXCH (SMC   , NSOIL, 5 , 5)
    CALL EXCH (CMC   , 1    , 5 , 5)
    CALL EXCH (STC   , NSOIL, 5 , 5)
    CALL EXCH (SH2O  , NSOIL, 5 , 5)
    CALL EXCH( ALBEDO, 1    , 5 , 5)
    CALL EXCH (PINT  , LM+1 , 5 , 5)
    CALL EXCH (Z     , LM+1 , 5 , 5)
    CALL EXCH (DWDT  , LM   , 5 , 5)
!
    DO J=MYJS_P4,MYJE_P4
        IVW(J) = IVWG(J + MY_JS_GLB - 1)
        IVE(J) = IVEG(J + MY_JS_GLB - 1)
        IHE(J) = IHEG(J + MY_JS_GLB - 1)
        IHW(J) = IHWG(J + MY_JS_GLB - 1)
    END DO
!
    CALL EXCH (ZEFFIJ, 4    , 5 , 5)
!
    RETURN
!
    END SUBROUTINE GOSSIP
