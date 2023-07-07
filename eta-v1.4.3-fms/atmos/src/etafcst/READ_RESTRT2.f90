    SUBROUTINE READ_RESTRT2
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE READ_RESTRT2
!>
!> SUBPROGRAM:READ_RESTRT2 - READ MULTIPLE SMALL RESTRT FILES
!> PROGRAMMER: BLACK
!> ORG: W/NP2
!> DATE: 99-09-01
!>
!> ABSTRACT:
!> READ_RESTRT2 READS IN QUANTITIES FROM THE SMALL RESTRT FILES WHICH WERE PREVIOUSLY WRITTEN BY
!> INDIVIDUAL NODES
!>
!> PROGRAM HISTORY LOG:
!> 99-09-01  BLACK      - REWRITTEN ROM READ_RESTRT
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
!>              BOCO
!>              CLDWTR
!>              CNVCLD
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPOT
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              OUTFIL
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              PHYS
!>              PRFHLD
!>              PVRBLS
!>              SOIL
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : INIT
!>              INITS
!>
!> CALLS      : GETENV
!>--------------------------------------------------------------------------------------------------
    USE ACMCLD
    USE ACMCLH
    USE ACMPRE
    USE ACMRDL
    USE ACMRDS
    USE ACMSFC
    USE BOCO
    USE CLDWTR
    USE CNVCLD
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPOT
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE OUTFIL
    USE PARMETA
    USE PARMSOIL
    USE PARM_TBL
    USE PHYS
    USE PRFHLD 
    USE PVRBLS
    USE SOIL
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: JMP1 = JM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: IMT  =  2 * IM - 1
!------------------  
! DECLARE VARIABLES
!------------------  
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & RUNB
!
    CHARACTER(LEN=4)                                                                            ::&
    & RESTHR
!
    CHARACTER(LEN=32)                                                                           ::&
    & LABEL
    CHARACTER(40)                                                                               ::&
    & CONTRL  , FILALL  , FILMST  , FILTMP  , FILTKE  , FILUNV  , FILCLD  , FILRAD  , FILSFC
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PSLP
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM, NSOIL)                                              ::&
    & TEMPSOIL
!
    INTEGER(KIND=I4KIND), DIMENSION(3)                                                          ::&
    & IDATB     
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!
#ifdef DP_REAL
!
    LOGICAL(KIND=L8KIND)                                                                        ::&
    & RUNX    , FIRSTX
!
    INTEGER(KIND=I8KIND), DIMENSION(3)                                                          ::&
    & IDATX
!
    INTEGER(KIND=I8KIND)                                                                        ::&
    & IHRSTX  , NTSDX   , IOUTX   , NSHDEX
!
    INTEGER(KIND=I8KIND), DIMENSION(IM, JM)                                                     ::&
    & ITEMPX  , ITEMP2X
!
#endif
!
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K       , N       , IER   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TSTART  , PDOMG   , RESOMG
!--------------------------------------------
! READ INITIAL CONDITIONS FROM RESTART FILES.
!--------------------------------------------
    IF (MYPE == 0) WRITE(LIST,*) 'INIT:  READ RESTART FILES'
!----------------------------- 
! CREATE NAME FOR RESTART FILE
!-----------------------------  
    TSTART = NINT(NSTART * DT / 3600.)
!-------------------- 
! NO LONGER HARDWIRED
!-------------------- 
    IF (NSTART == 0) THEN
        ITAG = 3
    ELSE
        ITAG = TSTART
    END IF
!
    IF (MYPE == 0) WRITE(6,*)'ITAG IN READ_RESTRT2 ', ITAG, TSTART, DT, NSTART, NTSD
!-------------------- 
! NO LONGER HARDWIRED
!-------------------- 
    CALL GETENV("TMMARKB", RESTHR)
!
    IF (RESTHR == '    ') THEN
        WRITE(RSTFIL,50) ITAG,MYPE
        50 FORMAT('RESTRT',I2.2,'.',I3.3)
    ELSE
        WRITE(RSTFIL,55)ITAG,MYPE,RESTHR
        55 FORMAT('RESTRT',I2.2,'.',I3.3,'.',A4)
    END IF
!---------------------------
! OPEN UNIT TO RESTART FILE.
!---------------------------
    LRSTRT = 8
!
    CLOSE(LRSTRT)
    OPEN (UNIT=LRSTRT, FILE=RSTFIL, FORM='UNFORMATTED', IOSTAT=IER)
    IF (IER /= 0) WRITE(LIST,*) ' LRSTRT OPEN UNIT ERROR IER=', IER
!
#ifdef DP_REAL
!
    READ(LRSTRT) RUNX, IDATX, IHRSTX, NTSDX, LABEL
    RUN   = RUNX
    IDAT  = IDATX
    IHRST = IHRSTX
    NTSD  = NTSDX
!
#else
!
    READ(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL
!
#endif
!
    NTSD = MAX(NTSD-1,0)
    READ(LRSTRT) PDOMG, RESOMG
!
    DO K=1,LM
        READ(LRSTRT) ((OMGALF(I,J,K),I=MYIS,MYIE),J=MYJS,MYJE)
    END DO
!
#ifdef DP_REAL
!
    READ(LRSTRT) RUNX, IDATX, IHRSTX, NTSDX, LABEL, FIRSTX, IOUTX, NSHDEX
    RUN   = RUNX
    IDAT  = IDATX
    IHRST = IHRSTX
    NTSD  = NTSDX
    FIRST = FIRSTX
    IOUT  = IOUTX
    NSHDE = NSHDEX
!
#else
!
    READ(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL, FIRST, IOUT, NSHDE
!
#endif
!
    NTSD  = MAX(NTSD-1,0)
    FIRST = .TRUE. 
!
    READ(LRSTRT) ((PD(I,J),I=MYIS,MYIE),J=MYJS,MYJE), ((RES(I,J),I=MYIS,MYIE),J=MYJS,MYJE),       &
    &           ((FIS(I,J),I=MYIS,MYIE),J=MYJS,MYJE)
!
    READ(LRSTRT)
!---------------------- 
! PRIMARY 3-D VARIABLES
!---------------------- 
    DO K=1,LM
        READ(LRSTRT) ((T    (I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT) ((Q    (I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT) ((U    (I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT) ((V    (I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT) ((Q2   (I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT)
        READ(LRSTRT) ((CWM  (I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT) ((TRAIN(I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT) ((TCUCN(I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
    END DO
!
#ifdef DP_REAL
!
    READ(LRSTRT) RUNX, IDATX, IHRSTX, NTSDX, LABEL,                                               &
    &            ((RSWIN (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((RSWOUT(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((TG    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((Z0    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((AKMS  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CZEN  (I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    RUN   = RUNX
    IDAT  = IDATX
    IHRST = IHRSTX
    NTSD  = NTSDX
!
#else
!
    READ(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL,                                                   &
    &            ((RSWIN (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((RSWOUT(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((TG    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((Z0    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((AKMS  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CZEN  (I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
#endif
!
    NTSD = MAX(NTSD-1,0)
!
    READ(LRSTRT) ((AKHS  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((THS   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((QS    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((TWBS  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((QWBS  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((HBOT  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CFRACL(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((THZ0  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((QZ0   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((UZ0   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((VZ0   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((USTAR (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((HTOP  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CFRACM(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((SNO   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SI    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CLDEFI(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((RF    (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((PSLP  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CUPPT (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CFRACH(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((SOILTB(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SFCEXC(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SMSTAV(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SMSTOT(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((GRNFLX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((PCTSNO(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((RLWIN (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
                 ((RADOT (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
                 ((CZMEAN(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
                 ((SIGT4 (I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!----------------------------------------------------------------
! ((U00(I,J),I=MYIS,MYIE),J=MYJS,MYJE) NOT USED IN NEW G/S SCHEME
! UL                                   NOT USED IN NEW G/S SCHEME
! ((LC(I,J),I=MYIS,MYIE),J=MYJS,MYJE)  NOT USED IN NEW G/S SCHEME 
!----------------------------------------------------------------
    READ(LRSTRT) ((U00(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                          &
    &            UL,                                                                              & 
    &            ((LC (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                          &
    &            ((SR (I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
#ifdef DP_REAL
!
    READ(LRSTRT) RUNX, IDATX, IHRSTX, NTSDX, LABEL,                                               & 
    &            ((PREC  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ACPREC(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ACCLIQ(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CUPREC(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    RUN   = RUNX
    IDAT  = IDATX
    IHRST = IHRSTX
    NTSD  = NTSDX
!
#else
!
    READ(LRSTRT) RUN, IDAT, IHRST, NTSD, LABEL,                                                   &
    &            ((PREC  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ACPREC(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ACCLIQ(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((CUPREC(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
#endif
!
    NTSD = MAX(NTSD-1,0)
!
    READ(LRSTRT) ((ACFRCV(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((NCFRCV(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ACFRST(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((NCFRST(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((ACSNOW(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ACSNOM(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SSROFF(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((BGROFF(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((SFCSHX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SFCLHX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SUBSHX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SNOPCX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SFCUVX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((SFCEVP(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((POTEVP(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ((ASWIN (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ASWOUT(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ASWTOA(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ALWIN (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ALWOUT(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((ALWTOA(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) ARDSW, ARDLW, ASRFC, AVRAIN, AVCNVC
!
    READ(LRSTRT) ((TH10  (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((Q10   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((U10   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((V10   (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((TSHLTR(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((QSHLTR(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
    &            ((PSHLTR(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                       &
!--------- 
! SM v100M
!--------- 
    &            ((TH100(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                        &
    &            ((Q100 (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                        &
    &            ((U100 (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                        &
    &            ((V100 (I,J), I=MYIS,MYIE), J=MYJS,MYJE)
!
    READ(LRSTRT) (((SMC  (I,J,N), I=MYIS,MYIE), J=MYJS,MYJE), N=1,NSOIL)
    READ(LRSTRT) ((CMC   (I,J)  , I=MYIS,MYIE), J=MYJS,MYJE)
    READ(LRSTRT) (((STC  (I,J,N), I=MYIS,MYIE), J=MYJS,MYJE), N=1,NSOIL)
    READ(LRSTRT) (((SH2O (I,J,N), I=MYIS,MYIE), J=MYJS,MYJE), N=1,NSOIL)
    READ(LRSTRT) ((ALBEDO(I,J)  , I=MYIS,MYIE), J=MYJS,MYJE)
!--------------------------------------------------------------------------------
! IF FORECAST IS NOT BEGINNING AT TIME 0 THEN WE MUST READ ADDITIONAL INFORMATION
!--------------------------------------------------------------------------------
    IF (NINT(TSTART) /= 0) THEN
        READ(LRSTRT) ((POTFLX(I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                   &
    &                ((TLMIN (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                   &
    &                ((TLMAX (I,J), I=MYIS,MYIE), J=MYJS,MYJE),                                   &
    &                ACUTIM, ARATIM, APHTIM
!    
        DO K=1,LM
            READ(LRSTRT)((RSWTT(I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
            READ(LRSTRT)((RLWTT(I,J,K), I=MYIS,MYIE), J=MYJS,MYJE)
        END DO
!
        READ(LRSTRT)((CNVBOT(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT)((CNVTOP(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT)((RSWTOA(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
        READ(LRSTRT)((RLWTOA(I,J), I=MYIS,MYIE), J=MYJS,MYJE)
    END IF
!-----------------------------------------------------------------------
! CALL RADIATION TO OBTAIN THE SHORT AND LONGWAVE TEMPERATURE TENDENCIES
!
! CALL RADTN
!
! DONE READING RESTART FILES.
!
! END OF SUBROUTINE READ_RESTRT2
!-----------------------------------------------------------------------
    IF (MYPE == 0) THEN
        WRITE(LIST,*) 'INIT:  EXIT READ_RESTRT2'
        WRITE(LIST,*) ' '
    END IF
!
    RETURN
!
    END SUBROUTINE READ_RESTRT2
