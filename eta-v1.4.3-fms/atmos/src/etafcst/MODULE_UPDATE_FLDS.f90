    MODULE UPDT
!>--------------------------------------------------------------------------------------------------
!> MODULE UPDT
!>
!> USE MODULES: F77KINDS
!>              PARMETA
!>
!> DRIVER     : READ_NHB
!>              SSTCH
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE PARMETA , ONLY : IM, JM ,IDIM1, IDIM2, JDIM1, JDIM2
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & LSST, LSSTMNTHLY  , LCO2
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & YR      , MON     , DAY     , UTC
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & YRSST      , MONSST     , DAYSST
!
    CHARACTER(LEN=10)                                                                           ::&
    & CHARDATE
!
    REAL   (KIND=R4KIND), DIMENSION(IM, JM)                                                     ::&
    & VARUPDATED
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, 600)                              ::&
    & SST_UP
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, 12)                               ::&
    & VEGFRM
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & INITSST  , INITCO2
!
    END MODULE UPDT


MODULE UPDATE_FLDS
PUBLIC FLDS_UPDATE_DRIVER
PUBLIC GETDATE
PUBLIC WEGHT_DAY_DIST
PUBLIC READ_SST
PUBLIC SST_UPDATE
PUBLIC SST_MNTHLY2DAY_UPDATE
PUBLIC READ_VGREEN
PUBLIC VEGUPDT

CONTAINS

  SUBROUTINE FLDS_UPDATE_DRIVER(SUB_ACT,FLD,CALC_SST,RECSST_REF)
    USE F77KINDS
    USE UPDT
    IMPLICIT NONE
    CHARACTER(*)        ,   INTENT(IN)                                                         ::&
    & SUB_ACT
    CHARACTER(*)        ,   INTENT(IN)                                                         ::&
    & FLD
    INTEGER(KIND=I4KIND),   INTENT(IN), OPTIONAL                                               ::&
    & CALC_SST
    INTEGER(KIND=I4KIND),   INTENT(IN), OPTIONAL                                               ::&
    & RECSST_REF
 
! NSST_UP=31
 SELECT CASE (TRIM(SUB_ACT))
   CASE ('READ')
     SELECT CASE (FLD)
       CASE ('SST')
         CALL READ_SST
       CASE ('VGREEN')
         CALL READ_VGREEN
       CASE ('CO2')
         CALL READ_CO2 	 
      END SELECT
   CASE('UPDATE')
     SELECT CASE (TRIM(FLD))
       CASE ('SST')
         CALL SST_UPDATE(CALC_SST,RECSST_REF)
       CASE ('SSTM2D')
         CALL SST_MNTHLY2DAY_UPDATE
       CASE ('VGREEN')
         CALL VGREEN_UPDATE
       CASE ('CO2CONCENTRATION')
         CALL CO2_CONC_UPDATE  	   
       CASE ('CO2COEFTRANSMISS')
         CALL CO2_COEF_UPDATE  	   
     END SELECT
  END SELECT  

  END SUBROUTINE FLDS_UPDATE_DRIVER
!
!
!
    SUBROUTINE GETDATE(IYR, IMO, IDY, IUTC, CYR, CMON, CDAY, CUTC)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE GETDATE
!>
!> SUBPROGRAM: GETDATE -
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 87-08-??  ?????     - ORIGINATOR
!> 18-01-15  LUCCI     - MODERNIZATION OF THE CODE, INCLUDING:
!>                       * F77 TO F90/F95
!>                       * INDENTATION & UNIFORMIZATION CODE
!>                       * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                       * DOCUMENTATION WITH DOXYGEN
!>                       * OPENMP FUNCTIONALITY
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
!> USE MODULES: F77KINDS
!>
!>
!> DRIVER     : SOLARD
!>              SSTCH
!>              SSTCH_GREGORIAN_CALENDAR
!>              SSTCH_360_DAY_CALENDAR
!>              VEGUPDT
!>              VEGUPDT_GREGORIAN_CALENDAR
!>              VEGUPDT_360_DAY_CALENDAR
!>              ZENITH
!>              ZENITH_GREGORIAN_CALENDAR
!>              ZENITH_360_DAY_CALENDAR
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & IYR     , IMO     , IDY     , IUTC
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(INOUT)       ::&
    & CYR     , CMON    , CDAY    , CUTC
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & YR      , MON     , DAY     , UTC
!
    INTEGER(KIND=I4KIND), DIMENSION(12)                                                         ::&
    & MONL
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DYSTEP
!
    YR  = IYR
    MON = IMO
    DAY = IDY
    UTC = IUTC
!
    DATA MONL /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
    DYSTEP = 0
!
    IF (MOD(YR,4) == 0) THEN
        MONL(2) = 29
    ELSE
        MONL(2) = 28
    END IF
!
    IF (UTC >= 24) THEN
        DYSTEP = DYSTEP + (UTC / 24)
        UTC    = (MOD(UTC, 24))
    END IF
!
    DAY = DAY + DYSTEP
!
    DO
        IF (DAY > MONL(MON)) THEN
            DAY = DAY - MONL(MON)
            MON = MON + 1
!
            IF (MON == 13) THEN
                MON = 1
                YR  = YR + 1
!
                IF (MOD(YR, 4) == 0) THEN
                    MONL(2) = 29
                ELSE
                    MONL(2) = 28
                END IF
!
            END IF
!
        ELSE
!
            EXIT
!
        END IF
!
    END DO
!
    CYR  = YR
    CMON = MON
    CDAY = DAY
    CUTC = UTC
!
      RETURN
!
      END SUBROUTINE GETDATE
!
!
!  SUBROUTINE WEGHT_DAY_DIST(YDAY,MA,MB,FA,FB)
  SUBROUTINE WEGHT_DAY_DIST(YDAY,MA,MB,RECFA,RECFB,FA,FB)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE WEGHT_DAY_DIST
!>
!> SUBPROGRAM: WEGHT_DAY_DIST - ?????
!> PROGRAMMER: ?????
!> ORG:?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-10-28  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> NTSD  -
!> DT    -
!>
!> OUTPUT ARGUMENT LIST:
!> FA    -
!> FB    -
!> MA    -
!> MB    -
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CTLBLK
!>              F77KINDS
!>
!> DRIVER     : 
!>
!> CALLS      : GETDATE
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE CTLBLK, ONLY: DT,IDAT,NTSD
    USE MPPCOM
    USE UPDT, ONLY: YRSST, MONSST, DAYSST
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & DAY     , MON     , YR      , UTC     , MA      , MB      , RECFA   , RECFB
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IYR     , IMO     , IDY    , IUTC   
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & YDAY    , ADD     , FA      , FB
!
    INTEGER(KIND=I4KIND), DIMENSION(12)                                                         ::&
    & MONL
!
    IYR  = IDAT(3)
    IMO  = IDAT(1)
    IDY  = IDAT(2)
   IF (MYPE == 0) PRINT*,"IN WEGHT_DAY_DIST",YRSST, MONSST, DAYSST
   IF (MYPE == 0) THEN
      PRINT*,"SST REFERENCE DATE ",YRSST, MONSST, DAYSST
      PRINT*,"IC DAY ",IDAT 
   ENDIF 
!
    IUTC = NTSD * DT / 3600
!
    CALL GETDATE(IYR, IMO, IDY, IUTC, YR, MON, DAY, UTC)
    IF (MYPE == 0) THEN
      PRINT*,"SST FCT DATE ",YR, MON, DAY 
    ENDIF 
   
!
    DATA MONL /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
    IF (MOD(YR,4) == 0) THEN
        MONL(2) = 29
    END IF

    YDAY = DAY
    MA   = MON - 1
    RECFA=((YR-YRSST)*12+(MON-MONSST))

    IF (MYPE == 0) THEN
      PRINT*,"RECFA: (",YR,"-",YRSST,")*12+(",MON,"-",MONSST,")= ",RECFA
    ENDIF 
!
    IF (MYPE == 0) PRINT*,"FCST DAY/MONTH ",YDAY,INT(MONL(MON)/2.0)
    IF ((YDAY >= (INT(MONL(MON)/2.0))) .OR. ((NTSD == 1) .AND. (YDAY > 15))) THEN
        MA = MON
        RECFA=RECFA+1
    ENDIF
    MB = MA + 1
    RECFB=RECFA+1
    IF (MYPE == 0) PRINT*,"REFERENCE MONTHS, RECFA and RECFB",MA,MB,RECFA,RECFB
!
    IF (MA <  1) MA = 12
    IF (MB > 12) MB =  1
!
    ADD = FLOAT(MONL(MA)) / 2.0 - 1.0
!
    IF (MA == MON) ADD = -ADD - 2.0
!
    FB = 2.0 * (YDAY + ADD) / FLOAT(MONL(MA) + MONL(MB))
    FA = 1.0 - FB
!
!
!
    RETURN
!
    END SUBROUTINE  WEGHT_DAY_DIST
!
!
    SUBROUTINE READ_SST
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE READ_SST
!>
!> SUBROUTINE: READ_SST- ?????
!> PROGRAMMER: ?????
!> ORG: ?????
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
!> USE MODULES: 
!>              CTLBLK
!>              F77KINDS
!>              MPPCOM
!>              PHYS
!>              PVRBLS
!>              SOIL
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : INIT
!>              INITS
!>
!> CALLS      : -----
!--------------------------------------------------------------------------------------------------
    USE CTLBLK
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA , ONLY : IM, JM
    USE UPDT, ONLY: SST_UP, YRSST, MONSST, DAYSST
    USE TEMPCOM, ONLY: TEMP1 
!
    INCLUDE "mpif.h"
!
!------------------
! DECLARE VARIABLES
!-------------------
!
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NSST_UP 
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ISST_UP , N
!
    ISST_UP = 90

!-------------------------------------------------
! ADAPTADO PARA A LEITURA DO ARQUIVO DE SST DIARIA
!-------------------------------------------------
!    OPEN(ISST_UP, FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
!
    READ(53,*) NSST_UP
    READ(53,*) YRSST, MONSST, DAYSST
    WRITE(6,*) "SST INIT: ",YRSST, MONSST, DAYSST
    CLOSE(53)

!---------------
! DISTRIBUTE SST_UP
!---------------
    DO N=1,NSST_UP
      IF (MYPE == 0) THEN
!        WRITE(6,*) "N, NSST_UP", N, NSST_UP
        READ(ISST_UP) TEMP1
      END IF
      CALL DSTRB(TEMP1, SST_UP, 1, 1, N)
    END DO
!
!------------------------------
! END OF SUBROUTINE READ_SST
!------------------------------
    RETURN
!
    END SUBROUTINE READ_SST
!
!
    SUBROUTINE READ_VGREEN
!>--------------------------------------------------------------------------------------------------
!> USE MODULES: USE CTLBLK
!>              USE F77KINDS
!>              USE GLB_TABLE
!>              USE MAPPINGS
!>              USE MPPCOM
!>		USE PARMETA , ONLY : IM, JM
!>		USE UPDT, ONLY: VEGFRM
!>		USE TEMPCOM, ONLY: TEMP1 
!>
!> DRIVER     : INIT
!>              INITS
!>
!> CALLS      : -----
!--------------------------------------------------------------------------------------------------
!
    USE CTLBLK
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA , ONLY : IM, JM
    USE UPDT, ONLY: VEGFRM
    USE TEMPCOM, ONLY: TEMP1 
!
    IMPLICIT NONE
!
    INTEGER(KIND=I8KIND)                                                                        ::&
    & N
    INCLUDE "mpif.h"
!------------------------------------------
! CHOU 02-02-2008 READ MONTHLY VEG GREENESS
!------------------------------------------
    OPEN(UNIT=9, FILE='VGREEN_12MO.dat', FORM='UNFORMATTED')
!
    DO N=1,12
        IF (MYPE == 0) THEN
            READ(9) TEMP1
        END IF
!
        CALL DSTRB(TEMP1, VEGFRM, 1, 12, N)
!
    END DO
!
    END SUBROUTINE READ_VGREEN
!
!
!
    SUBROUTINE READ_CO2
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE READ_CO2
!>
!> SUBROUTINE: READ_CO2 - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-09  JFP        - JOSE FERNNDO PESQUERO
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
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
!> USE MODULES: USE CTLBLK
!>              USE F77KIND
!>              USE MPPCOM
!>              USE PARMETA
!>              USE RDFSAV
!>
!> DRIVER     : INIT
!>              INITS
!>
!> CALLS      : -----
!--------------------------------------------------------------------------------------------------
!
    USE CTLBLK
    USE F77KINDS
    USE MPPCOM
    USE PARMETA
    USE RDFSAV
    USE UPDT
!
    IMPLICIT NONE
!    
!------------------------ 
! INCLUDE/SET PARAMETERS.
!------------------------ 
      INCLUDE "mpif.h"
!      
!------------------ 
! DECLARE VARIABLES
!------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & M      , IRTN      
!
    INTEGER ISTAT(MPI_STATUS_SIZE)
!
!------------------------------------------------ 
! ADAPTADO PARA A LEITURA DO ARQUIVO DE CO2 ANUAL 
! A DATA INICIAL DO ARQUIVO E 1875 
!------------------------------------------------- 
    OPEN (1, FORM='UNFORMATTED', FILE='CO2.bin', ACCESS='SEQUENTIAL') 
!---------------------------------
! DISTRIBUTE CO2 entre 1859 e 2100
!---------------------------------  
    DO M=1,242
!      
        IF (MYPE == 0) THEN
            READ(1) RACO2(M)
        END IF
!      
    END DO
!        
    CALL MPI_BCAST(RACO2(1), 242, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!-----------------------------------------     
! LE O INDICE DO VETOR PARA O TEMPO INCIAL
!----------------------------------------- 
    OPEN (8, FORM='FORMATTED', FILE='TEMPO_INIT_CO2.txt')
!
    READ(8,*) INITCO2
!---------------------------
! END OF SUBROUTINE READ_CO2     
!---------------------------
    RETURN
!
    END SUBROUTINE READ_CO2
!
!
!
  SUBROUTINE SST_UPDATE(CALC_SST,RECSST_REF)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SSTSTP
!>
!> SUBPROGRAM: SSTSTP - ?????
!> PROGRAMMER: ?????
!> ORG:?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-10-28  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> CALC_SSTUP (UNIT: DAY)
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CTLBLK
!>              F77KINDS
!>              PARMETA
!>              PHYS
!>              UPDT
!>
!> DRIVER     : EBU
!>
!>--------------------------------------------------------------------------------------------------
    USE CTLBLK
    USE F77KINDS
    USE MPPCOM
    USE PARMETA
    USE PHYS, ONLY: SST
    USE UPDT, ONLY: SST_UP
!
   INCLUDE "mpif.h"
!

    INTEGER(KIND=I4KIND)                                                                  :: NREC_SST_UP
    INTEGER(KIND=I4KIND),   INTENT(IN)                                                    :: CALC_SST
    INTEGER(KIND=I4KIND),   INTENT(IN)                                                    :: RECSST_REF
!
    NREC_SST_UP=(NTSD/CALC_SST)+RECSST_REF
    IF (MYPE == 0) write(0,*)"SUM SST BEFORE",sum(SST)
    IF (MYPE == 0) write(0,*)"SST UPDATE",NTSD,CALC_SST,NREC_SST_UP
    SST(:,:) =  SST_UP(:,:,NREC_SST_UP)
    IF (MYPE == 0) write(0,*)"SUM SST AFTER",sum(SST)
!
    RETURN
!
    END SUBROUTINE SST_UPDATE
!
!
!
 SUBROUTINE SST_MNTHLY2DAY_UPDATE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SST_MNTHLY2DAY_UPDATE
!>
!> SUBPROGRAM: SST - ?????
!> PROGRAMMER: ?????
!> ORG:?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-10-28  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> CALC_SSTUP (UNIT: DAY)
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CTLBLK
!>              F77KINDS
!>              PARMETA
!>              PHYS
!>              UPDT
!>
!> DRIVER     : EBU
!>
!>--------------------------------------------------------------------------------------------------
    USE CTLBLK, ONLY: DT,NTSD
    USE MPPCOM
    USE PARMETA
    USE PHYS, ONLY: SST
    USE UPDT, ONLY: SST_UP
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & MA      , MB      , RECFA      , RECFB
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & YDAY    , FA      , FB
!
    CALL WEGHT_DAY_DIST(YDAY,MA,MB,RECFA,RECFB,FA,FB)

    IF (MYPE == 0) PRINT*,"SST", NTSD, RECFA, RECFB, FA, FB
    IF (MYPE == 0) PRINT*,"SUM SST BEFORE",sum(SST)
!
    IF ((NTSD == 1) .AND. (INT(YDAY) == 15)) THEN
        SST(:,:) = SST_UP(:,:,RECFA)
    ELSE
        SST(:,:) = FA * SST_UP(:,:,RECFA) + FB * SST_UP(:,:,RECFB)
    END IF
!
    IF (MYPE == 0) PRINT*,"SUM SST AFTER",sum(SST)
    RETURN
!
    END SUBROUTINE SST_MNTHLY2DAY_UPDATE!
!
!
!
    SUBROUTINE VGREEN_UPDATE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VEG_UPDATE
!>
!> SUBPROGRAM: VEG_UPDATE - ?????
!> PROGRAMMER: ?????
!> ORG:?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-10-28  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SSTUPSTEP -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CTLBLK
!>              F77KINDS
!>              LOOPS
!>              MASKS
!>              PARMETA
!>              PARMSOIL
!>              PARM_TBL
!>              SOIL
!>
!> DRIVER     : EBU
!>
!> CALLS      : GETDATE
!>--------------------------------------------------------------------------------------------------
!
    USE F77KINDS
    USE CTLBLK, ONLY: DT,NTSD
    USE MPPCOM
    USE PARMETA , ONLY : IM, JM
    USE SOIL, ONLY: VEGFRC
    USE UPDT, ONLY: VEGFRM  
    INTEGER(KIND=I4KIND)                                                                        ::&
    & MA      , MB      , RECFA      , RECFB
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & YDAY    , FA      , FB
!
    CALL WEGHT_DAY_DIST(YDAY,MA,MB,RECFA,RECFB,FA,FB)

    IF (MYPE == 0) PRINT*,"VEGFRM", NTSD, MA, MB, FA, FB
!
    IF ((NTSD == 1) .AND. (INT(YDAY) == 15)) THEN
        VEGFRC(:,:) = VEGFRM(:,:,MA)
    ELSE
        VEGFRC(:,:) = (FA * VEGFRM(:,:,MA) + FB * VEGFRM(:,:,MB))
    END IF
!
         WHERE(VEGFRC>1.) VEGFRC = 1.
    RETURN
!
    END SUBROUTINE VGREEN_UPDATE
!
!
!
    SUBROUTINE CO2_CONC_UPDATE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE CO2_UPDATE
!>
!> SUBPROGRAM: CO2_UPDATE - ?????
!> PROGRAMMER: ?????
!> ORG:?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-10-28  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SSTUPSTEP -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CTLBLK
!>              F77KINDS
!>              LOOPS
!>              MASKS
!>              PARMETA
!>
!> DRIVER     : EBU
!>
!> CALLS      : CONRAD2
!>--------------------------------------------------------------------------------------------------
!
    USE F77KINDS
    USE CTLBLK, ONLY: DT,NTSD
    USE MPPCOM
    USE PARMETA , ONLY : IM, JM 
    USE RDFSAV
    USE UPDT
!
    IMPLICIT NONE
!     
!------------------------ 
! INCLUDE/SET PARAMETERS.
!------------------------ 
      INCLUDE "mpif.h"
!      
!------------------ 
! DECLARE VARIABLES
!------------------

!-------------------------
! CO2 CONCENTRATION UPDATE
!-------------------------
!
        INITCO2 = INITCO2 + 1
        RCO2    = RACO2(INITCO2)
        IF (MYPE == 0) print*, "INITCO2 ", INITCO2, RCO2, NTSD

    END SUBROUTINE CO2_CONC_UPDATE
!
!    
!
    SUBROUTINE CO2_COEF_UPDATE
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE CO2_UPDATE
!>
!> SUBPROGRAM: CO2_UPDATE - ?????
!> PROGRAMMER: ?????
!> ORG:?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> 93-10-28  ?????      - ORIGINATOR
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> SSTUPSTEP -
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: CTLBLK
!>              F77KINDS
!>              LOOPS
!>              MASKS
!>              PARMETA
!>
!> DRIVER     : EBU
!>
!> CALLS      : CONRAD2
!>--------------------------------------------------------------------------------------------------
!     
    USE F77KINDS
    USE CTLBLK, ONLY: DT,NTSD
    USE MPPCOM
    USE PARMETA , ONLY : IM, JM 
    USE RDFSAV
    USE UPDT
!
    IMPLICIT NONE
!     
!------------------------ 
! INCLUDE/SET PARAMETERS.
!------------------------ 
      INCLUDE "mpif.h"
!      
!------------------ 
! DECLARE VARIABLES
!------------------

!-------------------------------------------------------------
!CO2 TRANSMISSIVITY COEFF. UPDATE FOR CLIMATE CHANGE SCENARIOS 
!-------------------------------------------------------------
!
        CALL CONRAD2(67)
	IF (MYPE .eq. 0) THEN
	    OPEN(UNIT=44,FILE='LE_CO2.txt',FORM='formatted',STATUS='unknown',POSITION='append')
	    WRITE (44,*) "Enter CONRAD2, NTSD= ", NTSD, "Year= ", YR, "MYPE= ", MYPE   
	    PRINT*, ' time to update CO2 transm fcns again. year ', YR 
	    CLOSE (44) 
        ENDIF

    END SUBROUTINE CO2_COEF_UPDATE
!
!        
END MODULE UPDATE_FLDS
