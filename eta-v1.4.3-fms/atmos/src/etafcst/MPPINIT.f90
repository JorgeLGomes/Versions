    SUBROUTINE MPPINIT
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE MPPINIT
!>
!> SUBROUTINE: MPPINIT - ?????
!> PROGRAMMER: BLACK
!> ORG: W/NP22
!> DATE: 98-10-28
!>
!> ABSTRACT:
!> MPPINIT DETERMINES ALL RELEVANT VALUES FOR DIMENSIONS OF THE DISTRIBUTED SUBDOMAINS AND THEIR
!> HALOES.
!>
!> PROGRAM HISTORY LOG:
!> 97-??-??  MEYS     - ORIGINATOR
!> 97-??-??  BLACK    - CHANGES MADE FOR CLARITY
!> 98-10-29  BLACK    - REWRITTEN FOR CLARITY
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : EBU
!>
!> CALLS      : INDTABLE
!>--------------------------------------------------------------------------------------------------
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
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), DIMENSION(MPI_STATUS_SIZE)                                            ::&
    & ISTAT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IICHUNK           , IPE               , JNCHUNKS          ,  NCHUNKS   ,                    &
    & ICHUNK            , JCHUNK            ,                                                     &
    & MY_IE_CALC        , MY_JE_CALC        ,                                                     &
    & MY_IS_CALC        , MY_JS_CALC        ,                                                     &
    & I                 , J                 ,                                                     &
    & ICHUNK_CALC       , JCHUNK_CALC       ,                                                     &
    & IS_INC1_BND       , IS_INC2_BND       ,                                                     &
    & IE_INC1_BND       , IE_INC2_BND       ,                                                     &
    & JS_INC1_BND       , JS_INC2_BND       , JS_INC3_BND       , JS_INC4_BND       ,             &
    & JS_INC5_BND       ,                                                                         &
    & JE_INC1_BND       , JE_INC2_BND       , JE_INC3_BND       , JE_INC4_BND       ,             &
    & JE_INC5_BND       ,                                                                         &
    & ILOC              , JLOC              ,                                                     &
    & MYI               , MYJ               ,                                                     &
    & MYJS1_P5          , MYJS2_P5          ,                                                     &
    & MYJE1_P5          , MYJE2_P5          ,                                                     &
    & MYJE3_P1          , MYJS3_P1          ,                                                     &
    & MYJE3_P5          , MYJE4_P5
!--------------------------------------------------------------------------------------------------  
! INPES AND JNPES ARE THE NUMBER OF PEs REQUESTED IN X AND Y.
! ICHUNK AND JCHUNK ARE THE FIRST GUESS OF THE NUMBER OF I's AND J's IN EACH SUBDOMAIN OBTAINED BY
! SIMPLY DIVIDING THE GLOBAL DIMENSIONS BY THE NUMBER OF PEs REQUESTED IN EACH DIRECTION.
!-------------------------------------------------------------------------------------------------- 
    ICHUNK  = IM / INPES
    JCHUNK  = JM / JNPES
    IICHUNK = ICHUNK + 1
!-------------------------------------------------------------------------------------------------- 
! COMPUTE THE GLOBAL START AND END INDEX VALUES FOR I (MY_IS_GLB,MY_IE_GLB) AND 
! J (MY_JS_GLB,MY_JE_GLB) ON EACH PE.
! IN GENERAL, THE NUMBER OF POINTS IN EACH DIRECTION WILL NOT DIVIDE EVENLY WITH INPES AND JNPES. 
! THE LOGIC BELOW GIVES ONE EXTRA POINT TO AS MANY OF THE EARLIEST PEs IN EACH DIRECTION AS IT 
! TAKES TO USE UP THE REMAINDER POINTS (ITAIL AND JTAIL, WHICH ARE COMPUTED IN PARMETA).
!-------------------------------------------------------------------------------------------------- 
    IPE = 0
    MY_JS_CALC = 1
    JNCHUNKS = 0
!
    DO J=1,JNPES
        JCHUNK_CALC = JCHUNK
!
        IF (J <= JTAIL) JCHUNK_CALC = JCHUNK + 1
        JNCHUNKS   = JNCHUNKS + JCHUNK_CALC
        MY_JE_CALC = JNCHUNKS
        MY_IS_CALC = 1
        NCHUNKS    = 0
!    
        DO I=1,INPES
            ICHUNK_CALC = ICHUNK
!
            IF (I <= ITAIL) ICHUNK_CALC = ICHUNK + 1
            NCHUNKS    = NCHUNKS + ICHUNK_CALC
            MY_IE_CALC = NCHUNKS
!
            IF (MYPE == IPE) THEN
                MY_IS_GLB = MY_IS_CALC
                MY_IE_GLB = MY_IE_CALC
                MY_JS_GLB = MY_JS_CALC
                MY_JE_GLB = MY_JE_CALC
            END IF
!
            MY_IS_CALC = MY_IE_CALC + 1
           IPE = IPE + 1
 !
        END DO
!    
        MY_JS_CALC = MY_JE_CALC + 1
!
    END DO
!------------------------------------------------------------------------------------------
! ILPADx IS THE INCREMENT INTO THE LEFT HALO OF A SUBDOMAIN.
! IRPADx IS THE INCREMENT INTO THE RIGHT HALO OF A SUBDOMAIN.
! ILPADx IS ALWAYS 0 FOR SUBDOMAINS ALONG THE WEST GLOBAL BOUNDARY.
! IRPADx IS ALWAYS 0 FOR SUBDOMAINS ALONG THE EAST GLOBAL BOUNDARY.
!
! ILCOL IS A FLAG TELLING WHETHER OR NOT A SUBDOMAIN IS ON THE WEST (LEFT) GLOBAL BOUNDARY.
!
! IS_INCx_BND AND IE_INCx_BND ARE INCREMENTS FROM THE LOCAL STARTING OR ENDING 
! I VALUE AWAY FROM THE LOCAL BOUNDARY INTO THE SURBDOMAIN.
! THEY ARE NONZERO ONLY FOR SUBDOMAINS ON THE WESTERN AND EASTERN GLOBAL BOUNDARIES.
!------------------------------------------------------------------------------------------
    ILPAD1 = 1
    ILPAD2 = 2
    ILPAD3 = 3
    ILPAD4 = 4
    ILPAD5 = 5
!
    IRPAD1 = 1
    IRPAD2 = 2
    IRPAD3 = 3
    IRPAD4 = 4
    IRPAD5 = 5
!
    ILCOL  = 0
    IRCOL  = 0
!
    IS_INC1_BND = 0
    IS_INC2_BND = 0
    IE_INC1_BND = 0
    IE_INC2_BND = 0
!-----------------------  
! WESTERNMOST SUBDOMAINS
!-----------------------  
    IF (MOD(MYPE,INPES) == 0) THEN             
        ILPAD1 = 0
        ILPAD2 = 0
        ILPAD3 = 0
        ILPAD4 = 0
        ILPAD5 = 0
        ILCOL  = 1
        IS_INC1_BND = 1
        IS_INC2_BND = 2
    END IF
!-----------------------
! EASTERNMOST SUBDOMAINS
!-----------------------
    IF (MOD(MYPE,INPES) == INPES-1) THEN       
        IRPAD1 = 0
        IRPAD2 = 0
        IRPAD3 = 0
        IRPAD4 = 0
        IRPAD5 = 0
        IRCOL  = 1
        IE_INC1_BND = 1
        IE_INC2_BND = 2
        MY_IE_GLB   = IM
    END IF
!------------------------------------ 
! NOW DO THE SAME FOR THE J DIRECTION
!------------------------------------ 
    JBPAD1 = 1
    JBPAD2 = 2
    JBPAD3 = 3
    JBPAD4 = 4
    JBPAD5 = 5
!
    JTPAD1 = 1
    JTPAD2 = 2
    JTPAD3 = 3
    JTPAD4 = 4
    JTPAD5 = 5
!
    IBROW  = 0
    ITROW  = 0
!
    JS_INC1_BND = 0
    JS_INC2_BND = 0
    JS_INC3_BND = 0
    JS_INC4_BND = 0
    JS_INC5_BND = 0
!
    JE_INC1_BND = 0
    JE_INC2_BND = 0
    JE_INC3_BND = 0
    JE_INC4_BND = 0
    JE_INC5_BND = 0
!------------------------  
! SOUTHERNMOST SUBDOMAINS
!------------------------ 
    IF (MYPE/INPES == 0) THEN              
        JBPAD1 = 0
        JBPAD2 = 0
        JBPAD3 = 0
        JBPAD4 = 0
        JBPAD5 = 0
        IBROW  = 1
!
        JS_INC1_BND = 1
        JS_INC2_BND = 2
        JS_INC3_BND = 3
        JS_INC4_BND = 4
        JS_INC5_BND = 5
    END IF
!------------------------  
! NORTHERNMOST SUBDOMAINS
!------------------------ 
    IF (MYPE/INPES == JNPES-1) THEN        
        JTPAD1 = 0
        JTPAD2 = 0
        JTPAD3 = 0
        JTPAD4 = 0
        JTPAD5 = 0
        ITROW  = 1
!
        JE_INC1_BND = 1
        JE_INC2_BND = 2
        JE_INC3_BND = 3
        JE_INC4_BND = 4
        JE_INC5_BND = 5
        MY_JE_GLB   = JM
    END IF
!---------------------------------------------------------------- 
! THE FOLLOWING ARE THE LOCAL LIMITS OF I AND J IN EACH SUBDOMAIN
!---------------------------------------------------------------- 
    MY_IS_LOC = 1
    MY_IE_LOC = MY_IE_GLB - MY_IS_GLB + 1
    MY_JS_LOC = 1
    MY_JE_LOC = MY_JE_GLB - MY_JS_GLB + 1
!-------------------------------------------------------------------------------------------------- 
! EACH PE WILL NOW FILL ITS OWN SECTIONS OF THE GLOBAL-TO-LOCAL TRANSLATION ARRAYS 
! (DIMENSIONED GLOBALLY) AND LOCAL-TO-GLOBAL TRANSLATION ARRAYS (DIMENSIONED LOCALLY)
!-------------------------------------------------------------------------------------------------- 
    ILOC = 0
!
    DO I=MY_IS_GLB-1,MY_IE_GLB+1
        G2LI(I)    = ILOC
        L2GI(ILOC) = I
        ILOC       = ILOC + 1
    END DO
!
    JLOC = 0
!
    DO J=MY_JS_GLB-1,MY_JE_GLB+1
        G2LJ(J)    = JLOC
        L2GJ(JLOC) = J
        JLOC       = JLOC + 1
    END DO
!-------------------------------------------------------------------------------------------------- 
! EACH PE WILL NOW FILL THE ARRAY CALLED MY_NEB WHICH HOLDS THE NUMBER OF THE 8 PEs THAT ARE ITS 
! NEIGHBORS: 
! NORTH(1), EAST(2), SOUTH(3), WEST(4), NORTHEAST(5), SOUTHEAST(6), SOUTHWEST(7), AND NORTHWEST(8).
! THE VALUE IN THE ARRAY WILL BE -1 FOR THOSE NEIGHBORS THAT DO NOT EXIST BECAUSE THEY ARE BEYOND 
! THE GLOBAL DOMAIN BOUNDARY.
!-------------------------------------------------------------------------------------------------- 
    IPE = 0
!
    DO J=1,JNPES
        DO I=1,INPES
            ITEMP(I,J) = IPE
            IF (IPE == MYPE) THEN
                MYI = I
                MYJ = J
            END IF
            IPE = IPE + 1
        END DO
    END DO
!
    MY_N = -1
    IF (MYJ+1 <= JNPES)                                MY_N=ITEMP(MYI  , MYJ+1)
!
    MY_E = -1
    IF (MYI+1 <= INPES)                                MY_E=ITEMP(MYI+1, MYJ  )
!
    MY_S = -1
    IF (MYJ-1 >= 1)                                    MY_S=ITEMP(MYI  , MYJ-1)
!
    MY_W = -1
    IF (MYI-1 >= 1)                                    MY_W=ITEMP(MYI-1, MYJ  )
!
    MY_NE = -1
    IF ((MYI+1 <= INPES) .AND. (MYJ+1 <= JNPES))      MY_NE=ITEMP(MYI+1, MYJ+1)
!
    MY_SE = -1
    IF ((MYI+1 <= INPES) .AND. (MYJ-1 >= 1))          MY_SE=ITEMP(MYI+1, MYJ-1)
!
    MY_SW = -1
    IF ((MYI-1 >= 1) .AND. (MYJ-1 >= 1))              MY_SW=ITEMP(MYI-1, MYJ-1)
!
    MY_NW = -1
    IF ((MYI-1 >= 1) .AND. (MYJ+1 <= JNPES))          MY_NW=ITEMP(MYI-1, MYJ+1)
!
    MY_NEB(1) = MY_N
    MY_NEB(2) = MY_E
    MY_NEB(3) = MY_S
    MY_NEB(4) = MY_W
    MY_NEB(5) = MY_NE
    MY_NEB(6) = MY_SE
    MY_NEB(7) = MY_SW
    MY_NEB(8) = MY_NW
!-------------------------------------------------------------------------------------------------- 
! GENERATE THE TABLES (DIMENSIONED INPES*JNPES) THAT HOLD THE STARTING AND ENDING VALUES OF I AND J
! FOR EACH PE IN TERMS OF BOTH THE GLOBAL AND THE LOCAL DOMAINS.
!-------------------------------------------------------------------------------------------------- 
    CALL INDTABLE
!------------------------------------------ 
! CREATE ABBREVIATED NAMES FOR LOOP LIMITS.
!------------------------------------------
    MYIS     = MY_IS_LOC
    MYIS_P1  = MY_IS_LOC               - ILPAD1
    MYIS_P2  = MY_IS_LOC               - ILPAD2
    MYIS_P3  = MY_IS_LOC               - ILPAD3
    MYIS_P4  = MY_IS_LOC               - ILPAD4
    MYIS_P5  = MY_IS_LOC               - ILPAD5
!
    MYIS1    = MY_IS_LOC + IS_INC1_BND
    MYIS1_P1 = MY_IS_LOC + IS_INC1_BND - ILPAD1
    MYIS1_P2 = MY_IS_LOC + IS_INC1_BND - ILPAD2
    MYIS1_P3 = MY_IS_LOC + IS_INC1_BND - ILPAD3
    MYIS1_P4 = MY_IS_LOC + IS_INC1_BND - ILPAD4
!
    MYIS2    = MY_IS_LOC + IS_INC2_BND
!
    MYIE    = MY_IE_LOC
    MYIE_P1 = MY_IE_LOC + IRPAD1
    MYIE_P2 = MY_IE_LOC + IRPAD2
    MYIE_P3 = MY_IE_LOC + IRPAD3
    MYIE_P4 = MY_IE_LOC + IRPAD4
    MYIE_P5 = MY_IE_LOC + IRPAD5
!-------------------------------------------------------------------------------------------------- 
! THE SIZE OF THESE INCREMENTS IS ZERO UNLESS THE SUBDOMAIN LIES ALONG A GLOBAL BOUNDARY IN WHICH 
! CASE THE INCREMENT IS INDICATED BY THE NUMBER FOLLOWING 'INC'
!-------------------------------------------------------------------------------------------------- 
    MYIE1    = MY_IE_LOC - IE_INC1_BND         
    MYIE1_P1 = MY_IE_LOC - IE_INC1_BND + IRPAD1
    MYIE1_P2 = MY_IE_LOC - IE_INC1_BND + IRPAD2
    MYIE1_P3 = MY_IE_LOC - IE_INC1_BND + IRPAD3 
    MYIE1_P4 = MY_IE_LOC - IE_INC1_BND + IRPAD4  
!                                            
    MYIE2    = MY_IE_LOC - IE_INC2_BND
    MYIE2_P1 = MY_IE_LOC - IE_INC2_BND + IRPAD1
!
    MYJS     = MY_JS_LOC
    MYJS_P1  = MY_JS_LOC - JBPAD1
    MYJS_P2  = MY_JS_LOC - JBPAD2
    MYJS_P3  = MY_JS_LOC - JBPAD3
    MYJS_P4  = MY_JS_LOC - JBPAD4
    MYJS_P5  = MY_JS_LOC - JBPAD5
!
    MYJS1    = MY_JS_LOC + JS_INC1_BND
    MYJS1_P1 = MY_JS_LOC + JS_INC1_BND - JBPAD1
    MYJS1_P2 = MY_JS_LOC + JS_INC1_BND - JBPAD2
    MYJS1_P3 = MY_JS_LOC + JS_INC1_BND - JBPAD3
    MYJS1_P4 = MY_JS_LOC + JS_INC1_BND - JBPAD4
    MYJS1_P5 = MY_JS_LOC + JS_INC1_BND - JBPAD5
!
    MYJS2    = MY_JS_LOC + JS_INC2_BND
    MYJS2_P1 = MY_JS_LOC + JS_INC2_BND - JBPAD1
    MYJS2_P2 = MY_JS_LOC + JS_INC2_BND - JBPAD2
    MYJS2_P3 = MY_JS_LOC + JS_INC2_BND - JBPAD3
    MYJS2_P4 = MY_JS_LOC + JS_INC2_BND - JBPAD4
    MYJS2_P5 = MY_JS_LOC + JS_INC2_BND - JBPAD5
!
    MYJS3    = MY_JS_LOC + JS_INC3_BND
    MYJS3_P1 = MY_JS_LOC + JS_INC3_BND - JBPAD1
    MYJS3_P4 = MY_JS_LOC + JS_INC3_BND - JBPAD4
!
    MYJS4    = MY_JS_LOC + JS_INC4_BND
    MYJS4_P1 = MY_JS_LOC + JS_INC4_BND - JBPAD1
    MYJS4_P4 = MY_JS_LOC + JS_INC4_BND - JBPAD4
!
    MYJS5    = MY_JS_LOC + JS_INC5_BND
    MYJS5_P1 = MY_JS_LOC + JS_INC5_BND - JBPAD1
    MYJS5_P2 = MY_JS_LOC + JS_INC5_BND - JBPAD2
!
    MYJE     = MY_JE_LOC
    MYJE_P1  = MY_JE_LOC + JTPAD1
    MYJE_P2  = MY_JE_LOC + JTPAD2
    MYJE_P3  = MY_JE_LOC + JTPAD3
    MYJE_P4  = MY_JE_LOC + JTPAD4
    MYJE_P5  = MY_JE_LOC + JTPAD5
!
    MYJE1    = MY_JE_LOC - JE_INC1_BND
    MYJE1_P1 = MY_JE_LOC - JE_INC1_BND + JTPAD1
    MYJE1_P2 = MY_JE_LOC - JE_INC1_BND + JTPAD2
    MYJE1_P3 = MY_JE_LOC - JE_INC1_BND + JTPAD3
    MYJE1_P4 = MY_JE_LOC - JE_INC1_BND + JTPAD4
    MYJE1_P5 = MY_JE_LOC - JE_INC1_BND + JTPAD5
!
    MYJE2    = MY_JE_LOC - JE_INC2_BND
    MYJE2_P1 = MY_JE_LOC - JE_INC2_BND + JTPAD1
    MYJE2_P2 = MY_JE_LOC - JE_INC2_BND + JTPAD2
    MYJE2_P3 = MY_JE_LOC - JE_INC2_BND + JTPAD3
    MYJE2_P4 = MY_JE_LOC - JE_INC2_BND + JTPAD4
    MYJE2_P5 = MY_JE_LOC - JE_INC2_BND + JTPAD5
!
    MYJE3    = MY_JE_LOC - JE_INC3_BND
    MYJE3_P1 = MY_JE_LOC - JE_INC3_BND + JTPAD1
    MYJE3_P4 = MY_JE_LOC - JE_INC3_BND + JTPAD4
    MYJE3_P5 = MY_JE_LOC - JE_INC3_BND + JTPAD5
!
    MYJE4    = MY_JE_LOC - JE_INC4_BND
    MYJE4_P1 = MY_JE_LOC - JE_INC4_BND + JTPAD1
    MYJE4_P4 = MY_JE_LOC - JE_INC4_BND + JTPAD4
    MYJE4_P5 = MY_JE_LOC - JE_INC4_BND + JTPAD5
!
    MYJE5    = MY_JE_LOC - JE_INC5_BND
    MYJE5_P1 = MY_JE_LOC - JE_INC5_BND + JTPAD1
    MYJE5_P2 = MY_JE_LOC - JE_INC5_BND + JTPAD2
!
    RETURN
!
    END SUBROUTINE MPPINIT
!
!
!
    SUBROUTINE INDTABLE
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE INDTABLE
!>
!> SUBROUTINE: INDTABLE - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> THIS ROUTINE GENERATES THE TABLES THAT WILL GIVE THE STARTING AND ENDING VALUES OF I AND J FOR 
!> EACH PE ON THE GLOBAL AND LOCAL DOMAINS. 
!> EACH PE WILL HAVE A COPY OF THE FULL TABLES. 
!> THE ARGUMENT USED IS SIMPLY THE NUMBER OF THE PE FOR WHICH THESE VALUES ARE DESIRED.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????    - ORIGINATOR
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT  ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
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
!> DRIVER     : MPPINIT
!>
!> CALLS      : MPI_BARRIER
!>              MPI_BCAST
!>-------------------------------------------------------------------------------------------------- 
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
    INCLUDE "mpif.h"
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IPE     ,                                                                                   &
    & IRTN    ,                                                                                   &
    & IGLB    , JGLB    ,                                                                         &
    & IRECV
  
!
    IS_LOC_TABLE(MYPE) = MY_IS_LOC
    JS_LOC_TABLE(MYPE) = MY_JS_LOC
    IE_LOC_TABLE(MYPE) = MY_IE_LOC
    JE_LOC_TABLE(MYPE) = MY_JE_LOC
!
    IS_GLB_TABLE(MYPE) = MY_IS_GLB
    IE_GLB_TABLE(MYPE) = MY_IE_GLB
    JS_GLB_TABLE(MYPE) = MY_JS_GLB
    JE_GLB_TABLE(MYPE) = MY_JE_GLB
!
    DO IPE=0,NPES-1
        CALL MPI_BCAST(IS_LOC_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(JS_LOC_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(IE_LOC_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(JE_LOC_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
!    
        CALL MPI_BCAST(IS_GLB_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(JS_GLB_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(IE_GLB_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(JE_GLB_TABLE(IPE), 1, MPI_INTEGER, IPE, MPI_COMM_COMP, IRTN)
    END DO
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!-------------------------------------------------------------------------------------------------- 
! ALL OF THE PEs CAN NOW GENERATE A COMPLETE TABLE OF THE NUMBER OF GRID POINTS IN THE I DIRECTION 
! THAT ARE ON ALL OTHER PEs.
! THIS WILL BE USED IN THE MESINGER MSLP REDUCTION AS WELL AS IN THE BROADCAST BELOW.
!-------------------------------------------------------------------------------------------------- 
    DO IPE=0,NPES-1
        ICHUNKTAB(IPE) = IE_LOC_TABLE(IPE) - IS_LOC_TABLE(IPE) + 1
    END DO
!-------------------------------------------------------------------------------------------------- 
! SET UP A MAP OF THE GLOBAL DOMAIN THAT GIVES THE PE THAT OWNS EACH POINT. 
! (THIS APPEARS TO BE VESTIGIAL)
!-------------------------------------------------------------------------------------------------- 
!
!------------------------------------------------ 
! FIRST EACH PE FILLS IN ITS SECTION OF THE ARRAY
!------------------------------------------------
    DO JGLB=JS_GLB_TABLE(MYPE),JE_GLB_TABLE(MYPE)
        DO IGLB=IS_GLB_TABLE(MYPE),IE_GLB_TABLE(MYPE)
            ITEMP(IGLB,JGLB) = MYPE
        END DO
    END DO
!----------------------------------------------------------------- 
! NEXT, ALL PEs EXCHANGE THEIR SECTIONS SO EVERYONE HAS A FULL MAP
!----------------------------------------------------------------- 
    DO IPE=0,NPES-1
        DO JGLB=JS_GLB_TABLE(IPE),JE_GLB_TABLE(IPE)
            CALL MPI_BCAST(ITEMP(IS_GLB_TABLE(IPE), JGLB),                                        &
    &                               ICHUNKTAB(IPE), MPI_INTEGER, IPE, MPI_COMM_COMP, IRECV)
        END DO
    END DO
!
    CALL MPI_BARRIER(MPI_COMM_COMP, IRTN)
!
    RETURN
!
    END SUBROUTINE INDTABLE

