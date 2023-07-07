    SUBROUTINE BOCOH 
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE BOCOH
!>   
!> SUBROUTINE: BOCOH - UPDATE MASS POINTS ON BOUNDARY
!> PROGRAMMER: JANJIC          
!> ORG: W/NP22 
!> DATE: 94-03-08
!>     
!> ABSTRACT:
!> TEMPERATURE, SPECIFIC HUMIDITY, AND SURFACE PRESSURE ARE UPDATED ON THE DOMAIN BOUNDARY BY 
!> APPLYING THE PRE-COMPUTED TENDENCIES AT EACH TIME STEP.
!>     
!> PROGRAM HISTORY LOG:
!> 87-??-??  MESINGER   - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 96-12-13  BLACK      - FINAL MODIFICATION FOR NESTED RUNS
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 00-01-06  BLACK      - MODIFIED FOR JANJIC NONHYDROSTATIC CODE
!> 00-09-14  BLACK      - MODIFIED FOR DIRECT ACCESS READ
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
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
!> USE MODULES: BOCO
!>              CLDWTR
!>              CTLBLK
!>              DYNAM
!>              GLB_TABLE
!>              MAPOT
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              PARMETA
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : MPI_BCAST
!>--------------------------------------------------------------------------------------------------                 
    USE BOCO
    USE CLDWTR
    USE CTLBLK
    USE DYNAM
    USE GLB_TABLE
    USE MAPOT
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE NHYDRO
    USE PARMETA
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INCLUDE "mpif.h"
!
#include "sp.h"
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: ISIZ1 = 2  * LB
    INTEGER(KIND=I4KIND), PARAMETER :: ISIZ2 = 2  * LB * LM
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & NREC    , JK      , JN      , JL      , IRTN    , IIM     , JJM     , N       , I       ,   &
    & II      , J       , JJ
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & BCHR    , SHTM    , RHTM
!--------------------------------------
! READ FRESH BOUNDARY DATA IF NECESSARY
!--------------------------------------
    IF (NTSD-1 == NBOCO) THEN
        IF (MYPE == 0 .AND. NEST) THEN
            NREC = NINT((NTSD - 1) * DT / 3600.) + 2
            READ (NBC,REC=NREC) BCHR                          ,                                   &
    &             ((PDB(JK,JN)   , JK=1,LB)          , JN=1,2),                                   &
    &             (((TB(JK,JL,JN), JK=1,LB), JL=1,LM), JN=1,2),                                   &
    &             (((QB(JK,JL,JN), JK=1,LB), JL=1,LM), JN=1,2),                                   &
    &             (((UB(JK,JL,JN), JK=1,LB), JL=1,LM), JN=1,2),                                   &
    &             (((VB(JK,JL,JN), JK=1,LB), JL=1,LM), JN=1,2),                                   &
    &            (((Q2B(JK,JL,JN), JK=1,LB), JL=1,LM), JN=1,2),                                   &
    &           (((CWMB(JK,JL,JN), JK=1,LB), JL=1,LM), JN=1,2)
        END IF
!
        IF (MYPE == 0 .AND. .NOT. NEST) THEN
            READ(NBC) PDB
            READ(NBC) TB
            READ(NBC) QB
            READ(NBC) UB
            READ(NBC) VB
            READ(NBC) Q2B
            READ(NBC) CWMB
        END IF
!
        CALL MPI_BCAST(PDB,  ISIZ1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(TB,   ISIZ2, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(QB,   ISIZ2, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(UB,   ISIZ2, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(VB,   ISIZ2, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(Q2B,  ISIZ2, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
        CALL MPI_BCAST(CWMB, ISIZ2, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!----------------------------------
! FIND NEXT BOUNDARY CONDITION READ
!----------------------------------
        IF (NTSD < NTSTM) THEN
            IF (MYPE == 0 .AND. NEST) BCHR = BCHR + 1    ! THIS ASSUMES 1-HRLY BCs
            IF (MYPE == 0 .AND. .NOT. NEST) READ(NBC) BCHR
!
            CALL MPI_BCAST(BCHR, 1, MPI_REAL, 0, MPI_COMM_COMP, IRTN)
!
            NBOCO = INT(BCHR * TSPH + 0.5)
            IF (MYPE == 0) THEN
                PRINT*,'NEXT read BCHR, NBOCO = ', BCHR, NBOCO
            END IF
        END IF
!
    END IF
!
    IIM = IM - MY_IS_GLB + 1
    JJM = JM - MY_JS_GLB + 1
!----------------------------
! UPDATE THE SURFACE PRESSURE
!----------------------------
    N = 1
    DO 101 I=1,IM
        PDB(N,1) = PDB(N,1) + PDB(N,2) * DT
!
        IF (MY_JS_GLB == 1 .AND. I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
            II = I - MY_IS_GLB + 1
            PD(II,1) = PDB(N,1)
        END IF
!
        N = N + 1
        CONTINUE
101 END DO
!
    DO 102 I=1,IM
        PDB(N,1) = PDB(N,1) + PDB(N,2) * DT
!
        IF (MY_JE_GLB == JM .AND. I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
            II = I - MY_IS_GLB + 1
            PD(II,JJM) = PDB(N,1)
            IF (NTSD.lt.5760.and.mype.eq.539) then
               write(2000,*)'II:',II,PDB(N,1),PDB(N,2),PD(II,JJM),DT,NTSD,mype
            endif
        END IF
!
        N = N + 1
        CONTINUE
102 END DO
!
    DO 103 J=3,JM-2,2
        PDB(N,1) = PDB(N,1) + PDB(N,2) * DT
!
        IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
            JJ = J - MY_JS_GLB + 1
            PD(1,JJ) = PDB(N,1)
        END IF
!
        N = N + 1
        CONTINUE
103 END DO
!
    DO 104 J=3,JM-2,2
        PDB(N,1) = PDB(N,1) + PDB(N,2) * DT
!
        IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
            JJ = J - MY_JS_GLB + 1
            PD(IIM,JJ) = PDB(N,1)
        END IF
!
        N = N + 1
        CONTINUE
104 END DO
!------------------------------
! UPDATE THE 3-D MASS VARIABLES
!------------------------------
    DO 115 JL=1,LM
!
        N = 1
        DO 111 I=1,IM
              TB(N,JL,1) =   TB(N,JL,1) +   TB(N,JL,2) * DT
              QB(N,JL,1) =   QB(N,JL,1) +   QB(N,JL,2) * DT
             Q2B(N,JL,1) =  Q2B(N,JL,1) +  Q2B(N,JL,2) * DT
            CWMB(N,JL,1) = CWMB(N,JL,1) + CWMB(N,JL,2) * DT    
!
            IF (MY_JS_GLB == 1 .AND. I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
                II = I - MY_IS_GLB + 1
                   T(II,1,JL  ) =   TB(N,JL,1)
                   Q(II,1,JL  ) =   QB(N,JL,1)
                  Q2(II,1,JL  ) =  Q2B(N,JL,1)
                 CWM(II,1,JL  ) = CWMB(N,JL,1)                                                 
                PINT(II,1,JL+1) =   PD(II,1)   * ETA(JL+1) + PT
            END IF
!
            N = N + 1
            CONTINUE
    111 END DO
!
        DO 112 I=1,IM
              TB(N,JL,1) =   TB(N,JL,1) +   TB(N,JL,2) * DT
              QB(N,JL,1) =   QB(N,JL,1) +   QB(N,JL,2) * DT
             Q2B(N,JL,1) =  Q2B(N,JL,1) +  Q2B(N,JL,2) * DT
            CWMB(N,JL,1) = CWMB(N,JL,1) + CWMB(N,JL,2) * DT
!
            IF (MY_JE_GLB == JM .AND. I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
                II = I - MY_IS_GLB + 1
                   T(II,JJM,JL  )  =   TB(N,JL,1)
                   Q(II,JJM,JL  )  =   QB(N,JL,1)
                  Q2(II,JJM,JL  )  =  Q2B(N,JL,1)
                 CWM(II,JJM,JL  )  = CWMB(N,JL,1)     
                PINT(II,JJM,JL+1)  =   PD(II,JJM) * ETA(JL+1) + PT
            END IF
!
            N = N + 1
            CONTINUE
    112 END DO
!
        DO 113 J=3,JM-2,2
              TB(N,JL,1) =   TB(N,JL,1) +   TB(N,JL,2) * DT
              QB(N,JL,1) =   QB(N,JL,1) +   QB(N,JL,2) * DT
             Q2B(N,JL,1) =  Q2B(N,JL,1) +  Q2B(N,JL,2) * DT
            CWMB(N,JL,1) = CWMB(N,JL,1) + CWMB(N,JL,2) * DT
!
            IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JJ = J - MY_JS_GLB + 1
                     T(1,JJ,JL  ) =   TB(N,JL,1)
                     Q(1,JJ,JL  ) =   QB(N,JL,1)
                    Q2(1,JJ,JL  ) =  Q2B(N,JL,1)
                   CWM(1,JJ,JL  ) = CWMB(N,JL,1)  
                  PINT(1,JJ,JL+1) =   PD(1,JJ)   * ETA(JL+1) + PT
            END IF
!
            N = N + 1
            CONTINUE
    113 END DO
!
        DO 114 J=3,JM-2,2
              TB(N,JL,1) =   TB(N,JL,1) +   TB(N,JL,2) * DT
              QB(N,JL,1) =   QB(N,JL,1) +   QB(N,JL,2) * DT
             Q2B(N,JL,1) =  Q2B(N,JL,1) +  Q2B(N,JL,2) * DT
            CWMB(N,JL,1) = CWMB(N,JL,1) + CWMB(N,JL,2) * DT   
!
            IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JJ = J - MY_JS_GLB + 1
                   T(IIM,JJ,JL  ) =   TB(N,JL,1)
                   Q(IIM,JJ,JL  ) =   QB(N,JL,1)
                  Q2(IIM,JJ,JL  ) =  Q2B(N,JL,1)
                 CWM(IIM,JJ,JL  ) = CWMB(N,JL,1)
                PINT(IIM,JJ,JL+1) =   PD(IIM,JJ) * ETA(JL+1) + PT
            END IF
!
            N = N + 1
            CONTINUE
    114 END DO
!
        CONTINUE
115 END DO
!------------------------------------------------------
! SPACE INTERPOLATION OF PD AND T AT THE INNER BOUNDARY
!------------------------------------------------------
    IF (IBROW == 1) THEN
        DO 121 I=MYIS,MYIE1
            SHTM = HTM(I,1,LM) + HTM(I+1,1,LM) + HTM(I,3,LM) + HTM(I+1,3,LM)
!
            PD(I,2) = (PD(I  ,1) * HTM(I  ,1,LM)                                                  &
    &               +  PD(I+1,1) * HTM(I+1,1,LM)                                                  &
    &               +  PD(I  ,3) * HTM(I  ,3,LM)                                                  &
    &               +  PD(I+1,3) * HTM(I+1,3,LM))                                                 &
    &               /  SHTM
            CONTINUE
            IF (NTSD.lt.5760.and.mype.eq.539) then
               write(2001,*)'I:',I,PD(II,2),NTSD,SHTM,mype
               write(2001,*)'I:',I,PD(I ,1),PD(I+1,1),PD(I,3),PD(I+1,3)                                              
            endif
    121 END DO
    END IF
!
    IF (ITROW == 1) THEN
        DO 122 I=MYIS,MYIE1
            SHTM = HTM(I,JJM-2,LM) + HTM(I+1,JJM-2,LM) + HTM(I,JJM,LM) + HTM(I+1,JJM,LM)
!
            PD(I,JJM-1) = (PD(I  ,JJM-2) * HTM(I  ,JJM-2,LM)                                      &
    &                   +  PD(I+1,JJM-2) * HTM(I+1,JJM-2,LM)                                      &
    &                   +  PD(I  ,JJM  ) * HTM(I  ,JJM,  LM)                                      &
    &                   +  PD(I+1,JJM  ) * HTM(I+1,JJM  ,LM))                                     &
    &                   /  SHTM

            CONTINUE
    122 END DO
    END IF
!
    IF (ILCOL == 1) THEN
        DO 123 J=4,JM-3,2
            IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JJ = J - MY_JS_GLB + 1
                SHTM = HTM(1,JJ-1,LM) + HTM(2,JJ-1,LM) + HTM(1,JJ+1,LM) + HTM(2,JJ+1,LM)
!
                PD(1,JJ) = (PD(1,JJ-1) * HTM(1,JJ-1,LM) + PD(2,JJ-1) * HTM(2,JJ-1,LM)             &
    &                    +  PD(1,JJ+1) * HTM(1,JJ+1,LM) + PD(2,JJ+1) * HTM(2,JJ+1,LM))            &
    &                    /  SHTM
            END IF
            CONTINUE
    123 END DO
    END IF
!
    IF (IRCOL == 1) THEN
        DO 124 J=4,JM-3,2
            IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JJ = J - MY_JS_GLB + 1
                SHTM = HTM(IIM-1,JJ-1,LM) + HTM(IIM  ,JJ-1,LM)                                    &
    &                + HTM(IIM-1,JJ+1,LM) + HTM(IIM  ,JJ+1,LM)
!
                PD(IIM-1,JJ) = (PD(IIM-1,JJ-1) * HTM(IIM-1,JJ-1,LM)                               &
    &                        +  PD(IIM  ,JJ-1) * HTM(IIM  ,JJ-1,LM)                               &
    &                        +  PD(IIM-1,JJ+1) * HTM(IIM-1,JJ+1,LM)                               &
    &                        +  PD(IIM  ,JJ+1) * HTM(IIM  ,JJ+1,LM))                              &
    &                        /  SHTM
            END IF
            CONTINUE
    124 END DO
    END IF
!
    DO 135 JL=1,LM
!
        IF (IBROW == 1) THEN
            DO 131 I=MYIS,MYIE1
                RHTM = 1. / (HTM(I,1,JL) + HTM(I+1,1,JL) + HTM(I,3,JL) + HTM(I+1,3,JL))
!
                T(I,2,JL ) =    (T(I,1,JL  ) * HTM(I,1,JL) +    T(I+1,1,JL  ) * HTM(I+1,1,JL)     &
    &                      +     T(I,3,JL  ) * HTM(I,3,JL) +    T(I+1,3,JL  ) * HTM(I+1,3,JL))    &
    &                      *     RHTM
!
                Q(I,2,JL ) =    (Q(I,1,JL  ) * HTM(I,1,JL) +    Q(I+1,1,JL  ) * HTM(I+1,1,JL)     &
    &                      +     Q(I,3,JL  ) * HTM(I,3,JL) +    Q(I+1,3,JL  ) * HTM(I+1,3,JL))    &
    &                      *     RHTM
!
               Q2(I,2,JL ) =   (Q2(I,1,JL  ) * HTM(I,1,JL) +   Q2(I+1,1,JL  ) * HTM(I+1,1,JL)     &
    &                      +    Q2(I,3,JL  ) * HTM(I,3,JL) +   Q2(I+1,3,JL  ) * HTM(I+1,3,JL))    &
    &                      *    RHTM
!
              CWM(I,2,JL ) =  (CWM(I,1,JL  ) * HTM(I,1,JL) +  CWM(I+1,1,JL  ) * HTM(I+1,1,JL)     & 
    &                      +   CWM(I,3,JL  ) * HTM(I,3,JL) +  CWM(I+1,3,JL  ) * HTM(I+1,3,JL))    &
    &                      *   RHTM
!
             PINT(I,2,L+1) = (PINT(I,1,JL+1) * HTM(I,1,JL) + PINT(I+1,1,JL+1) * HTM(I+1,1,JL)     &
    &                      +  PINT(I,3,JL+1) * HTM(I,3,JL) + PINT(I+1,3,JL+1) * HTM(I+1,3,JL))    &
    &                      *  RHTM
                CONTINUE
        131 END DO
        END IF
!
        IF (ITROW == 1) THEN
            DO 132 I=MYIS,MYIE1
                RHTM = 1. / (HTM(I,JJM-2,JL) + HTM(I+1,JJM-2,JL) + HTM(I,JJM,JL) + HTM(I+1,JJM,JL))
!
                   T(I,JJM-1,JL ) =    (T(I  ,JJM-2,JL  ) * HTM(I  ,JJM-2,JL)                     &
    &                             +     T(I+1,JJM-2,JL  ) * HTM(I+1,JJM-2,JL)                     &
    &                             +     T(I  ,JJM  ,JL  ) * HTM(I  ,JJM  ,JL)                     &
    &                             +     T(I+1,JJM  ,JL  ) * HTM(I+1,JJM  ,JL))                    &
    &                             *     RHTM
!
                   Q(I,JJM-1,JL ) =    (Q(I  ,JJM-2,JL  ) * HTM(I  ,JJM-2,JL)                     &
    &                             +     Q(I+1,JJM-2,JL  ) * HTM(I+1,JJM-2,JL)                     &
    &                             +     Q(I  ,JJM  ,JL  ) * HTM(I  ,JJM  ,JL)                     &
    &                             +     Q(I+1,JJM  ,JL  ) * HTM(I+1,JJM  ,JL))                    &
    &                             *     RHTM
!
                  Q2(I,JJM-1,JL ) =   (Q2(I  ,JJM-2,JL  ) * HTM(I  ,JJM-2,JL)                     &
    &                             +    Q2(I+1,JJM-2,JL  ) * HTM(I+1,JJM-2,JL)                     &
    &                             +    Q2(I  ,JJM  ,JL  ) * HTM(I  ,JJM  ,JL)                     &
    &                             +    Q2(I+1,JJM  ,JL  ) * HTM(I+1,JJM  ,JL))                    &
    &                             *    RHTM
!
                 CWM(I,JJM-1,JL ) =  (CWM(I  ,JJM-2,JL  ) * HTM(I  ,JJM-2,JL)                     &
    &                             +   CWM(I+1,JJM-2,JL  ) * HTM(I+1,JJM-2,JL)                     &
    &                             +   CWM(I  ,JJM  ,JL  ) * HTM(I  ,JJM  ,JL)                     &
    &                             +   CWM(I+1,JJM  ,JL  ) * HTM(I+1,JJM  ,JL))                    &
    &                             *   RHTM
!
               PINT(I,JJM-1,JL+1) = (PINT(I  ,JJM-2,JL+1) * HTM(I  ,JJM-2,JL)                     &
    &                             +  PINT(I+1,JJM-2,JL+1) * HTM(I+1,JJM-2,JL)                     &
    &                             +  PINT(I  ,JJM  ,JL+1) * HTM(I  ,JJM  ,JL)                     &
    &                             +  PINT(I+1,JJM  ,JL+1) * HTM(I+1,JJM  ,JL))                    &
    &                             *  RHTM
                CONTINUE
        132 END DO
        END IF
!
        IF (ILCOL == 1) THEN
            DO 133 J=4,JM-3,2
                IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                    JJ = J - MY_JS_GLB + 1
                    RHTM = 1. / (HTM(1,JJ-1,JL) + HTM(2,JJ-1,JL) + HTM(1,JJ+1,JL) + HTM(2,JJ+1,JL))
!
                       T(1,JJ,JL ) =    (T(1,JJ-1,JL  ) * HTM(1,JJ-1,JL)                          & 
    &                              +     T(2,JJ-1,JL  ) * HTM(2,JJ-1,JL)                          &
    &                              +     T(1,JJ+1,JL  ) * HTM(1,JJ+1,JL)                          &
    &                              +     T(2,JJ+1,JL  ) * HTM(2,JJ+1,JL))                         &
    &                              *     RHTM
!
                       Q(1,JJ,JL ) =    (Q(1,JJ-1,JL  ) * HTM(1,JJ-1,JL)                          &
    &                              +     Q(2,JJ-1,JL  ) * HTM(2,JJ-1,JL)                          &
    &                              +     Q(1,JJ+1,JL  ) * HTM(1,JJ+1,JL)                          & 
    &                              +     Q(2,JJ+1,JL  ) * HTM(2,JJ+1,JL))                         &
    &                              *     RHTM
!
                      Q2(1,JJ,JL ) =   (Q2(1,JJ-1,JL  ) * HTM(1,JJ-1,JL)                          &
    &                              +    Q2(2,JJ-1,JL  ) * HTM(2,JJ-1,JL)                          &
    &                              +    Q2(1,JJ+1,JL  ) * HTM(1,JJ+1,JL)                          &
    &                              +    Q2(2,JJ+1,JL  ) * HTM(2,JJ+1,JL))                         &
    &                              *    RHTM
!
                     CWM(1,JJ,JL ) =  (CWM(1,JJ-1,JL )  * HTM(1,JJ-1,JL)                          &
    &                              +   CWM(2,JJ-1,JL )  * HTM(2,JJ-1,JL)                          &
    &                              +   CWM(1,JJ+1,JL )  * HTM(1,JJ+1,JL)                          &
    &                              +   CWM(2,JJ+1,JL )  * HTM(2,JJ+1,JL))                         &
    &                              *   RHTM
!
                   PINT(1,JJ,JL+1) = (PINT(1,JJ-1,JL+1) * HTM(1,JJ-1,JL)                          &
    &                              +  PINT(2,JJ-1,JL+1) * HTM(2,JJ-1,JL)                          &
    &                              +  PINT(1,JJ+1,JL+1) * HTM(1,JJ+1,JL)                          & 
    &                              +  PINT(2,JJ+1,JL+1) * HTM(2,JJ+1,JL))                         &
    &                              *  RHTM
                END IF
                CONTINUE
        133 END DO
        END IF
!
        IF (IRCOL == 1) THEN
            DO 134 J=4,JM-3,2
                IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                    JJ = J - MY_JS_GLB + 1
                    RHTM = 1. / (HTM(IIM-1,JJ-1,JL) + HTM(IIM  ,JJ-1,JL) + HTM(IIM-1,JJ+1,JL)     &
    &                    +       HTM(IIM  ,JJ+1,JL))
!
                       T(IIM-1,JJ,JL ) =    (T(IIM-1,JJ-1,JL  ) * HTM(IIM-1,JJ-1,JL)              &
    &                                  +     T(IIM  ,JJ-1,JL  ) * HTM(IIM  ,JJ-1,JL)              &
    &                                  +     T(IIM-1,JJ+1,JL  ) * HTM(IIM-1,JJ+1,JL)              &
    &                                  +     T(IIM  ,JJ+1,JL  ) * HTM(IIM  ,JJ+1,JL))             &
    &                                  *     RHTM
!
                       Q(IIM-1,JJ,JL ) =    (Q(IIM-1,JJ-1,JL  ) * HTM(IIM-1,JJ-1,JL)              &
    &                                  +     Q(IIM  ,JJ-1,JL  ) * HTM(IIM  ,JJ-1,JL)              &
    &                                  +     Q(IIM-1,JJ+1,JL  ) * HTM(IIM-1,JJ+1,JL)              &
    &                                  +     Q(IIM  ,JJ+1,JL  ) * HTM(IIM  ,JJ+1,JL))             &
    &                                  *     RHTM
!
                      Q2(IIM-1,JJ,JL ) =   (Q2(IIM-1,JJ-1,JL  ) * HTM(IIM-1,JJ-1,JL)              &
    &                                  +    Q2(IIM  ,JJ-1,JL  ) * HTM(IIM  ,JJ-1,JL)              &
    &                                  +    Q2(IIM-1,JJ+1,JL  ) * HTM(IIM-1,JJ+1,JL)              &
    &                                  +    Q2(IIM  ,JJ+1,JL  ) * HTM(IIM  ,JJ+1,JL))             &
    &                                  *    RHTM
!
                     CWM(IIM-1,JJ,JL ) =  (CWM(IIM-1,JJ-1,JL  ) * HTM(IIM-1,JJ-1,JL)              &
    &                                  +   CWM(IIM  ,JJ-1,JL  ) * HTM(IIM  ,JJ-1,JL)              &
    &                                  +   CWM(IIM-1,JJ+1,JL  ) * HTM(IIM-1,JJ+1,JL)              &
    &                                  +   CWM(IIM  ,JJ+1,JL  ) * HTM(IIM  ,JJ+1,JL))             &
    &                                  *   RHTM
!
                   PINT(IIM-1,JJ,JL+1) = (PINT(IIM-1,JJ-1,JL+1) * HTM(IIM-1,JJ-1,JL)              &
    &                                  +  PINT(IIM  ,JJ-1,JL+1) * HTM(IIM  ,JJ-1,JL)              &
    &                                  +  PINT(IIM-1,JJ+1,JL+1) * HTM(IIM-1,JJ+1,JL)              &
    &                                  +  PINT(IIM  ,JJ+1,JL+1) * HTM(IIM  ,JJ+1,JL))             &
    &                                  *  RHTM           
                END IF
                CONTINUE
        134 END DO
        END IF
        CONTINUE
135 END DO
!
    RETURN
!
    END SUBROUTINE BOCOH
