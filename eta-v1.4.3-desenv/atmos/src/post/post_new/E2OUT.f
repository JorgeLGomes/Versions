      SUBROUTINE E2OUT(ITAG1,ITAG2,EGRID1,EGRID2,
     1     GRID1,GRID2,IMOUT,JMOUT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    E2OUT       INTRP E-GRID TO OUTPUT GRID
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-21       
C     
C ABSTRACT:  
C     THIS SUBROUTINE PERFORMS ALL ETA POST PROCESSOR 
C     INTERPOLATION/CONVERSION OF DATA ON THE E-GRID TO
C     THE OUTPUT GRID.  ADDITIONALLY, THE ROUTINE WILL
C     SMOOTH OR FILTER THE FIELD(S) AT ANY OF THREE 
C     STEPS IN THE ROUTINE.  DATA MAY BE SMOOTHED ON 
C     THE INPUT E-GRID, FILTERED ON A FILLED E-GRID, OR
C     FILTERED ON THE OUTPUT GRID.  VORTICITY FIELDS ARE
C     GIVEN AN ADDITIONAL HEAVY HANDED SMOOTHING TO PRODUCE
C     A PLEASING PRODUCT.  CONTROL OF SMOOTHING/FILTERING
C     IS VIA SWITCHES SET IN THE CONTROL FILE.  
C   .     
C     
C PROGRAM HISTORY LOG:
C   92-12-21  RUSS TREADON
C   93-06-13  RUSS TREADON - ADDED INTERPOLATION TO LAT-LON GRID.
C   98-06-01  BLACK - CONVERSION FROM 1-D TO 2-D
C   00-01-05  JIM TUCCILLO - MPI VERSION
C     
C USAGE:    CALL E2OUT(ITAG1,ITAG2,EGRID1,EGRID2,
C                      GRID1,GRID2,IMOUT,JMOUT)
C   INPUT ARGUMENT LIST:
C     ITAG1    - INTEGER ID FOR DATA IN EGRID1
C     ITAG2    - INTEGER ID FOR DATA IN EGRID2
C     EGRID1   - FIRST FIELD ON E-GRID
C     EGRID2   - SECOND FIELD ON E-GRID
C     IMOUT    - FIRST DIMENSION OF OUTPUT GRID
C     JMOUT    - SECOND DIMENSION OF OUTPUT GRID
C
C   OUTPUT ARGUMENT LIST: 
C     GRID1    - FIRST FIELD ON OUTPUT GRID
C     GRID2    - SECOND FIELD ON OUTPUT GRID
C     
C   OUTPUT FILES:
C     STDOUT  - RUN TIME STANDARD OUT.
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       LOAD     - LOAD FILLED E-GRID INTO 2-D OUTPUT GRID
C       EFILL    - REPLACES MISSING VALUES ON E-GRID WITH FIELD MEAN
C       FILLH    - FILL MASS POINTS FOR VELOCITY POINT E-GRID
C       FILLV    - FILL VELOCITY POINTS FOR MASS POINT E-GRID
C       CETLIH4  - INTERPOLATE E-GRID TO OUPUT GRID CONSERVING
C                  THE AREA INTEGRAL OF THE INPUT FIELD OVER
C                  SPECIFIED SUB-GRIDS.
C       EUVGUV   - ROTATE ETA (U,V) TO OUTPUT GRID (U,V)
C       INTERP3  - BILINEAR INTERPOLATION TO OUTPUT GRID
C       P2FILTF  - SMOOTH MASS POINT DATA ON E-GRID
C       P2FLTVF  - SMOOTH VELOCITY POINT DATA ON E-GRID
C       EFILT    - HEAVY-HANDED MASS POINT SMOOTHER.  CURRENTLY
C                  HARDWIRED FOR USE ON VORITCITY FIELDS.
C       FILTER   - 25 POINT BLECK FILTER ON A REGULAR GRID.
C     LIBRARY:
C       COMMON   - OPTIONS
C                  LLGRDS
C                  RQSTFLD
C                  MASKS
C                  IOUNIT
C                  OUTGRD
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C     
C     INCLUDE ETA MODEL DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
C     
      INCLUDE "parmeta"
      INCLUDE "parmout"
      PARAMETER (IMJM=IM*JM-JM/2,IMT=2*IM-1,JMT=JM,LP1=LM+1)
      PARAMETER (D00=0.0)
C     
C     DECLARE VARIABLES.
      CHARACTER*6  PROJ
      LOGICAL NORTH
      REAL EGRID1(IM,JM),EGRID2(IM,JM)
      REAL GRID1(IMOUT,JMOUT),GRID2(IMOUT,JMOUT)
      REAL HFUL(IMT,JMT),UFUL(IMT,JMT),VFUL(IMT,JMT)
C     
C     INCLUDE COMMONS.
      INCLUDE "OPTIONS.comm"
      INCLUDE "LLGRDS.comm"
      INCLUDE "RQSTFLD.comm"
      INCLUDE "MASKS.comm"
      INCLUDE "BITMAP.comm"
      INCLUDE "IOUNIT.comm"
      INCLUDE "OUTGRD.comm"
      INCLUDE "CTLBLK.comm"
C
      common/jjt/time_output, time_e2out
      real(8) ist, rtc, time_output, time_e2out
C*****************************************************************************
C     START SUBROUTINE E2OUT.
C     
C
C     GATHER EGRID1 AND EGRID2
C
C      ist = rtc()
      CALL COLLECT(EGRID1)
c      print*,'writing egrid1 in E2OUT after calling COLLECT'
c      do i=1,im
c      do j=jsta,jend
c      print*,i,j,egrid1(i,j)
c      end do
c      end do
      CALL COLLECT(EGRID2)
C
      IF ( ME .EQ. 0 ) THEN
C
C     ZERO OUTPUT GRIDS.
C
      GRID1(1:IMOUT,1:JMOUT)=D00
      GRID2(1:IMOUT,1:JMOUT)=D00
C     
C     GRID TYPE KGTYPE=90,92,94,96 IS OUTPUT ON THE STAGGERED MODEL
C     GRIDS. SIMPLY LOAD DATA INTO OUTPUT GRID ARRAYS.  IF WANTED,
C     SMOOTH DATA ON OUTPUT GRID.  (YOU NEVER KNOW!)
C    
   
      IF (KGTYPE.EQ.90.OR.KGTYPE.EQ.92.OR.
     X      KGTYPE.EQ.94.OR.KGTYPE.EQ.96.OR.
     X      KGTYPE.EQ.99.OR.KGTYPE.EQ.190.OR.KGTYPE.EQ.192
Cmp
     X      .OR.KGTYPE.EQ.194.OR.KGTYPE.EQ.196.or.KGTYPE.eq.255) THEN
cmp
c          call ps7
         CALL EFILL(EGRID1,IM,JM)
c         print*,'writing egrid1 in E2OUT after calling EFILL'
c         do i=1,im
c         do j=jsta,jend
c         print*,i,j,egrid1(i,j)
c         end do
c         end do
c           call ps8
         DO J=1,JM
         DO I=1,IM
           GRID1(I,J)=EGRID1(I,J)
           IBMAP(I,J)=1
         ENDDO
         ENDDO
cc          call ps9
         IF (ITAG2.GT.0) THEN
c              call psa
            CALL EFILL(EGRID2,IM,JM)
c
            DO J=1,JM
            DO I=1,IM
              GRID2(I,J)=EGRID2(I,J)
            ENDDO
            ENDDO
c              call psb
         ENDIF
         GOTO 400
      ENDIF
      
C     
C     IF SELECTED, SMOOTH DATA ON STAGGERED E-GRID.
C     
      ISMTH = ISMSTG(IGET(ITAG1))
      IF (IFILV(ITAG1).GT.0) THEN
         IF (ISMTH.GT.0) THEN
              CALL P2FILTF(ISMTH,HBM2,EGRID1)
         END IF
         IF ((ISMTH.GT.0).AND.((ITAG1.EQ.10).OR.(ITAG1.EQ.21))) THEN
              CALL EFILT(EGRID1)
         END IF
         IF(ITAG2.GT.0)THEN
           ISMTH = ISMSTG(IGET(ITAG2))
           IF(ISMTH.GT.0) THEN 
              CALL P2FILTF(ISMTH,HBM2,EGRID2)
           END IF
         ENDIF
      ELSE
         IF (ISMTH.GT.0) THEN
              CALL P2FLTVF(ISMTH,VBM2,EGRID1)
         END IF
         IF(ITAG2.GT.0)THEN
           ISMTH = ISMSTG(IGET(ITAG2))
           IF(ISMTH.GT.0) THEN
              CALL P2FLTVF(ISMTH,VBM2,EGRID2)
           END IF
         ENDIF
      ENDIF
C     
C     INTERPOLATE H-POINT FIELDS FROM STAGGERED E-GRID TO OUTPUT GRID.
C     
      IF ( IFILV(ITAG1).GT.0) THEN
C     
C        FILL H-POINT FIELD AT V POINTS.
c         call psd
         CALL FILLV(EGRID1,HFUL,IFLAG,IMT,JMT)
c          call pse
C     
C        IF REQUESTED, SMOOTH DATA ON FILLED E-GRID.
         ISMTH = ISMFUL(IGET(ITAG1))
         IF (ISMTH.GT.0) CALL FILTER(IMT,JMT,HFUL,ISMTH)
C     
C        INTERPOLATE TO OUTPUT GRID.
         IF ((KGTYPE.LT.90).OR.(KGTYPE.GT.97
     1        .AND.KGTYPE.NE.99.AND.KGTYPE.NE.190.AND.KGTYPE.NE.192
Cmp
     1 .AND.KGTYPE.NE.194.AND.KGTYPE.NE.196.and.kgtype.ne.255))THEN
C     
C           INTERPOLATION FOR GENERIC FIELDS.
c       call ps1
          if (itag1.eq.050) then
            CALL OUT_MASKS(HFUL,GRID1,IMOUT,JMOUT)
C             CALL INTERP3(HFUL,GRID1,IMOUT,JMOUT)
C               where(GRID1> 0.5)
C                 GRID1=1.      
C               else where
C                 GRID1=0 
C               end where
          else
            CALL INTERP3(HFUL,GRID1,IMOUT,JMOUT)
          endif
C
C           U-V FIELDS ORGINALLY AT H POINTS NEED TO BE ROTATED.
            IF ( (ITAG1.EQ.056).OR.(ITAG1.EQ.057).OR.
     X           (ITAG1.EQ.060).OR.(ITAG1.EQ.061).OR.
     X           (ITAG1.EQ.064).OR.(ITAG1.EQ.065).OR.
     X           (ITAG1.EQ.073).OR.(ITAG1.EQ.074).OR.
     X           (ITAG1.EQ.095).OR.(ITAG1.EQ.096) ) THEN
               CALL FILLV(EGRID2,VFUL,IFLAG,IMT,JMT)
               IF(ITAG2.GT.0)THEN
                 ISMTH = ISMFUL(IGET(ITAG2))
                 IF(ISMTH.GT.0)CALL FILTER(IMT,JMT,VFUL,ISMTH)
               ENDIF
               CALL EUVGUV(HFUL,VFUL,FVTLON,IMT,JMT,EVLAT,
     X              EVLON,ALATVT,ALONVT,NORTH,PROJ)
c               call ps3
               CALL INTERP3(HFUL,GRID1,IMOUT,JMOUT)
               CALL INTERP3(VFUL,GRID2,IMOUT,JMOUT)
            ENDIF

C     
C           PRECIPITATION FIELDS USE AREA CONSERVING INTERPOLATION.
c            call ps4
            IF ( (ITAG1.EQ.033).OR.(ITAG1.EQ.034).OR.
     X           (ITAG1.EQ.087) )
     X           CALL CETLIH4(EGRID1,GRID1,IMOUT,JMOUT,KSB,IOFFS)
C     
C        OUTPUT ON FILLED E-GRID REQUIRES NO INTERPOLATION.
         ELSE
            DO J=1,JMT
            DO I=1,IMT
              GRID1(I,J)=HFUL(I,J)
              IBMAP(I,J)=1
            ENDDO
            ENDDO
C
            IF ( (ITAG1.EQ.056).OR.(ITAG1.EQ.057).OR.
     X           (ITAG1.EQ.060).OR.(ITAG1.EQ.061).OR.
     X           (ITAG1.EQ.064).OR.(ITAG1.EQ.065).OR.
     X           (ITAG1.EQ.073).OR.(ITAG1.EQ.074).OR.
     X           (ITAG1.EQ.095).OR.(ITAG1.EQ.096) ) THEN
               CALL FILLV(EGRID2,VFUL,IFLAG,IMT,JMT)
               IF(ITAG2.GT.0)THEN
                 ISMTH = ISMFUL(IGET(ITAG2))
                 IF(ISMTH.GT.0)CALL FILTER(IMT,JMT,VFUL,ISMTH)
               ENDIF
C
               DO J=1,JMT
               DO I=1,IMT
                 GRID2(I,J)=VFUL(I,J)
               ENDDO
               ENDDO
C
            ENDIF
         ENDIF
C     
C     NOW HANDLE FIELDS AT V-POINTS.
C
      ELSE
         IF (ITAG1.NE.053.AND.ITAG1.NE.162) THEN 
            CALL FILLH(EGRID1,UFUL,IMT,JMT)
            CALL FILLH(EGRID2,VFUL,IMT,JMT)
            ISMTH = ISMFUL(IGET(ITAG1))
            IF (ISMTH.GT.0) CALL FILTER(IMT,JMT,UFUL,ISMTH)
            IF(ITAG2.GT.0)THEN
              ISMTH = ISMFUL(IGET(ITAG2))
              IF(ISMTH.GT.0)
     X           CALL FILTER(IMT,JMT,VFUL,ISMTH)
            ENDIF
C
            IF ((KGTYPE.LT.90).OR.(KGTYPE.GT.97
     1          .AND.KGTYPE.NE.99.AND.KGTYPE.NE.190.AND.KGTYPE.NE.192
Cmp
     1    .AND.KGTYPE.NE.194.AND.KGTYPE.NE.196.and.KGTYPE.NE.255))THEN
               CALL EUVGUV(UFUL,VFUL,FVTLON,IMT,JMT,
     X              EVLAT,EVLON,ALATVT,ALONVT,NORTH,PROJ)
               CALL INTERP3(UFUL,GRID1,IMOUT,JMOUT)
               CALL INTERP3(VFUL,GRID2,IMOUT,JMOUT)
            ELSE
C
               DO J=1,JMT
               DO I=1,IMT
                 GRID1(I,J)=UFUL(I,J)
                 IBMAP(I,J)=1
               ENDDO
               ENDDO
C
               DO J=1,JMT
               DO I=1,IMT
                 GRID2(I,J)=VFUL(I,J)
               ENDDO
               ENDDO
C
            ENDIF
         ELSE
            CALL FILLH(EGRID1,HFUL,IMT,JMT)
            ISMTH = ISMFUL(ITAG1)
            IF (ISMTH.GT.0) CALL FILTER(IMT,JMT,HFUL,ISMTH)
            IF ((KGTYPE.LT.90).OR.(KGTYPE.GT.97
     1          .AND.KGTYPE.NE.99.AND.KGTYPE.NE.190.AND.KGTYPE.NE.192
Cmp
     1   .AND.KGTYPE.NE.194.AND.KGTYPE.NE.196.and.KGTYPE.ne.255))THEN
c        call ps6
               CALL INTERP3(HFUL,GRID1,IMOUT,JMOUT)
            ELSE
               DO J=1,JMT
               DO I=1,IMT
                 GRID1(I,J)=HFUL(I,J)
                 IBMAP(I,J)=1
               ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF
C     
C     ZERO FILLED EGRID ARRAYS.
C
      HFUL(1:IMT,1:JMT) = 0.0
      UFUL(1:IMT,1:JMT) = 0.0
      VFUL(1:IMT,1:JMT) = 0.0
C     
C     IF SELECTED, APPLY SMOOTHER TO OUTPUT GRID(S).
C     
 400  CONTINUE
C
      IF(ITAG1.GT.0)THEN
        IF(IGET(ITAG1).GT.0)THEN
          ISMTH=ISMOUT(IGET(ITAG1))
          IF(ISMTH.GT.0)CALL FILTER(IMOUT,JMOUT,GRID1,ISMTH)
        ENDIF
      ENDIF
C
      IF(ITAG2.GT.0)THEN
        IF(IGET(ITAG2).GT.0)THEN
          ISMTH=ISMOUT(IGET(ITAG2))
          IF(ISMTH.GT.0)CALL FILTER(IMOUT,JMOUT,GRID2,ISMTH)
        ENDIF
      ENDIF
C     
      END IF
C
C     SCATTER EGRID1 AND EGRID2
C
c      if(itag1.eq.23)then
c      print*,'printing grid1 before calling DIST'
c      do i=1,im
c      do j=jsta,jend
c       print*,i,j,grid1(i,j)
c      end do
c      end do 
c      end if
c      CALL DIST(GRID1)
c      CALL DIST(GRID2)
c      if(itag1.eq.23)then
c      print*,'printing grid1 after calling DIST'
c      do i=1,im
c      do j=jsta,jend
c       print*,i,j,grid1(i,j)
c      end do
c      end do
c      end if
      CALL DIST(EGRID1)
      CALL DIST(EGRID2)
C
C      time_e2out = time_e2out + rtc() - ist
C     END OF ROUTINE.
C     
      RETURN
      END
