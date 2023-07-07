     SUBROUTINE SSTG2EGRD(GLATR,GLONR,INCR_LAT,INCR_LON,NX,NY,LAT_P1,LON_P1,LAT_PNY,LON_PNX         &
                       & ,AR1,AR2,AR3,INH,JNH)
    USE PARMCONF
    USE CONSTANTS
    USE DIAGNOSTIC
!
    IMPLICIT NONE
!
!
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(IN)              :: GLATR
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(IN)              :: GLONR  
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: INCR_LAT
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: INCR_LON
    INTEGER(KIND=I4KIND)                          , INTENT(IN)              :: NX
    INTEGER(KIND=I4KIND)                          , INTENT(IN)              :: NY
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: LAT_P1
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: LON_P1
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: LAT_PNY
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: LON_PNX
    INTEGER(KIND=I4KIND)   , DIMENSION(4,IM  ,JM) , INTENT(INOUT)           :: INH
    INTEGER(KIND=I4KIND)   , DIMENSION(4,IM  ,JM) , INTENT(INOUT)           :: JNH
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(INOUT)           :: AR1
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(INOUT)           :: AR2
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(INOUT)           :: AR3
!
    REAL   (KIND=R4KIND)   , DIMENSION(3,IM  ,JM)                           :: COH
    INTEGER(KIND=R4KIND)                                                    :: I
    INTEGER(KIND=R4KIND)                                                    :: J 

!----  INTERPOLATE 0.25-DEG GLOBAL SST TO ETA GRID  -------
!
!-CP NOTE:  THIS SUBROUTINE AND INTERPOLATION ALGORITHM ASSUME
!-CP A 0.25-DEG GLOBAL SST FIELD IN THE FOLLOWING FORMAT:  
!
!	I=1 at 0.125 E, I=1440 at 359.875E, I=1441 at 0.125 E
!	J=1 at 89.875 S, J=720 at 89.875 N 
!
!	Old 1 degree data
!	I=1 AT 0.5 E,  I=2 AT 1.5 E, ... , I=360 at 0.5W
!	J=1 AT 89.5S, J=2 AT 88.5 S, ..., J=180 at 89.5N
!
!	Old 0.5 degree data
!	I=1 at 0.25 E, I=720 at 359.75E, I=721 at 0.25 E
!	J=1 at 89.75 S, J=360 at 89.75 N 
!
!-CP  
!-CP In the interpolation algorithm below, glon is positive westward,
!-CP from 0 to 360, with 0 at the greenwich meridian.  Elon is positive 
!-CP eastward, thus the need to subtract glon from 360 to get the index
!-CP of the correct oisst point.  If your input 1 deg SST field is in
!-CP a different indexing scheme, you will need to change the algorithm
!-CP below - see "grdeta.oldoi"

      DO J=1,JM
        DO I=1,IM
 	  CALL CED_IJ(GLATR(I,J),GLONR(I,J),INH(1,I,J),JNH(1,I,J),                                     &
     &         COH(1,I,J),COH(2,I,J),COH(3,I,J),                                     &
     &         NX,NY,LAT_P1,LON_P1,LAT_PNY,LON_PNX,INCR_LON,INCR_LAT)
          if (inh(1,I,J) .eq. 0 ) then
            inh(1,I,J)=NX
          end if
          inh(2,I,J)=inh(1,I,J)+1
          inh(3,I,J)=inh(1,I,J)
          inh(4,I,J)=inh(1,I,J)+1
!
          if(INH(1,I,J).eq.NX) then
            INH(2,I,J)=1
            INH(4,I,J)=1
          endif
          JNH(2,I,J)=JNH(1,I,J)
          JNH(3,I,J)=JNH(1,I,J)+1
          JNH(4,I,J)=JNH(1,I,J)+1

          AR1(I,J)=COH(2,I,J)
          AR2(I,J)=COH(3,I,J)
          AR3(I,J)=COH(1,I,J)       
         
          IF (DIAG) THEN
	    if ( (mod(I,(IM/2)) .eq. 0) .and. (mod(J,(JM/2)) .eq. 0) ) then
	      write(6,*) 'I, J: ',I, J
	      write(6,*) 'weights: ',AR1(I,J),AR2(I,J),AR3(I,J)
	      write(6,*) 'LAT,LON: ', GLATR(I,J),GLONR(I,J)
	      write(6,*) '------------------------------------------'
	     endif
          ENDIF
        END DO
      END DO
!
      END SUBROUTINE SSTG2EGRD





