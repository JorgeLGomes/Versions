    SUBROUTINE DRIVE_SST (SST,SST2,SM,IM,JM,DLMD,DPHD,TLM0D,TPH0D,SSTG,NX,NXP1,NY,INCR_LAT,INCR_LON) 
!
    USE F77KINDS
    IMPLICIT NONE
!
!
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: IM
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: JM  
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: NX
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: NXP1
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: NY
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: DLMD
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: DPHD
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: TLM0D
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: TPH0D
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: INCR_LAT
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: INCR_LON
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(IN)              :: SM
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(INOUT)           :: SST
    REAL   (KIND=R4KIND)   , DIMENSION(IM  ,JM)   , INTENT(INOUT)           :: SST2
    REAL   (KIND=R4KIND)   , DIMENSION(NXP1,NY)   , INTENT(IN)              :: SSTG

    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)                               :: GLAT
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)                               :: GLON
    INTEGER(KIND=R4KIND)   , DIMENSION(IM,JM)                               :: LON1INDX 
    INTEGER(KIND=R4KIND)   , DIMENSION(IM,JM)                               :: LON2INDX
    INTEGER(KIND=R4KIND)   , DIMENSION(IM,JM)                               :: LAT1INDX
    INTEGER(KIND=R4KIND)   , DIMENSION(IM,JM)                               :: LAT2INDX 
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)                               :: AR1
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)                               :: AR2
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)                               :: AR3
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)                               :: AR4 

!
    INTEGER(KIND=R4KIND)                                                    :: ILAT1
    INTEGER(KIND=R4KIND)                                                    :: ILON1
    INTEGER(KIND=R4KIND)                                                    :: ILAT2
    INTEGER(KIND=R4KIND)                                                    :: ILON2 
    REAL   (KIND=R4KIND)                                                    :: ELAT
    REAL   (KIND=R4KIND)                                                    :: ELON
    REAL   (KIND=R4KIND)                                                    :: DIF
    REAL   (KIND=R4KIND)                                                    :: W1 
    REAL   (KIND=R4KIND)                                                    :: W2
    INTEGER(KIND=R4KIND)                                                    :: I
    INTEGER(KIND=R4KIND)                                                    :: J
    INTEGER(KIND=R4KIND)                                                    :: INSST
    INTEGER(KIND=R4KIND)                                                    :: INDXST
 
!
    CALL Geographic(GLAT, GLON,IM,JM,DLMD,DPHD,TLM0D,TPH0D)
! 
    CALL SSTG2EGRD                                                                                &       
    &(GLAT,GLON,INCR_LAT,INCR_LON,NX,AR1,AR2,AR3,AR4,LON1INDX,LON2INDX,LAT1INDX,LAT2INDX)

    DO J=1,JM
      DO I=1,IM
        SST2(I,J)=AR1(I,J)*SSTG(LON2INDX(I,J),LAT2INDX(I,J))+                                    &
        &         AR2(I,J)*SSTG(LON2INDX(I,J),LAT1INDX(I,J))+                                    &
        &         AR3(I,J)*SSTG(LON1INDX(I,J),LAT1INDX(I,J))+                                    &
        &         AR4(I,J)*SSTG(LON1INDX(I,J),LAT2INDX(I,J))
        SST(I,J)=SSTG(LON2INDX(I,J),LAT2INDX(I,J))
        IF ((SST(I,J).LT.100).AND.(SM(I,J).GT.0.5)) THEN
          WRITE(6,*)'SST(',I,',',J,')=',SST(I,J),'SSTG(',LON2INDX(I,J),',',LAT2INDX(I,J),')=',SSTG(LON2INDX(I,J),LAT2INDX(I,J))
        END IF
        IF ((SST2(I,J).LT.100).AND.(SM(I,J).GT.0.5)) THEN
          WRITE(6,*)'SST2(',I,',',J,')=',SST2(I,J),'SM(',I,',',J,')=',SM(I,J)
        END IF
      ENDDO
    ENDDO
!
    END  SUBROUTINE DRIVE_SST





