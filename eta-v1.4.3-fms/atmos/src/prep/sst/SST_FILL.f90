   SUBROUTINE SST_FILL(GSST,NX,NY,Undef)
   USE DIAGNOSTIC
!
    IMPLICIT NONE
!
!
    REAL   (KIND=R4KIND)                                                    :: Undef
    INTEGER(KIND=R4KIND)                                                    :: I
    INTEGER(KIND=R4KIND)                                                    :: J
    INTEGER(KIND=R4KIND)                                                    :: NP
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: NX
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: NY
    REAL   (KIND=R4KIND)   , DIMENSION(NX,NY)     , INTENT(INOUT)           :: GSST
    INTEGER(KIND=R4KIND)   , DIMENSION(400)                                 :: NFILL

!------Over the land, the value of SST data are equal -999000000.000000 leaving shadows at the edges,
!------so, to solve the problem we fill with neighbors value
     WRITE(6,*)'Inside SST_FILL GSST MINVAL ',MINVAL(GSST)
     WRITE(6,*)'Inside SST_FILL_GSST MAXVAL ',MAXVAL(GSST)
!-     WRITE(6,*)'Inside SS_FILL_GSST LOC MINVAL ',MINLOC(GSST,DIM=1,MASK=GSST.EQ.MINVAL(GSST))
!-     WRITE(6,*)'Inside SS_FILL_GSST LOC MAXVAL ',MAXLOC(GSST,DIM=1,MASK=GSST.EQ.MAXVAL(GSST))
     WRITE(6,*)'Inside SST_FILL_GSST Undef ',Undef

     WHERE((GSST>400.).OR.(GSST<250.).AND.(GSST/=Undef))GSST=Undef 
       WRITE(6,*)'nside SST_FILL_GSST GSST MAXVAL 2 ',MAXVAL(GSST)
       IF (MAXVAL(GSST).NE.Undef) THEN
         WRITE(6,*)'nside SST_FILL_GSST NO Undef FOUND'
       ELSE
       NFILL=0  
       DO NP=1,400
          DO I=2,NX-1
           DO J=1,NY
             IF ( GSST(I,J).EQ.Undef.AND.GSST(I+1,J).ne.Undef) THEN
                GSST(I,J)=GSST(I+1,J)
                NFILL(NP)=NFILL(NP)+1
               ENDIF
           END DO
         END DO
         DO I=NX-1,2,-1
           DO J=1,NY
             IF ( GSST(I,J).EQ.Undef.AND.GSST(I-1,J).NE.Undef) THEN
                GSST(I,J)=GSST(I-1,J)    
                NFILL(NP)=NFILL(NP)+1
             END IF
           END DO
         END DO
         DO I=1,NX
           DO J=2,NY-1
             IF ( GSST(I,J).EQ.Undef.AND.GSST(I,J+1).NE.Undef) THEN
               GSST(I,J)=GSST(I,J+1)
               NFILL(NP)=NFILL(NP)+1
            END IF
           END DO
         END DO
         DO I=1,NX
           DO J=NY-1,2,-1
             IF ( GSST(I,J).EQ.Undef.AND.GSST(I,J-1).NE.Undef) THEN
               GSST(I,J)=GSST(I,J-1)
               NFILL(NP)=NFILL(NP)+1
            END IF
           END DO
         END DO
       END DO     
       WRITE(6,*)'End of Fill Data',NFILL(1),NFILL(400)
       END IF
  
     WRITE(6,*)'End of Fill Data'
     WRITE(6,*)'After FILL GSST MINVAL ',MINVAL(GSST)
     WRITE(6,*)'After FILL GSST MAXVAL ',MAXVAL(GSST)
     WRITE(6,*)'-------------------------------------------------------------------'
!---------End of Fill Data -----------!
     END SUBROUTINE SST_FILL
