!$$$  sub  PROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! sub  PROGRAM:  gridst   READ GLOBAL Grib SST
!                            
!
! ABSTRACT: SUBSTITUTE THE GRIBST SUB PROGRAM
!           READ THE SST FIELD IN A BINARY FORMAT 
!           INSTED A GRIB FORMAT
!           IN THIS CASE THE BINARY FILE IS GENERATE
!           FROM SST GRIB2 FILE WITH WGRIB2 PROGRAM: 
!           wgrib2 -d 1 sst2dvar_grb_0.5.grib2 -v \
!           -no_header big_endian -bin ssthiresgrd.bin
!
! PROGRAM HISTORY LOG:
!   08-11-03  JORGE LUIS GOMES
!
!
    PROGRAM READ_SST
!  
    USE F77KINDS
    USE DIAGNOSTIC
!
    IMPLICIT NONE 
!------------------------------------------------------- 
! SET PRIMARY GRID DIMENSIONS
!-------------------------------------------------------
    INTEGER(KIND=I4KIND)						    :: IM
    INTEGER(KIND=I4KIND)						    :: JM
    REAL   (KIND=R4KIND)						    :: DLMD
    REAL   (KIND=R4KIND)						    :: DPHD
    REAL   (KIND=R4KIND) 						    :: TLM0D
    REAL   (KIND=R4KIND)                                                    :: TPH0D
    
 !-------------------------------------------------------    
    INTEGER(KIND=R4KIND)                                                    :: NX
    INTEGER(KIND=R4KIND)                                                    :: NXP1
    INTEGER(KIND=R4KIND)                                                    :: NY
    INTEGER(KIND=I4KIND)                                                    :: jj
    INTEGER(KIND=I4KIND)                                                    :: kk
    INTEGER(KIND=I4KIND)                                                    :: I
    INTEGER(KIND=I4KIND)                                                    :: J
    INTEGER(KIND=I4KIND)                                                    :: N
    INTEGER(KIND=I4KIND)                                                    :: INSST
    INTEGER(KIND=I4KIND)                                                    :: OUTSST
    INTEGER(KIND=I4KIND)                                                    :: OUTNREC
    INTEGER(KIND=I4KIND)                                                    :: INLST
    INTEGER(KIND=I4KIND)                                                    :: IOstatus
    INTEGER(KIND=I4KIND)                                                    :: NREC
    INTEGER(KIND=I4KIND)                                                    :: NSMOOT
    INTEGER(KIND=I4KIND)                                                    :: i00
    INTEGER(KIND=I4KIND)                                                    :: i10
    INTEGER(KIND=I4KIND)                                                    :: i01
    INTEGER(KIND=I4KIND)                                                    :: i11
    INTEGER(KIND=I4KIND)                                                    :: j00
    INTEGER(KIND=I4KIND)                                                    :: j10
    INTEGER(KIND=I4KIND)                                                    :: j01
    INTEGER(KIND=I4KIND)                                                    :: j11
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: SSTIN
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: GSST
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: GLATR
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: GLONR   
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: SST
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: SST_R8
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: HGT
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: SM
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: SST_MSK
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: AR1
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: AR2
    REAL   (KIND=R4KIND), DIMENSION(:,:)       , ALLOCATABLE                :: AR3
    INTEGER(KIND=R4KIND), DIMENSION(:,:,:)     , ALLOCATABLE                :: INH 
    INTEGER(KIND=R4KIND), DIMENSION(:,:,:)     , ALLOCATABLE                :: JNH 
    REAL   (KIND=R4KIND)                                                    :: INCR_LAT
    REAL   (KIND=R4KIND)                                                    :: INCR_LON 
    REAL   (KIND=R4KIND)                                                    :: LAT_1
    REAL   (KIND=R4KIND)                                                    :: LON_1
    REAL   (KIND=R4KIND)                                                    :: LAT_NY
    REAL   (KIND=R4KIND)                                                    :: LON_NX
    REAL   (KIND=R4KIND)                                                    :: UNDEF 
    CHARACTER(LEN=20)                                                       :: FILE_ACCESS
    LOGICAL                                                                 :: SMOOT      
    LOGICAL                                                                 :: YREV      

!   
    REAL   (KIND=R4KIND)                                                    :: SST_MINVAL
    REAL   (KIND=R4KIND)                                                    :: SST_MAXVAL
!
    NAMELIST /PARMCONF/                                                                         & 
    & IM               ,                                                                        &
    & JM               ,                                                                        &
    & DLMD             ,                                                                        &
    & DPHD             ,                                                                        &
    & TLM0D            ,                                                                        &
    & TPH0D 
!
    NAMELIST /SSTGRD_IN/                                                                        &
    & INCR_LAT          ,                                                                       &
    & INCR_LON          ,                                                                       &
    & NX                ,                                                                       &
    & NY                ,                                                                       & 
    & LAT_1             ,                                                                       &
    & LON_1             ,                                                                       & 
    & LAT_NY            ,                                                                       &
    & LON_NX            ,                                                                       & 
    & SMOOT             ,                                                                       & 
    & NSMOOT            ,                                                                       & 
    & FILE_ACCESS       ,                                                                       &
    & YREV              ,                                                                       &
    & UNDEF

!---------------------------------------------------------------------------------------
! READ NAMELIST SSTPARMCONF WHICH DEFINE SET PRIMARY GRID DIMENSIONS
!---------------------------------------------------------------------------------------
    INLST = 11
    OPEN(INLST, FILE='NAMELIST_PARMCONF', STATUS='OLD')
    REWIND INLST
    READ(INLST,PARMCONF)
    CLOSE(11)

!---------------------------------------------------------------------------------------
! READ NAMELIST SSTGRD WHICH DEFINE THE SST INPUT PARAMETERS
!---------------------------------------------------------------------------------------
    INLST = 11
    OPEN(INLST, FILE='NAMELIST_SSTGR', STATUS='OLD')
    REWIND INLST
    READ(INLST, SSTGRD_IN)
    CLOSE(11)
    NXP1=NX+1
    
    INSST =59
    OUTSST=57
    OUTNREC=53
    write(6,*)LAT_1,LON_1,LAT_NY,LON_NX
    ALLOCATE(SSTIN(NX,NY))
    ALLOCATE(GSST(NXP1,NY))
    ALLOCATE(SST(IM,JM))
    ALLOCATE(SST_R8(IM,JM))
    ALLOCATE(GLATR(IM,JM))
    ALLOCATE(GLONR(IM,JM))
    ALLOCATE(AR1(IM,JM))
    ALLOCATE(AR2(IM,JM))
    ALLOCATE(AR3(IM,JM))
    ALLOCATE(INH(4,IM,JM))
    ALLOCATE(JNH(4,IM,JM))
    ALLOCATE(HGT(IM,JM))
    ALLOCATE(SM(IM,JM))
    ALLOCATE(SST_MSK(IM,JM))


!    READ SEA-LAND MASK ETA-GRID
    write(6,*)IM,JM
    open(10,file='etatopo.dat',status='unknown',form='unformatted')
    READ(10) HGT,SM

 !            OPEN UNIT FOR READING GRID FILE
 !            OPEN UNIT FOR WRITING EGRID FILE
     OPEN (95,file='sstin.bin',form='unformatted',convert='big_endian')
     OPEN (96,file='sstinfilled.bin',form='unformatted',convert='big_endian')
  
    CALL Geographic(GLATR, GLONR,IM,JM,DLMD,DPHD,TLM0D,TPH0D)
! 
    CALL SSTG2EGRD                                                                       &
    &(GLATR,GLONR,INCR_LAT,INCR_LON,NX,NY,LAT_1,LON_1,LAT_NY,LON_NX                      &
                       & ,AR1,AR2,AR3,INH,JNH)

       WRITE(6,*)'FILE_ACCESS ',TRIM(FILE_ACCESS)
      SELECT CASE (TRIM(FILE_ACCESS))
        CASE ('DIRECT')
          OPEN (INSST,form='unformatted',access='direct',recl=nx*ny*4)
        CASE ('SEQUENTIAL')
          OPEN (INSST,form='unformatted',convert='big_endian')
          WRITE(6,*)'FILE_ACCESS ',TRIM(FILE_ACCESS)
      END SELECT
    NREC=0
    DO
!      READ(INSST,IOSTAT=IOstatus,rec=NREC+1)((sstin(I,ny+1-J),I=1,nx),J=1,ny)
      SELECT CASE (TRIM(FILE_ACCESS))
       CASE ('DIRECT')
         IF (YREV) THEN
           READ(INSST,IOSTAT=IOstatus,REC=NREC+1)((sstin(I,J),I=1,nx),J=ny,1,-1)	 
	 ELSE  
           READ(INSST,IOSTAT=IOstatus,REC=NREC+1)((sstin(I,J),I=1,nx),J=1,ny)
	 ENDIF
       CASE ('SEQUENTIAL')
         IF (YREV) THEN    
	   READ(INSST,IOSTAT=IOstatus)((sstin(I,J),I=1,nx),J=ny,1,-1) 
	 ELSE
           READ(INSST,IOSTAT=IOstatus)((sstin(I,J),I=1,nx),J=1,ny)
	 ENDIF
         WRITE(6,*)'FILE_ACCESS ',TRIM(FILE_ACCESS),NREC
      END SELECT

      WRITE(6,*)'IOSTAT=',IOstatus
      SST_MINVAL=MINVAL(sstin(1:NX,1:NY)) 
      SST_MAXVAL=MAXVAL(sstin(1:NX,1:NY))
       
      WRITE(6,*)'SSTIN MINVAL: ',SST_MINVAL
      WRITE(6,*)'SSTIN MAXVAL: ',SST_MAXVAL
      IF (IOstatus/=0) EXIT
      NREC=NREC+1
      WRITE(6,*)'READ SST REC ',NREC, 'NX:',NX,'NY:',NY           
      GSST(1:NX,1:NY) = SSTIN(1:NX,1:NY)
!...   add greenich to right side of grid
!
      GSST(NXP1,1:NY) = GSST(1,1:NY)
      SST_MINVAL=MINVAL(GSST(1:NX,1:NY)) 
      SST_MAXVAL=MAXVAL(GSST(1:NX,1:NY))
       
      WRITE(6,*)'SSTIN MINVAL: ',SST_MINVAL
      WRITE(6,*)'SSTIN MAXVAL: ',SST_MAXVAL
      WRITE(95)GSST
      CALL SST_FILL(GSST,NXP1,NY,UNDEF)
      WRITE(96)GSST
      SST=-9.99E+20
      SST_MSK=0
!      DO CONCURRENT (I=1:IM, J=1:JM, SM(I,J) > .5)
      DO J=1,JM
        DO I=1,IM
           IF (SM(I,J) > .5) THEN
              i00=inh(1,i,j)
              i10=inh(2,i,j)
              i01=inh(3,i,j)
              i11=inh(4,i,j)
!
              j00=jnh(1,i,j)
              j10=jnh(2,i,j)
              j01=jnh(3,i,j)
              j11=jnh(4,i,j)
!
              SST_R8(I,J) =            (GSST(i00,j00))                                 &
     &                      +AR1(I,J)*((GSST(i10,j10))- (GSST(i00,j00)))           &
     &                      +AR2(I,J)*((GSST(i01,j01))- (GSST(i00,j00)))           &
     &                      +AR3(I,J)*((GSST(i00,j00)) - (GSST(i10,j10))           &
     &                                -(GSST(i01,j01)) + (GSST(i11,j11)) )
                            SST_MSK(I,J)=1
           ENDIF
         ENDDO
       ENDDO
      WRITE(6,*)"Done!"
!     SST_R8 to SST_R4
      SST(1:IM,1:JM) = SST_R8(1:IM,1:JM)
      SST_MINVAL=MINVAL(SST(1:IM,1:JM),SM>0.5) 
      SST_MAXVAL=MAXVAL(SST(1:IM,1:JM),SM>0.5)
      IF(SMOOT)  THEN
       DO N=1,NSMOOT
         DO J=2,JM-1
           DO I=2,IM-1
             IF (SM(I,J) > .5) THEN
                SST(I,J)=     (   SST(I,J)*SM(I,J)                                                 &
                &              + SST(I+1,J  )*SM(I+1,J  )                                          &
                &              + SST(I-1,J  )*SM(I-1,J  )                                          &
                &              + SST(I  ,J+1)*SM(I  ,J+1)                                          &
                &              + SST(I  ,J-1)*SM(I  ,J-1)                                          &
                &              + SST(I+1,J+1)*SM(I+1,J+1)                                          &
                &              + SST(I+1,J-1)*SM(I+1,J-1)                                          &
                &              + SST(I-1,J+1)*SM(I-1,J+1)                                          &
                &              + SST(I-1,J-1)*SM(I-1,J-1))                                         &
                &            /(  SM(I,J)+SM(I+1,J  )+SM(I-1,J  )+SM(I  ,J+1)                       &
                &              + SM(I  ,J-1)+SM(I+1,J+1)+SM(I+1,J-1)                               &
                &              + SM(I-1,J+1)+SM(I-1,J-1))     
             ENDIF
           ENDDO
         ENDDO
       ENDDO
      ENDIF
     
      WRITE(OUTSST)SST
       
      SST_MINVAL=MINVAL(SST(1:IM,1:JM),SM>0.5) 
      SST_MAXVAL=MAXVAL(SST(1:IM,1:JM),SM>0.5)
       
      WRITE(6,*)'SST EGRD MINVAL: ',SST_MINVAL
      WRITE(6,*)'SST EGRD MAXVAL: ',SST_MAXVAL
      WRITE(6,*)'-------------------------------------------------------------------'
      WRITE(6,*)'-------------------------------------------------------------------'

    END DO       
 
    WRITE(OUTNREC,"(I4)")NREC
    WRITE(6,*)'NUMBER OF SST REC ',NREC           
    WRITE(6,*) 'leaving GRIDST'
   END PROGRAM READ_SST
