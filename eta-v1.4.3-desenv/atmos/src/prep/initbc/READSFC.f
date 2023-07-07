      SUBROUTINE READSFC(HLAT,HLON,SM,SICE,SMC,STC,TG)
C     ******************************************************************
C
C	Needed inputs:
C
C	HLAT(IM,JM), HLON(IM,JM)
C	SM(IM,JM), SICE(IM,JM)
C
C	OUTPUTS:
C
C	SMC(IM,JM,NSOIL),STC(IM,JM,NSOIL),TG(IM,JM)
C
C----------------------------------------------------------------------
      integer,parameter::real_32=selected_real_kind(6,30)
      INCLUDE "parmeta"
      INCLUDE "parmsoil"
C----------------------------------------------------------------------
       REAL SICE (IM,JM), SM(IM,JM),SMC(IM,JM,NSOIL),STC(IM,JM,NSOIL)
       REAL HLAT(IM,JM),HLON(IM,JM)
                             D I M E N S I O N
     &  IGDATE(4)
C
       REAL(REAL_32),ALLOCATABLE,DIMENSION(:)::GLATS
       REAL(REAL_32),ALLOCATABLE,DIMENSION(:,:)::GGTG,GGSLI,GGRID1
       REAL(REAL_32),ALLOCATABLE,DIMENSION(:,:,:)::GGSMC2,GGSTC2

C	real*8 etime
C
       CHARACTER*32 GLABEL
C-----------------------------------------------------------------------
      DPR=180./3.141592654
C
C	btime=etime()
	open(unit=33,file='sfcanl',form='unformatted',
     +		access='sequential')	
Csomelinux	open(unit=33,file='sfcanl',form='unformatted',
Csomelinux     +		access='sequential',convert='big_endian')	
      READ (33) 
      READ (33) GFHR,IGDATE,NGLONS,NGLATS,NUM
      print *,' READSFC: nglats,nglons = ',gfhr,igdate,nglats,nglons,num
      ALLOCATE(GGRID1(NGLONS,NGLATS))
      ALLOCATE(GLATS(NGLATS))
      ALLOCATE(GGSMC2(NGLONS,NGLATS,2))
      ALLOCATE(GGSTC2(NGLONS,NGLATS,2))
      ALLOCATE(GGTG(NGLONS,NGLATS))
      ALLOCATE(GGSLI(NGLONS,NGLATS))
C

	CALL GAULAT(glats,nglats)

C	write(6,*) 'back from GAULAT'
	do J=1,nglats,4
C	write(6,*) 'J,glats(J)= ', J,glats(J)
	enddo


C
C NOW READ IN THE FIELDS, NOTE WE NEED THE GLOBAL SEA,LAND ICE MASK
C TO DO A PROPER INTERPOLATION.
C FIRST TWO RECORDS NOW READ IN ECONVK
c     OPEN(33,FILE='fort.33',
c    &     FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IRET)
C
      print *,'ok after gfhr in sfcanl'
C  BYPASS THE SURFACE TEMPERATURE RECORD
      READ (33) GGRID1
      print *,'ok after rec 1 in sfcanl'
C  READ AND SAVE THE TWO-LAYER SOIL MOISTURE RECORDS
      READ (33) GGSMC2
      print *,'ok after rec 2 in sfcanl'
C  BYPASS THE SNOW DEPTH RECORD
      READ (33) GGRID1
      print *,'ok after rec 3 in sfcanl'
C  READ AND STORE THE THE SOIL TEMPERATURES (TWO-LAYERS)
      READ (33) GGSTC2
      print *,'ok after rec 4 in sfcanl'
C  READ AND STORE THE DEEP LAYER SOIL TEMPERATURE
      READ (33) GGTG
      print *,'ok after rec 5 in sfcanl'
C  BYPASS 5 MORE RECORDS AND THEN READ AND STORE THE SEA,LAND ICE RECORD
C  BYPASSING Z0,CONV. CLD FRACT.,CONV. CLD BASE,CONV CLD TOP, AND ALBEDO
      DO KRD = 1,6
      READ (33) GGSLI
      END DO
      print *,'ok after rec 6-12 in sfcanl'

C	read_tim=etime()-btim

C		write(6,*) 'read time ', read_tim

C	btim=etime()
C	write(6,*) 'some gautoeta input arguments'
C	write(6,*) 'GLATS(25) ', glats(25)
C	write(6,*) 'nglats,nglons= ', nglats,nglons
C	write(6,*) 'HLAT,HLON: ', HLAT(25,25),HLON(25,25)
C	write(6,*) 'SM,SICE : ', SM(25,25),sice(25,25)
C	write(6,*) 'GGSLI,GGSMC2 ', GGSLI(25,25),GGSMC2(25,25,1)
      CALL GAUTOETA(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGSMC2(1,1,1),SMC(1,1,1))
      CALL GAUTOETA(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGSTC2(1,1,1),STC(1,1,1))
      CALL GAUTOICE(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGSTC2(1,1,1),STC(1,1,1))
Cmp
C	write(6,*) 'gaussian soil level 1'
	do J=nglats,1,-30
C	write(6,666) (GGSTC2(I,J,1),I=1,nglons,40)
	enddo
C	write(6,*) 'gaussian soil level 2'
	do J=nglats,1,-30
C	write(6,666) (GGSTC2(I,J,2),I=1,nglons,40)
	enddo
  666	format(10(f5.1,1x))
Cmp
      DO KK=2,NSOIL
       CALL GAUTOETA(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGSMC2(1,1,2),SMC(1,1,KK))
       CALL GAUTOETA(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGSTC2(1,1,2),STC(1,1,KK))
       CALL GAUTOICE(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGSTC2(1,1,2),STC(1,1,KK))
      ENDDO
Cmp
C	write(6,*) 'layer 1 STC '
	do J=JM,1,-4
C	write(6,633) (STC(I,J,1),I=1,IM,IM/18)
	enddo
C	write(6,*) 'layer 4 STC '
	do J=JM,1,-4
C	write(6,633) (STC(I,J,4),I=1,IM,IM/18)
	enddo
  633	format(20(f4.0,1x))
C
C NOW SET LOWER BOUNDARY CONDITION ON SUB-SURFACE
C TEMPERATURE, TREAT UNDER SEA-ICE AFTER CALL
      CALL GAUTOETA(GLATS,NGLATS,NGLONS,HLAT,HLON,SM,SICE,GGSLI,
     &              GGTG,TG)
C
      DEALLOCATE(GGRID1)
      DEALLOCATE(GLATS)
      DEALLOCATE(GGSMC2)
      DEALLOCATE(GGSTC2)
      DEALLOCATE(GGTG)
      DEALLOCATE(GGSLI)
C-----------------------------------------------------------------------
                             RETURN
                             END
