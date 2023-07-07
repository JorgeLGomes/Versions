      program read_degrib
!
      implicit none
!  
      parameter(nz=21,nxny=22608)
      integer nxny,nz,i,datsav,n,l
!
      real variab(21,10,22608),dummy,pr(nz),tmp(nxny,1),tmp2(nxny)
!
      integer kgds(200),kpds(50),len,kerr,jpds(200),jgds(200) &
     &         ,lenpds,lenkgds,nwords,kpdsl &
     &        ,j,k,nx,IRETO,KNUM,IRET1
!

      character(LEN=1):: pds(50)
      character(LEN=12):: atime
      character(LEN=2):: gproj
      logical bitmap(nxny)
      common /gdsinfo/datsav(11)
!_______________________________________________________________________________
!
!
! *** Check that the input dimensions match grib file dimensions. 
!
!	write(6,*) 'inside read_degrib '

!        write(*,*)nxny,nz
!	write(*,*)pr(1:nz)
!  Should be able to derive all header info in the DEGRIB file
! from the GDS/PDS of the GRIB file.  
!
!
!
! *** Fill time stamp (assume all records are for same time).
!
! *** Fill a missing value flag into first space of each variable.
!
!       do n=1, 10
!        do l=1, nz
!	 do k=1, nxny
!           write(*,*)n,l,k,variab(l,n,k)   
!	 enddo
!	enddo
!       enddo	   
      stop
!
      end
!
