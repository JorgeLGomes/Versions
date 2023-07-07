      subroutine rh_fix(nx,ny,nz,rh,pr)
      implicit none
!
      integer         :: nx,ny,nz,i,j,k,kk
!
      real, dimension(nx,ny,nz)           :: rh
      real, dimension(nz)                 :: pr
!_______________________________________________________________________________
!

! *** Fix bottom levels if necessary.
!
      if (rh(1,1,1) .eq. -99999.) then
	 do k=2,nz
	    if (rh(1,1,k) .ne. -99999.) then

	       DO kk=k-1,1,-1
	       print *,'Copying',nint(pr(kk+1)),' mb to' , nint(pr(kk)),' mb.'
	       do j=1,ny
	       do i=1,nx
		  rh(i,j,kk)=rh(i,j,kk+1)
               enddo
               enddo

               ENDDO

	       goto 10
            endif
         enddo
	 print *,'RH patch did not work.'
	 stop
      endif
!
! *** Fix upper levels if necessary.
!
10    continue
      if (rh(1,1,nz) .eq. -99999.) then
	 do k=nz-1,1,-1
	write(6,*) 'checking RH at lev= ', pr(k)
	    if (rh(1,1,k) .ne. -99999.) then
	       do kk=k+1,nz
	       print *,'Copying',nint(pr(kk-1)),' mb to', nint(pr(kk)),' mb.'
	       do j=1,ny
	       do i=1,nx
		  rh(i,j,kk)=rh(i,j,kk-1)
               enddo
               enddo
               enddo
	       goto 20
            endif
         enddo      
      endif
!
20    continue
      do k=1,nz
	 if (rh(1,1,k) .eq. -99999.) then
	    print *,'RH patch did not work.'
	    stop
	 endif
      enddo
!
      return
      end