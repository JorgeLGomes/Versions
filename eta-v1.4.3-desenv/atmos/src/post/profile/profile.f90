      program sfchist
!
!  generalized for use for any vertical resolution (up to 50 layers)
!  and for class 1 data
!
!----------------------------------------------------------------------
!---read number of stations of the file staids ------------------------
!----------------------------------------------------------------------
!
      implicit none
      integer :: nsta,npnt  
!
 333  read(5,*) nsta,npnt
!      write(6,*)"nsta ",nsta," npnt ", npnt
      call get_dat(nsta,npnt)
!
!      stop	
      end program sfchist
!
      subroutine get_dat(nsta,npnt)
      use parmeta
      
!      LCL0ML - NUMBER OF MULTI-LAYER VARIABLES OUTPUT FOR CLASS 0
!      LCL1ML - NUMBER OF MULTI-LAYER VARIABLES OUTPUT FOR CLASS 1
!      LCL0SL - NUMBER OF SINGLE LAYER VARIABLES OUTPUT FOR CLASS 0
!      LCL1SL - NUMBER OF SINGLE LAYER VARIABLES OUTPUT FOR CLASS 1
!
      implicit none
      integer :: nsta,ierr,ier,nhr,ns,nrec,ihrst,ifcst,id,n,icls,lmh,i !isr:declaracoes incluidas
      integer, parameter            :: LCL0ML=6            
      integer, parameter            :: LCL1ML=15
      integer, parameter            :: LCL1SL=50
      integer, parameter            :: LRECPR=4*(8+9+LCL1ML*LM+LCL1SL)
      integer, parameter            :: nword=999
      integer                       :: npnt
      real   , parameter            :: DEG=180./3.14159265
      real   , parameter            :: zero=0.
      real   , parameter            :: t0=273.15  
      real   , dimension(nword)     :: fpack
      real   , dimension(:,:), allocatable :: slp, pds, alt, tsf, tmi, tma, ppt, cup, sm
      real   , dimension(:,:), allocatable :: u10m, v10m, u, v, t2ms, rh2ms, q2ms
      real   , dimension(:,:), allocatable :: cfracl, cfracm, cfrach, aswin, alwin
      integer, dimension(3)         :: idate
      character(len=8)              :: CISTAT
!
      print*,"LRECPR ",LRECPR, "LCL1ML ", LCL1ML,"LM ", LM,"LCL1SL ",LCL1SL
      allocate(slp(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "slp : allocation failed"
      allocate(pds(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "pds : allocation failed"
      allocate(alt(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "alt : allocation failed"
      allocate(tsf(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "tsf : allocation failed"
      allocate(tmi(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "tmi : allocation failed"
      allocate(tma(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "tma : allocation failed"
      allocate(ppt(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "ppt : allocation failed"
      allocate(cup(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "cup : allocation failed"
      allocate(sm(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "sm : allocation failed"
      allocate(u10m(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "u10m : allocation failed"
      allocate(v10m(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "v10m : allocation failed"
      allocate(u(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "u : allocation failed"
      allocate(v(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "v : allocation failed"
      allocate(t2ms(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "t2ms : allocation failed"
      allocate(aswin(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "aswin : allocation failed"
      allocate(alwin(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "alwin : allocation failed"
      allocate(rh2ms(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "rh2ms : allocation failed"
      allocate(q2ms(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "q2ms : allocation failed"
      allocate(cfracl(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "cfracl : allocation failed"
      allocate(cfracm(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "cfracm : allocation failed"
      allocate(cfrach(npnt,nsta), stat=ierr)
      if (ierr /= 0) print*, "cfrach : allocation failed"
!
      open(66,file='profilm.c1.t00s',ACCESS='DIRECT',RECL=LRECPR,IOSTAT=IER,status='old')
      open(30,status='unknown',form='unformatted')
!      open(30,status='unknown',form='unformatted',iostat=istat)
  
      DO nhr=1,npnt
      DO NS=1,nsta     
        NREC=(nhr-1)*nsta+ns
!
!        print*,"NREC ",NREC,LRECPR,nhr,NS
        read(66,REC=NREC) ihrst,idate,ifcst,id,CISTAT,(fpack(n),n=1,9), &
                          (fpack(n),n=10,LCL1ML*fpack(4)+LCL1SL+9)
!isr    id=ns !isr: nao comentar essa linha se o numero de estacoes for diferente do valor maximo,.e.g.,147.       
!       print*,"Passei!!",NREC,ihrst,idate,ifcst,id,CISTAT,LCL1ML*fpack(4)+LCL1SL+9, LRECPR
!       nhr=nint(ifcst/3600.)+1
        icls=nint(fpack(8))
        lmh=nint(fpack(4))            
        alt(nhr,id)=nint(fpack(3)) 
!	  if (lmh>=36) then
!          u(nhr,id)=(fpack(1+9+lmh*2)+fpack(2+9+lmh*2))*0.5  !(?+9+lmh*2) ? => nivel
!          v(nhr,id)=(fpack(1+9+lmh*3)+fpack(2+9+lmh*3))*0.5
!	  else
!          u(nhr,id)=fpack(1+9+lmh*2)           !(?+9+lmh*2) ? => nivel
!          v(nhr,id)=fpack(1+9+lmh*3)
!       endif
        u(nhr,id)=fpack(1+9+lmh*2)           !(?+9+lmh*2) ? => nivel
        v(nhr,id)=fpack(1+9+lmh*3)
        slp(nhr,id)=fpack(LCL1ML*lmh+10)    ! MSLP
        pds(nhr,id)=fpack(LCL1ML*lmh+11)    ! Sfc Pressure
        tsf(nhr,id)=fpack(LCL1ML*lmh+12)    ! Sfc Temperature    
        tmi(nhr,id)=fpack(LCL1ML*lmh+13)
        tma(nhr,id)=fpack(LCL1ML*lmh+14)
        ppt(nhr,id)=fpack(LCL1ML*lmh+16)    ! Accum Total Precip.
        cup(nhr,id)=fpack(LCL1ML*lmh+17)    ! Accum Convc Precip.
        u10m(nhr,id)=fpack(LCL1ML*lmh+37)
        v10m(nhr,id)=fpack(LCL1ML*lmh+38)
        t2ms(nhr,id)=fpack(LCL1ML*lmh+41)
        q2ms(nhr,id)=fpack(LCL1ML*lmh+42)
        sm(nhr,id)=fpack(LCL1ML*lmh+54)
        cfracl(nhr,id)=fpack(LCL1ML*lmh+55)
        cfracm(nhr,id)=fpack(LCL1ML*lmh+56)
        cfrach(nhr,id)=fpack(LCL1ML*lmh+57)
	aswin(nhr,id)=fpack(LCL1ML*lmh+23)
	alwin(nhr,id)=fpack(LCL1ML*lmh+25)
!       print*,nhr,id,alt(nhr,id),slp(nhr,id),tsf(nhr,id),ppt(nhr,id),t2ms(nhr,id)
!       print*,nhr,id,alt(nhr,1),slp(nhr,1),tsf(nhr,1),ppt(nhr,1),t2ms(nhr,1)
        if (ppt(nhr,id).lt.0.)ppt(nhr,id)=0.
!isr    if (id.eq.20) write(110,*)nhr,id,alt(nhr,id),slp(nhr,id),tsf(nhr,id),ppt(nhr,id),t2ms(nhr,id)
      enddo
      enddo      
!	 
!-----------------------------------------------------------------------
!
!----- Print alt and lsm
!
!-----------------------------------------------------------------------
      write(33,*)'      alt       sm'
      do i=1,nsta
      write(33,'(i4,2f9.2)') i,alt(1,i),sm(1,i)
      enddo
!
!-----------------------------------------------------------------------
!----- Calculating relative humidity -----------------------------------
!-----------------------------------------------------------------------
! *** using sfc temperature instead of 2m temp. ***
!  **  until 2m temp is fixed   ***

      call calrh(pds,tsf,q2ms,rh2ms,npnt,nsta)
!
!-----------------------------------------------------------------------
!----- time to write out       -----------------------------------------
!-----------------------------------------------------------------------
      do i=1,npnt
      write(30) (slp(i,id)/100,id=1,nsta)
!isr  write(112,*)i,(id,slp(i,id)/100,id=1,nsta)
      write(30) (tsf(i,id),id=1,nsta)
      write(30) (ppt(i,id),id=1,nsta)
      write(30) (t2ms(i,id),id=1,nsta)
      write(30) (aswin(i,id),id=1,nsta)
      write(30) (alwin(i,id),id=1,nsta)  
      write(30) (q2ms(i,id),id=1,nsta)
      write(30) (rh2ms(i,id),id=1,nsta)
      write(30) (u10m(i,id),id=1,nsta)
      write(30) (u10m(i,id),id=1,nsta)
      write(30) (v10m(i,id),id=1,nsta)
      write(30) (v10m(i,id),id=1,nsta)
      write(30) (cfracl(i,id),id=1,nsta)
      write(30) (cfracm(i,id),id=1,nsta)
      write(30) (cfrach(i,id),id=1,nsta)
      write(30) (u(i,id),id=1,nsta)
      write(30) (u(i,id),id=1,nsta)
      write(30) (v(i,id),id=1,nsta)
      write(30) (v(i,id),id=1,nsta)
      enddo
!	 	
      close(30)
!
      return
      end subroutine get_dat
!
      subroutine calrh(Ps,T,Q,RH,npnt,nsta)
      real, parameter            :: To= 273.16
      real, parameter            :: Cpd= 1005.46
      real, parameter            :: VTMPC1= .60776868
      real, parameter            :: C2ES= 379.892958
      real, parameter            :: C3IES= 21.875
      real, parameter            :: C3LES= 17.269
      real, parameter            :: C4IES= 7.66 
      real, parameter            :: C4LES= 35.86
      real, parameter            :: C5IES= 5807.8125
      real, parameter            :: C5LES= 4097.9337
      real, parameter            :: als= 2834500.
      real, parameter            :: alv= 2500800. 
      real, dimension(npnt,nsta) :: Ps, T, RH, Q  

!      write(7,*) ' starting to calculate rel humidity '
!      print*, nsta,npnt
      do i=1,npnt
        do j=1,nsta
          do it=1,2
            if (T(i,j)-To.lt.0.) then
              c3=C3IES
              c4=C4IES
              c5=C5IES*als/cpd
              aL=als
            else
              c3=C3LES
              c4=C4LES
              c5=C5LES*alv/cpd
              aL=alv
            endif
            rs=(c2es/Ps(i,j))*exp(c3*(T(i,j)-To)/(T(i,j)-c4))
            cor=1./(1.- vtmpc1 * rs)
            qs=rs*cor
            zcond=(q(i,j)-qs)/(1.+c5*qs*cor*(1./(T(i,j)-c4))**2)
            zcond=max(zcond,0.)
!
!  correct T and Q is supersaturation occurs
!
            T(i,j)=T(i,j)+(aL/Cpd)*zcond
            Q(i,j)=Q(i,j)-zcond
          enddo      
          RH(i,j)=(Q(i,j)/qs)*100.
!          RH(i,j)=(Q(i,j)/rs)*100.
         enddo
        enddo
      return
      end subroutine calrh

