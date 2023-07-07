!*********************************************************************
      subroutine thetae(im,jm,kmax,pr,temp,q,thte)
      use constants
!
!    Calculo  de theta_E
!
!    im,jm unica dimensao na horizontal, e kt=kq=kv=ku=kw na vertical
!
      integer, intent(in)                           :: im
      integer, intent(in)                           :: jm
      integer, intent(in)                           :: kmax
      real, dimension(im,jm,kmax), intent(in)       :: temp
      real, dimension(im,jm,kmax), intent(in)       :: q
      real, dimension(kmax), intent(inout)          :: pr
      real, dimension(im,jm,kmax), intent(out)      :: thte
      real, dimension(im,jm,kmax)                   :: tht
      real, dimension(im,jm,kmax)                   :: ur

      do k=1,kmax
      pr(k) = pr(k)*100.
      if(pr(k).eq.50000.) k5=k
      if(pr(k).eq.100000.) k10=k
      enddo
!
!   Calculo de theta, ummidade relativa e theta_e
!

      do k=1,kmax         !felizmente, todas as variaveis tem klyr =
      do j= 1,jm
      do i= 1,im

        tht(i,j,k) = Temp(i,j,k)*(p0/pr(k))**RDC
        qs       = q(i,j,k)
        icall = 0
        call TQadj (temp(i,j,k),qs,pr(k),icall)
        if ((Temp(i,j,k)-To).lt.0.) then
           aL=als
        else
           aL=alv
        endif
        rsat    = q(i,j,k)/qs
        ur(i,j,k) = rsat*100.
        TL= (1./((1./(temp(i,j,k)-55.))-(alog(rsat)/2840.))) + 55.
        thte(i,j,k)= tht(i,j,k)*exp((aL*(q(i,j,k)))/(Cpd*TL))

      enddo
      enddo
      enddo
!
      return 
      end subroutine thetae
