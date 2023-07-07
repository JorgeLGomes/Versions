!-------------------------------------------------------------------
!     TQadj
!-------------------------------------------------------------------
      subroutine TQadj(T,Q,P,icall)
      use constants
      real, intent(inout)  :: T
      real, intent(inout)  :: Q
      real, intent(in)     :: P
      integer, intent(in)  :: icall

      itime=0
      do itime=1,4
       if (T-To.lt.0.) then
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
       rs= (c2es/p)*exp(c3*(t-To)/(t-c4))
       rs= min(0.5,rs)
       cor = 1./(1.- vtmpc1 * rs)
       qs = rs*cor
       if (icall.eq.0) then
         zcond = 0.
         Q = qs
         exit
       elseif (icall.eq.1) then
         zcond = (q-qs)/(1.+c5*qs*cor*(1./(T-c4))**2)
         zcond = max(zcond,0.)
         T = T + (aL/Cpd) * zcond
         Q = Q - zcond
         if (zcond.eq.0.) exit
       endif
      enddo
      return
      end subroutine TQadj
