
!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE DDFF2UV                    ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
     SUBROUTINE ddff2uv(dd,ff,u,v)
!
!#######################################################################
!
!     PURPOSE:
!
!     Calculate u and v wind components from direction and speed.
!
!#######################################################################
!
!
!     AUTHOR: Keith Brewster
!     3/11/1996
!
!     09/feb/2004 (ERNANI L. NASCIMENTO)
!       THIS SUBROUTINE IS EXACTLY LIKE SUBROUTINE  ddff2uv IN
!     /src/adas/thermo3d.f  I´VE JUST CHANGED THE NAME.
!
!#######################################################################
!
  implicit none
  integer j,k
  real :: u,v,dd,ff
!   dimension  u(100,35),v(100,35),dd(100,35),ff(100,35)
  real :: arg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!   arg = (dd(j,k) * (3.141592654/180.))
  arg = (dd * (3.141592654/180.))
!   if (m==1) then
!    write(*,*) 'Dentro do ddff2uv: ',dd,ff,arg
!   endif
!   u(j,k) = -ff(j,k) * sin(arg)
  u = -ff * sin(arg)
!   v(j,k) = -ff(j,k) * cos(arg)
  v = -ff * cos(arg)
  RETURN
END SUBROUTINE ddff2uv
