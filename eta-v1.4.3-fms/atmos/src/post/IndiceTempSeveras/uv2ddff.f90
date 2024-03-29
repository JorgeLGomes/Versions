!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UV2DDFF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uv2ddff(u,v,dd,ff)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate direction and speed from u and v.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  3/11/1996
!
!-----------------------------------------------------------------------
!
 IMPLICIT NONE
 REAL :: u,v,dd,ff
 REAL :: dlon
 REAL :: r2deg
 PARAMETER (r2deg=180./3.141592654)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
 ff = SQRT(u*u + v*v)

 IF(v > 0.) THEN
   dlon=r2deg*ATAN(u/v)
 ELSE IF(v < 0.) THEN
   dlon=180. + r2deg*ATAN(u/v)
 ELSE IF(u >= 0.) THEN
   dlon=90.
 ELSE
   dlon=-90.
 END IF

 dd= dlon + 180.
 dd= dd-360.*(nint(dd)/360)
 RETURN
END SUBROUTINE uv2ddff
