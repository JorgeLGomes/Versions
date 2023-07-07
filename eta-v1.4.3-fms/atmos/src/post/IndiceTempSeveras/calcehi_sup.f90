!
!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE CALCEHI_SUP                ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######        and by Ernani L. Nascimento at SIMEPAR        ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
SUBROUTINE calcehi_sup(nrad,lm,nlev,presshpa,cape,heli,dnrv,ehi,sup)
!
!#######################################################################
!
!     PURPOSE:
!
!     Calculate the energy helicity-index (EHI) (DAVIES 1993; RASMUSSEN AND
!     BLANCHARD 1998), and the supercell composite index (THOMPSON et al
!     2003).
!
!#######################################################################
!
!     AUTHOR: Ernani L. Nascimento, adapting code from arpspltderive.f
!     09 February, 2004
!
!     MODIFICATION HISTORY:
!     THIS IS A NEW SUBROUTINE, NOT AVAILABLE IN THE ORIGINAL ARPS CODE
!
!     (04/19/2004) (Ernani L. Nascimento)
!     Code was modified to fit in índices_severos3.f90
!
!     (08/02/2004) (Ernani L. Nascimento)
!     Code was modified to include the computation of the supercell
!     composite index.
!
!#######################################################################
!
!     Calculates EHI and SUP from cape, helicity and bulk Richardson
!     shear.
!
!
!     nrad          Number of soundings for which EHI will be computed
!     cape          Convective available potential energy (J/kg)
!     heli          Helicity, storm relative (m²/s²)
!     dnrv          Bulk Richardson number shear (m²/s²)
!     ehi           Energy-helicity index (non-dimensional)
!     sup           Supercell composite index
!
!#######################################################################
!
   implicit none
!
!#######################################################################
!
!     Input variables
!
!#######################################################################
!
   integer :: nrad,lm,k,j
   integer :: nlev(nrad)
   real :: cape(nrad), heli(nrad), dnrv(nrad)
   real :: presshpa(lm,nrad)
!
!#######################################################################
!
!     Output variables
!
!#######################################################################
!
   real :: ehi(nrad)      ! Energy-helicity index
   real :: sup(nrad)      ! Supercell composite index

!#######################################################################
!     Work variable (constant)
!#######################################################################

   real :: denomi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   denomi=160000.0
   DO  k=1,nrad
    if (nlev(k)==0) then
     ehi(k)=-999.9
     sup(k)=-999.9
    else
     j=nlev(k)
     if (presshpa(j,k)>300.0) then
      cycle
     endif
     ehi(k) = ( cape(k)*heli(k) )/denomi
     sup(k) = ( (cape(k)/1000.)*(heli(k)/150.)*(dnrv(k)/40.) )
    endif
   ENDDO

 RETURN
END SUBROUTINE calcehi_sup
