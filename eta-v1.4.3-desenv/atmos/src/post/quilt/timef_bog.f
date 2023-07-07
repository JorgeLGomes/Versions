C
C  	function written early Dec. 1999 by M. Pyle to support workstation
C	Eta for users with etime but not timef functionality (like certain
Cmp	HPs)  Designed to duplicate timef (elapsed time in milliseconds)
C	
	function timef()
	real et(2)
	real*8 timef
	timef=etime(et)
	timef=timef*1.e3
	end
