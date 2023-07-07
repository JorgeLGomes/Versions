      FUNCTION WDIR(X,Y)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    WDIR        CMPTE VEC DIR FROM U-V
C   PRGRMMR: TREADON         ORG: W/NMC2     DATE: 93-05-12       
C     
C ABSTRACT:  
C     GIVEN U AND V WIND COMPONENTS, THIS FUNCTION CALCULATES
C     THE VECTOR WIND DIRECTION (IN DEGREES)
C   .     
C     
C PROGRAM HISTORY LOG:
C   ??-??-??  CHRIS PETERS ??
C   93-05-12  RUSS TREADON - ADDED DOCBLOC
C     
C USAGE:   WIND = WDIR(X,Y)
C   INPUT ARGUMENT LIST:
C     X        - U WIND COMPONENT
C     Y        - V WIND COMPONENT
C
C   OUTPUT ARGUMENT LIST: 
C     WDIR     - VECTOR WIND DIRECTION IN DEGREES
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       DIR    - CONPUTE VECTOR WIND DIRECTION FROM U-V WIND.
C     LIBRARY:
C       NONE
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY Y-MP
C$$$  
C     
C**********************************************************************
C     START WDIR HERE.
C
      RDIR = 270.-DIR(X,Y)
      IF (RDIR.LT.0.) RDIR = RDIR+360.
      IF (RDIR.GE.360.) RDIR = RDIR-360.
      WDIR = RDIR
      RETURN
      END
