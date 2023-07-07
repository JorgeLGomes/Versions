      SUBROUTINE CHR2INT(CHR,ILEN,IVAL)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    CHR2INT     CONVERTS CHR STRING TO INTEGER
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-31
C     
C ABSTRACT:  
C     GIVEN A CHARACTER STRING, CHR, OF LENGTH LEN CONTAINING
C     INTEGERS OR BLANKS, THIS ROUTINE RETURNS THE INTEGER
C     VALUE IN IVAL.
C   .     
C     
C PROGRAM HISTORY LOG:
C   92-12-31  RUSS TREADON
C     
C USAGE:    CALL CHR2INT(CHR,ILEN,IVAL)
C   INPUT ARGUMENT LIST:
C     CHR      - CHARACTER STRING CONTAINING INTEGERS OR BLANKS.
C     LEN      - LENGTH OF CHR.
C
C   OUTPUT ARGUMENT LIST: 
C     IVAL     - INTEGER EQUIVALENT TO CHR.
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       NONE
C     LIBRARY:
C       NONE
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C     
C     SET PARAMETER IOFF.  IOFF IS THE ASCII CODE FOR
C     THE ZERO CHARACTER.  CHARACTERS 48 THROUGH 57
C     ARE INTEGER 0 THROUGH 9, RESPECTIVELY.
C
      PARAMETER (IOFF=48)
C
C     DECLARE VARIABLES.
C     
      CHARACTER*10 CHR
C
C     
C***********************************************************************
C     START CHR2INT HERE.
C     
C     LOOP OVER LENGTH OF CHARACTER STRING. CONSTRUCT
C     THE INTEGER EQUIVALENT TO THE INTEGER IN THE 
C     CHARACTER STRING.
C     
      ISUM = 0
      DO 10 IPOS = ILEN,1,-1
         ICVAL = ICHAR(CHR(IPOS:IPOS))
         IPOWER = ILEN-IPOS
         IF ((ICVAL.GE.48).AND.(ICVAL.LE.57)) THEN
            INTVAL = ICVAL-IOFF
            ISUM = ISUM + INTVAL*10**IPOWER
         ENDIF
 10   CONTINUE
C     
C     SET INTEGER VALUE.
C
      IVAL = ISUM
C     
C     END OF ROUTINE.
C     
      RETURN
      END
