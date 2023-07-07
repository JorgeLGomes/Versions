C----------------------------------------------------------------------
C  REVERSE A WIERD INTEGER
C----------------------------------------------------------------------
      FUNCTION IREV(N)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CINT,DINT
      EQUIVALENCE(CINT,INT)
      EQUIVALENCE(DINT,JNT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NREV.EQ.0) THEN
         IREV = N
      ELSE
         INT = N
         DO I=1,NBYTW
         DINT(I:I) = CINT(IORD(I):IORD(I))
         ENDDO
         IREV = JNT
      ENDIF
 
      RETURN
      END
