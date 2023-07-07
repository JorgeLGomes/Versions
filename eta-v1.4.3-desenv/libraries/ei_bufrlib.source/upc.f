C----------------------------------------------------------------------
C  COPY CHARACTERS FROM A BIT ARRAY
C----------------------------------------------------------------------
      SUBROUTINE UPC(CHR,NCHR,IBAY,IBIT)
 
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*(*) CHR
      CHARACTER*8   CVAL
      DIMENSION     IBAY(*)
      EQUIVALENCE   (CVAL,IVAL)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      LB = IORD(NBYTW)
      DO I=1,NCHR
      CALL UPB(IVAL,8,IBAY,IBIT)
      CHR(I:I) = CVAL(LB:LB)
      IF(IASCII.EQ.0) CALL IPKM(CHR(I:I),1,IATOE(IUPM(CHR(I:I),8)))
      ENDDO
 
      RETURN
      END
