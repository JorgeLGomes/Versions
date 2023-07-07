C----------------------------------------------------------------------
C  COPY CHARACTERS INTO A BIT ARRAY
C----------------------------------------------------------------------
      SUBROUTINE PKC(CHR,NCHR,IBAY,IBIT)
 
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*(*) CHR
      CHARACTER*1   CVAL(8)
      DIMENSION     IBAY(*)
      EQUIVALENCE   (CVAL,IVAL)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NCHR.GT.LEN(CHR)) CALL BORT('PKC - CHR < NCHR')
      LB = IORD(NBYTW)
      IVAL = 0
      NBIT = 8
 
      DO I=1,NCHR
      CVAL(LB) = CHR(I:I)
      IF(IASCII.EQ.0) CALL IPKM(CVAL(LB),1,IETOA(IUPM(CVAL(LB),8)))
      NWD  = IBIT/NBITW + 1
      NBT  = MOD(IBIT,NBITW)
      INT = ISHFT(IVAL,NBITW-NBIT)
      INT = ISHFT(INT,-NBT)
      MSK = ISHFT(  -1,NBITW-NBIT)
      MSK = ISHFT(MSK,-NBT)
      IBAY(NWD) = IREV(IOR(IAND(IREV(IBAY(NWD)),NOT(MSK)),INT))
      IF(NBT+NBIT.GT.NBITW) THEN
         INT = ISHFT(IVAL,2*NBITW-(NBT+NBIT))
         MSK = ISHFT(  -1,2*NBITW-(NBT+NBIT))
         IBAY(NWD+1) = IREV(IOR(IAND(IREV(IBAY(NWD+1)),NOT(MSK)),INT))
      ENDIF
      IBIT = IBIT + NBIT
      ENDDO
 
      RETURN
      END
