C----------------------------------------------------------------------
C  PACK UP A NUMBER ACCORDING TO SPECS
C----------------------------------------------------------------------
      SUBROUTINE PKB(NVAL,NBITS,IBAY,IBIT)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      DIMENSION IBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      NWD  = IBIT/NBITW + 1
      NBT  = MOD(IBIT,NBITW)
      IVAL = NVAL
      IF(ISHFT(IVAL,-NBITS).GT.0) IVAL = -1
      INT = ISHFT(IVAL,NBITW-NBITS)
      INT = ISHFT(INT,-NBT)
      MSK = ISHFT(  -1,NBITW-NBITS)
      MSK = ISHFT(MSK,-NBT)
      IBAY(NWD) = IREV(IOR(IAND(IREV(IBAY(NWD)),NOT(MSK)),INT))
      IF(NBT+NBITS.GT.NBITW) THEN
         INT = ISHFT(IVAL,2*NBITW-(NBT+NBITS))
         MSK = ISHFT(  -1,2*NBITW-(NBT+NBITS))
         IBAY(NWD+1) = IREV(IOR(IAND(IREV(IBAY(NWD+1)),NOT(MSK)),INT))
      ENDIF
 
      IBIT = IBIT + NBITS
 
      RETURN
      END
