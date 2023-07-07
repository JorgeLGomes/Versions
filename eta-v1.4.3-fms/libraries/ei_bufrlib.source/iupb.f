C----------------------------------------------------------------------
C  UNPACK UP AN INTEGER FROM A PACKED INTEGER ARRAY (FUNCTION)
C----------------------------------------------------------------------
      FUNCTION IUPB(MBAY,NBYT,NBIT)
 
      DIMENSION MBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      MBIT = (NBYT-1)*8
      CALL UPB(IRET,NBIT,MBAY,MBIT)
      IUPB = IRET
      RETURN
      END
