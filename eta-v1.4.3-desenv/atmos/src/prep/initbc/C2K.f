C-CP
      SUBROUTINE C2K(M,N,A)
C-CP
C-CP CONVERTS AN ARRAY OF TEMPERATURES FROM CENTIGRADE TO KELVIN
C-CP
      DIMENSION A(M,N)
C
      IMAX=M
      JMAX=N
C
      doout50: DO J =1,JMAX
      doin50: DO  I =1,IMAX
        A(I,J) = A(I,J) + 273.16
      END DO doin50
      END DO doout50
C
      RETURN
      END
