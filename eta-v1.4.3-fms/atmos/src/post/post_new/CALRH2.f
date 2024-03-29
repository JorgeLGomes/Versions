      SUBROUTINE CALRH2(P1,T1,Q1,ICE1,RH,IM,JM)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    CALRH2      COMPUTES RELATIVE HUMIDITY
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
C     
C ABSTRACT:  
C     THIS ROUTINE COMPUTES RELATIVE HUMIDITY GIVEN PRESSURE, 
C     TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND CLOUD 
C     ICE/WATER FLAG.  THE CODE IS BASED ON SUBROUTINE GSCOND
C     IN THE ETA MODEL.  AN UPPER AND LOWER BOUND
C     OF 100 AND 1 PERCENT RELATIVE HUMIDITY IS ENFORCED.  WHEN
C     THESE BOUNDS ARE APPLIED THE PASSED SPECIFIC HUMIDITY 
C     ARRAY IS ADJUSTED AS NECESSARY TO PRODUCE THE SET RELATIVE
C     HUMIDITY.
C   .     
C     
C PROGRAM HISTORY LOG:
C   ??-??-??  DENNIS DEAVEN
C   92-12-22  RUSS TREADON - MODIFIED AS DESCRIBED ABOVE.
C   98-06-08  T BLACK      - CONVERSION FROM 1-D TO 2-D
C   98-08-18  MIKE BALDWIN - MODIFY TO COMPUTE RH OVER ICE AS IN MODEL
C   98-12-16  GEOFF MANIKIN - UNDO RH COMPUTATION OVER ICE
C   00-01-04  JIM TUCCILLO - MPI VERSION
C     
C USAGE:    CALL CALRH2(P1,T1,Q1,ICE1,RH,IM,JM)
C   INPUT ARGUMENT LIST:
C     P1     - PRESSURE (PA)
C     T1     - TEMPERATURE (K)
C     Q1     - SPECIFIC HUMIDITY (KG/KG)
C     ICE1   - CLOUD ICE (KG/KG)
C     IM,JM  - ARRAY DIMENSIONS
C
C   OUTPUT ARGUMENT LIST: 
C     RH     - RELATIVE HUMIDITY  (DECIMAL FORM)
C     Q1     - ADJUSTED SPECIFIC HUMIDITY (KG/KG)
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C     LIBRARY:
C       NONE
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C
      INCLUDE "CTLBLK.comm"
C     
C     SET PARAMETER.
C
      PARAMETER (A2=17.2693882,A3=273.16,A4=35.86,
     &  PQ0=379.90516)
C
C     DECLARE VARIABLES.
C     
C      REAL QI,QINT,QC,P1(IM,JM),T1(IM,JM),Q1(IM,JM),RH(IM,JM)
      REAL QC,P1(IM,JM),T1(IM,JM),Q1(IM,JM),RH(IM,JM) 
C      REAL ICE1(IM,JM)
C***************************************************************
C**  THE CODE WRITTEN TO ADD IN THE ICE COMPUTATION HAS BEEN
C**     COMMENTED OUT.   SIMPLY UNCOMMENT THE SECTIONS BELOW AND
C**     COMMENT OUT THE 8 LINES USED 17 LINES BELOW THIS ONE
C**     IF YOU WISH TO ADD THE ICE COMPUTATION BACK.
C
C     START CALRH2.
C
C      DO J=JSTA,JEND
C        DO I=1,IM
C        IF (ABS(P1(I,J)).GT.1) THEN
C           TMT0=T1(I,J)-273.16
C           TMT15=AMIN1(TMT0,-15.)
C           AI=0.008855
C           BI=1.
C           IF(TMT0.LT.-20.)THEN
C             AI=0.007225
C             BI=0.9674
C           ENDIF
C           QW=PQ0/P1(I,J)
C     1          *EXP(A2*(T1(I,J)-A3)/(T1(I,J)-A4))
C           QI=QW*(BI+AI*AMIN1(TMT0,0.))
C           QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
C
      DO J=JSTA,JEND
        DO I=1,IM
        IF (ABS(P1(I,J)).GT.1) THEN
           QC=PQ0/P1(I,J)
     1          *EXP(A2*(T1(I,J)-A3)/(T1(I,J)-A4))
C
C
           if(qc.lt.1.0e-20)print*,'0 qc',i,j,pq0,P1(I,J),A2
     1,    T1(I,J),A3,A4,qc
           RH(I,J)=Q1(I,J)/QC

C   IF TEMP IS BELOW -15 C OR IF THERE IS ANY CLOUD ICE
C   AND TEMP IS BETWEEN -15 AND 0 C, THEN DO RH OVER ICE
C
C           IF(TMT0.LT.-15.)THEN
C               QC=QI
C           ELSEIF(TMT0.GE.0.)THEN
C               QC=QINT
C           ELSE
C             IF(ICE1(I,J).GT.0.0) THEN
C               QC=QI
C             ELSE
C               QC=QINT
C             ENDIF
C           ENDIF
C
C           RH(I,J)=Q1(I,J)/QC
C
C   BOUNDS CHECK
C
           IF (RH(I,J).GT.1.0) THEN
            RH(I,J)=1.0
            Q1(I,J)=RH(I,J)*QC
           ENDIF
           IF (RH(I,J).LT.0.01) THEN
            RH(I,J)=0.01
            Q1(I,J)=RH(I,J)*QC
           ENDIF
C
        ENDIF
        ENDDO
      ENDDO

      RETURN
      END
