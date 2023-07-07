      PROGRAM GETGRD
C*
      INTEGER LGRB,LREC,LENF,IN,MD(12),C,
     *        I,T,YI,MI,DI,HI,Y,M,D,H,FF
      CHARACTER FNAME*256
      LOGICAL EX
      DATA MD /31,28,31,30,31,30,31,31,30,31,30,31/
C*
      CALL GETARG(1,FNAME)
      IN=INDEX(FNAME//' ',' ')-1
      T=IN-24 !identifica o inicio da data da condicao inicial do arquivo .grb
C*
      READ(FNAME(T:T+3),'(I4)')YI
      READ(FNAME(T+4:T+5),'(I2)')MI
      READ(FNAME(T+6:T+7),'(I2)')DI
      READ(FNAME(T+8:T+9),'(I2)')HI
      READ(FNAME(T+11:T+20),'(I11)')FF    !chou  pula o character '+'
      Y=YI
      M=MI
      D=DI
      H=HI
      C=0
      DOWHILE((YI*1000000+MI*10000+DI*100+HI).LT.FF)
      C=C+1
      HI=HI+1
      IF(HI.GT.23) THEN
      HI=0
      DI=DI+1
      IF (MOD(YI,4).EQ.0) THEN
      MD(2)=29
      ELSE
      MD(2)=28
      ENDIF
      IF(DI.GT.MD(MI)) THEN
      DI=1
      MI=MI+1
      IF(MI.GT.12)THEN
      MI=1
      YI=YI+1
      ENDIF
      ENDIF
      ENDIF
      ENDDO
C*
      OPEN(12,FILE=FNAME(1:IN),STATUS='OLD',
     *        FORM='UNFORMATTED',
     *        ACCESS='STREAM')
      INQUIRE(FILE=FNAME(1:IN)//'.new',EXIST=EX)
      IF (EX) CALL SYSTEM('rm -f '//FNAME(1:IN)//'.new'//char(0))
      OPEN(13,FILE=FNAME(1:IN)//'.new',STATUS='NEW',
     *        FORM='UNFORMATTED',
     *        ACCESS='STREAM')
C*
   10 CALL HDRGRB(LGRB,LREC,LENF,FNAME(1:IN)) 
      CALL UPKGRB(LGRB-8,LREC,LENF,FNAME(1:IN),Y,M,D,H,C)
      GOTO 10
C*
      END
C*
      SUBROUTINE HDRGRB(LGRB,LREC,LENF,FNAME)
C*
      INTEGER I,LGRB,LREC
      CHARACTER*4 CH,GRIB,GRBL
      CHARACTER*(*) FNAME
      INTEGER ISTAT,STAT,LENF,NSTAT(64)
      INTEGER LG(4)
      CHARACTER*1 GRBI(4)
      EQUIVALENCE (GRBI,GRBL)
C*
      ISTAT=STAT(FNAME,NSTAT)
C*
      READ(12,END=99)GRIB
      IF (GRIB .NE. 'GRIB') 
     *    STOP' *** ERROR FROM HDRGRB: IT IS NOT A GRIB FILE ***'
      READ(12)GRBL
      DO I=1,4
      WRITE(CH,'(Z2)')GRBI(I)
      READ(CH,'(Z2)')LG(I)
      ENDDO
      LGRB=LG(1)*256*256+LG(2)*256+LG(3)
      LREC=LGRB/4
      LENF=NSTAT(8)/LGRB
C*
Chou      WRITE(*,'(A,A4,A,I8,A,I3,2(A,I8))')' GRIB = ',GRIB,
Chou     *        ' LGRB = ',LGRB,' EDITION = ',LG(4),
Chou     *        ' LREC = ',LREC,' LENF = ',LENF
C*
      WRITE(13)GRIB
      WRITE(13)GRBL
      RETURN
   99 CLOSE(12)
      CLOSE(13)
      CALL SYSTEM('mv -f '//FNAME//'.new '//FNAME//char(0))
      
      STOP'*** END GRIB FILE ***'
      END
C*
      SUBROUTINE UPKGRB(LGRB,LREC,LENF,FNAME,Y,M,D,H,C)
C*
      INTEGER I,J,K,LGRB,LREC,LENF,ITYLEV,IPRM,
     *        Y,M,D,H,C,JMAX,JH,JA,JB,IC
      CHARACTER*1 GRD(LGRB)*1
      CHARACTER*4 CH
      CHARACTER*(*) FNAME
      PARAMETER S=8
C*
      READ(12)GRD
      WRITE(CH,'(Z2)')GRD(17-S)
      READ(CH,'(Z2)')IPRM
C*    GRD(12) -> PDS(4) -> NUMERO DA VERSAO DA TABELA DE PARAMETROS
      IF (IPRM .GE. 128) THEN
      WRITE(GRD(12-S),'(A1)')254
      ELSE
      WRITE(GRD(12-S),'(A1)')1
      ENDIF
C*    GRD(13) -> PDS(5) -> IDENTIFICACAO DO CENTRO
      WRITE(GRD(13-S),'(A1)')46
C*    GRD(14) -> PDS(6) -> NUMERO (ID) DO PROCESSO GERADOR
      WRITE(GRD(14-S),'(A1)')100          ! modelo Eta de 100 a 199
C*    GRD(18) -> PDS(10) -> INDICADOR DO TIPO DE NIVEL OU CAMADA (1-SUPERFICIE)
      WRITE(CH,'(Z2)')GRD(18-S)
      READ(CH,'(Z2)')ITYLEV
      IF (ITYLEV .EQ. 99)WRITE(GRD(18-S),'(A1)')1
C*    GRD(21) -> PDS(13) -> ANO DO SECULO
C*    GRD(22) -> PDS(14) -> MES DO ANO
C*    GRD(23) -> PDS(15) -> DIA DO MES
C*    GRD(24) -> PDS(16) -> HORA DO DIA
      WRITE(GRD(21-S),'(A1)') MOD(Y,100)
       IF (MOD(Y,100).EQ.0) WRITE(GRD(21-S),'(A1)')100
C*    GRD(21) -> PDS(25) -> SECULO
       IC=INT(((Y-1)*.01)+1)
       WRITE(GRD(33-S),'(A1)')IC
      WRITE(GRD(22-S),'(A1)')M
      WRITE(GRD(23-S),'(A1)')D
      WRITE(GRD(24-S),'(A1)')H
C*    GRD(26) -> PDS(18) -> UNIDADE DO TEMPO DE PREVISAO (1-HORA)
      WRITE(GRD(26-S),'(A1)')1
C*    GRD(27) -> PDS(19) -> P1: PERIODO DE TEMPO (0-ANALISE)
      WRITE(GRD(27-S),'(A1)')C
C*    GRD(34) -> PDS(26) -> IDENTIFICADOR DO SUB-CENTRO
      WRITE(GRD(34-S),'(A1)')0
C*    GRD(62) -> GDS(26)
C*    GRD(63) -> GDS(27)
C*    GDS(26)-GDS(27) -> NUMERO DE LATITUDES EM UM HEM. NA GR. GAUSSIANA
Chou      CALL GETJMX(FNAME,JMAX)
Chou      JH=JMAX/2
Chou      JA=JH/256
Chou      JB=JH-JA*256
Chou      WRITE(GRD(62),'(A1)')JA
Chou      WRITE(GRD(63),'(A1)')JB
C*    GRD(64) -> GDS(28) -> INDICADOR DO MODO DE VARREDURA
      WRITE(GRD(64-S),'(A1)')64
C*    GRD(65) -> GDS(29) -> RESERVADO (=0)
      WRITE(GRD(65-S),'(A1)')0
      WRITE(13)GRD
C*
      RETURN
      END
C*
      SUBROUTINE GETJMX(FNAME,JMAX)
C*
      INTEGER JMAX
      CHARACTER*(*) FNAME
C*
      INTEGER IT,I
      INTEGER JR(3)
      CHARACTER JMC(3)*2
      DATA JMC/'20','40','80'/
      DATA JR /157 ,157 , 87 /
C*
      IT=INDEX(FNAME,'ll')+2
C*
      DO I=1,3
      IF (FNAME(IT:IT+1) .EQ. JMC(I)) THEN
      JMAX=JR(I)
      RETURN
      ENDIF
      ENDDO
C*
      WRITE(*,'(A)')
     *' DISASTER IN GETJMX: DID NOT FIND -> '//FNAME(IT:IT+1)
      STOP999
      END
