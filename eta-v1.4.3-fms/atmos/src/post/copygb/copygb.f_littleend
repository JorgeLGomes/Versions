      PROGRAM COPYGB
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: COPYGB       COMPARE GRIB FILES
C   PRGMMR: IREDELL          ORG: NP23        DATE: 97-03-05  
C
C ABSTRACT: The command copygb copies all or part of one GRIB file
C   to another GRIB file, interpolating if necessary.  Unless
C   otherwise directed (-x option), the GRIB index file is also used
C   to speed the reading. The fields are interpolated to an output grid
C   if specified (-g option). The interpolation type defaults to
C   bilinear but may be specified directly (-i option).  The copying
C   may be limited to specific fields (-k option). It may also be
C   limited to a specified subgrid of the output grid or to a subrange
C   of the input fields (-B and -b, -A, and -K options). Fields can be
C   identified as scalars or vectors (-v option), which are interpolated
C   differently.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C   97-03-05  IREDELL  CORRECTED THE COPYING OF THE V-WIND FIELD
C                      WHEN NO INTERPOLATION IS DONE
C   98-07-03  EBISUZAKI Linux port
C
C COMMAND LINE OPTIONS:
C   -A "<> mapthreshold"
C      Inequality and threshold used in determining
C      where on the map statistics will be computed.
C      Statistics are computed only where the given 
C      map field is on the correct side of the threshold.
C      The mapthresh defaults to '>-1.e30'; in this case,
C      only the map field's bitmap will limit the domain.
C
C   -b mapindex   
C      Optional index file used to get the map field.
C
C   -B mapgrib    
C      GRIB file used to get the map field.  The map field
C      is read from the GRIB file and compared to the
C      map threshold to determine for which region on the map
C      the statistics will be computed.  mapgrib can be the
C      name of an actual GRIB file (in which case the index
C      file may be specified with the -b option) or it can
C      be '-1'.  If mapgrib is '-1', then GRIB file 1
C      (first positional argument) is used.
C      The -K option specifies which field to read from
C      the mapgrib GRIB file.  If mapgrib is an actual file,
C      then the first field is taken if -K is not specified.
C      On the other hand, if mapgrib is '-1', then if the
C      if -K is not specified, the current field is taken
C      as the map field.  A special exception is if -K '-1'
C      is specified, in which case the current field is
C      taken as the map field and it is applied before any
C      interpolation; otherwise the map field is always
C      applied after interpolation.
C
C   -g "grid [kgds]"
C      Verification grid identification.  If grid=-1
C      (the default), then the grid is taken from the first
C      GRIB field in GRIB file 1.  If grid=-4,
C      then the grid is taken from the first GRIB field
C      in the mapgrib file.
C      If 0<grid<255, then grid designates an NCEP grid.
C      If grid=255, then the grid must be specified by the
C      full set of kgds parameters determining a GRIB GDS
C      (grid description section) in the W3FI63 format.
C
C   -i "ip [ipopts]"
C      Interpolation options.  The default is bilinear
C      interpolation (ip=0).  Other interpolation options
C      are bicubic (ip=1), neighbor (ip=2), budget (ip=3),
C      and spectral (ip=4).  See the documentation for iplib
C      for further details about interpolation.
C
C   -k "kpds"
C      Full set of kpds parameters determing a GRIB PDS
C      (product definition section) in the W3FI63 format
C      determining the field(s) from GRIB file 1
C      for which statistics are computed.  Note that
C      kpds(5) is the parameter indicator (PDS octet 9).
C      A wildcard is specified by -1 (the defaults).
C      If the -k is not specified, then copygb will attempt
C      to copy every field in GRIB file 1.
C
C   -K "mapkpds"
C      Full set of kpds parameters determing a GRIB PDS
C      (product definition section) in the W3FI63 format
C      determining the map field to be used to determine
C      where on the map statistics will be computed.  
C      A wildcard is specified by -1 (the defaults).
C
C   -M "mask"
C      Mask used to fill out bitmapped areas of the map.
C      If specified, there will be no bitmap in the output.
C      The mask must be in the format '#value' where value
C      is the real number used to fill out the field.
C
C   -s "ids,ibs"
C      Decimal scaling and optionally binary scaling.
C      Each scaling determines the number of places
C      to the right of the decimal point that will be saved.
C      Decimal and binary scalings are combined to give the
C      final total scaling (e.g. a scaling of "2,-3" means
C      that the values are saved to the nearest 0.08).
C      If binary scaling is greater than 100, then it
C      actually represents 100 plus the number of bits
C      to pack each value in, overriding other scalings.
C      Binary scaling defaults to zero if not specified.
C      Decimal scaling defaults to that of the input field
C      if not specified or if it is less than -100.
C
C   -v "uparms"
C      Parameter indicator(s) for the u-component of vectors.
C      The parameter indicator for the v-component is assumed
C      to be one more than that of the u-component.
C      If the -v option is not specified, then the wind
C      components (parameters 33 and 34) are the only fields
C      assumed to be vector components in the GRIB file.
C
C   -x
C      Turns off the use of an index file.  The index records
C      are then extracted from the GRIB file, which
C      will increase the time taken by copygb.
C
C   -X
C      Turns on verbose printout.  This option is incompatible
C      with GRIB output to standard output (designated as '-').
C
C INPUT FILES:
C   UNIT   11    GRIB FILE 1
C   UNIT   14    MAP GRIB FILE
C   UNIT   31    GRIB INDEX FILE 1
C   UNIT   34    MAP GRIB INDEX FILE
C
C OUTPUT FILES:
C   UNIT   51    GRIB FILE 2
C
C SUBPROGRAMS CALLED:
C   FILENV
C   IARGC
C   ISTRLEN
C   GETARG
C   ERRMSG
C   EUSAGE
C   EXIT
C   FPARSEI
C   FPARSER
C   BAOPENR
C   BAOPEN
C   CPGB
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CHARACTER*256 CARG,CG1,CX1,CGB,CXB,CG2
      INTEGER KARG(100)
      INTEGER KGDSI(200),IPOPT(20),JPDS1(200),JPDSB(200),IUV(100)
      REAL RARG(100)
      INTEGER ISS(2)
      CHARACTER*400 GDS
      DATA IGI/-1/,KGDSI/19*0,255,180*0/
      DATA IP/0/,IPOPT/20*-1/
      DATA JPDS1/200*-1/,JPDSB/200*-1/,IUV/33,99*0/,NUV/1/
      DATA LXX/0/,LX/1/,KZ1/-1/,KZ2/-2/
      DATA ISS/-999,0/
      DATA JB/0/,JBK/0/,LAB/1/,AB/-1.E30/,LAM/0/,AM/0./
      DATA CGB/' '/,CXB/' '/
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PARSE COMMAND LINE OPTIONS
Cmp      CALL FILENV
      NARG=IARGC()
      IARG=1
      LSTOPT=0
      DOWHILE(IARG.LE.NARG.AND.LSTOPT.EQ.0)
        call GETARG(IARG,CARG)
        LARG=istrlen(CARG)
        IARG=IARG+1
        IF(CARG(1:1).NE.'-') THEN
          LSTOPT=1
          IARG=IARG-1
        ELSEIF(LARG.EQ.1) THEN
          CALL ERRMSG('copygb: invalid option -')
          CALL EUSAGE
          CALL EXIT(1)
        ELSE
          L=2
          DOWHILE(L.LE.LARG)
            IF(CARG(L:L).EQ.'-') THEN
              LSTOPT=1
            ELSEIF(CARG(L:L).EQ.'A') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              IF(CARG(L+1:L+1).EQ.'>') THEN
                LAB=1
                L=L+1
              ELSEIF(CARG(L+1:L+1).EQ.'<') THEN
                LAB=-1
                L=L+1
              ELSE
                CALL ERRMSG('copygb: invalid threshold '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              CALL FPARSER(CARG(L+1:LARG),1,AB)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'B') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              LCGB=LARG-L
              CGB=CARG(L+1:LARG)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'b') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              LCXB=LARG-L
              CXB=CARG(L+1:LARG)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'g') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              KARG(1)=IGI
c             KARG(2:100)=KGDSI(1:99)
              do ii = 1, 99
                 KARG(ii+1)=KGDSI(ii)
              enddo
              CALL FPARSEI(CARG(L+1:LARG),100,KARG)
              IGI=KARG(1)
              IF(IGI.GT.0.AND.IGI.LT.255) THEN
                CALL MAKGDS(IGI,KGDSI,GDS,LGDS,IRET)
                IF(IRET.NE.0) IGI=-1
              ELSEIF(IGI.EQ.255) THEN
c               KGDSI(1:99)=KARG(2:100)
                do ii = 1, 99
                   KGDSI(ii)=KARG(ii+1)
                enddo
              ENDIF
              IF(IGI.LT.-4.OR.IGI.EQ.0.OR.IGI.GT.255) THEN
                CALL ERRMSG('copygb: invalid output grid '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              MI=LENGDS(KGDSI)
C	write(6,*) 'LENGDS= ', MI
              IF(MI.LE.0) THEN
                CALL ERRMSG('copygb: unsupported output grid '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'i') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              KARG(1)=IP
c             KARG(2:21)=IPOPT
              do ii = 1, 20
                 KARG(ii+1)=IPOPT(ii)
              enddo
              CALL FPARSEI(CARG(L+1:LARG),21,KARG)
              IP=KARG(1)
c             IPOPT=KARG(2:21)
              do ii = 1, 20
                  IPOPT(ii)=KARG(ii+1)
              enddo
              L=LARG
            ELSEIF(CARG(L:L).EQ.'K') THEN
              IF(L.EQ.LARG) THEN
                L=0
                call GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              JBK=1
              CALL FPARSEI(CARG(L+1:LARG),100,JPDSB)
              IF(JPDSB(5).EQ.0) THEN
                CALL ERRMSG('copygb: invalid PDS parms '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'k') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              CALL FPARSEI(CARG(L+1:LARG),100,JPDS1)
              IF(JPDS1(5).EQ.0) THEN
                CALL ERRMSG('copygb: invalid PDS parms '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'M') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              IF(CARG(L+1:L+1).EQ.'#') THEN
                L=L+1
                CALL FPARSER(CARG(L+1:LARG),1,AM)
                LAM=1
              ELSE
                CALL ERRMSG('copygb: invalid merge parameter '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'s') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              CALL FPARSEI(CARG(L+1:LARG),2,ISS)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'v') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=istrlen(CARG)
                IARG=IARG+1
              ENDIF
              CALL FPARSEI(CARG(L+1:LARG),100,IUV)
              NUV=1
              DO JUV=2,100
                IF(IUV(JUV).NE.0) NUV=JUV
              ENDDO
              L=LARG
            ELSEIF(CARG(L:L).EQ.'x') THEN
              LX=0
            ELSEIF(CARG(L:L).EQ.'X') THEN
              LXX=1
            ELSE
              CALL ERRMSG('copygb: invalid option '//CARG(L:L))
              CALL EUSAGE
              CALL EXIT(1)
            ENDIF
            L=L+1
          ENDDO
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PARSE COMMAND LINE POSITIONAL ARGUMENTS
      NXARG=LX+2
      IF(NARG-IARG+1.NE.NXARG) THEN
        CALL ERRMSG('copygb: incorrect number of arguments')
        CALL EUSAGE
        CALL EXIT(NXARG)
      ENDIF
      CALL GETARG(IARG,CG1)
      LCG1=istrlen(CG1)
      IARG=IARG+1
      LG1=11
      CALL BAOPENR(LG1,CG1(1:LCG1),IRETBA)
      IF(IRETBA.NE.0) THEN
        CALL ERRMSG('copygb:  error accessing file '//CG1(1:LCG1))
        CALL EXIT(8)
      ENDIF
      IF(LX.GT.0) THEN
        CALL GETARG(IARG,CX1)
        LCX1=istrlen(CX1)
        IARG=IARG+1
        LX1=31
        CALL BAOPENR(LX1,CX1(1:LCX1),IRETBA)
        IF(IRETBA.NE.0) THEN
          CALL ERRMSG('copygb:  error accessing file '//CX1(1:LCX1))
          CALL EXIT(8)
        ENDIF
      ELSE
        LX1=0
      ENDIF
      CALL GETARG(IARG,CG2)
      LCG2=istrlen(CG2)
      IARG=IARG+1
      IF(CG2(1:LCG2).EQ.'-') THEN
        IF(LXX.GT.0) THEN
          CALL ERRMSG('copygb:  piping incompatible with the X option')
          CALL EXIT(1)
        ENDIF
        LG2=6
      ELSE
        LG2=51
        CALL BAOPEN(LG2,CG2(1:LCG2),IRETBA)
        IF(IRETBA.NE.0) THEN
          CALL ERRMSG('copygb:  error accessing file '//CG2(1:LCG2))
          CALL EXIT(8)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  OPEN MAP FILE
      IF(CGB.NE.' ') THEN
        IF(CGB(1:2).EQ.'-1') THEN
          IF(JPDSB(5).EQ.-1) THEN
            JB=1
          ELSE
            JB=4
            LGB=LG1
            LXB=LX1
          ENDIF
        ELSE
          JB=4
          LGB=14
          CALL BAOPENR(LGB,CGB(1:LCGB),IRETBA)
          IF(IRETBA.NE.0) THEN
            CALL ERRMSG('copygb:  error accessing file '//CGB(1:LCGB))
            CALL EXIT(8)
          ENDIF
          IF(CXB(1:1).NE.' ') THEN
            LXB=34
            CALL BAOPENR(LXB,CXB(1:LCXB),IRETBA)
            IF(IRETBA.NE.0) THEN
              CALL ERRMSG('copygb:  error accessing file '//CXB(1:LCXB))
              CALL EXIT(8)
            ENDIF
          ELSE
            LXB=0
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GO
      IF(LXX.GT.0) THEN

Cmp      CALL W3LOG('$S97064.78','COPYGB  ')
Cmp      CALL W3TAGB('COPYGB  ',0097,0064,0078,'NP23   ')                   


      ENDIF
      CALL CPGB(LG1,LX1,LGB,LXB,LG2,
     &          IGI,KGDSI,IP,IPOPT,JPDS1,ISS,NUV,IUV,
     &          JPDSB,JB,JBK,LAB,AB,LAM,AM,LXX)
      IF(LXX.GT.0) THEN

Cmp        CALL W3LOG('$E')
Cmp      CALL W3TAGE('COPYGB  ') 

      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE EUSAGE
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    EUSAGE      PRINT PROPER USAGE TO STDERR
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT PROPER USAGE TO STDERR.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL EUSAGE
C
C SUBPROGRAMS CALLED:
C   ERRMSG
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL ERRMSG('Usage: copygb'//
     & ' [-g "grid [kgds]"] [-i "ip [ipopts]"]'//
     & ' [-k "kpds"] [-v "uparms"]')
      CALL ERRMSG('             '//
     & ' [-B mapgrib [-b mapindex] [-A "<> mapthreshold"]'//
     & ' [-K "mapkpds"]]')
      CALL ERRMSG('             '//
     & ' [-M "mask"] [-X] [-s "ids,ibs"]')
      CALL ERRMSG('       then either:')
      CALL ERRMSG('             '//
     & ' grib1 index1 grib2')
      CALL ERRMSG('            or:')
      CALL ERRMSG('             '//
     & ' -x grib1 grib2')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE CPGB(LG1,LX1,LGB,LXB,LG2,
     &                IGI,KGDSI,IP,IPOPT,JPDS1,ISS,NUV,IUV,
     &                JPDSB,JB,JBK,LAB,AB,LAM,AM,LXX)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    CPGB        COPY GRIB FILES
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: COPY GRIB FILES.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL CPGB(LG1,LX1,LGB,LXB,LG2,
C    &                IGI,KGDSI,IP,IPOPT,JPDS1,ISS,NUV,IUV,
C    &                JPDSB,JB,JBK,LAB,AB,LAM,AM,LXX)
C   INPUT ARGUMENTS:
C     LG1          INTEGER UNIT NUMBER FOR GRIB FILE 1
C     LX1          INTEGER UNIT NUMBER FOR GRIB INDEX FILE 1
C     LGB          INTEGER UNIT NUMBER FOR GRIB FILE MAP
C     LXB          INTEGER UNIT NUMBER FOR GRIB INDEX FILE MAP
C     LG2          INTEGER UNIT NUMBER FOR GRIB FILE 2
C     IGI          INTEGER OUTPUT GRID IDENTIFICATION
C     KGDSI        INTEGER (200) OUTPUT GRID PARAMETERS
C     IP           INTEGER INTERPOLATION TYPE
C     IPOPT        INTEGER (20) INTERPOLATION OPTIONS
C     JPDS1        INTEGER (100) KPDS SEARCH OPTIONS
C     ISS          INTEGER (2) DECIMAL AND BINARY SCALINGS
C     NUV          INTEGER NUMBER OF VECTOR PARAMETER IDS
C     IUV          INTEGER (100) VECTOR PARAMETER IDS
C     JPDSB        INTEGER (100) KPDS SEARCH OPTIONS (MAP)
C     JB           INTEGER FLAG FOR MAP OPTION
C     JBK          INTEGER FLAG FOR MAP OPTION
C     LAB          INTEGER FLAG FOR MAP THRESHOLD INEQUALITY
C     AB           REAL MAP THRESHOLD
C     LAM          INTEGER FLAG FOR MASK VALUE
C     AM           REAL MASK VALUE
C     LXX          INTEGER FLAG FOR VERBOSE OUTPUT
C
C SUBPROGRAMS CALLED:
C   GETGBMH
C   CPGB1  
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      PARAMETER(MBUF=256*1024)
      CHARACTER*1 CBUF1(MBUF),CBUFB(MBUF)
      INTEGER JPDS1(100),JPDSB(100),IUV(100)
      INTEGER ISS(2)
      INTEGER KGDSI(200)
      INTEGER IPOPT(20)
      INTEGER JPDS(200),JGDS(200),KPDS1(200),KGDS1(200)
      INTEGER KPDSB(200),KGDSB(200)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  READ GRIB HEADERS
      IF(JB.EQ.4) THEN
c       JPDS=-1
c       JGDS=-1
c       KPDSB=0
c       KGDSB=0
        do ii = 1, 200
            JPDS(ii)=-1
            JGDS(ii)=-1
            KPDSB(ii)=0
            KGDSB(ii)=0
        enddo
        KRB=-1
        CALL GETGBMH(LGB,LXB,KRB,JPDSB,JGDS,
     &               MBUF,CBUFB,NLENB,NNUMB,MNUMB,
     &               KB,MB,KRB,KPDSB,KGDSB,IRET)
        IF(IRET.NE.0) THEN
          CALL ERRMSG('copygb: error retrieving bitmap')
          CALL EXIT(IRET)
        ENDIF
      ENDIF
      KR1=-1
c     KPDS1=0
c     KGDS1=0
      do ii = 1, 200
          KPDS1(ii)=0
          KGDS1(ii)=0
      enddo
      CALL GETGBMH(LG1,LX1,KR1,JPDS1,JGDS,
     &             MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &             K1,M1,KR1,KPDS1,KGDS1,IRET)
      IF(IRET.NE.0) THEN
        CALL ERRMSG('copygb: error retrieving requested fields')
        CALL EXIT(IRET)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  LOOP UNTIL DONE
      NO=0
      DOWHILE(IRET.EQ.0)
        IF(IGI.EQ.-1) THEN
          IGI=KPDS1(3)
c         KGDSI=KGDS1
          do ii = 1, 200
             KGDSI(ii)=KGDS1(ii)
          enddo
          MI=M1
        ELSEIF(IGI.EQ.-4.AND.JB.EQ.4) THEN
          IGI=KPDSB(3)
c         KGDSI=KGDSB
          do ii = 1, 200
             KGDSI(ii)=KGDSB(ii)
          enddo
          MI=MB
        ELSE
          MI=LENGDS(KGDSI)
        ENDIF
        IF(IGI.GT.0.AND.IGI.LE.255) THEN
          MF=MAX(M1,MB)
C	write(6,*) 'calling cpgb1 with KS1= ', KR1-1
          CALL CPGB1(LG1,LX1,M1,CBUF1,NLEN1,NNUM1,MNUM1,
     &               MBUF,MF,MI,
     &               IGI,KGDSI,IP,IPOPT,JPDS1,ISS,NUV,IUV,
     &               JPDSB,JB,JBK,LAB,AB,LAM,AM,
     &               LGB,LXB,MB,CBUFB,NLENB,NNUMB,MNUMB,
     &               LG2,LXX,KR1-1,NO,IRET1)
        ENDIF
c       KPDS1=0
c       KGDS1=0
        do ii = 1, 200
            KPDS1(ii)=0
            KGDS1(ii)=0
        enddo
        CALL GETGBMH(LG1,LX1,KR1,JPDS1,JGDS,
     &               MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &               K1,M1,KR1,KPDS1,KGDS1,IRET)
        IF(LXX.GT.0) THEN
          IF(IRET.NE.0.AND.IRET.NE.99) THEN
            PRINT *,'copygb GRIB unpacking error code ',IRET
          ENDIF
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LXX.GT.0) THEN
        PRINT *,'copygb wrote ',NO,' total records'
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE CPGB1(LG1,LX1,M1,CBUF1,NLEN1,NNUM1,MNUM1,
     &                 MBUF,MF,MI,
     &                 IGI,KGDSI,IP,IPOPT,JPDS1,ISS,NUV,IUV,
     &                 JPDSB,JB,JBK,LAB,AB,LAM,AM,
     &                 LGB,LXB,MB,CBUFB,NLENB,NNUMB,MNUMB,
     &                 LG2,LXX,KS1,NO,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    CPGB1       COPY ONE GRIB FIELD
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: COPY ONE GRIB FIELD.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL CPGB1(LG1,LX1,M1,CBUF1,NLEN1,NNUM1,MNUM1,
C    &                 MBUF,MF,MI,
C    &                 IGI,KGDSI,IP,IPOPT,JPDS1,ISS,NUV,IUV,
C    &                 JPDSB,JB,JBK,LAB,AB,LAM,AM,
C    &                 LGB,LXB,MB,CBUFB,NLENB,NNUMB,MNUMB,
C    &                 LG2,LXX,KS1,NO,IRET)
C   INPUT ARGUMENTS:
C     LG1          INTEGER UNIT NUMBER FOR GRIB FILE 1
C     LX1          INTEGER UNIT NUMBER FOR GRIB INDEX FILE 1
C     M1           INTEGER DIMENSION OF GRIB FIELD 1
C     CBUF1        CHARACTER (MBUF) INDEX BUFFER 1
C     NLEN1        INTEGER RECORD LENGTH OF INDEX BUFFER 1
C     NNUM1        INTEGER NUMBER OF RECORDS IN INDEX BUFFER 1
C     NLEN1        INTEGER LENGTH OF EACH INDEX RECORD 1
C     NNUM1        INTEGER NUMBER OF INDEX RECORDS 1
C     MNUM1        INTEGER NUMBER OF INDEX RECORDS 1 SKIPPED
C     MBUF         INTEGER DIMENSION OF INDEX BUFFERS
C     MF           INTEGER DIMENSION OF FIELD
C     MI           INTEGER DIMENSION OF OUTPUT GRID
C     IGI          INTEGER OUTPUT GRID IDENTIFICATION
C     KGDSI        INTEGER (200) OUTPUT GRID PARAMETERS
C     IP           INTEGER INTERPOLATION TYPE
C     IPOPT        INTEGER (20) INTERPOLATION OPTIONS
C     JPDS1        INTEGER (100) KPDS SEARCH OPTIONS
C     ISS          INTEGER (2) DECIMAL AND BINARY SCALINGS
C     NUV          INTEGER NUMBER OF VECTOR PARAMETER IDS
C     IUV          INTEGER (100) VECTOR PARAMETER IDS
C     JPDSB        INTEGER (100) KPDS SEARCH OPTIONS (MAP)
C     JB           INTEGER FLAG FOR MAP OPTION
C     JBK          INTEGER FLAG FOR MAP OPTION
C     LAB          INTEGER FLAG FOR MAP THRESHOLD INEQUALITY
C     AB           REAL MAP THRESHOLD
C     LAM          INTEGER FLAG FOR MASK VALUE
C     AM           REAL MASK VALUE
C     LGB          INTEGER UNIT NUMBER FOR GRIB FILE MAP
C     LXB          INTEGER UNIT NUMBER FOR GRIB INDEX FILE MAP
C     MB           INTEGER DIMENSION OF GRIB FIELD MAP
C     CBUFB        CHARACTER (MBUF) INDEX BUFFER MAP
C     NLENB        INTEGER RECORD LENGTH OF INDEX BUFFER MAP
C     NNUMB        INTEGER NUMBER OF RECORDS IN INDEX BUFFER MAP
C     NLENB        INTEGER LENGTH OF EACH INDEX RECORD MAP
C     NNUMB        INTEGER NUMBER OF INDEX RECORDS MAP
C     MNUMB        INTEGER NUMBER OF INDEX RECORDS MAP SKIPPED
C     LG2          INTEGER UNIT NUMBER FOR GRIB FILE 2
C     LXX          INTEGER FLAG FOR VERBOSE OUTPUT
C     KS1          INTEGER INPUT RECORD COUNTER
C     NO           INTEGER OUTPUT RECORD COUNTER
C   OUTPUT ARGUMENTS:
C     NO           INTEGER OUTPUT RECORD COUNTER
C     IRET         INTEGER RETURN CODE
C
C SUBPROGRAMS CALLED:
C   GETGBM
C   INTGRIB
C   PUTGB
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CHARACTER*1 CBUF1(MBUF),CBUFB(MBUF)
      INTEGER JPDS1(100),JPDSB(100),IUV(100)
      INTEGER ISS(2)
      INTEGER KGDSI(200)
      INTEGER IPOPT(20)
      INTEGER JPDS(200),JGDS(200),KPDS1(200),KGDS1(200)
      INTEGER KPDSB(200),KGDSB(200)
      LOGICAL LR(MF),L1I(MI),LBI(MI)
      REAL FR(MF),F1I(MI),FBI(MI)
      REAL GR(MF),G1I(MI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET FIELD FROM FILE 1
c     JGDS=-1
c     KPDS1=0
c     KGDS1=0
      do ii = 1, 200
         JGDS(ii)=-1
         KPDS1(ii)=0
         KGDS1(ii)=0
      enddo
Cmp
C	write(6,*) ' '
C	write(6,*) 'first getgbm call...KS1= ', ks1
Cmp
      CALL GETGBM(LG1,LX1,M1,KS1,JPDS1,JGDS,
     &            MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &            K1,KR1,KPDS1,KGDS1,LR,FR,IRET)
      IV=0
      KRV=0
      IF(IRET.EQ.0) THEN
        JUV=1
        DOWHILE(JUV.LE.NUV.AND.KPDS1(5).NE.IUV(JUV).AND.
     &          KPDS1(5).NE.IUV(JUV)+1)
          JUV=JUV+1
        ENDDO
Cmp
C	write(6,*) 'KPDS1(5),IUV(JUV) ', KPDS1(5),IUV(JUV)
Cmp
        IF(JUV.LE.NUV.AND.KPDS1(5).EQ.IUV(JUV)) THEN
C	write(6,*) 'in vector part'
Cmp
Cmp	look for JPDS that is same as KPDS coming out from
Cmp	above, but with PDS(5) one greater (v component)
          IV=1
c         JPDS=-1
          do ii = 1, 200
              JPDS(ii)=-1
          enddo
c         JPDS(1:21)=KPDS1(1:21)
          do ii = 1, 21
              JPDS(ii)=KPDS1(ii)
          enddo
          JPDS(5)=KPDS1(5)+1
c         JGDS=KGDS1
          do ii = 1, 200
              JGDS(ii)=KGDS1(ii)
          enddo
C	write(6,*) 'calling 2nd getgbm...looking for following'
C	write(6,*) 'JPDS(1:12)' ,(JPDS(I),I=1,12)

Cmp
	NLEN1=-1	
	NNUM1=-1	
	MNUM1=-1	

Cmp
          CALL GETGBM(LG1,LX1,M1,KRV,JPDS,JGDS,
     &                MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &                K1,KRV,KPDS1,KGDS1,LR,GR,IRET)

C	write(6,*) 'leaving 2nd getgbm...IRET= ', IRET
C	write(6,*) 'KPDS1(5),JPDS(5)-1 = ', KPDS1(5),JPDS(5)-1
C
          KPDS1(5)=JPDS(5)-1
        ELSEIF(JUV.LE.NUV.AND.KPDS1(5).EQ.IUV(JUV)+1) THEN
          IRET=-1
        ENDIF
      ENDIF
      IF(LXX.GT.0) THEN
        IF(IRET.EQ.-1) THEN
          PRINT *,'copygb skipping 2nd vector component field'
	ENDIF
        IF(KRV.EQ.0) THEN
          PRINT *,'copygb read scalar field from record ',KR1
          PRINT *,'       ...KPDS(1:16)=',(KPDS1(I),I=1,16)
        ELSE
          PRINT *,'copygb read vector field from records ',KR1,KRV
          PRINT *,'       ...KPDS(1:16)=',(KPDS1(I),I=1,16)
          PRINT *,'       ...KPDS(1:16)=',(KPDS1(I),I=1,4),
     &            KPDS1(5)+1,(KPDS1(I),I=6,16)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INVOKE MAP MASK BEFORE INTERPOLATION
      IF(IRET.EQ.0.AND.JBK.EQ.1.AND.JB.EQ.1) THEN
        DO I=1,K1
          IF(LR(I)) THEN
            IF((LAB.EQ.1.AND.FR(I).LE.AB).OR.
     &         (LAB.EQ.-1.AND.FR(I).GE.AB)) THEN
              IB1=1
              LR(I)=.FALSE.
            ENDIF
          ENDIF
        ENDDO
        IF(LXX.GT.0) THEN
          PRINT *,'       applied pre-interpolation map mask'
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE FIELD 1
      IF(IRET.EQ.0) THEN
        IB1=MOD(KPDS1(4)/64,2)
        CALL INTGRIB(IV,IP,IPOPT,KGDS1,K1,IB1,LR,FR,GR,KGDSI,MI,
     &               IB1I,L1I,F1I,G1I,IRET)
        IF(LXX.GT.0) THEN
          IF(IRET.EQ.0) THEN
            PRINT *,'       interpolated to grid ',IGI
          ELSEIF(IRET.GT.0) THEN
            PRINT *,'       interpolation error code ',IRET
          ENDIF
        ENDIF
        IF(IRET.EQ.-1) IRET=0
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET MAP FIELD
      IF(IRET.EQ.0.AND.JB.EQ.4) THEN
        KRB=0
c       JGDS=-1
        do ii = 1, 200
           JGDS(ii)=-1
        enddo
        CALL GETGBM(LGB,LXB,MB,KRB,JPDSB,JGDS,
     &              MBUF,CBUFB,NLENB,NNUMB,MNUMB,
     &              KB,KRB,KPDSB,KGDSB,LR,FR,IRET)
        IF(LXX.GT.0) THEN
          IF(IRET.EQ.0) THEN
            PRINT *,'       map field retrieved'
            PRINT *,'       ...KPDS(1:16)=',(KPDSB(I),I=1,16)
          ELSEIF(IRET.EQ.99) THEN
            PRINT *,'       map field not found'
          ELSE
            PRINT *,'       map field retrieval error code ',IRET
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE MAP FIELD
        IF(IRET.EQ.0) THEN
          IBB=MOD(KPDSB(4)/64,2)
          CALL INTGRIB(0,IP,IPOPT,KGDSB,KB,IBB,LR,FR,GR,KGDSI,MI,
     &                 IBBI,LBI,FBI,GBI,IRET)
          IF(LXX.GT.0) THEN
            IF(IRET.EQ.0) THEN
              PRINT *,'       interpolated to grid ',IGI
            ELSEIF(IRET.GT.0) THEN
              PRINT *,'       interpolation error code ',IRET
            ENDIF
          ENDIF
          IF(IRET.EQ.-1) IRET=0
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INVOKE MAP MASK
      IF(IRET.EQ.0) THEN
        IF(JBK.EQ.0.AND.JB.EQ.1) THEN
          DO I=1,MI
            IF(L1I(I)) THEN
              IF((LAB.EQ.1.AND.F1I(I).LE.AB).OR.
     &           (LAB.EQ.-1.AND.F1I(I).GE.AB)) THEN
                IB1I=1
                L1I(I)=.FALSE.
              ENDIF
            ENDIF
          ENDDO
          IF(LXX.GT.0) THEN
            PRINT *,'       applied post-interpolation map mask'
          ENDIF
        ELSEIF(JB.EQ.4) THEN
          DO I=1,MI
            IF(LBI(I)) THEN
              IF((LAB.EQ.1.AND.FBI(I).LE.AB).OR.
     &           (LAB.EQ.-1.AND.FBI(I).GE.AB)) THEN
                IB1I=1
                L1I(I)=.FALSE.
              ENDIF
            ELSE
              IB1I=1
              L1I(I)=.FALSE.
            ENDIF
          ENDDO
          IF(LXX.GT.0) THEN
            PRINT *,'       applied fixed map mask'
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  MASK VALUES
        IF(LAM.EQ.1) THEN
          IB1I=0
          DO I=1,MI
            IF(.NOT.L1I(I)) THEN
              L1I(I)=.TRUE.
              F1I(I)=AM
              IF(KRV.GT.0) G1I(I)=AM
            ENDIF
          ENDDO
          IF(LXX.GT.0) THEN
            PRINT *,'       substituted mask fill value'
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  WRITE OUTPUT FIELD
      IF(IRET.EQ.0) THEN
        KPDS1(3)=IGI
        KPDS1(4)=128+64*IB1I
	IF(ISS(1).GT.-100) KPDS1(22)=ISS(1)
        IF(ISS(2).LE.100) THEN
          NBIT=0
          IBS=ISS(2)
        ELSE
          NBIT=ISS(2)-100
          IBS=0
        ENDIF
Cmp        CALL PUTGB(LG2,MI,KPDS1,KGDSI,L1I,F1I,IRET)
	CALL PUTGBN(LG2,MI,KPDS1,KGDSI,IBS,NBIT,L1I,F1I,IRET)
        IF(IRET.EQ.0) NO=NO+1
        IF(IRET.EQ.0.AND.KRV.GT.0) THEN
          KPDS1(5)=KPDS1(5)+1
          CALL PUTGB(LG2,MI,KPDS1,KGDSI,L1I,G1I,IRET)
          IF(IRET.EQ.0) NO=NO+1
          KPDS1(5)=KPDS1(5)-1
        ENDIF
        IF(LXX.GT.0) THEN
          IF(IRET.EQ.0) THEN
            IF(KRV.EQ.0) THEN
              PRINT *,'       wrote scalar field to record ',NO
              PRINT *,'       ...KPDS(1:16)=',(KPDS1(I),I=1,16)
            ELSE
              PRINT *,'       wrote vector field to records ',NO-1,NO
              PRINT *,'       ...KPDS(1:16)=',(KPDS1(I),I=1,16)
              PRINT *,'       ...KPDS(1:16)=',(KPDS1(I),I=1,4),
     &                KPDS1(5)+1,(KPDS1(I),I=6,16)
            ENDIF
          ELSE
            PRINT *,'       GRIB packing error code ',IRET
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE INTGRIB(IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,
     &                   IB2,L2,F2,G2,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    INTGRIB     INTERPOLATE FIELD
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: INTERPOLATE FIELD.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL INTGRIB(IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,
C    &                   IB2,L2,F2,G2,IRET)
C   INPUT ARGUMENTS:
C     IV           INTEGER VECTOR FLAG
C     IP           INTEGER INTERPOLATION TYPE
C     IPOPT        INTEGER (20) INTERPOLATION OPTIONS
C     KGDS1        INTEGER (200) INPUT GRID PARAMETERS
C     K1           INTEGER INPUT DIMENSION
C     IB1          INTEGER INPUT BITMAP FLAG
C     L1           LOGICAL (K1) INPUT BITMAP IF IB1=1
C     F1           REAL (K1) INPUT FIELD
C     G1           REAL (K1) INPUT Y-COMPONENT IF IV=1
C     KGDS2        INTEGER (200) OUTPUT GRID PARAMETERS
C     K2           INTEGER OUTPUT DIMENSION
C     IB2          INTEGER OUTPUT BITMAP FLAG
C     L2           LOGICAL (K2) OUTPUT BITMAP
C     F2           REAL (K2) OUTPUT FIELD
C     G2           REAL (K2) OUTPUT Y-COMPONENT IF IV=1
C
C SUBPROGRAMS CALLED:
C   LENGDSF
C   INTGRIB1
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER IPOPT(20)
      INTEGER KGDS1(200),KGDS2(200)
      LOGICAL L1(K1),L2(K2)
      REAL F1(K1),F2(K2)
      REAL G1(K1),G2(K2)
      INTEGER KGDS1F(200),KGDS2F(200)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  DETERMINE WHETHER INTERPOLATION IS NECESSARY
      IF(IP.EQ.4) THEN
        INT=1
      ELSE
        INT=0
        DO I=1,200
          INT=MAX(INT,ABS(KGDS1(I)-KGDS2(I)))
        ENDDO
        INT=MIN(INT,1)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COPY FIELD
      IF(INT.EQ.0) THEN
        IB2=IB1
        DO I=1,K1
          L2(I)=L1(I)
          F2(I)=F1(I)
          IF(IV.NE.0) G2(I)=G1(I)
        ENDDO
        IRET=-1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE REGULARIZED GRIDS AND INTERPOLATE FIELD
      ELSE
        K1F=LENGDSF(KGDS1,KGDS1F)
        IF(K1F.EQ.K1) K1F=1
        K2F=LENGDSF(KGDS2,KGDS2F)
        IF(K2F.EQ.K2) K2F=1
        MRL=MAX(K2,K2F)
        IF(IV.EQ.0) THEN
          MRO=1
        ELSE
          MRO=MRL
        ENDIF
        IF(K1F.GT.0.AND.K2F.GT.0) THEN
          CALL INTGRIB1(K1F,KGDS1F,K2F,KGDS2F,MRL,MRO,
     &                  IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,
     &                  IB2,L2,F2,G2,IRET)
        ELSE
          IRET=101
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE INTGRIB1(K1F,KGDS1F,K2F,KGDS2F,MRL,MRO,
     &                    IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,
     &                    IB2,L2,F2,G2,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    INTGRIB1    INTERPOLATE FIELD
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: INTERPOLATE FIELD.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL INTGRIB1(K1F,KGDS1F,K2F,KGDS2F,MRL,MRO,
C    &                    IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,
C    &                    IB2,L2,F2,G2,IRET)
C   INPUT ARGUMENTS:
C     K1F          INTEGER REGULARIZED INPUT DIMENSION
C     KGDS1F       INTEGER (200) REGULARIZED INPUT GRID PARAMETERS
C     K2F          INTEGER REGULARIZED OUTPUT DIMENSION
C     KGDS2F       INTEGER (200) REGULARIZED OUTPUT GRID PARAMETERS
C     MRL          INTEGER DIMENSION OF RLAT AND RLON
C     MRO          INTEGER DIMENSION OF CROT AND SROT
C     IV           INTEGER VECTOR FLAG
C     IP           INTEGER INTERPOLATION TYPE
C     IPOPT        INTEGER (20) INTERPOLATION OPTIONS
C     KGDS1        INTEGER (200) INPUT GRID PARAMETERS
C     K1           INTEGER INPUT DIMENSION
C     IB1          INTEGER INPUT BITMAP FLAG
C     L1           LOGICAL (K1) INPUT BITMAP IF IB1=1
C     F1           REAL (K1) INPUT FIELD
C     G1           REAL (K1) INPUT Y-COMPONENT IF IV=1
C     KGDS2        INTEGER (200) OUTPUT GRID PARAMETERS
C     K2           INTEGER OUTPUT DIMENSION
C     IB2          INTEGER OUTPUT BITMAP FLAG
C     L2           LOGICAL (K2) OUTPUT BITMAP
C     F2           REAL (K2) OUTPUT FIELD
C     G2           REAL (K2) OUTPUT Y-COMPONENT IF IV=1
C
C SUBPROGRAMS CALLED:
C   LENGDSF
C   INTGRIB1
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER IPOPT(20)
      INTEGER KGDS1(200),KGDS2(200)
      LOGICAL L1(K1),L2(K2)
      REAL F1(K1),F2(K2),G1(K1),G2(K2)
      INTEGER KGDS1F(200),KGDS2F(200)
      LOGICAL L1F(K1F),L2F(K2F)
      REAL F1F(K1F),F2F(K2F),G1F(K1F),G2F(K2F)
      REAL RLAT(MRL),RLON(MRL),CROT(MRO),SROT(MRO)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  REGLR TO REGLR SCALAR
      IF(K1F.EQ.1.AND.K2F.EQ.1.AND.IV.EQ.0) THEN
        CALL IPOLATES(IP,IPOPT,KGDS1,KGDS2,K1,K2,1,IB1,L1,F1,
     &                KI,RLAT,RLON,IB2,L2,F2,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  IRREG TO REGLR SCALAR
      ELSEIF(K1F.NE.1.AND.K2F.EQ.1.AND.IV.EQ.0) THEN
        CALL IPXWAFS(1,K1,K1F,1,KGDS1,F1,KGDS1F,F1F,IRET)
        IF(IRET.EQ.0) THEN
          CALL IPOLATES(IP,IPOPT,KGDS1F,KGDS2,K1F,K2,1,IB1,L1F,F1F,
     &                  KI,RLAT,RLON,IB2,L2,F2,IRET)
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  REGLR TO IRREG SCALAR
      ELSEIF(K1F.EQ.1.AND.K2F.NE.1.AND.IV.EQ.0) THEN
        CALL IPOLATES(IP,IPOPT,KGDS1,KGDS2F,K1,K2F,1,IB1,L1,F1,
     &                KI,RLAT,RLON,IB2,L2F,F2F,IRET)
        IF(IRET.EQ.0) THEN
          CALL IPXWAFS(-1,K2,K2F,1,KGDS2,F2,KGDS2F,F2F,IRET)
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  IRREG TO IRREG SCALAR
      ELSEIF(K1F.NE.1.AND.K2F.NE.1.AND.IV.EQ.0) THEN
        CALL IPXWAFS(1,K1,K1F,1,KGDS1,F1,KGDS1F,F1F,IRET)
        IF(IRET.EQ.0) THEN
          CALL IPOLATES(IP,IPOPT,KGDS1F,KGDS2F,K1F,K2F,1,IB1,L1F,F1F,
     &                  KI,RLAT,RLON,IB2,L2F,F2F,IRET)
          IF(IRET.EQ.0) THEN
            CALL IPXWAFS(-1,K2,K2F,1,KGDS2,F2,KGDS2F,F2F,IRET)
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  REGLR TO REGLR VECTOR
      ELSEIF(K1F.EQ.1.AND.K2F.EQ.1.AND.IV.NE.0) THEN
        CALL IPOLATEV(IP,IPOPT,KGDS1,KGDS2,K1,K2,1,IB1,L1,F1,G1,
     &                KI,RLAT,RLON,CROT,SROT,IB2,L2,F2,G2,IRET)
        IF(IRET.EQ.0.AND.KI.EQ.K2-1) THEN
          F2(K2)=0
          G2(K2)=0
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  IRREG TO REGLR VECTOR
      ELSEIF(K1F.NE.1.AND.K2F.EQ.1.AND.IV.NE.0) THEN
        CALL IPXWAFS(1,K1,K1F,1,KGDS1,F1,KGDS1F,F1F,IRET)
        CALL IPXWAFS(1,K1,K1F,1,KGDS1,G1,KGDS1F,G1F,IRET)
        IF(IRET.EQ.0) THEN
          CALL IPOLATEV(IP,IPOPT,KGDS1F,KGDS2,K1F,K2,1,IB1,L1F,F1F,G1F,
     &                  KI,RLAT,RLON,CROT,SROT,IB2,L2,F2,G2,IRET)
          IF(IRET.EQ.0.AND.KI.EQ.K2-1) THEN
            F2(K2)=0
            G2(K2)=0
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  REGLR TO IRREG VECTOR
      ELSEIF(K1F.EQ.1.AND.K2F.NE.1.AND.IV.NE.0) THEN
        CALL IPOLATEV(IP,IPOPT,KGDS1,KGDS2F,K1,K2F,1,IB1,L1,F1,G1,
     &                KI,RLAT,RLON,CROT,SROT,IB2,L2F,F2F,G2F,IRET)
        IF(IRET.EQ.0) THEN
          CALL IPXWAFS(-1,K2,K2F,1,KGDS2,F2,KGDS2F,F2F,IRET)
          CALL IPXWAFS(-1,K2,K2F,1,KGDS2,G2,KGDS2F,G2F,IRET)
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  IRREG TO IRREG VECTOR
      ELSEIF(K1F.NE.1.AND.K2F.NE.1.AND.IV.NE.0) THEN
        CALL IPXWAFS(1,K1,K1F,1,KGDS1,F1,KGDS1F,F1F,IRET)
        CALL IPXWAFS(1,K1,K1F,1,KGDS1,G1,KGDS1F,G1F,IRET)
        IF(IRET.EQ.0) THEN
         CALL IPOLATEV(IP,IPOPT,KGDS1F,KGDS2F,K1F,K2F,1,IB1,L1F,F1F,G1F,
     &                  KI,RLAT,RLON,CROT,SROT,IB2,L2F,F2F,G2F,IRET)
          IF(IRET.EQ.0) THEN
            CALL IPXWAFS(-1,K2,K2F,1,KGDS2,F2,KGDS2F,F2F,IRET)
            CALL IPXWAFS(-1,K2,K2F,1,KGDS2,G2,KGDS2F,G2F,IRET)
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      FUNCTION LENGDSF(KGDS,KGDSF)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    LENGDSF     RETURN THE LENGTH OF A FILLED GRID
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: GIVEN A GRID DESCRIPTION SECTION (IN W3FI63 FORMAT),
C   RETURN THE GRID DESCRIPTION SECTION AND SIZE OF ITS REGULARIZED
C   COUNTERPART.  THAT IS, IF THE INPUT GRID IS REGULAR, THEN ITSELF
C   IS RETURNED ALONG WITH ITS GRID SIZE; HOWEVER IF THE INPUT GRID IS
C   ONLY QUASI-REGULAR (SUCH AS THE WAFS GRIDS), THEN ITS FILLED REGULAR
C   VERSION IS RETURNED ALONG WITH ITS FILLED GRID SIZE.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL LENGDSF(KGDS,KGDSF)
C   INPUT ARGUMENTS:
C     KGDS         INTEGER (200) GDS PARAMETERS IN W3FI63 FORMAT
C   OUTPUT ARGUMENTS:
C     KGDSF        INTEGER (200) REGULAR GDS PARAMETERS IN W3FI63 FORMAT
C     LENGDSF      INTEGER SIZE OF REGULARIZED GRID
C
C SUBPROGRAMS CALLED:
C   IPXWAFS
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER KGDS(200),KGDSF(200)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(KGDS(1).EQ.201) THEN
c       KGDSF=KGDS
        do ii = 1, 200
            KGDSF(ii)=KGDS(ii)
        enddo
        LENGDSF=KGDS(7)*KGDS(8)-KGDS(8)/2
      ELSEIF(KGDS(1).EQ.202) THEN
c       KGDSF=KGDS
        do ii = 1, 200
            KGDSF(ii)=KGDS(ii)
        enddo
        LENGDSF=KGDS(7)*KGDS(8)
      ELSEIF(KGDS(20).NE.255) THEN
        CALL IPXWAFS(1,1,1,0,KGDS,DUM,KGDSF,DUMF,IRET)
        IF(IRET.EQ.0) THEN
          LENGDSF=KGDSF(2)*KGDSF(3)
        ELSE
          LENGDSF=0
        ENDIF
      ELSE
c       KGDSF=KGDS
        do ii = 1, 200
            KGDSF(ii)=KGDS(ii)
        enddo
        LENGDSF=KGDS(2)*KGDS(3)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END

	INTEGER FUNCTION ISTRLEN(C)
	CHARACTER*(*) C
	L = LEN(C)
	DO ISTRLEN = L, 1, -1
	   IF (C(ISTRLEN:ISTRLEN).NE.' ') RETURN
	ENDDO
	ISTRLEN = 0
	RETURN
	END

C-----------------------------------------------------------------------
      SUBROUTINE PUTGBN(LUGB,KF,KPDS,KGDS,IBS,NBITS,LB,F,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: PUTGBN         PACKS AND WRITES A GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 94-04-01
C
C ABSTRACT: PACK AND WRITE A GRIB MESSAGE.
C   THIS SUBPROGRAM IS NEARLY THE INVERSE OF GETGB.
C
C PROGRAM HISTORY LOG:
C   94-04-01  IREDELL
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C
C USAGE:    CALL PUTGBN(LUGB,KF,KPDS,KGDS,NBITS,LB,F,IRET)
C   INPUT ARGUMENTS:
C     LUGB         INTEGER UNIT OF THE UNBLOCKED GRIB DATA FILE
C     KF           INTEGER NUMBER OF DATA POINTS
C     KPDS         INTEGER (200) PDS PARAMETERS
C          (1)   - ID OF CENTER
C          (2)   - GENERATING PROCESS ID NUMBER
C          (3)   - GRID DEFINITION
C          (4)   - GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)
C          (5)   - INDICATOR OF PARAMETER
C          (6)   - TYPE OF LEVEL
C          (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
C          (8)   - YEAR INCLUDING (CENTURY-1)
C          (9)   - MONTH OF YEAR
C          (10)  - DAY OF MONTH
C          (11)  - HOUR OF DAY
C          (12)  - MINUTE OF HOUR
C          (13)  - INDICATOR OF FORECAST TIME UNIT
C          (14)  - TIME RANGE 1
C          (15)  - TIME RANGE 2
C          (16)  - TIME RANGE FLAG
C          (17)  - NUMBER INCLUDED IN AVERAGE
C          (18)  - VERSION NR OF GRIB SPECIFICATION
C          (19)  - VERSION NR OF PARAMETER TABLE
C          (20)  - NR MISSING FROM AVERAGE/ACCUMULATION
C          (21)  - CENTURY OF REFERENCE TIME OF DATA
C          (22)  - UNITS DECIMAL SCALE FACTOR
C          (23)  - SUBCENTER NUMBER
C          (24)  - PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
C                  128 IF FORECAST FIELD ERROR
C                   64 IF BIAS CORRECTED FCST FIELD
C                   32 IF SMOOTHED FIELD
C                  WARNING: CAN BE COMBINATION OF MORE THAN 1
C          (25)  - PDS BYTE 30, NOT USED
C     KGDS         INTEGER (200) GDS PARAMETERS
C          (1)   - DATA REPRESENTATION TYPE
C          (19)  - NUMBER OF VERTICAL COORDINATE PARAMETERS
C          (20)  - OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE
C                  PARAMETERS
C                  OR
C                  OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS
C                  IN EACH ROW
C                  OR
C                  255 IF NEITHER ARE PRESENT
C          (21)  - FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID
C          (22)  - NUMBER OF WORDS IN EACH ROW
C       LATITUDE/LONGITUDE GRIDS
C          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
C          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LA(2) LATITUDE OF EXTREME POINT
C          (8)   - LO(2) LONGITUDE OF EXTREME POINT
C          (9)   - DI LATITUDINAL DIRECTION OF INCREMENT
C          (10)  - DJ LONGITUDINAL DIRECTION INCREMENT
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C       GAUSSIAN  GRIDS
C          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
C          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LA(2) LATITUDE OF EXTREME POINT
C          (8)   - LO(2) LONGITUDE OF EXTREME POINT
C          (9)   - DI LATITUDINAL DIRECTION OF INCREMENT
C          (10)  - N - NR OF CIRCLES POLE TO EQUATOR
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - NV - NR OF VERT COORD PARAMETERS
C          (13)  - PV - OCTET NR OF LIST OF VERT COORD PARAMETERS
C                             OR
C                  PL - LOCATION OF THE LIST OF NUMBERS OF POINTS IN
C                       EACH ROW (IF NO VERT COORD PARAMETERS
C                       ARE PRESENT
C                             OR
C                  255 IF NEITHER ARE PRESENT
C       POLAR STEREOGRAPHIC GRIDS
C          (2)   - N(I) NR POINTS ALONG LAT CIRCLE
C          (3)   - N(J) NR POINTS ALONG LON CIRCLE
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LOV GRID ORIENTATION
C          (8)   - DX - X DIRECTION INCREMENT
C          (9)   - DY - Y DIRECTION INCREMENT
C          (10)  - PROJECTION CENTER FLAG
C          (11)  - SCANNING MODE (RIGHT ADJ COPY OF OCTET 28)
C       SPHERICAL HARMONIC COEFFICIENTS
C          (2)   - J PENTAGONAL RESOLUTION PARAMETER
C          (3)   - K      "          "         "
C          (4)   - M      "          "         "
C          (5)   - REPRESENTATION TYPE
C          (6)   - COEFFICIENT STORAGE MODE
C       MERCATOR GRIDS
C          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
C          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LA(2) LATITUDE OF LAST GRID POINT
C          (8)   - LO(2) LONGITUDE OF LAST GRID POINT
C          (9)   - LATIT - LATITUDE OF PROJECTION INTERSECTION
C          (10)  - RESERVED
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - LONGITUDINAL DIR GRID LENGTH
C          (13)  - LATITUDINAL DIR GRID LENGTH
C       LAMBERT CONFORMAL GRIDS
C          (2)   - NX NR POINTS ALONG X-AXIS
C          (3)   - NY NR POINTS ALONG Y-AXIS
C          (4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
C          (5)   - LO1 LON OF ORIGIN (LOWER LEFT)
C          (6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LOV - ORIENTATION OF GRID
C          (8)   - DX - X-DIR INCREMENT
C          (9)   - DY - Y-DIR INCREMENT
C          (10)  - PROJECTION CENTER FLAG
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - LATIN 1 - FIRST LAT FROM POLE OF SECANT CONE INTER
C          (13)  - LATIN 2 - SECOND LAT FROM POLE OF SECANT CONE INTER
C     IBS          INTEGER BINARY SCALE FACTOR (0 TO IGNORE)
C     NBITS        INTEGER NUMBER OF BITS IN WHICH TO PACK (0 TO IGNORE)
C     LB           LOGICAL*1 (KF) BITMAP IF PRESENT
C     F            REAL (KF) DATA
C   OUTPUT ARGUMENTS:
C     IRET         INTEGER RETURN CODE
C                    0      ALL OK
C                    OTHER  W3FI72 GRIB PACKER RETURN CODE
C
C SUBPROGRAMS CALLED:
C   R63W72         MAP W3FI63 PARAMETERS ONTO W3FI72 PARAMETERS
C   GETBIT         GET NUMBER OF BITS AND ROUND DATA
C   W3FI72         PACK GRIB
C   WRYTE          WRITE DATA
C
C REMARKS: SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C   DO NOT ENGAGE THE SAME LOGICAL UNIT FROM MORE THAN ONE PROCESSOR.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  CRAY, WORKSTATIONS
C
C$$$
      INTEGER KPDS(200),KGDS(200)
      LOGICAL LB(KF)
      REAL F(KF)
      PARAMETER(MAXBIT=16)
      INTEGER IBM(KF),IPDS(200),IGDS(200),IBDS(200)
      REAL FR(KF)
      CHARACTER PDS(400),GRIB(1000+KF*(MAXBIT+1)/8)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET W3FI72 PARAMETERS
      CALL R63W72(KPDS,KGDS,IPDS,IGDS)
      IBDS=0
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COUNT VALID DATA
      KBM=KF
      IF(IPDS(7).NE.0) THEN
        KBM=0
        DO I=1,KF
          IF(LB(I)) THEN
            IBM(I)=1
            KBM=KBM+1
          ELSE
            IBM(I)=0
          ENDIF
        ENDDO
        IF(KBM.EQ.KF) IPDS(7)=0
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET NUMBER OF BITS AND ROUND DATA
      IF(KBM.EQ.0) THEN
        DO I=1,KF
          FR(I)=0
        ENDDO
        NBIT=0
      ELSEIF(NBITS.GT.0) THEN
        IF(IPDS(7).EQ.0) THEN
          DO I=1,KF
            FR(I)=F(I)
          ENDDO
        ELSE
          DO I=1,KF
            IF(LB(I)) THEN
              FR(I)=F(I)
            ELSE
              FR(I)=0
            ENDIF
          ENDDO
        ENDIF
        NBIT=NBITS
      ELSE
        CALL GETBIT(IPDS(7),IBS,IPDS(25),KF,IBM,F,FR,FMIN,FMAX,NBIT)
        NBIT=MIN(NBIT,MAXBIT)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PACK AND WRITE GRIB DATA
      CALL W3FI72(0,FR,0,NBIT,0,IPDS,PDS,
     &            1,255,IGDS,0,0,IBM,KF,IBDS,
     &            KFO,GRIB,LGRIB,IRET)
      IF(IRET.EQ.0) CALL WRYTE(LUGB,LGRIB,GRIB)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

