!
!           PROGRAMA DE LEITURA de arquivos tipo ETAGRDhh.tmhh
!
      program leit

      use variables

      integer               :: dt,kgtype,imdlty,igout,jgout
      real                  :: dx,lonw,alatvt,lats,dy
      character(len=20)     :: field
      character(len=6)      :: nufile,outype,proj,readco,readll
      character(len=12)     :: compl
      character(len=6)      :: itag,datset,dataset
      logical               :: ltsoil, lqsoil,north
      logical               :: infloop
      logical               :: split
      integer               :: fcthr
      integer               :: FreqOut

      integer, parameter    :: NSOIL=8
      real,    parameter    :: UNDEF=1.E+20
!
! Dimension 2D variables
!
      real,allocatable,dimension(:,:)      :: aghf      !
      real,allocatable,dimension(:,:)      :: alb       !SFC MIDDAY ALBEDO
      real,allocatable,dimension(:,:)      :: ashf      !AVE GROUND HEAT FX
      real,allocatable,dimension(:,:)      :: bost      !
      real,allocatable,dimension(:,:)      :: cape      !CNVCT AVBL POT ENRGY
      real,allocatable,dimension(:,:)      :: cbh       !
      real,allocatable,dimension(:,:)      :: cbps      !CLOUD BOT PRESSURE
      real,allocatable,dimension(:,:)      :: cine      !CNVCT INHIBITION
      real,allocatable,dimension(:,:)      :: cld       !AVG TOTAL CLD FRAC
      real,allocatable,dimension(:,:)      :: clsf      !AVE SFC LATHEAT FX
      real,allocatable,dimension(:,:)      :: cssf      !AVE SFC SENHEAT FX
      real,allocatable,dimension(:,:)      :: cth       !
      real,allocatable,dimension(:,:)      :: ctps      !CLOUD TOP PRESSURE
      real,allocatable,dimension(:,:)      :: cttp      !
      real,allocatable,dimension(:,:)      :: cvprec    !ACM CONVCTIVE PRECIP
      real,allocatable,dimension(:,:)      :: evpp      !
      real,allocatable,dimension(:,:)      :: ffcs      !
      real,allocatable,dimension(:,:)      :: frvel     !
      real,allocatable,dimension(:,:)      :: gvco      !
      real,allocatable,dimension(:,:)      :: hatt      !
      real,allocatable,dimension(:,:)      :: higcf     !HIGH CLOUD FRACTION
      real,allocatable,dimension(:,:)      :: htfrz     !HEIGHT OF FRZ LVL
      real,allocatable,dimension(:,:)      :: icld      !
      real,allocatable,dimension(:,:)      :: ighf      !
      real,allocatable,dimension(:,:)      :: ilwr      !
      real,allocatable,dimension(:,:)      :: ipty      !
      real,allocatable,dimension(:,:)      :: iswr      !
      real,allocatable,dimension(:,:)      :: li        !LIFTED INDEX--BEST
      real,allocatable,dimension(:,:)      :: lowcf     !LOW CLOUD FRACTION
      real,allocatable,dimension(:,:)      :: lsmk      !LAND/SEA MASK
      real,allocatable,dimension(:,:)      :: lsprec    !ACM GRD SCALE PRECIP
      real,allocatable,dimension(:,:)      :: midcf     !MID CLOUD FRACTION
      real,allocatable,dimension(:,:)      :: mwhl      !
      real,allocatable,dimension(:,:)      :: mwpl      !MAX WIND PRESS LEVEL
      real,allocatable,dimension(:,:)      :: oces      !AVE OUTGO SFC SW RAD
      real,allocatable,dimension(:,:)      :: ocis      !AVE INCMG SFC SW RAD
      real,allocatable,dimension(:,:)      :: oles      !AVE OUTGO SFC LW RAD
      real,allocatable,dimension(:,:)      :: olis      !AVE INCMG SFC LW RAD
      real,allocatable,dimension(:,:)      :: pnmm      !MESINGER MEAN SLP
      real,allocatable,dimension(:,:)      :: pnms      !
      real,allocatable,dimension(:,:)      :: potevap   !ACC POT EVAPORATION
      real,allocatable,dimension(:,:)      :: ps        !SURFACE PRESSURE
      real,allocatable,dimension(:,:)      :: ptrop     !PRESS AT TROPOPAUSE
      real,allocatable,dimension(:,:)      :: pw        !PRECIPITABLE WATER
      real,allocatable,dimension(:,:)      :: qs        !SURFACE SPEC HUMID
      real,allocatable,dimension(:,:)      :: qsoil01   !SOIL MOISTURE
      real,allocatable,dimension(:,:)      :: qsoil20   !SOIL MOISTURE
      real,allocatable,dimension(:,:)      :: rctt      !
      real,allocatable,dimension(:,:)      :: rhfrz     !REL HUMID AT FRZ LVL
      real,allocatable,dimension(:,:)      :: rnof      !ACM STORM SFC RNOFF
      real,allocatable,dimension(:,:)      :: rnofsg    !ACM BSFL-GDWR RNOFF
      real,allocatable,dimension(:,:)      :: roce      !AVE OUTGO TOA SW RAD
      real,allocatable,dimension(:,:)      :: role      !AVE OUTGO TOA LW RAD
      real,allocatable,dimension(:,:)      :: sexc      !
      real,allocatable,dimension(:,:)      :: smal      !
      real,allocatable,dimension(:,:)      :: smav      !SOIL MOISTURE AVAIL
      real,allocatable,dimension(:,:)      :: smof      !
      real,allocatable,dimension(:,:)      :: snofal    !ACM SNOWFALL
      real,allocatable,dimension(:,:)      :: srhel     !
      real,allocatable,dimension(:,:)      :: strop     !
      real,allocatable,dimension(:,:)      :: tdsh      !SHELTER DEWPOINT
      real,allocatable,dimension(:,:)      :: us2m      !SHELTER SPEC HUMID
      real,allocatable,dimension(:,:)      :: ts        !SFC (SKIN) TEMPRATUR
      real,allocatable,dimension(:,:)      :: tsh       !SHELTER TEMPERATURE
      real,allocatable,dimension(:,:)      :: tmax      !MAX TEMPERATURE
      real,allocatable,dimension(:,:)      :: tmin      !MIN TEMPERATURE
      real,allocatable,dimension(:,:)      :: tsoil01   !SOIL TEMPERATURE
      real,allocatable,dimension(:,:)      :: tsoil20   !SOIL TEMPERATURE
      real,allocatable,dimension(:,:)      :: tsom      !
      real,allocatable,dimension(:,:)      :: ttprec    !ACM TOTAL PRECIP
      real,allocatable,dimension(:,:)      :: ttrop     !
      real,allocatable,dimension(:,:)      :: uanem     !U WIND AT ANEMOM HT
      real,allocatable,dimension(:,:)      :: maguv     !MAG U10 V10
      real,allocatable,dimension(:,:)      :: u100m     !U WIND AT 100M
      real,allocatable,dimension(:,:)      :: v100m     !V WIND AT 100M
      real,allocatable,dimension(:,:)      :: ucmw      !U COMP MAX WIND
      real,allocatable,dimension(:,:)      :: usst      !SFC U WIND STRESS
      real,allocatable,dimension(:,:)      :: ustrm     !U COMP STORM MOTION
      real,allocatable,dimension(:,:)      :: utrop     !
      real,allocatable,dimension(:,:)      :: vanem     !V WIND AT ANEMOM HT
      real,allocatable,dimension(:,:)      :: vcmw      !V COMP MAX WIND
!      real,allocatable,dimension(:,:)      :: vegtyp    !VEGETATION TYPE 
      real,allocatable,dimension(:,:)      :: vsst      !SFC V WIND STRESS
      real,allocatable,dimension(:,:)      :: vstrm     !V COMP STORM MOTION
      real,allocatable,dimension(:,:)      :: vtrop     !
      real,allocatable,dimension(:,:)      :: zorl      !ROUGHNESS LENGTH
      real,allocatable,dimension(:,:)      :: zs        !SURFACE HEIGHT
      real,allocatable,dimension(:,:)      :: GRIB
!GSM
      real,allocatable,dimension(:,:)      :: ur2m1
      real,allocatable,dimension(:,:)      :: ur2m      !Relative Humidity 2m
!GSM
      real,allocatable,dimension(:,:)      :: QUint
      real,allocatable,dimension(:,:)      :: QVint      ! Vert Int Cloud Water
      real,allocatable,dimension(:,:)      :: CWint
      real,allocatable,dimension(:,:)      :: CIint      ! Vert Int Cloud Ice
      real,allocatable,dimension(:,:)      :: snodp
      real,allocatable,dimension(:,:)      :: snofr      ! Snow Depth, Snow Fraction
      real,allocatable,dimension(:,:)      :: snoeq      ! Snow Water Equivalent
      real,allocatable,dimension(:,:)      :: Hpbl       ! Boundary Lyr Height
      real,allocatable,dimension(:,:)      :: pmsl       ! Shuell msl Pressure
!
      !CHOU  3D SOIL
      real,allocatable,dimension(:,:,:) :: tsoil8    !Soil Temperature (8 layers)
      real,allocatable,dimension(:,:,:) :: qsoil8    !Soil Moisture    (8 layers)
!
! Dimension 3D variables
!
      real,allocatable,dimension(:,:,:) :: cips  !CLOUD ICE ON P SFCS
      real,allocatable,dimension(:,:,:) :: clwt  !CLOUD WATR ON P SFCS
      real,allocatable,dimension(:,:,:) :: ome   !OMEGA ON PRESS SFCS
      real,allocatable,dimension(:,:,:) :: phi   !HEIGHT OF PRESS SFCS
      real,allocatable,dimension(:,:,:) :: q      !SPEC HUM ON P SFCS
      real,allocatable,dimension(:,:,:) :: rh    !REL HUMID ON P SFCS
      real,allocatable,dimension(:,:,:) :: temp  !TEMP ON PRESS SFCS
      real,allocatable,dimension(:,:,:) :: thte  !THETAE
      real,allocatable,dimension(:,:,:) :: u     !U WIND ON PRESS SFCS
      real,allocatable,dimension(:,:,:) :: v     !V WIND ON PRESS SFCS

      infloop=.true.
      read(5,*)ITAG,NRSTRT,NPINCR,dataset,FreqOut,split
      write(6,*)ITAG,NRSTRT,NPINCR,dataset,FreqOut,split
      
      open (unit=15,file='namelist_'//itag,STATUS='old',form='formatted')
      read(15,1500) fcthr
1500  format(I6)
      print*, 'fcthr= ',fcthr         

!     READ OUTPUT GRID SPECIFICATIONS.
!
      OPEN(18,FILE='cntrl.parm_'//dataset,status='old')
      rewind(18)
      read(18,1000) kgtype
      read(18,1000) imdlty
      read(18,1030) datset
      read(18,1030) outype
      read(18,1030) nufile
      read(18,1030) proj
      read(18,1010) north
      read(18,1000) im
      read(18,1000) jm
      read(18,1020) dx
      read(18,1020) lonw
      read(18,1020) alatvt
      read(18,1020) lats
      read(18,1020) dy
      read(18,1030) readll
      read(18,1030) readco
 1000 format(T28,I5)
 1010 format(T28,L1)
 1020 format(T28,F11.6)
 1030 format(T28,A6)
      print*,"im: ",im," jm: ",jm



! Allocate variables
      allocate(GRIB(im,jm))
      allocate(pnmm(im,jm))
      allocate(ps(im,jm))
      allocate(zs(im,jm))
      allocate(ts(im,jm))
      allocate(qs(im,jm))
      allocate(tsh(im,jm))
      allocate(tmax(im,jm))
      allocate(tmin(im,jm))
      allocate(tdsh(im,jm))
      allocate(us2m(im,jm))
      allocate(uanem(im,jm))
      allocate(maguv(im,jm))
      allocate(vanem(im,jm))
      allocate(u100m(im,jm))
      allocate(v100m(im,jm))
      allocate(ttprec(im,jm))
      allocate(cvprec(im,jm))
      allocate(lsprec(im,jm))
      allocate(snofal(im,jm))
      allocate(clsf(im,jm))
      allocate(cssf(im,jm))
      allocate(cape(im,jm))
      allocate(cine(im,jm))
      allocate(pw(im,jm))
      allocate(cld(im,jm))
      allocate(role(im,jm))
      allocate(ustrm(im,jm))
      allocate(vstrm(im,jm))
      allocate(ptrop(im,jm))
      allocate(htfrz(im,jm))
      allocate(rhfrz(im,jm))
      allocate(li(im,jm))
      allocate(zorl(im,jm))
      allocate(usst(im,jm))
      allocate(vsst(im,jm))
      allocate(ocis(im,jm))
      allocate(olis(im,jm))
      allocate(oces(im,jm))
      allocate(oles(im,jm))
      allocate(roce(im,jm))
      allocate(potevap(im,jm))
      allocate(tsoil01(im,jm))
      allocate(tsoil20(im,jm))
      allocate(qsoil01(im,jm))
      allocate(qsoil20(im,jm))
      allocate(rnof(im,jm))
      allocate(rnofsg(im,jm))
      allocate(alb(im,jm))
      allocate(lsmk(im,jm))
!       Variaveis novas introduzidas pelo ETA 2d
!
      allocate(smav(im,jm))
      allocate(aghf(im,jm))
      allocate(mwpl(im,jm))
      allocate(ucmw(im,jm))
      allocate(vcmw(im,jm))
      allocate(lowcf(im,jm))
      allocate(midcf(im,jm))
      allocate(higcf(im,jm))
      allocate(cbps(im,jm))
      allocate(ctps(im,jm))
!GSM
      allocate(ur2m1(im,jm))
      allocate(ur2m(im,jm))
!GSM
      allocate(quint(im,jm))
      allocate(qvint(im,jm))
!Chou
      allocate(Cwint(im,jm))
      allocate(Ciint(im,jm))
      allocate(Snodp(im,jm))
      allocate(snofr(im,jm))
      allocate(snoeq(im,jm))
      allocate(hpbl(im,jm))
      allocate(pmsl(im,jm))

!CHOU 8-layer  soil
      allocate(tsoil8(im,jm,nsoil))
      allocate(qsoil8(im,jm,nsoil))
!      
!       Variaveis 3D
      allocate(ome(im,jm,kmax))
      allocate(phi(im,jm,kmax))
      allocate(temp(im,jm,kmax))
      allocate(q(im,jm,kmax))
      allocate(rh(im,jm,kmax))
      allocate(u(im,jm,kmax))
      allocate(v(im,jm,kmax))
      allocate(clwt(im,jm,kmax))
      allocate(cips(im,jm,kmax))
      allocate(thte(im,jm,kmax))

      OPEN(unit=51,file=datset//itag//'.t00s' &
     &     ,status='old',form='unformatted')

! get I and J indexes for the smaller domain
!      IW= ABS((-83) - (-83))/ 0.40  + 1
!      IE= ABS((-83) - (-25.8))/ 0.40 + 1
!      JS= ABS((-50.2) - (-50.2))/ 0.40 + 1
!      JN= ABS((-50.2) - (12.2))/ 0.40 + 1
      IW= 1
      IE= im
      JS= 1
      JN= jm

!
! get the nx and ny  number of points in x and y direction
      nx = ie - iw + 1
      ny = jn - js + 1
!
! get the exactly lat/lon of the smaller domain
!
      clonll= -lonw + ((iw-1) * dx)
      clonur= -lonw + ((ie-1) * dx)
      clatll= lats + ((js-1) * dy)
      clatur= lats + ((jn-1) * dy)

      print*, iw, ie, js, jn
      print*, clonll, clonur, clatll, clatur

      READ(51)   IHRST,IMM,IDD,IYY,IHH
      print*, ' reading header'
      READ(51)   KGTYP,PROJ,NORTH,IMOUT,JMOUT,POLEI,POLEJ,ALATVT,ALONVT,XMESHL
      kqs=0
      kts=0
      kt=0
      kq=0
      kr=0
      kp=0
      kw=0
      ku=0
      kv=0
      kc=0
      ki=0
!Chou 20221227 8-layer soil
      ks=0
      kx=0
     
      ltsoil=.false.
      lqsoil=.false.
      do while (infloop)
        READ(51,END=999) FIELD,SFC
        print*, field, sfc
        READ(51) GRIB
        SELECT CASE (field)
        CASE ('MESINGER MEAN SLP   ') ; pnmm=grib/100
        CASE ('SHUELL MEAN SLP     ') ; pmsl=grib/100
        CASE ('SURFACE PRESSURE    ') ; ps=grib/100
        CASE ('SURFACE HEIGHT      ') ; zs=grib
        CASE ('SFC (SKIN) TEMPRATUR') ; ts=grib
        CASE ('SURFACE SPEC HUMID  ') ; qs=grib
        CASE ('SHELTER TEMPERATURE ') ; tsh=grib
        CASE ('MAX TEMPERATURE     ') ; tmax=grib
        CASE ('MIN TEMPERATURE     ') ; tmin=grib
        CASE ('SHELTER DEWPOINT    ') ; tdsh=grib
        CASE ('SHELTER SPEC HUMID  ') ; us2m=grib
        CASE ('U WIND AT ANEMOM HT ') ; uanem=grib
        CASE ('V WIND AT ANEMOM HT ') ; vanem=grib
        CASE ('MAG U10 V10         ') ; maguv=grib
        CASE ('U WIND 100M         ') ; u100m=grib
        CASE ('V WIND 100M         ') ; v100m=grib
        CASE ('ACM TOTAL PRECIP    ') ; ttprec=grib/1000
        CASE ('ACM CONVCTIVE PRECIP') ; cvprec=grib/1000
        CASE ('ACM GRD SCALE PRECIP') ; lsprec=grib/1000
        CASE ('ACM SNOWFALL        ') ; snofal=grib/1000
        CASE ('AVE SFC LATHEAT FX  ') ; clsf=-1*grib
        CASE ('AVE SFC SENHEAT FX  ') ; cssf=-1*grib
        CASE ('CNVCT AVBL POT ENRGY') ; cape=grib
        CASE ('CNVCT INHIBITION    ') ; cine=grib
        CASE ('PRECIPITABLE WATER  ') ; pw=grib
        CASE ('AVG TOTAL CLD FRAC  ') ; cld=grib/100
        CASE ('AVE OUTGO TOA LW RAD') ; role=grib
        CASE ('U COMP STORM MOTION ') ; ustrm=grib
        CASE ('V COMP STORM MOTION ') ; vstrm=grib
        CASE ('PRESS AT TROPOPAUSE ') ; ptrop=grib/100
        CASE ('HEIGHT OF FRZ LVL   ') ; htfrz=grib
        CASE ('REL HUMID AT FRZ LVL') ; rhfrz=grib
        CASE ('LIFTED INDEX--BEST  ') ; li=grib
        CASE ('ROUGHNESS LENGTH    ') ; zorl=grib
        CASE ('SFC U WIND STRESS   ') ; usst=grib
        CASE ('SFC V WIND STRESS   ') ; vsst=grib
        CASE ('AVE INCMG SFC SW RAD') ; ocis=grib
        CASE ('AVE INCMG SFC LW RAD') ; olis=grib
        CASE ('AVE OUTGO SFC SW RAD') ; oces=-1*grib
        CASE ('AVE OUTGO SFC LW RAD') ; oles=-1*grib
        CASE ('AVE OUTGO TOA SW RAD') ; roce=grib
        CASE ('ACC POT EVAPORATION ') ; potevap=grib/1000
!
!CHOU 20221227   ADD 8-layer soil T, Q
        CASE ('SOIL TEMPERATURE    ') ;ks=ks+1; tsoil8(1:im,1:jm,ks)=grib
        CASE ('SOIL MOISTURE       ') ;kx=kx+1; qsoil8(1:im,1:jm,kx)=grib;print*,"qsoil",kx,minval(grib),maxval(grib)

        CASE ('ACM STORM SFC RNOFF ') ; rnof=grib/1000
        CASE ('ACM BSFL-GDWR RNOFF ') ; rnofsg=grib/1000
        CASE ('SFC MIDDAY ALBEDO   ') ; alb=grib
        CASE ('LAND/SEA MASK       ') ; lsmk=grib
!
!  retrieving 3d fields
!
        CASE ('OMEGA ON PRESS SFCS ') ; kw=kw+1; ome(1:im,1:jm,kw)=grib
        CASE ('HEIGHT OF PRESS SFCS') ; kp=kp+1; phi(1:im,1:jm,kp)=grib
        CASE ('TEMP ON PRESS SFCS  ') ; kt=kt+1; temp(1:im,1:jm,kt)=grib
        CASE ('SPEC HUM ON P SFCS  ') ; kq=kq+1; q(1:im,1:jm,kq)=grib
        CASE ('REL HUMID ON P SFCS ') ; kr=kr+1; rh(1:im,1:jm,kr)=grib
        CASE ('U WIND ON PRESS SFCS') ; ku=ku+1; u(1:im,1:jm,ku)=grib
        CASE ('V WIND ON PRESS SFCS') ; kv=kv+1; v(1:im,1:jm,kv)=grib
        CASE ('CLOUD WATR ON P SFCS') ; kc=kc+1; clwt(1:im,1:jm,kc)=grib
        CASE ('CLOUD ICE ON P SFCS ') ; ki=ki+1; cips(1:im,1:jm,ki)=grib

!
!       Variaveis novas introduzidas pelo ETA 2d
!
        CASE ('SOIL MOISTURE AVAIL ') ; smav=grib/100
        CASE ('AVE GROUND HEAT FX  ') ; aghf=grib
        CASE ('MAX WIND PRESS LEVEL') ; mwpl=grib/100
        CASE ('U COMP MAX WIND     ') ; ucmw=grib
        CASE ('V COMP MAX WIND     ') ; vcmw=grib
        CASE ('LOW CLOUD FRACTION  ') ; lowcf=grib/100
        CASE ('MID CLOUD FRACTION  ') ; midcf=grib/100
        CASE ('HIGH CLOUD FRACTION ') ; higcf=grib/100
        CASE ('CLOUD BOT PRESSURE  ') ; cbps=grib/100
        CASE ('CLOUD TOP PRESSURE  ') ; ctps=grib/100
        CASE ('VERT INT MOISTxU    ') ; quint=grib
        CASE ('VERT INT MOISTxV    ') ; qvint=grib
        CASE ('VERT INT CLOUD WATER') ; CWint=grib
        CASE ('VERT INT CLOUD ICE  ') ; Ciint=grib
        CASE ('SNOW DEPTH          ') ; snodp=grib*1000 !convert to mm
        CASE ('PERCENT SNOW COVER  ') ; snofr=grib    
        CASE ('SNOW WATER EQUIVALNT') ; snoeq=grib  
        CASE ('BNDRY LYR HEIGHT    ') ; hpbl=grib
        END SELECT
      enddo
 999  continue
!
! Calculate Theta_e
!
      call thetae(im,jm,kmax,pr,temp,q,thte)
!
!GSM Calculate UR2M
          ur2m1(:,:)=(7.5*(tdsh(:,:)-273.15)/(237.7+(tdsh(:,:)-273.15)))-(7.5*(tsh(:,:)-273.15)/(237.7+(tsh(:,:)-273.15)))
          ur2m(:,:)=(10**ur2m1(:,:))*100
!GSM Calculate UR2M
!DCR Calculate MAG U10 V10
          MAGUV(:,:)=((UANEM(:,:)**2)+(VANEM(:,:)**2))**0.5
!DCR Calculate MAG U10 V10
!CHOU standardize to LS4P Protocolo FIx values in the undef area
          do  I=iw,ie
          DO  J=js,jn
           RNOFSG(I,J)   = RNOFSG(I,J) + RNOF(I,J)    !RNOFSG  becomes total runoff
           TSOIL01(I,J)  = (TSOIL8(I,J,1)*0.1 + TSOIL8(I,J,2)*0.3    &
                           & + TSOIL8(I,J,3)*0.6 + TSOIL8(I,J,4)*1.0 &
                           & + TSOIL8(I,J,5)*1.5)/3.5
           QSOIL01(I,J)  = (QSOIL8(I,J,1)*0.1 + QSOIL8(I,J,2)*0.3    &
                           & + QSOIL8(I,J,3)*0.6 + QSOIL8(I,J,4)*1.0 &
                           & + QSOIL8(I,J,5)*1.5)/3.5
           if (SMAV(I,J).GT. UNDEF*0.1) THEN          !points outside integration domain set undef
              PNMM(I,J)  = UNDEF
              PMSL(I,J)  = UNDEF
              PS(I,J)    = UNDEF
              UR2M(I,J)  = UNDEF
              CSSF(I,J)  = UNDEF
              CLSF(I,J)  = UNDEF
              AGHF(I,J)  = UNDEF
              ROLE(I,J)  = UNDEF
              CLD(i,j)   = UNDEF
              ALB(i,j)   = UNDEF
              ROLE(I,J)  = UNDEF
              SNODP(I,J) = UNDEF
              RNOFSG(I,J)= UNDEF
              RNOF(I,J)  = UNDEF
              TSOIL01(I,J)=UNDEF
              QSOIL01(I,J)=UNDEF
             DO K=1,KMAX
              Q(I,J,K)   = UNDEF
             END DO
           ENDIF
          END DO
          END DO
!
!  save fields for Grads plotting. First 2D, then 3D fields
!
      IF ( .true. .AND. itag == '000000') THEN
        OPEN(unit=53,file='topo_lsmsk.bin',STATUS='UNKNOWN', &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny) 
        irec=1
        write(53,rec=irec)((zs(i,j),i=iw,ie),j=js,jn)
        irec=irec+1
        write(53,rec=irec)((lsmk(i,j),i=iw,ie),j=js,jn)
        CLOSE(53) 
      ENDIF     
!
      OPEN(unit=53,file='soil_'//itag,STATUS='UNKNOWN',      &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny)
      irec=1
      do k=1,8
        write(53,rec=irec)((tsoil8(i,j,k),i=iw,ie),j=js,jn)
        irec=irec+1
        write(53,rec=irec)((qsoil8(i,j,k),i=iw,ie),j=js,jn)
        irec=irec+1
      enddo
      CLOSE(53)
!
!
      IF ( split ) THEN
      OPEN(unit=53,file='latlon_2d_'//itag,STATUS='UNKNOWN', &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny)
      ELSE
      OPEN(unit=53,file='latlon_'//itag,STATUS='UNKNOWN', &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny)
      ENDIF
       irec=1      !1 Mesinger M S L Pressure  (hPa)
       write(53,rec=irec)((pnmm(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !2 Surface Pressure         (hPa)
       write(53,rec=irec)((ps(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !3 Shelter Temperature      (K)
       write(53,rec=irec)((tsh(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !4 Max Temperature          (K)
       write(53,rec=irec)((tmax(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !5 Min Temperature          (K)
       write(53,rec=irec)((tmin(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !6 Shelter Dew Temp.        (K)
       write(53,rec=irec)((tdsh(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !7 Shelter specific Humidity 2m     (kg/kg)
       write(53,rec=irec)((us2m(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !8 Relative Humidity 2m     (%)
       write(53,rec=irec)((ur2m(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !9 U 10m                    (m/s)
       write(53,rec=irec)((uanem(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !10 V 10m                    (m/s)
       write(53,rec=irec)((vanem(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !11  
       write(53,rec=irec)((maguv(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !12 U 100m                   (m/s)
       write(53,rec=irec)((u100m(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !13 V 100m                   (m/s)
       write(53,rec=irec)((v100m(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !14 Total  6h Precip.        (m)
       where(ttprec<=0.) ttprec=0.
       write(53,rec=irec)((ttprec(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !15 Conv   6h Precip.        (m)
       where(cvprec<=0.) cvprec=0.
       write(53,rec=irec)((cvprec(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !16 Large Scl  6h Prec.      (m)
       where(lsprec<=0.) lsprec=0.
       write(53,rec=irec)((lsprec(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !17 Snowfall 6h              (m)
       write(53,rec=irec)((snofal(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !18 Time Ave Sfc Lat Ht Flx  ()
       write(53,rec=irec)((clsf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !19 Time Ave Sfc Sen Ht Flx  ()
       write(53,rec=irec)((cssf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !20 Time Ave Ground  Ht Flx  ()
       write(53,rec=irec)((aghf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !21 Sfc (skin) Temperature   (K)
       write(53,rec=irec)((ts(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !22 Soil Temperature 0.1m    (K)
       write(53,rec=irec)((tsoil01(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !23 Soil Moisture Cont. 0.1m (0-1)
       write(53,rec=irec)((qsoil01(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !24 Soil Moisture Avail      ()
       write(53,rec=irec)((smav(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !25 Storm Sfc Rnoff 6h       ()
       write(53,rec=irec)((rnof(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !26 Storm Sfc Rnoff SG 6h    ()
       write(53,rec=irec)((rnofsg(i,j),i=iw,ie),j=js,jn)       
       irec=irec+1 !27 Sfc  U wind Stress       ()
       write(53,rec=irec)((usst(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !28 Sfc  V wind Stress       ()
       write(53,rec=irec)((vsst(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !29 Low cloud fraction       ()
       write(53,rec=irec)((lowcf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !30 Mid cloud fraction       ()
       write(53,rec=irec)((midcf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !31 High cloud fraction      ()
       write(53,rec=irec)((higcf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !32 
       write(53,rec=irec)((cld(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !33 Ave Incmg Sfc SW Rad     (W/m2)
       write(53,rec=irec)((ocis(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !34 Ave Incmg Sfc LW Rad     (W/m2)
       write(53,rec=irec)((olis(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !35 Ave Outgo Sfc SW Rad     (W/m2)
       write(53,rec=irec)((oces(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !36 Ave Outgo Sfc LW Rad     (W/m2)
       write(53,rec=irec)((oles(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !37 Ave Outgo TOA SW Rad     (W/m2)
       write(53,rec=irec)((roce(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !38 Ave Outgo TOA LW Rad     (W/m2)
       write(53,rec=irec)((role(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !39
       write(53,rec=irec)((alb(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !40 CAPE                     (W/m2)
       write(53,rec=irec)((cape(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !41 PRECIPITABLE WATER       (W/m2)
       write(53,rec=irec)((pw(i,j),i=iw,ie),j=js,jn)
!CHOU
       irec=irec+1 !42 VERT INTGRTD Q x U       (W/m2)
       write(53,rec=irec)((quint(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !43 VERT INTGRTD Q x V       (W/m2)
       write(53,rec=irec)((qvint(i,j),i=iw,ie),j=js,jn)
!CHOU  CLOUD integrated water/ice
       irec=irec+1 !44 VERT INTG CLOUD WATER    (kg/kg)
       write(53,rec=irec)((cwint(i,j),i=iw,ie),j=js,jn)
       irec=irec+1 !45 VERT INTG CLOUD ICE      (kg/kg)
       write(53,rec=irec)((ciint(i,j),i=iw,ie),j=js,jn)
!CHOU  Boundary layer height
       irec=irec+1 !46 Plan Bound Lyr Height    (m)
       write(53,rec=irec)((hpbl(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
!
      IF ( mod(fcthr,FreqOut).eq.0 .and. split ) THEN
      close(53)

      OPEN(unit=53,file='latlon_3d_'//itag,STATUS='UNKNOWN', &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny)
      irec=1
      ENDIF 
  
      do k=kmax,1,-1
       write(53,rec=irec)((phi(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((u(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((v(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((temp(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((rh(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((ome(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((q(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((thte(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((clwt(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo

      do k=kmax,1,-1
       write(53,rec=irec)((cips(i,j,k),i=iw,ie),j=js,jn)
       irec=irec+1
      enddo
      
      close(53)    
      stop 
      end program leit 

