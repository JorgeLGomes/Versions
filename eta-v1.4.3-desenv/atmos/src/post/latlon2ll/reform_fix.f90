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
      integer               :: fcthr

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
      real,allocatable,dimension(:,:)      :: ur2m1,ur2m !Relative Humidity 2m
!GSM
!
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
      read(5,*)ITAG,NRSTRT,NPINCR,dataset
      write(6,*)ITAG,NRSTRT,NPINCR,dataset
      
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
      allocate(uanem(im,jm))
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

      OPEN(unit=51,file=datset//itag//'.t00s',status='old',form='unformatted')

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
      
      ltsoil=.false.
      lqsoil=.false.
      do while (infloop)
        READ(51,END=999) FIELD,SFC
        print*, field, sfc
        READ(51) GRIB
        elseif (field.eq.'SURFACE HEIGHT      ') then; zs=grib
        elseif (field.eq.'LAND/SEA MASK       ') then; lsmk=grib
!
      enddo
 999  continue
!
      OPEN(unit=53,file='latlon_fix2d_'//itag,STATUS='UNKNOWN', &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny)

       irec=1
       write(53,rec=irec)((zs(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((lsmk(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
      stop 
      end program leit 

