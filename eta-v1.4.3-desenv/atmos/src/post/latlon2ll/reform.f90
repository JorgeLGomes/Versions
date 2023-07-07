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
      real,allocatable,dimension(:,:)      :: vegtyp    !VEGETATION TYPE
      real,allocatable,dimension(:,:)      :: vsst      !SFC V WIND STRESS
      real,allocatable,dimension(:,:)      :: vstrm     !V COMP STORM MOTION
      real,allocatable,dimension(:,:)      :: vtrop     !
      real,allocatable,dimension(:,:)      :: zorl      !ROUGHNESS LENGTH
      real,allocatable,dimension(:,:)      :: zs        !SURFACE HEIGHT
      real,allocatable,dimension(:,:)      :: GRIB
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
        if (field.eq.    'MESINGER MEAN SLP   ') then; pnmm=grib/100
        elseif (field.eq.'SURFACE PRESSURE    ') then; ps=grib/100
        elseif (field.eq.'SURFACE HEIGHT      ') then; zs=grib
        elseif (field.eq.'SFC (SKIN) TEMPRATUR') then; ts=grib
        elseif (field.eq.'SURFACE SPEC HUMID  ') then; qs=grib
        elseif (field.eq.'SHELTER TEMPERATURE ') then; tsh=grib
        elseif (field.eq.'MAX TEMPERATURE     ') then; tmax=grib
	elseif (field.eq.'MIN TEMPERATURE     ') then; tmin=grib
        elseif (field.eq.'SHELTER DEWPOINT    ') then; tdsh=grib
        elseif (field.eq.'U WIND AT ANEMOM HT ') then; uanem=grib
        elseif (field.eq.'V WIND AT ANEMOM HT ') then; vanem=grib
        elseif (field.eq.'U WIND 100M         ') then; u100m=grib
        elseif (field.eq.'V WIND 100M         ') then; v100m=grib
        elseif (field.eq.'ACM TOTAL PRECIP    ') then; ttprec=grib/1000
        elseif (field.eq.'ACM CONVCTIVE PRECIP') then; cvprec=grib/1000
        elseif (field.eq.'ACM GRD SCALE PRECIP') then; lsprec=grib/1000
        elseif (field.eq.'ACM SNOWFALL        ') then; snofal=grib/1000
        elseif (field.eq.'AVE SFC LATHEAT FX  ') then; clsf=grib
        elseif (field.eq.'AVE SFC SENHEAT FX  ') then; cssf=grib
        elseif (field.eq.'CNVCT AVBL POT ENRGY') then; cape=grib
        elseif (field.eq.'CNVCT INHIBITION    ') then; cine=grib
        elseif (field.eq.'PRECIPITABLE WATER  ') then; pw=grib
        elseif (field.eq.'AVG TOTAL CLD FRAC  ') then; cld=grib/100
        elseif (field.eq.'AVE OUTGO TOA LW RAD') then; role=grib
        elseif (field.eq.'U COMP STORM MOTION ') then; ustrm=grib
        elseif (field.eq.'V COMP STORM MOTION ') then; vstrm=grib
        elseif (field.eq.'PRESS AT TROPOPAUSE ') then; ptrop=grib/100
        elseif (field.eq.'HEIGHT OF FRZ LVL   ') then; htfrz=grib
        elseif (field.eq.'REL HUMID AT FRZ LVL') then; rhfrz=grib
        elseif (field.eq.'LIFTED INDEX--BEST  ') then; li=grib
        elseif (field.eq.'ROUGHNESS LENGTH    ') then; zorl=grib
        elseif (field.eq.'SFC U WIND STRESS   ') then; usst=grib
        elseif (field.eq.'SFC V WIND STRESS   ') then; vsst=grib
        elseif (field.eq.'AVE INCMG SFC SW RAD') then; ocis=grib
        elseif (field.eq.'AVE INCMG SFC LW RAD') then; olis=grib
        elseif (field.eq.'AVE OUTGO SFC SW RAD') then; oces=-1*grib
        elseif (field.eq.'AVE OUTGO SFC LW RAD') then; oles=-1*grib
        elseif (field.eq.'AVE OUTGO TOA SW RAD') then; roce=grib
        elseif (field.eq.'ACC POT EVAPORATION ') then; potevap=grib/1000
        elseif (field.eq.'SOIL TEMPERATURE    ') then
          if (.not.ltsoil) then  ! it is the first soil layer
            ltsoil=.true.
            tsoil01=grib
           else
             tsoil20=grib
           endif
        elseif (field.eq.'SOIL MOISTURE       ') then
          if (.not.lqsoil) then  ! it is the first soil layer
            lqsoil=.true.
            qsoil01=grib
          else
            qsoil20=grib
          endif
        elseif (field.eq.'ACM STORM SFC RNOFF ') then; rnof=grib/1000
        elseif (field.eq.'SFC MIDDAY ALBEDO   ') then; alb=grib
        elseif (field.eq.'LAND/SEA MASK       ') then; lsmk=grib
!
!  retrieving 3d fields
!
        elseif (field.eq.'OMEGA ON PRESS SFCS ') then; kw=kw+1; ome(1:im,1:jm,kw)=grib
        elseif (field.eq.'HEIGHT OF PRESS SFCS') then; kp=kp+1; phi(1:im,1:jm,kp)=grib
        elseif (field.eq.'TEMP ON PRESS SFCS  ') then; kt=kt+1; temp(1:im,1:jm,kt)=grib
        elseif (field.eq.'SPEC HUM ON P SFCS  ') then; kq=kq+1; q(1:im,1:jm,kq)=grib
        elseif (field.eq.'REL HUMID ON P SFCS ') then; kr=kr+1; rh(1:im,1:jm,kr)=grib
        elseif (field.eq.'U WIND ON PRESS SFCS') then; ku=ku+1; u(1:im,1:jm,ku)=grib
        elseif (field.eq.'V WIND ON PRESS SFCS') then; kv=kv+1; v(1:im,1:jm,kv)=grib
        elseif (field.eq.'CLOUD WATR ON P SFCS') then; kc=kc+1; clwt(1:im,1:jm,kc)=grib
        elseif (field.eq.'CLOUD ICE ON P SFCS ') then; ki=ki+1; cips(1:im,1:jm,ki)=grib
!
!       Variaveis novas introduzidas pelo ETA 2d
!
        elseif (field.eq.'SOIL MOISTURE AVAIL ') then; smav=grib/100
        elseif (field.eq.'AVE GROUND HEAT FX  ') then; aghf=grib
        elseif (field.eq.'MAX WIND PRESS LEVEL') then; mwpl=grib/100
        elseif (field.eq.'U COMP MAX WIND     ') then; ucmw=grib
        elseif (field.eq.'V COMP MAX WIND     ') then; vcmw=grib
        elseif (field.eq.'LOW CLOUD FRACTION  ') then; lowcf=grib/100
        elseif (field.eq.'MID CLOUD FRACTION  ') then; midcf=grib/100
        elseif (field.eq.'HIGH CLOUD FRACTION ') then; higcf=grib/100
        elseif (field.eq.'CLOUD BOT PRESSURE  ') then; cbps=grib/100
        elseif (field.eq.'CLOUD TOP PRESSURE  ') then; ctps=grib/100
        endif
      enddo
 999  continue
!
! Calculate Theta_e
!
      call thetae(im,jm,kmax,pr,temp,q,thte)
!
!
!  save fields for Grads plotting. First 2D, then 3D fields
!
      OPEN(unit=53,file='latlon_'//itag,STATUS='UNKNOWN', &
     &        form='unformatted',access='direct' ,           &
     &        recl=4*nx*ny)

       irec=1
       write(53,rec=irec)((pnmm(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((ps(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((zs(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((lsmk(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((tsh(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((tmax(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((tmin(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((tdsh(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((uanem(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((vanem(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((u100m(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((v100m(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       where(ttprec<=0.) ttprec=0.
       write(53,rec=irec)((ttprec(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       where(cvprec<=0.) cvprec=0.
       write(53,rec=irec)((cvprec(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       where(lsprec<=0.) lsprec=0.
       write(53,rec=irec)((lsprec(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((snofal(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((clsf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((cssf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((aghf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((ts(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((qs(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((tsoil01(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((tsoil20(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((qsoil01(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((qsoil20(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((smav(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((rnof(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((zorl(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((potevap(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((usst(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((vsst(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((lowcf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((midcf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((higcf(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((cld(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((ocis(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((olis(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((oces(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((oles(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((roce(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((role(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((alb(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((cape(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((cine(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((li(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((pw(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((ptrop(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((htfrz(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((rhfrz(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((mwpl(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((ucmw(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((vcmw(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((cbps(i,j),i=iw,ie),j=js,jn)
       irec=irec+1
       write(53,rec=irec)((ctps(i,j),i=iw,ie),j=js,jn)


      print*, ' 2D irec=',irec-1
      irec=irec+1
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


!      stop
      end program leit
