    program etatopo

    include "ecommons.h"
    call eta_commons
    call mainprog

    end program etatopo

! CCCCCCCC

    subroutine mainprog

    implicit none

    include "ecommons.h"
    
    integer*4 :: nsorc2,nexvlf,isrch,zeffl,zefflc,check1,check2,check3
    integer*4 :: nssbsp,nssb,nssb2m,nssbp,nssbpp,kne,knw,kse,ksw
    integer*4 :: khl00,khh00,khl22,khh22,ninc
    real*4 :: dtr,pihf,a,cd,rvk2,rz0,dmin,hmin,phl,dlmdl,dphdl,rssb
    real*4 :: sumtxp,sumty,sumtyp,vk,wtamin,area,rarea,htdir,ndir
    real*4 :: sl,saob,shhob,hi,hob,wtarea,hghmn,hodir
    real*4 :: rtsubi,z0,hhdir,alah,sumtx,psbm,pav,psm,shs,hwelp
    integer :: ifill,jx,ix,mthdz0,imin,imax,is,jmin,jmax,llsb
    integer :: js,kssb,hmin1,hmax2,hmin3,nexvlk,nlob,ilob,nhghmn
    integer :: nlim,imn,nhmm,iex,ndirp,kndirt,kldirt,kwdirt,km,kndirs
    integer :: kldirs,kwdirs,ktdirs,nhused,m,ksdb,kspb,kssdbf,kssdb
    integer :: ndbmin,nsbmin,nsbhd,NINTC,nphd
    integer :: mssdbf,msdb,mspb,mssdb,kstart,kend
    integer*4 :: LONNW,JTSUBH,NTSUBH,NTSUBI,ITSUBI, &
    NHREAD,issb,jssb,iqsb,jqsb,kqsb

! p
    real :: lata,latb,lona,lonb,tlm,tph,stph,ctph,sinphi, &
    fact,coslam,wb,sb,tph0,dlm &
    ,tdlm,dph
! p
    real :: tlat,tlon,r30sh,tslat,rtsubh,alm,aph,assb,r30s

   
    REAL,ALLOCATABLE:: GLAT(:),GLON(:)

    INTEGER,ALLOCATABLE::NDB(:,:),NPB(:,:)
    REAL, ALLOCATABLE::HSB(:,:),PSB(:,:),AGTK(:)

    integer :: NDPSB(4)
    real ::  HS(4)
    real*4 :: dlo,dla
    integer*4 :: ntiles,ISTART,IEND,LL


!	PARAMETER (NSSB=04)
    PARAMETER (NSSB=04)
!      PARAMETER (NTSUBH=01)
    PARAMETER (NTSUBH=02)
    parameter (nsorc2=4)
    parameter (nssbsp=nssb*nssb+1)
    parameter (nssb2m=nssb*2-1)

    parameter (nssbp=nssb+1)
    parameter (nssbpp=nssbp+1)

    PARAMETER (PIHF=1.570796,DTR=0.01745329)
    parameter (ndbmin=1)
    parameter (nsbmin=1)
    parameter (dlo=1./120.,dla=1./120.)
    parameter (ntiles=18*36)


! arrays for subbox calculations

    REAL,ALLOCATABLE:: hssb(:,:,:),hssa(:)
    INTEGER,ALLOCATABLE:: ndsa(:),ndsb(:,:,:)

! arrays for effective z0

    REAL,ALLOCATABLE:: z0eff(:,:), z0eff2(:,:,:)
    REAL,ALLOCATABLE:: zha(:,:),zas(:,:)
    INTEGER,ALLOCATABLE:: isorc(:)

    real :: ao(4)
    real :: ha(4)
    real :: ds(4)
    integer :: msq(2,2)
! arrays for statistics
    integer :: nssdbf(nssbsp)

    INTEGER,ALLOCATABLE:: nsdb(:),nspb(:),nssdb(:)
    INTEGER,ALLOCATABLE:: kndir(:,:,:),kldir(:,:,:),kwdir(:,:,:)

    integer :: kndirx(4),kldirx(4),kwdirx(4)
    real :: andirx(4),aldirx(4),awdirx(4)

    INTEGER,ALLOCATABLE::kisrch(:,:,:),klsrch(:,:,:),kwsrch(:,:,:)
    integer :: kountm(nsorc2)
    integer :: kountl(nsorc2)
    integer :: kountw(nsorc2)
    real :: zemean(4,nsorc2),zlmean(4,nsorc2),zwmean(4,nsorc2)
    real :: hamean(4,nsorc2),hlmean(4,nsorc2),hwmean(4,nsorc2)
    real :: asmean(4,nsorc2),almean(4,nsorc2),awmean(4,nsorc2)
    real :: zemin(4,nsorc2),zemax(4,nsorc2)
    real :: hamin(4,nsorc2),hamax(4,nsorc2)
    real :: asmin(4,nsorc2),asmax(4,nsorc2)
!      integer izemin(4,nsorc2),izemax(4,nsorc2)
    integer :: ihamin(4,nsorc2),ihamax(4,nsorc2)
    integer :: iasmin(4,nsorc2),iasmax(4,nsorc2)
    integer :: jzemin(4,nsorc2),jzemax(4,nsorc2)
    integer :: jhamin(4,nsorc2),jhamax(4,nsorc2)
!      integer jasmin(4,nsorc2),jasmax(4,nsorc2)
    integer :: kzemin(4,nsorc2),kzemax(4,nsorc2)
    integer :: khamin(4,nsorc2),khamax(4,nsorc2)
    integer :: khl2(jm),khh2(jm)
    integer :: kasmin(4,nsorc2),kasmax(4,nsorc2)
! arrays for exteme values
    real :: exvl(nssbpp)
    real :: exvlc(nssb)
! arrays for accumulating obstacle heights
    integer :: nadir(4)
    real :: htadir(4)
    real :: hhadir(4)

    REAL,ALLOCATABLE:: hgtk(:),dum2d(:,:),tmp(:,:),smk(:)
    REAL,ALLOCATABLE:: hgt(:,:),sm(:,:),hbm2(:)

    REAL,ALLOCATABLE:: tfn(:),tfx(:),tln(:),tlx(:)

    real :: &
    almx,almn,apmx,apmn &
    ,awlo,aelo,asla,anla &
    ,alt0d,almd,aphd,p,q &
    ,ald1,ald2,ald3,ald4,ald5,ald6,ald7,ald8 &
    ,apd1,apd2,apd3,apd4,apd5,apd6,apd7,apd8


    INTEGER*2,ALLOCATABLE:: iht(:,:)

    integer*4 :: i,j,k,kk,l,n,lat,lon,khl,khh &
    ,nlon,nlat,irecl &
    ,nwd,ned,nsd,nnd &
    ,nlo1,nlo2,nla1,nla2 &
    ,nrrd,ilast,ihw(jm),ihe(jm),iter
    character(255) :: fname
    character(8) ::   ctile(ntiles)
    character(3) ::   alon
    character(2) ::   alat
    character(1) ::   lndir,ltdir,akm

!	NAMELIST /SEARESO/ seares
    data hwelp/99999./
    DATA msq/4,1,2,3/

    logical :: last,inside


! *** Create topo for eta grid.
!     Uses 30 s global data obtained from ftp site edcftp.cr.usgs.gov

! *** Code obtained from U. of Athens and modified at FSL.

! ______________________________________________________________________________

! *** Fill eta model common blocks.

    call eta_commons

    write(6,*) 'back from eta_commons' , seares
! *** Create eta topo file.

!===============================================================================

! *** Create topo tile filenames and lat, lon extents.
!     At FSL topo tiles are 10x10 degrees with a file name indicating
!     the SW corner of the tile (e.g. U40N110W).


! ew (formerly parameter statements in oper)


! p
! ph0=cos(tph0d*dtr)
    stph0=sin(tph0d*dtr)
    KNE=IM
    KNW=IM-1
    KSE=1-IM
    KSW=-IM
    KHL00=1
    KHH00=IMJM
    KHL22=2*IM+1
    KHH22=IMJM-2*IM
    NINC=2*IM-1

! ew

    ALLOCATE(tfn(ntiles),tfx(ntiles),tln(ntiles),tlx(ntiles))
    ALLOCATE(IHT(1200,1200))
    n=0
    do lat=-90,80,10
        do lon=-180,170,10
            n=n+1
            if (lon < 0) then
                lndir='W'
            else
                lndir='E'
            endif
            if (lat < 0) then
                ltdir='S'
            else
                ltdir='N'
            endif
            write(alon,'(i3.3)') abs(lon)
            write(alat,'(i2.2)') abs(lat)
            ctile(n)='U'//alat//ltdir//alon//lndir
            tfn(n)=float(lat)
            tfx(n)=float(lat+10)
            tln(n)=float(lon)
            tlx(n)=float(lon+10)
        enddo
    enddo

! *** Find lat-lon extremes for eta grid.

    write(6,*) 'wbd, sbd ', wbd,sbd


    call corners(im,jm,imjm,tph0d,tlm0d,dlmd,dphd,apmn,almn,apmx,almx)


    print *,'ETA grid window:'
    print *,'   Northwest =',apmx,almn
    print *,'   Southeast =',apmn,almx
    print *,' '

!	stop

! *** Boundary point heights are set to zero (where hbm2=0).

    do j=1,jm
        khl2(j)=im*(j-1)-(j-1)/2+2
        khh2(j)=im*j-j/2-1
    enddo


    ALLOCATE(HBM2(IMJM))
    do k=1,imjm
        hbm2(k)=0.
    enddo

    do j=3,jm2
        khl=khl2(j)
        khh=khh2(j)
    !	write(6,*) 'khl, khh ', khl,khh
        do k=khl,khh
            hbm2(k)=1.
        enddo
    enddo

! *** Initialize height (hsb), ocean pts (psb), and total pts (ndb) within
!        each eta grid diamond to zero.


    ALLOCATE(NDB(IMJM,4))
    ALLOCATE(HSB(IMJM,4))
    ALLOCATE(PSB(IMJM,4))

    do m=1,4
        do k=1,imjm
            hsb(k,m)=0.
            psb(k,m)=0.
            ndb(k,m)=0
        enddo
    enddo

! *** Process tiles that cover eta grid domain.

    ALLOCATE( HSSB(nssb,nssb,imjm),NDSB(nssb,nssb,imjm) )

    HSSB=0.
    NDSB=0

    do n=1,ntiles
        check1=0
        check2=0
    ! eck3=0


        if (apmn <= tfx(n) .AND. apmx >= tfn(n)) then


            IF (almn <= tlx(n) .AND. almx >= tln(n)) check1=1

        ! p	West extending into East
            if (almn < -180 .AND. tln(n) > 0 .AND. &
            (almn+360) <= tlx(n)) check2=1

        ! p	East extending into West
            if (almx > 180 .AND. tlx(n) < 0 .AND. &
            almn >= tln(n) .AND. (almx-360) >= tln(n) ) check3=1

        endif


    !	write(6,*) 'lon bounds: this tile: ', almn,almx,tln(n),tlx(n)

        if (check1 == 1 .OR. check2 == 1 .OR. check3 == 1) then

        ! p	write(6,*) 'inside for ', tlx(n),tfx(n),check1,check2,check3
        
            fname=topo_in(1:index(topo_in,' ')-1)//ctile(n)
        
        ! ********* Check that tile parameters match expected 10x10 degree size.
        
            awlo=tln(n)
            aelo=tlx(n)
            asla=tfn(n)
            anla=tfx(n)
            nlon=int((aelo-awlo)/dlo)+1
            nlat=int((anla-asla)/dla)+1
            write(6,*)"n:",n,"nlat:",nlat,"nlon:",nlon
        
        !	conditional compilation for word length on DEC being different

        ! ifdef DEC
        !        irecl=nlon*0.5
        ! else
            irecl=nlon*2
        ! endif

        ! st            nwd=(awlo-almn)/dlo
            if (almn < -180 .AND. awlo > 0) then
            !	write(6,*) 'here (1)'
                nwd=(awlo-(almn+360.))/dlo
            elseif (almn > 0 .AND. awlo < 0) then
                nwd=(awlo-(almn-360.))/dlo
            else
                nwd=(awlo-almn)/dlo
            endif

            if (aelo > 0 .AND. almx < 0)  then
            !	write(6,*) 'here (2)'
                ned=(aelo-(almx+360.))/dlo
            elseif (aelo < 0 .AND. almx > 180) then
            !	write(6,*) 'here (2a)'
                ned=(aelo-(almx-360.))/dlo
            !	write(6,*) 'aelo, almx-360,ned ', aelo, almx-360, ned
            else
                ned=(aelo-almx)/dlo
            endif

            nsd=(asla-apmn)/dla
            nnd=(anla-apmx)/dla
        
            if (nwd <= 0) then
                nlo1=-nwd-4
            else
                nlo1=1
            endif
        
            if (ned >= 0) then
                nlo2=nlon-ned+4
            else
                nlo2=nlon
            endif
        
            if (nnd >= 0) then
                nla1=nnd-4
            else
                nla1=1
            endif
        
            if (nsd <= 0) then
                nla2=nlat+nsd+4
            else
                nla2=nlat
            endif
        
            if (nlat /= 1200 .OR. nlon /= 1200) then
                print *,'Unexpected tile size:',nlon,nlat
                print *,'Should be 1200 x 1200'
                print *,'Abort...'
                stop
            endif
        !           print*,' ','(',nlo1,'-',nlo2,'),(',nla1,'-',nla2,')'
        
        ! ********* Open input topo file.
        
            open(unit=1,file=fname,status='old' &
            ,form='unformatted',access='direct',recl=irecl &
            ,err=900)
            l=index(fname//' ',' ')-1
        ! p            print *,'reading topo data ',fname(1:l)
            print *,'reading topo data ',fname(l-8:l)
        ! p	write(6,*) ' '
        
        ! ********* Process row by row.
        
        ! o            last=.false.
        ! o            nrrd=nla1

        ! CCCCCCCCCCCCCCCCCCCCCCCCCCCC

        ! p      alt0d=anla-dla/2.-(nrrd-1)*dla
            alt0d=tfx(n)
            assb=float(nssb)
            rssb=1.0/assb
            r30s=1.0/120.0
            r30sh=r30s*0.5
            RTSUBH=R30sh/(NTSUBh*1.0)

        ! ew            do while (.not. last)

            do 710 nrrd=1,nlat
                read(1,rec=nrrd)(iht(i,nrrd),i=1,nlon)

            ! p     swap the topo data element by element

                do I=1,NLON
                    call swapval(iht(I,nrrd))
                enddo

            !	if (nrrd.eq.(nlat/2)) write(6,*)nrrd, fname,(iht(I,NRRD),I=20,30)

            710 END DO
            close(1)

        ! p	at this point the 1200 X 1200 block of data has been read

            do 7109 nrrd=1,nlat
                tlat=(ALT0D-(2*nrrd-1)*r30sh)
                if (tlat > apmx) go to 7109
                if (tlat < apmn) go to 7109
                do 7108 i=nlo1,nlo2
                    tlon=(awlo*1.0+(2*i-1)*(1./240.))
                !      if ((tlon.gt.blonem).and.(tlon.lt.blonw ).or.
                !     1    (tlon.gt.blone ).and.(tlon.lt.blonwp).or.
                !     2    (tlon.gt.blonep).and.(tlon.lt.blonwq)) go to 7108
                ! ew
                    DO  jtsubh=1,ntsubh
                        tslat=tlat+(jtsubh*2-ntsubh)*rtsubh
                        ntsubi=ntsubh-ifix(ntsubh*(1.0-cos(dtr*tslat)))
                        rtsubi=r30sh/(ntsubi*1.0)

                        DO itsubi=1,ntsubi
                            ALM=(tlon+(ITSUBi*2-ntsubi)*RTSUBi)
                            APH=tslat
                            nhread=nhread+1
                            ALM=ALM*DTR
                            APH=APH*DTR
                            call PQK_NEW(alm,aph,inside,p,q,kk)
                             
                            IF ( .NOT. INSIDE) cycle

                        !  assign the subboxes based on the values of p and q
                            issb = 1 + ifix(p*assb)
                            jssb = 1 + ifix(q*assb)
                            iqsb = 1 + ifix(p*2.0)
                            jqsb = 1 + ifix(q*2.0)
                        
                            kqsb=msq(iqsb,jqsb)

                            if (NDB(kk,kqsb) >= 0) then
                                HSB(k,kqsb)=HSB(kk,kqsb)+iht(i,nrrd)
                                NDB(kk,kqsb)=NDB(kk,kqsb)+1
                            endif

                            if (ndsb(issb,jssb,kk) >= 0) then
                                hssb(issb,jssb,k)=hssb(issb,jssb,k)+iht(i,nrrd)
                                ndsb(issb,jssb,k)=ndsb(issb,jssb,k)+1
                                if ((issb < 1) .OR. (issb > nssb) .OR. &
                                (jssb < 1) .OR. (jssb > nssb) .OR. &
                                (kk   < 1) .OR. (kk   > imjm)) then
                                    write(6,'('' data written outside the box'',4i6,2f7.3)') &
                                    issb,jssb,kk   ,kk,p,q
                                endif
                            endif
                            nhused=nhused+1
                        
                        enddo
                    enddo
                7108 END DO
            7109 END DO

        ! p    endif on whole section?
        ENDIF
    ! p	enddo on number of tiles?
    enddo


    DEALLOCATE(tfx,tfn,tlx,tln)
    DEALLOCATE(IHT)

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    write(6,*) 'done reading data?'
    write(6,*) 'nssbsp= ', nssbsp
    write(6,*) 'imjm= ', imjm
! zero out counters for filling of grid boxes

    ALLOCATE(NSDB(102))
    ALLOCATE(NSPB(102))
    ALLOCATE(NSSDB(102))
    do 4115 k=1,102
        nsdb(k)=0
        nspb(k)=0
        nssdb(k)=0
    4115 END DO
    do 4117 k=1,nssbsp
        nssdbf(k)=0
    4117 END DO

! convert positive numbers of points, fill counters


    ALLOCATE(ISORC(IMJM))
    ISORC=0
    ALLOCATE(NPB(IMJM,4))
    do 5177 k=1,imjm
        do 5163 m=1,4
            if (ndb(k,m) > 0) then
                if (mod(isorc(k),2) == 0) then
                    isorc(k)=isorc(k)+1
                endif
            endif
        ! est
        !          if (ndb(k,m).ge.ndbmin) then
        !            ndb(k,m)=-ndb(k,m)
        !          endif
        ! est
            ksdb=min(1+iabs(ndb(k,m)),102)
            kspb=min(1+iabs(npb(k,m)),102)
            nsdb(ksdb)=nsdb(ksdb)+1
            nspb(kspb)=nspb(kspb)+1
        5163 END DO
        kssdbf=1
        do 5167 j=1,nssb
            do 5165 i=1,nssb

            ! est?
            !            if (ndsb(i,j,k).ge.nsbmin) then
            !              ndsb(i,j,k)=-ndsb(i,j,k)
            !            endif
            ! est?
                kssdb=min(1+iabs(ndsb(i,j,k)),102)
                nssdb(kssdb)=nssdb(kssdb)+1
                if (ndsb(i,j,k) < 0) kssdbf=kssdbf+1
            5165 END DO
        5167 END DO
        isorc(k)=isorc(k)*2
        nssdbf(kssdbf)=nssdbf(kssdbf)+1
    5177 END DO

! report the numbers of boxes and subboxes with each number of points

    write(6,'('' '')')

    mssdbf=0
    do 5403 k=1,nssbsp
        km=k-1
        write(6,'('' boxes with'',i8,'' subboxes filled:'',i10)') &
        km,nssdbf(k)
        mssdbf=mssdbf+nssdbf(k)
    5403 END DO
    write  (6,'(''           '',8x,''     total boxes:'',i10)') &
    mssdbf

    write(6,'('' '')')

    msdb=0
    mspb=0
    mssdb=0
    do 5405 k=1,102
        km=k-1
        if (km == 101) then
            akm='>'
            km=km-1
        else
            akm=' '
        endif
        write(6,'('' topo,land,subboxes with '',a1,i3, &
        '' points:'',3i8)') &
        akm,km,nsdb(k),nspb(k),nssdb(k)
        msdb=msdb+nsdb(k)
        mspb=mspb+nspb(k)
        mssdb=mssdb+nssdb(k)
    5405 END DO
    write  (6,'('' topo,land,subboxes      '',1x,3x, &
    '' totals:'',3i8)') &
    msdb,mspb,mssdb

    write(6,'('' '')')


! *** All rows of topo data have been processed.
!     Now calculate the heights (hgtk) and sea mask (smk).

    ALLOCATE(GLAT(IMJM))
    ALLOCATE(GLON(IMJM))

!	dtr=acos(-1.)/180.

    tph0=tph0d*dtr
    wb=wbd*dtr
    sb=sbd*dtr
    dlm=dlmd*dtr
    dph=dphd*dtr
    tph=sb-dph
    tdlm=dlm+dlm

!	write(6,*) 'original values '
!	write(6,*) tph0,wb,sb,dlm,dph,tph,dtr

    do j=1,jm
        khl=im*(j-1)-(j-1)/2+1
        khh=im*j-j/2

    
        tlm=wb-tdlm+mod(j+1,2)*dlm
        tph=tph+dph
        stph=sin(tph)
        ctph=cos(tph)
    
        do k=khl,khh
            tlm=tlm+tdlm
            sinphi=ctph0*stph+stph0*ctph*cos(tlm)
            glat(k)=asin(sinphi)
            coslam=ctph*cos(tlm)/(cos(glat(k))*ctph0) &
            -tan(glat(k))*tan(tph0)
            coslam=min(coslam,1.)
            fact=1.
            if (tlm > 0.0) fact=-1.
            glon(k)=-tlm0d*dtr+fact*acos(coslam)
            glat(k)=glat(k)/dtr
            glon(k)=-glon(k)/dtr
        !	write(6,*) glat(k),glon(k)
        enddo
    enddo

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    do  M=1,4
        do  k=1,imjm

            if (NDB(K,M) <= 0) then
                write(6,*) 'K,M,NDB(K,M)= ', K,M,NDB(K,M)
            endif
        !     IF (NDB(K,M).GT.0) GO TO 131
        !       NO DATA IN THE SUB-BOX: ASSUME PERCENTAGE OF WATER
        !       SURFACE OF THE SUB-BOX TO BE 50 PERCENT, AND LEAVE
        !       ELEVATION ZERO:
            if (npb(k,m) <= 0) then
                PSB(K,M)=50.
            else
                PSB(K,M)=PSB(K,M)/npb(K,M)
            endif

            IF (NDB(K,M) <= 0) THEN
                HSB(K,M)=0.0
            ELSE
                HSB(K,M)=HSB(K,M)/NDB(K,M)
            ! p	write(6,*) 'ave hgt= ', HSB(K,M),NDB(K,M)
            ENDIF

        enddo
    enddo

    DEALLOCATE(NPB)

    ALLOCATE(HGTK(IMJM))
    ALLOCATE(SMK(IMJM))

! EW

    HGTK(1:IMJM)=0.


    if ( .NOT. sigma) then
        kstart=KHL22
        kend=KHH22
    endif

    if (sigma) then
        kstart=1
        kend=IMJM
    endif

    write(6,*) 'looping on 109 from  ', kstart, ' to ', kend


!      DO 109 K=KHL22,KHH22
    DO 109 K=kstart,kend

    
        PSM=0.
        DO 111 M=1,4
            PSM=PSM+PSB(K,M)
        111 END DO
        PAV=PSM*0.25

        IF (PAV <= 50.) GO TO 121

    !       (IF NO DATA WAS AVAILABLE WITHIN ANY OF THE FOUR SUB-BOXES,
    !       OR THE PERCENTAGE OF WATER SURFACE IS EXACTLY 50, THE POINT
    !       IS DECLARED AS LAND)
    
    !       THE POINT IS DEFINED AS SEA/LAKE:
    !***  CALCULATE THE ELEVATION BY ASSIGNING TO THE FOUR SUB-BOX SQUARE
    !***  THE ELEVATION OF THE SUB-BOX WHICH HAD THE MAXIMUM PERCENTAGE OF
    !***  WATER COVERAGE, PSB(K,M). (NOTE: "LARGE LAKES OR INLAND SEAS
    !***  WILL NOT BE CODED AS 100" ON THE NAVY TAPE !!)
    !***
        PSBM=50.
        DO 411 M=1,4
            IF(NDB(K,M) == 0 .OR. PSB(K,M) <= PSBM) GO TO 411
            if ( .NOT. sigma) then
                HGTK(K)=HBM2(K)*HSB(K,M)
            else
                HGTK(K)=HSB(K,M)
            endif
            PSBM   =        PSB(K,M)
        411 END DO
        IF (HGTK(K) <= HWELP) GO TO 109
        WRITE(6,412) K,PSBM,HGTK(K)
        412 FORMAT('  K=',I5,'  PSBM=',F5.1,'  HGTK(K)=',F5.0)
    !       WATER POINT IS HIGHER THAN HWELP.  DECLARE IT LAND:
        SMK(K)=0.
        GO TO 109
    
    !       THE POINT IS DEFINED AS LAND:
    
        121 SMK(K)=0.
    !	write(6,*) 'ever here?????'
    
        HS   (1)=AMAX1(HSB(K,4),HSB(K,1))
        NDPSB(1)=      NDB(K,4)+NDB(K,1)
        HS   (2)=AMAX1(HSB(K,3),HSB(K,2))
        NDPSB(2)=      NDB(K,3)+NDB(K,2)
        HS   (3)=AMAX1(HSB(K,4),HSB(K,3))
        NDPSB(3)=      NDB(K,4)+NDB(K,3)
        HS   (4)=AMAX1(HSB(K,1),HSB(K,2))
        NDPSB(4)=      NDB(K,1)+NDB(K,2)
    
    !     NOW CALCULATE THE NUMBER OF PAIRS OF SUB-BOXES WHICH HAVE
    !     NAVY MOUNTAIN DATA IN AT LEAST ONE OF THE TWO SUB-BOXES
    
    !       HOW MANY PAIRS OF SUB-BOXES ARE THERE WITH HEIGHT DATA?

        NPHD=4
        DO 151 N=1,4
            IF(NDPSB(N) == 0) NPHD=NPHD-1
        151 END DO
        IF(NPHD > 0) GO TO 112

    
    !       THERE WAS NO DATA WITHIN THE GRID BOX AND ITS ELEVATION IS
    !       LEFT ZERO.  DON'T WORRY ABOUT IT: IF IT'S A LAND POINT ITS
    !       ELEVATION WILL BE RAISED IN 'PTETA' TO MATCH THAT UNDER
    !       THE LOWEST NEIGHBORING NON-ZERO WIND POINT.  HOWEVER:
    !       SHOULD IT BE A WATER POINT AFTER ALL?  DECLARE IT A WATER
    !       POINT IF THE AVERAGE OF THE FOUR SUB-BOXES NEIGHBORING TO
    !       K AND SOUTH OF IT IS MORE THAN 50 PERCENT WATER:

        PAV=0.25*(PSB(K+KSE,3)+PSB(K+KSE,2) &
        +PSB(K+KSW,1)+PSB(K+KSW,2))
    ! p      IF(PAV.GT.50.) SMK(K)=1.
        GO TO 109
    
        112 SHS=0.
        DO 110 N=1,4
            SHS=SHS+HS(N)
        110 END DO

        if ( .NOT. sigma) then
            HGTK(K)=HBM2(K)*SHS/NPHD
        else
            HGTK(K)=SHS/NPHD
        endif
    
    109 END DO

    DEALLOCATE(PSB)

!****SMK**smk************************

!	2700,1350 for 8 minute
!	5400,2700 for 4 minute

!	Use the eta_sea_us call if domain is within box extending from
!	10N,180W to 80N,30W.  This uses 2 minute data



!	PAUSE

    if     (seares == 2) then
    ! LL get_sea(smk,4501,2101)
    elseif (seares == 4) then
    ! LL get_sea(SMK,5400,2700)
    elseif (seares == 8) then
    ! LL get_sea(SMK,2700,1350)
    elseif (seares == 30) then
     get_sea30(SMK)
    else
        write(6,*) 'BAD VALUE FOR SEARES!!!!!! ' , seares
        STOP
    endif

!	PAUSE


!*************************************

!     FOR AESTHETIC REASONS, MAKE SURE POINTS ALONG THE BUFFER LINES
!     REMAIN DECLARED WATER
    if ( .NOT. sigma) then
        DO 431 K=KHL22,KHH22
            SMK(K)=1.+HBM2(K)*(SMK(K)-1.)
        431 END DO
    endif

    ALLOCATE(DUM2D(IMT,JMT))

    call conh12t(smk,imjm,1,dum2d,imt,jmt)
    write(6,*) 'sea mask  '
    do J=JMT,1,-2
    !        write(6,273)(dum2d(I,J),I=1,IMT,4)
    enddo

    DEALLOCATE(DUM2D)

!-----------------------------------------------------------
!***  CALCULATE AND SAVE AVERAGE ELEVATIONS
!***

    ALLOCATE(AGTK(IMJM))

    DO 301 K=KHL00,KHH00
        AGTK(K)=0.
    301 END DO
    DO 302 K=KHL00,KHH00
        NSBHD=4
        DO 303 N=1,4
            IF(NDB(K,N) == 0)NSBHD=NSBHD-1
        303 END DO
    !	if (NSBHD .ne. 0) write(6,*) 'NSBHD= ', NSBHD
        IF(NSBHD == 0)GO TO 302
        AGTK(K)=(HSB(K,1)+HSB(K,2)+HSB(K,3)+HSB(K,4))/ &
        NSBHD
    !	write(6,*) 'AGTK(K)= ', AGTK(K)
    302 END DO

    DEALLOCATE(HSB)
    DEALLOCATE(NDB)
!***
!***  THE silhouette/mean MOUNTAINS:  AT LAND POINTS WITH CONCAVE
!***  ACTUAL (AVERAGE) TOPOGRAPHY, REPLACE THE
!***  SILHOUETTE ELEVATION WITH THE AVERAGE ELEVATION
!***  (22 FEB 94)
    DO 311 K=KHL22,KHH22
        IF(SMK(K) > 0.5)GO TO 311
        ALAH=AGTK(K+1)+AGTK(K+NINC)+AGTK(K-1) &
        +AGTK(K-NINC)+2.*(AGTK(K+KNE) &
        +AGTK(K+KNW)+AGTK(K+KSW)+AGTK(K+KSE)) &
        -12.*AGTK(K)
        IF(ALAH >= 0.)HGTK(K)=AGTK(K)
    311 END DO

!     If K is a point seen as a "valley" when looking at three-point
!     averages in any one of the four directions, and it is at the
!     same time not higher than each of its four nearest neighbors,
!     and also not higher than each of its four second-nearest
!     neighbors, choose always the mean elevation (23 June 97)

!     (Saddle points are expected to be frequently declared mean as
!     a result)

    DO 316 K=KHL22,KHH22
        IF(SMK(K) > 0.5)GO TO 316
        if (agtk(k) > agtk(k+kne) .AND. agtk(k) > agtk(k+ knw) .AND. &
        agtk(k) > agtk(k+ksw) .AND. agtk(k) > agtk(k+ kse)) &
        go to 316
        if (agtk(k) > agtk(k+  1) .AND. agtk(k) > agtk(k+ninc) .AND. &
        agtk(k) > agtk(k-  1) .AND. agtk(k) > agtk(k-ninc)) &
        go to 316
        sumtx =agtk(k-   1)+agtk(k)+agtk(k+   1)
        sumtxp=agtk(k+ ksw)+agtk(k)+agtk(k+ kne)
        sumty =agtk(k-ninc)+agtk(k)+agtk(k+ninc)
        sumtyp=agtk(k+ kse)+agtk(k)+agtk(k+ knw)
        if (sumtx < agtk(k+ ksw)+agtk(k-ninc)+agtk(k+ kse) .AND. &
        sumtx < agtk(k+ knw)+agtk(k+ninc)+agtk(k+ kne)) &
        hgtk(k)=agtk(k)
        if (sumtxp < agtk(k-ninc)+agtk(k+ kse)+agtk(k+   1) .AND. &
        sumtxp < agtk(k-   1)+agtk(k+ knw)+agtk(k+ninc)) &
        hgtk(k)=agtk(k)
        if (sumty < agtk(k+ kse)+agtk(k+   1)+agtk(k+ kne) .AND. &
        sumty < agtk(k+ ksw)+agtk(k-   1)+agtk(k+ knw)) &
        hgtk(k)=agtk(k)
        if (sumtyp < agtk(k+   1)+agtk(k+ kne)+agtk(k+ninc) .AND. &
        sumtyp < agtk(k-ninc)+agtk(k+ ksw)+agtk(k-   1)) &
        hgtk(k)=agtk(k)
    316 END DO

!    the silhouette/mean topography is done
    write(6,'('' silhouette/mean topography is done'')')


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCC       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCC Z0EFF CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCC       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!    process the high-resolution arrays for roughness data

! zero out the arrays for the grid box averages

    ALLOCATE(HSSA(IMJM),NDSA(IMJM))

    do 4203 k=1,imjm
        hssa(k)=0.0
        ndsa(k)=0
    4203 END DO

! invert the negative numbers and find grid-box average height

    do 4303 k=1,imjm
        do 4302 jssb=1,nssb
            do 4301 issb=1,nssb
                if (ndsb(issb,jssb,k) < 0) then
                    ndsb(issb,jssb,k)=-ndsb(issb,jssb,k)
                endif
                if (ndsb(issb,jssb,k) > 0) then
                    hssa(k)=hssa(k)+hssb(issb,jssb,k)
                    ndsa(k)=ndsa(k)+ndsb(issb,jssb,k)
                    hssb(issb,jssb,k)=hssb(issb,jssb,k)/ &
                    ndsb(issb,jssb,k)
                !	write(6,*) 'found an hssb value! ',
                !     + issb,jssb,hssb(issb,jssb,k)
                endif
            4301 END DO
        4302 END DO
    4303 END DO
! fill in the blanks with the average over the whole box

    do 4343 k=1,imjm
        ifill=0
        if (ndsa(k) == 0) then
            write(6,'('' We have problems, no heights in k'',i8,2f6.1)') &
            k,glat(k),glon(k)
            ndsa(k)=1
            hssa(k)=0
        else
            hssa(k)=hssa(k)/ndsa(k)
        endif
        do 4342 jssb=1,nssb
            do 4341 issb=1,nssb
                if (ndsb(issb,jssb,k) == 0.0) then
                    hssb(issb,jssb,k)=hssa(k)
                    ndsb(issb,jssb,k)=1
                    ifill=ifill+1
                endif
            4341 END DO
        4342 END DO
        if (ifill > 0) then
        !          write(6,'('' '',i6,'' filled points for k='',i8)')
        !     1      ifill,k
        endif
    4343 END DO

    DEALLOCATE(GLAT,GLON)

    write(6,'('' empty subboxes have been filled'')')

    do 4439 k=1,imjm
    !         IMT=2*IM-1
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX <= IM)THEN
            I=2*IX-1
            J=2*JX-1
        ELSE
            I=2*(IX-IM)
            J=2*JX
        ENDIF
    !       write(6,'('' '',3i6,''=k,i,j'',f20.10,i10,''=h,n'')')
    !    1    k,i,j,hssa(k),ndsa(k)
    4439 END DO

    DEALLOCATE(HSSA,NDSA)

!  constants:

    a=6371229.0
    z0 = 0.15
!     cd = 0.3
!     cd = 0.8
    cd = 0.6
    vk = 0.4
    rvk2 = 1.0 / ( vk ** 2 )
    rz0 = 1.0 / z0
    dmin=1.0e-6
    hmin=2.0*z0+dmin
    wtamin=dmin

    write(6,'('' constants defined'')')

! zero out statistics


    ALLOCATE(KNDIR(NSSBSP,4,NSORC2))
    ALLOCATE(KLDIR(NSSBSP,4,NSORC2))
    ALLOCATE(KWDIR(NSSBSP,4,NSORC2))

    ALLOCATE(kisrch(3,4,nsorc2),klsrch(3,4,nsorc2))
    ALLOCATE(kwsrch(3,4,nsorc2))

    do 4494 l=1,nsorc2
        do 4493 m=1,4
            hamean(m,l)=0.0
            hlmean(m,l)=0.0
            hwmean(m,l)=0.0
            asmean(m,l)=0.0
            almean(m,l)=0.0
            awmean(m,l)=0.0
            zemean(m,l)=0.0
            zlmean(m,l)=0.0
            zwmean(m,l)=0.0
            hamin(m,l)= 1.0e30
            hamax(m,l)=-1.0e30
            asmin(m,l)= 1.0e30
            asmax(m,l)=-1.0e30
            zemin(m,l)= 1.0e30
            zemax(m,l)=-1.0e30
            ihamax(m,l)=0
            ihamin(m,l)=0
            iasmax(m,l)=0
            iasmin(m,l)=0
        !          izemax(m,l)=0
        !          izemin(m,l)=0
            jhamax(m,l)=0
            jhamin(m,l)=0
        !          jasmax(m,l)=0
        !          jasmin(m,l)=0
            jzemax(m,l)=0
            jzemin(m,l)=0
            khamax(m,l)=0
            khamin(m,l)=0
            kasmax(m,l)=0
            kasmin(m,l)=0
            kzemax(m,l)=0
            kzemin(m,l)=0
            do 4491 k=1,nssbsp
                kndir(k,m,l)=0
                kldir(k,m,l)=0
                kwdir(k,m,l)=0
            4491 END DO
            do 4492 k=1,3
                kisrch(k,m,l)=0
                kwsrch(k,m,l)=0
                klsrch(k,m,l)=0
            4492 END DO
        4493 END DO
    4494 END DO

    write(6,'('' statistics zeroed out'')')

! calculate 4 z0's for each grid box


    write(6,*) 'memory here?'

!	PAUSE

    ALLOCATE (zha(imjm,4),zas(imjm,4))
    ALLOCATE (z0eff(imjm,4))

    do 4557 k=1,imjm
    !	write(6,*) 'here 1 ', k
    ! p
    !	goto 4557
    ! p
    
        l=isorc(k)+1
    
    !       find the local mesh lengths
    
    !----------------------------------------------------------
    !***
    !***  CONVERT FROM ONE-DIMENSIONAL K VALUES ON THE E-GRID
    !***  TO I,J VALUES ON THE FILLED E-GRID.
    !***
    !----------------------------------------------------------
    
    !       IMT=2*IM-1
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX <= IM)THEN
            I=2*IX-1
            J=2*JX-1
        ELSE
            I=2*(IX-IM)
            J=2*JX
        ENDIF

        phl = ( j-1 ) * dphd + sbd
    
        dlmdl = dlmd * dtr * cos(phl*dtr)
        dphdl = dphd * dtr
    
    !       find the length of the sides of the obstacles
    
        ds(1)=a*rssb*dphdl
        ds(3)=a*rssb*dlmdl
        ds(2)=sqrt(ds(1)**2+ds(3)**2)
        ds(4)=ds(2)
        area=dphdl*dlmdl*a*a*2.0
        rarea=1.0/area
    
        if ((i == j) .OR. ((im-i) == (jm-j))) then
        ! p            write(6,'('' '',i6,''=i'',i6,''=j'')') i,j
        ! p            write(6,'('' ds  ew nesw ns nwse'',4f12.3)')
        ! p     1        ds
        endif
    
    
    !       test whether this is a boundary point
    
        if (hbm2(k) == 0.0) then
        
        !         set boundary points to z0
        
            do 4007 m=1,4
                z0eff(k,m)=z0
                zha(k,m)=0.0
                zas(k,m)=0.0
            4007 END DO
        !	write(6,*) 'setting boundary point'

        else

        
        !         calculate z0eff for interior points
        
            mthdz0=2
        
        !         mthdz0 1, find obstacles by searching for
        !         relative minima and maxima along lines in each direction
        !   Abandoned, as this was an older and an inferior option, having had
        !              no consolidation of obstacles.  FM, June 97
        
        !         mthdz0 2, find obstacles by searching for
        !         relative minima and maxima along lines in each direction
        !         and eliminating insignificant relative extrema
        
            if (mthdz0 == 2) then
            
            !         loop for directions
            
            !	write(6,*) 'hit loop 4287'

                do 4287 m=1,4
                    hhdir=0.0
                    htdir=0.0
                    ndir=0.0
                    imin=1
                    if (mod(m,2) == 1) then
                    !             cardinal direction, diagonal on grid
                        imax=nssb2m
                    else
                    !             diagonal direction, along grid
                        imax=nssb
                    endif
                
                !           loop to select a line of subboxes
                
                    do 4285 is=imin,imax
                    !	write(6,*) 'inside 4285 ',imin,is,imax
                        if (mod(m,2) == 1) then
                        !               cardinal direction, diagonal on grid
                            jmin=max(0,is-nssb)
                            jmax=min(nssbp,is+1)
                        else
                        !               diagonal direction, along grid
                            jmin=0
                            jmax=nssbp
                        endif
                        llsb=jmax-jmin+1
                        nexvlf=0
                    
                    !             loop to search for extrema along the line of subboxes
                    
                        do 4223 js=jmin,jmax
                        !	write(6,*) 'start 4223 ', jmin,js,jmax
                        !               select subbox indices by direction indicator
                            if (m == 1) then
                            !                 east to west
                                issb=1+is-js
                                jssb=js
                                kssb=k
                            endif
                            if (m == 2) then
                            !                 southwest to northeast
                                issb=js
                                jssb=is
                                kssb=k
                            endif
                            if (m == 3) then
                            !                 south to north
                                issb=js
                                jssb=nssb-is+js
                                kssb=k
                            endif
                            if (m == 4) then
                            !                 southeast to northwest
                                issb=is
                                jssb=js
                                kssb=k
                            endif
                        !               reindex to the interior of a box
                            if (issb < 1) then
                                issb=issb+nssb
                                kssb=kssb+ksw
                            endif
                            if (jssb < 1) then
                                jssb=jssb+nssb
                                kssb=kssb+kse
                            endif
                            if (issb > nssb) then
                                issb=issb-nssb
                                kssb=kssb+kne
                            endif
                            if (jssb > nssb) then
                                jssb=jssb-nssb
                                kssb=kssb+knw
                            endif
                        !               find minima and maxima
                            if (js == jmin) then
                            !                 first point, taken to be first minimum
                                hmin1=hssb(issb,jssb,kssb)
                                isrch=1
                            else
                            !                 any other point
                                if (isrch == 1) then
                                !                   looking for the first minimum
                                    if (hmin1 > hssb(issb,jssb,kssb)) then
                                    !                     lower the first minimum
                                        hmin1=hssb(issb,jssb,kssb)
                                    elseif (hmin1 < hssb(issb,jssb,kssb)) then
                                    !                     we have a minimum, now look for maximum
                                        nexvlf=nexvlf+1
                                        exvl(nexvlf)=hmin1
                                        hmax2=hssb(issb,jssb,kssb)
                                        isrch=2
                                    endif
                                elseif (isrch == 2) then
                                !                   looking for a maximum
                                    if (hmax2 < hssb(issb,jssb,kssb)) then
                                    !                     raise the maximum
                                        hmax2=hssb(issb,jssb,kssb)
                                    elseif (hmax2 > hssb(issb,jssb,kssb)) then
                                    !                     we have a maximum, now look for minimum
                                        hmin3=hssb(issb,jssb,kssb)
                                        isrch=3
                                    endif
                                elseif (isrch == 3) then
                                !                   looking for a minimum
                                    if (hmin3 > hssb(issb,jssb,kssb)) then
                                    !                     lower the minimum
                                        hmin3=hssb(issb,jssb,kssb)
                                    endif
                                    if ((hmin3 < hssb(issb,jssb,kssb)) .OR. &
                                    (js == jmax)) then
                                    !                     we either have another minimum or are done
                                    !                     add an obstacle to the totals
                                    !                     hodir=hmax2-0.5*(hmin1+hmin3)
                                    !                     htdir=htdir+hodir
                                    !                     hhdir=hhdir+hodir**2
                                    !                     ndir=ndir+1
                                    !                     add an obstacle to the array of extrema
                                        nexvlf=nexvlf+1
                                        exvl(nexvlf)=hmax2
                                        nexvlf=nexvlf+1
                                        exvl(nexvlf)=hmin3
                                    !                     move the minimum value to front of next mt
                                        hmin1=hmin3
                                    !                     now we are looking for another maximum
                                        hmax2=hssb(issb,jssb,kssb)
                                        isrch=2
                                    endif
                                endif
                            endif
                        !		write(6,*) 'end 4223 '
                        4223 END DO
                    !	write(6,*) 'BEYOND 4223'
                    !	write(6,*) 'isrch, m, l= ', isrch,m,l
                        kisrch(isrch,m,1)=kisrch(isrch,m,1)+1
                        kisrch(isrch,m,l)=kisrch(isrch,m,l)+1
                        if (smk(k) > 0.5) then
                            kwsrch(isrch,m,1)=kwsrch(isrch,m,1)+1
                            kwsrch(isrch,m,l)=kwsrch(isrch,m,l)+1
                        else
                            klsrch(isrch,m,1)=klsrch(isrch,m,1)+1
                            klsrch(isrch,m,l)=klsrch(isrch,m,l)+1
                        endif
                    !	write(6,*) 'BEYOND SRCH'
                    
                    !  consolidation of obstacles
                    
                    !	write(6,*) 'how many obstacles? ', nexvlf
                        if (nexvlf >= 5) then
                        
                        !  calculate z0eff of the line just processed,
                        !  assuming that the width of the obstacles is unity
                        
                        !  sl = length of the line along which obstacles are calculated
                        !  nlob = n of the peak of the last obstacle along the line
                        
                            if (m == 1) then
                                sl=llsb*ds(3)*2.0
                            elseif (m == 3) then
                                sl=llsb*ds(1)*2.0
                            else
                                sl=llsb*ds(2)
                            endif
                        
                            nlob=nexvlf-1
                            saob=0.0
                            shhob=0.0
                            do 4235 ilob=2,nlob,2
                                hi=exvl (ilob)-0.5*(exvl (ilob-1)+exvl (ilob+1))
                                saob=saob+hi
                                shhob=shhob+hi**2
                            4235 END DO
                            if (saob > hmin) then
                                hob=shhob/saob
                            else
                                hob=hmin
                                saob=hmin
                                shhob=hmin**2
                            endif
                        
                        !  calculate a tentative value of z0eff for this line
                        
                            if (0.5*rz0*hob < 0) write(6,*) 'neg alog', 0.5*rz0*hob
                            wtarea= &
                            &             0.5*cd*rvk2*saob/sl + &
                            &             1.0/((alog(0.5*rz0*hob))**2)

                            if (wtarea < wtamin) then
                                write(6,'('' nonpositive tentative wtarea'', &
                                f20.10,2i8)') &
                                wtarea,k,m
                                wtarea=wtamin
                            endif
                            zeffl =0.5*hob/exp(1.0/sqrt(wtarea))

                        ! p	write(6,*) 'in z0eff code, wtarea,zeffl = ', wtarea,zeffl
                        
                        !  of the one or more inside minima, find the highest
                        
                            4243 continue
                            hghmn=exvl(3)
                            nhghmn=3
                            nlim=3
                        
                        !  nlim = n of the last minimum
                        
                        !  nlim = 5 if there is only one inside minimum; it is the highest
                        
                            if(nexvlf > 5) then
                                nlim=nexvlf-2
                                do 4245 imn=5,nlim,2
                                    if (hghmn < exvl(imn)) then
                                        hghmn=exvl(imn)
                                        nhghmn=imn
                                    endif
                                4245 END DO
                            endif
                        
                        !  the highest inside minimum has been found.
                        !  does removing it and compressing exvl increase z0eff of the line?
                        
                            nhmm=nhghmn-1
                            do 4251 iex=1,nhmm
                                exvlc(iex)=exvl(iex)
                            4251 END DO
                            if (exvl(nhghmn+1) > exvl(nhmm)) then
                                exvlc(nhmm)=exvl(nhghmn+1)
                            endif
                        
                        !  the higher of the two maxima has been selected and stored
                        !  now transfer the remainder of exvl to exvlc
                        
                            do 4253 iex=nhghmn,nlim
                                exvlc(iex)=exvl(iex+2)
                            4253 END DO
                        
                        !  calculate z0eff of the compressed array, zefflc
                        
                            nlob=nlim-1
                            saob=0.0
                            shhob=0.0
                            do 4255 ilob=2,nlob,2
                                hi=exvlc(ilob)-0.5*(exvlc(ilob-1)+exvlc(ilob+1))
                                saob=saob+hi
                                shhob=shhob+hi**2
                            4255 END DO
                            if (saob > hmin) then
                                hob=shhob/saob
                            else
                                hob=hmin
                                saob=hmin
                                shhob=hmin**2
                            endif
                        
                        !  calculate a tentative value of z0eff for this line
                        
                            wtarea= &
                            &             0.5*cd*rvk2*saob/sl + &
                            &             1.0/((alog(0.5*rz0*hob))**2)
                            if (wtarea < wtamin) then
                                write(6,'('' nonpositive tentative wtarea'', &
                                f20.10,2i8)') &
                                wtarea,k,m
                                wtarea=wtamin
                            endif
                            zefflc=0.5*hob/exp(1.0/sqrt(wtarea))
                        
                            if (zeffl < zefflc) then
                            
                            !  switch to compressed exvl and either exit the consolidation code
                            !  or search for another minimum to remove
                            
                                do 4257 iex=nhmm,nlim
                                    exvl(iex)=exvlc(iex)
                                4257 END DO
                                nexvlf=nlim
                            
                                if (nexvlf < 5) go to 4263
                            
                                zeffl =zefflc
                                go to 4243
                            else
                                go to 4263
                            endif
                        
                        endif
                    
                        4263 continue
                    
                    !  There has been at most only one obstacle to start with,
                    !  or consolidation of obstacles is completed.
                    !  Add obstacles, if any, to the totals
                    
                        if (nexvlf >= 3) then
                            nlob=nexvlf-1
                            do 4267 ilob=2,nlob,2
                                hodir=exvl(ilob)-0.5*(exvl(ilob-1)+exvl(ilob+1))
                                if (hodir == 0) then
                                    write(6,*) 'increasing hodir!!! '
                                    hodir=hmin
                                endif
                                htdir=htdir+hodir
                                hhdir=hhdir+hodir**2
                                ndir=ndir+1
                            4267 END DO
                        endif
                    !	write(6,*) 'end 4285 ', is
                    4285 END DO
                !	write(6,*) 'past 4285, m= ', m
                !           in case we didn't find any obstacles
                    if (ndir == 0) then
                        ndir=1
                        htdir=hmin
                        hhdir=hmin**2
                    endif
                !           fill directional obstacle totals
                    nadir(m)=ndir
                    htadir(m)=htdir
                    hhadir(m)=hhdir
                    ndirp=ndir+1
                    kndir(ndirp,m,1)=kndir(ndirp,m,1)+1
                    kndir(ndirp,m,l)=kndir(ndirp,m,l)+1
                    if (smk(k) > 0.5) then
                        kwdir(ndirp,m,1)=kwdir(ndirp,m,1)+1
                        kwdir(ndirp,m,l)=kwdir(ndirp,m,l)+1
                    else
                        kldir(ndirp,m,1)=kldir(ndirp,m,1)+1
                        kldir(ndirp,m,l)=kldir(ndirp,m,l)+1
                    endif
                4287 END DO
            !	write(6,*) 'past 4287 ', K
            
            !         end of mthdz0 2
            
            else
            
            !         default mthdz0, all h and a are based on hmin
            
                do 4428 m=1,4
                    nadir(m)=1
                    htadir(m)=hmin
                    hhadir(m)=hmin**2
                4428 END DO
            
            !         end of default mthdz0
            
            endif
        
            if ((i == j) .OR. ((im-i) == (jm-j))) then
            ! p            write(6,'('' na  ew nesw ns nwse'',4i12)')
            !     1        nadir
            !            write(6,'('' hta ew nesw ns nwse'',4f12.4)')
            !     1        htadir
            !            write(6,'('' hha ew nesw ns nwse'',4f12.2)')
            ! p     1        hhadir
            endif
        
            do 4437 m=1,4
            
            !         find the total area of obstacles in each direction
            
                ao(m)=htadir(m)*ds(m)*rarea
            
            !         find the average heights of obstacles, use hmin to avoid zero
            !           the average is weighted by the height of the obstacle
            
                if (htadir(m) < 0) write(6,*) 'divide by zero!'
                ha(m)=hhadir(m)/htadir(m)
            
            4437 END DO
        
        
            if ((i == j) .OR. ((im-i) == (jm-j))) then
            ! p            write(6,'('' ao  ew nesw ns nwse'',4f12.5)')
            !     1        ao
            !            write(6,'('' ha  ew nesw ns nwse'',4f12.3)')
            ! p     1        ha
            endif
        
        !         calculate z0eff
        
            do 4552 m=1,4
                zha(k,m)=ha(m)
                zas(k,m)=ao(m)
            
            ! accumulate statistics on zha and zas
            
                if (hamin(m,1) > zha(k,m)) then
                    hamin(m,1)=zha(k,m)
                    ihamin(m,1)=i
                    jhamin(m,1)=j
                    khamin(m,1)=k
                endif
                if (hamin(m,l) > zha(k,m)) then
                    hamin(m,l)=zha(k,m)
                    ihamin(m,l)=i
                    jhamin(m,l)=j
                    khamin(m,l)=k
                endif
                if (hamax(m,1) < zha(k,m)) then
                    hamax(m,1)=zha(k,m)
                    ihamax(m,1)=i
                    jhamax(m,1)=j
                    khamax(m,1)=k
                endif
                if (hamax(m,l) < zha(k,m)) then
                    hamax(m,l)=zha(k,m)
                    ihamax(m,l)=i
                    jhamax(m,l)=j
                    khamax(m,l)=k
                endif
                hamean(m,1)=hamean(m,1)+zha(k,m)
                hamean(m,l)=hamean(m,l)+zha(k,m)
                if (smk(k) > 0.5) then
                    hwmean(m,1)=hwmean(m,1)+zha(k,m)
                    hwmean(m,l)=hwmean(m,l)+zha(k,m)
                else
                    hlmean(m,1)=hlmean(m,1)+zha(k,m)
                    hlmean(m,l)=hlmean(m,l)+zha(k,m)
                endif
            
                if (asmin(m,1) > zas(k,m)) then
                    asmin(m,1)=zas(k,m)
                    iasmin(m,1)=i
                !              jasmin(m,1)=j
                    kasmin(m,1)=k
                endif
                if (asmin(m,l) > zas(k,m)) then
                    asmin(m,l)=zas(k,m)
                    iasmin(m,l)=i
                !              jasmin(m,l)=j
                    kasmin(m,l)=k
                endif
                if (asmax(m,1) < zas(k,m)) then
                    asmax(m,1)=zas(k,m)
                    iasmax(m,1)=i
                !              jasmax(m,1)=j
                    kasmax(m,1)=k
                endif
                if (asmax(m,l) < zas(k,m)) then
                    asmax(m,l)=zas(k,m)
                    iasmax(m,l)=i
                !              jasmax(m,l)=j
                    kasmax(m,l)=k
                endif
                asmean(m,1)=asmean(m,1)+zas(k,m)
                asmean(m,l)=asmean(m,l)+zas(k,m)
                if (smk(k) > 0.5) then
                    awmean(m,1)=awmean(m,1)+zas(k,m)
                    awmean(m,l)=awmean(m,l)+zas(k,m)
                else
                    almean(m,1)=almean(m,1)+zas(k,m)
                    almean(m,l)=almean(m,l)+zas(k,m)
                endif
            
            ! test for out-of-range values
            
                if (zha(k,m) <= hmin) then
                !             write(6,'('' nonpositive zha(k,m,1)'',f20.10,2i8)')
                !    1          zha(k,m),k,m
                    zha(k,m)=hmin
                endif
            
                if (zas(k,m) < 0.0) then
                    write(6,'('' negative zha(k,m,1)'',f20.10,2i8)') &
                    zas(k,m),k,m
                    zas(k,m)=0.0
                endif
            
            ! calculate the effective z0
            ! from zas and zha
            
                wtarea= &
                &         0.5*cd*rvk2*zas(k,m)+ &
                &         1.0/((alog(0.5*rz0*zha(k,m)))**2)
                if (wtarea < wtamin) then
                    write(6,'('' nonpositive wtarea'',f20.10,2i8)') &
                    wtarea,k,m
                    wtarea=wtamin
                endif
                z0eff(k,m)=0.5*zha(k,m)/exp(1.0/sqrt(wtarea))
            !	if (z0eff(k,m).ne. 0.15) then
            !	write(6,*) 'zoeff= ', z0eff(k,m)
            !	endif
            
            ! accumulate statistics
            
                if (zemin(m,1) > z0eff(k,m)) then
                    zemin(m,1)=z0eff(k,m)
                !              izemin(m,1)=i
                    jzemin(m,1)=j
                    kzemin(m,1)=k
                endif
                if (zemin(m,l) > z0eff(k,m)) then
                    zemin(m,l)=z0eff(k,m)
                !              izemin(m,l)=i
                    jzemin(m,l)=j
                    kzemin(m,l)=k
                endif
                if (zemax(m,1) < z0eff(k,m)) then
                    zemax(m,1)=z0eff(k,m)
                !	izemax(m,1)=i
                    jzemax(m,1)=j
                    kzemax(m,1)=k
                endif
                if (zemax(m,l) < z0eff(k,m)) then
                    zemax(m,l)=z0eff(k,m)
                !              izemax(m,l)=i
                    jzemax(m,l)=j
                    kzemax(m,l)=k
                endif
                zemean(m,1)=zemean(m,1)+z0eff(k,m)
                zemean(m,l)=zemean(m,l)+z0eff(k,m)
                if (smk(k) > 0.5) then
                    zwmean(m,1)=zwmean(m,1)+z0eff(k,m)
                    zwmean(m,l)=zwmean(m,l)+z0eff(k,m)
                else
                    zlmean(m,1)=zlmean(m,1)+z0eff(k,m)
                    zlmean(m,l)=zlmean(m,l)+z0eff(k,m)
                endif

            
            4552 END DO
        
            kountm(1)=kountm(1)+1
            kountm(l)=kountm(l)+1
            if (smk(k) > 0.5) then
                kountw(1)=kountw(1)+1
                kountw(l)=kountw(l)+1
            else
                kountl(1)=kountl(1)+1
                kountl(l)=kountl(l)+1
            endif
        
        endif
    !	write(6,*) 'end 4557 ', K
    4557 END DO

    DEALLOCATE(ISORC)
    DEALLOCATE(HBM2)



    write(6,'('' z0eff calculation done'')')

! print statistics

    do 4876 l=nsorc2,1,-1
    
    ! report the numbers of boxes with each number of obstacles
        write(6,'('' '')')
        write(6,'('' statistics for topography source data'',i4)') &
        l
        write(6,'('' '')')
    
        kndirt=0
        kldirt=0
        kwdirt=0
        do 4561 m=1,4
            kndirx(m)=0
            kldirx(m)=0
            kwdirx(m)=0
            andirx(m)=0.0
            aldirx(m)=0.0
            awdirx(m)=0.0
        4561 END DO
        do 4563 k=1,nssbsp
            km=k-1
            kndirs=kndir(k,1,l)+kndir(k,2,l)+kndir(k,3,l)+kndir(k,4,l)
            kldirs=kldir(k,1,l)+kldir(k,2,l)+kldir(k,3,l)+kldir(k,4,l)
            kwdirs=kwdir(k,1,l)+kwdir(k,2,l)+kwdir(k,3,l)+kwdir(k,4,l)
            ktdirs=kndirs+kldirs+kwdirs
            if (ktdirs > 0) then
                write(6,'('' '',i3,''o'',4(1x,i6,i5,i6))') &
                km,(kndir(k,m,l),kldir(k,m,l),kwdir(k,m,l),m=1,4)
                do 4562 m=1,4
                    kndirx(m)=kndirx(m)+kndir(k,m,l)
                    kldirx(m)=kldirx(m)+kldir(k,m,l)
                    kwdirx(m)=kwdirx(m)+kwdir(k,m,l)
                    andirx(m)=andirx(m)+km*kndir(k,m,l)
                    aldirx(m)=aldirx(m)+km*kldir(k,m,l)
                    awdirx(m)=awdirx(m)+km*kwdir(k,m,l)
                4562 END DO
                kndirt=kndirt+kndirs
                kldirt=kldirt+kldirs
                kwdirt=kwdirt+kwdirs
            endif
        4563 END DO
        do 4564 m=1,4
            if (kndirx(m) > 0) then
                andirx(m)=andirx(m)/kndirx(m)
            endif
            if (kldirx(m) > 0) then
                aldirx(m)=aldirx(m)/kldirx(m)
            endif
            if (kwdirx(m) > 0) then
                awdirx(m)=awdirx(m)/kwdirx(m)
            endif
        4564 END DO

        do 4587 m=1,4
            if (kountm(l) > 0) then
                hamean(m,l)=hamean(m,l)/kountm(l)
                asmean(m,l)=asmean(m,l)/kountm(l)
                zemean(m,l)=zemean(m,l)/kountm(l)
            endif
            if (kountl(l) > 0) then
                hlmean(m,l)=hlmean(m,l)/kountl(l)
                almean(m,l)=almean(m,l)/kountl(l)
                zlmean(m,l)=zlmean(m,l)/kountl(l)
            endif
            if (kountw(l) > 0) then
                hwmean(m,l)=hwmean(m,l)/kountw(l)
                awmean(m,l)=awmean(m,l)/kountw(l)
                zwmean(m,l)=zwmean(m,l)/kountw(l)
            endif
        4587 END DO
    4876 END DO

    DEALLOCATE(KNDIR)
    DEALLOCATE(KLDIR)
    DEALLOCATE(KWDIR)
    DEALLOCATE(HSSB,NDSB)

    write(6,*) 'to conh12t'


!	call conh12t(AGTK,imjm,1,AGT,imt,jmt)

!	write(6,*) 'zoeff values'
    do I=1,IMJM,50
    !	write(6,297) (z0eff(I,M),M=1,4)
    enddo
    297	format(4(e10.4,x))

    ALLOCATE(Z0EFF2(IM,JM,4))
    do 8903 m=1,4
        call conh12(z0eff(1,m),imjm,1,z0eff2(1,1,m),im,jm)
    8903 END DO

    write(6,*) '2d zoeff values'
    do J=JM,1,-2
    !	write(6,299) (z0eff2(I,J,2),I=1,IM,3)
    enddo
    299	format(33(f4.2,x))


! CCC write out the ZEFF file CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    open(unit=27,file='ZEFF',form='unformatted',access='sequential' &
    ,status='unknown')


    ALLOCATE(TMP(IM,JM))
    do m=1,4
    !	write(6,*) 'BOGUS ZEFF!!!!!!!!!!!!!!!!!!!!'
    !	write(6,*) 'BOGUS ZEFF!!!!!!!!!!!!!!!!!!!!'
    !	write(6,*) 'BOGUS ZEFF!!!!!!!!!!!!!!!!!!!!'
        tmp=-999.
        do j=1,jm
            do i=1,im
                TMP(I,J)=z0eff2(I,J,M)
            ! og	TMP(I,J)=z0eff2(I,J,M)*5.
            enddo
        enddo
        write(27) TMP
    enddo

    DEALLOCATE(Z0EFF2)
    DEALLOCATE(Z0EFF)
    DEALLOCATE(ZHA)
    DEALLOCATE(ZAS)
    DEALLOCATE(TMP)
    


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! *** Reorder hgtk and smk array to old 2-d configuration.


    ALLOCATE (HGT(IM,JM))
    ALLOCATE (SM(IM,JM))

    call conh12(hgtk,imjm,1,hgt,im,jm)
    call conh12(smk,imjm,1,sm,im,jm)
!      open(1,file=topo_out,status='unknown',form='unformatted')
!      write(1) hgtk,smk
!      close(1)

    DEALLOCATE(SMK)
    DEALLOCATE(HGTK)


    write(6,*) 'sea mask at end of topo '
    do J=JM,1,-JM/30
        write(6,273)(sm(I,J),I=1,IM,IM/15)
    enddo
    273 format(80f2.0)


!	write(6,*) 'topo data before vertical interp (W,center,E) '
    ISTART=1

    do J=jm,1,-4
        write(6,275) (hgt(I,J),I=1,IM,IM/15)
        275 format(22(f5.0,x))
    enddo


! Check for bogus, elevated water points (water islands)

    do J=1,jm
        ihw(j)=-1+mod(j+1,2)
        ihe(j)=ihw(j)+1
    enddo


    do ITER=1,5
        write(6,*) ' START ITERATION ', ITER

        do J=3,jm-2
            do I=3,im-2
                
                if (hgt(i,j) > 0  .AND. sm(i,j) == 1) then
                ! elevated water
                    if ( ( hgt(I+IHE(J),J+1) == 0 .AND. sm(I+IHE(J),J+1) == 1) .OR. &
                    ( hgt(I+IHW(J),J+1) == 0 .AND. sm(I+IHW(J),J+1) == 1) .OR. &
                    ( hgt(I+IHW(J),J-1) == 0 .AND. sm(I+IHW(J),J-1) == 1) .OR. &
                    ( hgt(I+IHE(J),J-1) == 0 .AND. sm(I+IHE(J),J-1) == 1) .OR. &
                    ( hgt(I-1,J) == 0 .AND. sm(I-1,J) == 1 ) .OR. &
                    ( hgt(I+1,J) == 0 .AND. sm(I+1,J) == 1 ) .OR. &
                    ( hgt(I,J+2) == 0 .AND. sm(I,J+2) == 1 ) .OR. &
                    ( hgt(I,J-2) == 0 .AND. sm(I,J-2) == 1 )) then

                        if (hgt(i,J) < 100) then
                            write(6,*) 'making elevated water point zero elev at : ', i,j
                            hgt(i,j)=0.
                        else
                            write(6,*) 'making elevated water point a land point at: ', i,j
                            sm(i,j)=0.
                        endif

                    endif
                endif

            enddo
        enddo

    enddo
    

! *** Write eta grid topo.


    open(1,file=topo_out,status='unknown',form='unformatted')
    write(1) hgt,sm
    close(1)

    DEALLOCATE(HGT)
    DEALLOCATE(SM)
    stop
!     open(1,file='planine')
!     do j=jmt,1,-1
!        write(1,'(1000F6.0)')(hgt(i,j),i=1,imt)
!     enddo
!     close(1)

! tatic      return

! *** Error trapping.

    900 print *,'Error opening input topo data.**********************'
    print *,'   Filename = ',fname
    print*, '*************************************************'
    print*, '***********!!!!!!!!!!!!!!!!!********************'
    print*, '***********!!!!!!!!!!!!!!!!!********************'
    print*, '***********!!!!!!!!!!!!!!!!!********************'
    print*, ' M I S S I N G   S O M E   D A T A !!!!!!!!!!!!'
    print*, '***********!!!!!!!!!!!!!!!!!********************'
    print*, '***********!!!!!!!!!!!!!!!!!********************'
    print*, '***********!!!!!!!!!!!!!!!!!********************'
    print*, '*************************************************'
    stop

    return
    end subroutine mainprog


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!C************************************************************

    SUBROUTINE PQK_NEW &
    (ALM,APH,INSIDE,P,Q,K      )
!     ******************************************************************
!     *                                                                *
!     *  ROUTINE TO TRANSFORM ALAMBDA,APHI (LONGITUDE,LATITUDE) INTO   *
!     *  TRANSFORMED (ROTATED) LONGITUDE,LATITUDE, CHECK WHETHER THE   *
!     *  VALUES OBTAINED ARE WITHIN THE MODEL REGION, AND, IF THEY     *
!     *  ARE, CALCULATE VALUES OF P,Q (COORDINATES WITHIN A SQUARE     *
!     *  FORMED BY CONNECTING FOUR NEIGHBORING HEIGHT POINTS) AND K.   *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  GRID CONSTANTS:                                               *
!     *  WBD,SBD - WESTERN AND SOUTHERN BOUNDARIES ON LL OR TLL GRID   *
!     *            IN DEGREES                                          *
!     *  TLM0D,TPH0D - ANGLES OF ROTATION OF THE LL COORDINATE SYSTEM  *
!     *            IN THE DIRECTION OF LAMBDA AND PHI RESPECTIVELY     *
!     *            IN ORDER TO OBTAIN THE TLL COORDINATE SYSTEM        *
!     *  DLMD,DPHD - MESH SIDES IN DEGREES                             *
!     *                                                                *
!     ******************************************************************
!     *  5/95  Wobus  use whole domain, calculate the actual eta box k *
!     *               and internal coordinates p and q here            *
!     ******************************************************************
!     include "parmeta"
!      include "parmgrd"
    include 'ecommons.h'
!      include "parmsub"
!      include "parmloc"
    PARAMETER (DTR=0.01745329,D50=0.50,H1=1.)
    common /pqkcnt/nwest,neast,nsouth,nnorth
    common /pqkcnt/kwest,keast,ksouth,knorth
    common /pqkcnt/kpneg,kqneg,kklow,kkhigh
!----------------------------------------------------------------------
    LOGICAL ::  INSIDE

!---------------------------CONSTANTS-----------------------------------

    RIM1=IM1
    TLM0=TLM0D*DTR
    TPH0=TPH0D*DTR
    DLM=DLMD*DTR
    DPH=DPHD*DTR

    RDLM=H1/DLM
    RDPH=H1/DPH
! lwb
! test for whole domain
!     ALMWB=WBD*DTR +DLM
!     APHSB=SBD*DTR +DPH
!     ALMEB=ALMWB   +2*(IM-2)*DLM
!     APHNB=APHSB   +  (JM-3)*DPH
!     WB=ALMWB-DLM
!     SB=APHSB-DPH
    wb=wbd*dtr
    sb=sbd*dtr
    ALMWB=WB     - DLM
    APHSB=SB     - DPH
    ALMEB=ALMWB  +(2*IM  )*DLM
    APHNB=APHSB  +(  JM+1)*DPH
! lwe
    STPH0= SIN(TPH0)
    CTPH0= COS(TPH0)

!--------------ALM,APH IS LON,LAT IN RADIAN MEASURE;--------------------
!------------TRANSFORM INTO ROTATED LONGITUDE,LATITUDE------------------

    RLM=ALM-TLM0
    SRLM= SIN(RLM)
    CRLM= COS(RLM)
    SAPH= SIN(APH)
    CAPH= COS(APH)
    CC=CRLM*CAPH

    ANUM=SRLM*CAPH
    DENOM=CTPH0*CC+STPH0*SAPH

    ALM= ATAN2(ANUM,DENOM)
    APH= ASIN(CTPH0*SAPH-STPH0*CC)

    if (neast < 5) then
    !	write(6,*) 'bounds (w,s,e,n)'
    !	write(6,*) almwb,aphsb,almeb,aphnb
    !	write(6,*) 'we are at ',alm,aph
    endif

!-----IS ALM,APH INSIDE THE MODEL DOMAIN?-------------------------------

    kmiss=0
    if (alm < almwb) then
        nwest=nwest+1
        kmiss=kmiss-1
    endif
    if (alm > almeb) then
        neast=neast+1
        kmiss=kmiss-2
    endif
    if (aph < aphsb) then
        nsouth=nsouth+1
        kmiss=kmiss-4
    endif
    if (aph > aphnb) then
        nnorth=nnorth+1
        kmiss=kmiss-8
    endif
!	write(6,*) 'w,e,s,n ', nwest,neast,nsouth,nnorth
    if (kmiss /= 0) go to 102
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

!---------X1,Y1 IS A COORDINATE SYSTEM WITH DLM,DPH AS LENGTH UNITS-----
!            add y offset to index the centers of the boxes
!            instead of the southern corners
    X1=(ALM-WB)*RDLM
    Y1=(APH-SB)*RDPH+1.0
!---------X2,Y2 ROTATED FOR +45 DEG. & TRANSLATED FOR IM1---------------
    X2=D50*( X1+Y1)
    Y2=D50*(-X1+Y1)+RIM1
!---------I2,J2 ARE COORDINATES OF center POINTS OF GRID BOXES--------
    I2=  INT(X2)
    J2=  INT(Y2)
!---------REMAINING PARAMETERS NEEDED TO KNOW THE POSITION--------------
!         OF THE TRANSFORMED POINT WITHIN THE MODEL REGION
    P=X2-I2
    Q=Y2-J2
!-----------------INDEX K CORRESPONDS TO I2,J2--------------------------
    JR=J2-IM1
    I3=I2-JR
    J3=I2+JR
    K=J3*IM-J3/2+(I3+2)/2
!---------------see if the point is in the domain-----------------------
!  note, i3 runs from 0 thru imt-1
!        j3 runs from 0 thru jm-1
!        P and q must be nonnegative for correct 0,0 point
    kmiss=0
    if (i3 < 0) then
        kwest=kwest+1
        kmiss=kmiss+1
    endif
    if (i3 >= imt) then
        keast=keast+1
        kmiss=kmiss+2
    endif
    if (j3 < 0) then
        ksouth=ksouth+1
        kmiss=kmiss+4
    endif
    if (j3 >= jm) then
        knorth=knorth+1
        kmiss=kmiss+8
    endif
    if (p < 0) then
        kpneg=kpneg+1
        kmiss=kmiss+16
    endif
    if (q < 0.0) then
        kqneg=kqneg+1
        kmiss=kmiss+32
    endif
    if (k <= 0) then
        kklow=kklow+1
        kmiss=kmiss+64
    endif
    if (k > imjm) then
        kkhigh=kkhigh+1
        kmiss=kmiss+128
    endif
    if (kmiss > 0) go to 102

    INSIDE= .TRUE. 
    RETURN
!-----------------------------------------------------------------------

    102 continue
    if (kmiss == 64) then
        write(6,'('' unexpected k too low '',3i8,2f10.6)') &
        k,i3,j3,p,q
    endif
    if (kmiss == 128) then
        write(6,'('' unexpected k too high'',3i8,2f10.6)') &
        k,i3,j3,p,q
    endif
    INSIDE= .FALSE. 
    RETURN
    END SUBROUTINE PQK_NEW
                                 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCC  GET_SEA  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!	subroutine eta_sea(isea,idim,jdim)
    subroutine get_sea(isea,idim,jdim)

    implicit none

    include 'ecommons.h'

    INTEGER :: IRES
    integer :: nbx,nby,idim,jdim
    real :: cxsta,cxsto,cysta,cysto,factor

    real :: nbdnew,wbdnew,sbdnew,ebdnew &
    ,almin,almax,apmin,apmax,testlat,apmaxnew,apminnew &
    ,dummy,testlon,alminnew,almaxnew

    integer :: ista,isto,jsta,jsto &
    ,nx,ny,k,jx,ix,ji &
    ,i,ii,j,jj,l,istatus,khl22,khh22,GDS(200) &
    ,III,JJJ,icount,ilimit,jlimit

    real :: isea(imjm), &
    r2d,tph0,wb,sb,dlm,dph,tdlm,tdph,glatd,glond &
    ,tph,tlm,stph,ctph,sinphi,coslam,fact,dtr,x,y &
    ,dlat,frac,sum

    PARAMETER (DTR=0.01745329)

    REAL, ALLOCATABLE:: DUM2D(:,:)
    REAL, ALLOCATABLE:: SEABYTE(:,:),GLATR(:,:),GLONR(:,:)

    nbx=idim
    nby=jdim

    write(6,*) 'starting sea program, seares= ', seares
    write(6,*) 'sigma= ', sigma

    if (seares == 4) then

        write(6,*) 'opening 4 minute'
        open(unit=2,file='global_4m.ieee',form='unformatted', &
        access='sequential',err=920,status='old')

    elseif (seares == 8) then

        write(6,*) 'opening 8 minute'
        open(unit=2,file='global_8m.ieee',form='unformatted', &
        access='sequential',err=920,status='old')

    elseif (seares == 2) then

        open(2,file='US_2m_slm.ieee',form='unformatted', &
        access='sequential',status='old',ERR=920)
             
    else
            
        write(6,*) 'opening something with an unsupported IRES value...'
        stop

    endif

    ALLOCATE (SEABYTE(IDIM,JDIM))

    write(6,*) 'unit 2 open ', nbx,nby
    do J=1,nby
        read(2) (seabyte(I,J),I=1,nbx)
    enddo
    goto 921

    920 write(6,*) '!!!!!!!!!trouble opening the seamask data!!!!!!'
    921 write(6,*) 'data read successfully '
    close(2)

    GDS=0
    
    if (seares == 2) then
        GDS(1)=0
        GDS(2)=4501
        GDS(3)=2101
        GDS(4)=9983
        GDS(5)=179983
        GDS(6)=128
        GDS(7)=80017
        GDS(8)=-29983
        GDS(9)=033
        GDS(10)=033
    elseif (seares == 4) then
        GDS(1)=0
        GDS(2)=5400
        GDS(3)=2700
        GDS(4)=-89933
        GDS(5)=-179933
        GDS(6)=128
        GDS(7)=89933
        GDS(8)=179933
        GDS(9)=067
        GDS(10)=067
    elseif (seares == 8) then
        GDS(1)=0
        GDS(2)=2700
        GDS(3)=1350
        GDS(4)=-89967
        GDS(5)=-179967
        GDS(6)=128
        GDS(7)=89967
        GDS(8)=179967
        GDS(9)=133
        GDS(10)=133
    endif

    dlat=GDS(9)/1000.
    ilimit=0.5*(dlmd/dlat)
    jlimit=0.5*(dphd/dlat)

!	ilimit=max0(1,ilimit)
!	jlimit=max0(1,jlimit)

    write(6,*) 'if possible, will go +/- ', ilimit, 'in I and +/-', &
    jlimit, 'in J'


! p     values will be ones or twos (1=land, 2=water)

! p     for output (0=land, 1=water)

!	compute lat lon on e-grid
    r2d=57.2957795
    tph0=tph0d*dtr
    wb=wbd*dtr
    sb=sbd*dtr
    dlm=dlmd*dtr
    dph=dphd*dtr
    tdlm=dlm+dlm
    tdph=dph+dph
    tph=sb-dph
    ctph0=cos(tph0)
    stph0=sin(tph0)

    ALLOCATE(GLATR(IM,JM),GLONR(IM,JM))

    do j=1,jm
        tlm=wb-tdlm+mod(j+1,2)*dlm
        tph=tph+dph
        stph=sin(tph)
        ctph=cos(tph)
    
        do i=1,im
            tlm=tlm+tdlm
            sinphi=ctph0*stph+stph0*ctph*cos(tlm)
            glatr(i,j)=asin(sinphi)
            coslam=ctph*cos(tlm)/(cos(glatr(i,j))*ctph0) &
            -tan(glatr(i,j))*tan(tph0)
            coslam=min(coslam,1.)
            fact=1.
            if (tlm > 0.0) fact=-1.
            glonr(i,j)=-tlm0d*dtr+fact*acos(coslam)
        
        enddo
    enddo

!	Initialize as all sea...boundaries arent recomputed and
!	will remain sea points

    ISEA=1
   

    if ( .NOT. SIGMA) then
        KHL22=2*IM+1
        KHH22=IMJM-2*IM
    else
        KHL22=1
        KHH22=IMJM
    endif

    do K=KHL22,KHH22

    ! added 001003
    ! COMPUTE OUTPUT GRID lats AND lons

    !	convert K value to get I,J
        IMT=2*IM-1
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX <= IM)THEN
            I=IX
            J=2*JX-1
        ELSE
            I=IX-IM
            J=2*JX
        ENDIF

    ! DETERMINE LAT/LON OF TARGET (output) GRID POINT
    
        GLATD = GLATR(i,j)*r2d
        GLOND = GLONR(I,J)*r2d
        IF (GLOND < 0) GLOND=GLOND+360
        GLOND = 360. - GLOND
    
    ! DETERMINE NEAREST NEIGHBOR FROM INPUT GRID
    
    !	write(6,*) 'i,j,glat,glon: ', i,j,glatd,glond
        call ced_ij(GLATD,glond,x,y,gds)
    !	write(6,*) 'found x,y: ', x,y

        II=INT(X+0.5)
        JI=INT(Y+0.5)

    !	If input data is 1 --> land --> isea=0
    !	If input data is 0 --> water --> isea=1
    


        if (II >= ilimit .AND. II <= nbx-ilimit .AND. &
        JI >= jlimit .AND. JI <= nby-jlimit) then

            ICOUNT=0
            sum=0.

            do JJJ=JI-jlimit,JI+jlimit
                do III=II-ilimit,II+ilimit
                    sum=sum+seabyte(III,JJJ)
                    ICOUNT=ICOUNT+1
                enddo
            enddo
        

            frac=sum/float(icount)

            if (frac > FRACLK) then
                isea(k)=0
            else
                isea(k)=1
            endif



        ELSE
            
            write(6,*) 'too close to boundary...just use n.n.'

            if (seabyte(II,JI) == 1) then
                isea(k)=0
            elseif (seabyte(II,JI) == 0) then
                isea(k)=1
            else
                write(6,*) 'bad seabyte value!!!!! ', II,JI, seabyte(II,JI)
            endif
            

        ENDIF

    enddo

    DEALLOCATE(SEABYTE)
    DEALLOCATE(GLATR,GLONR)

    write(6,*) 'GOT TO THIS POINT'

    ALLOCATE (dum2d(im,jm))
    call conh12(isea,imjm,1,dum2d,im,jm)
    do J=JM-1,1,-JM/40
    !        write(6,273)(dum2d(I,J),I=1,IM,IM/53)
    enddo
    273 format(71f2.0)

    DEALLOCATE(DUM2D)

    return
    end subroutine get_sea

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    subroutine CORNERS(im,jm,IMJM,tph0d,tlm0d,dlmd,dphd,minlat,wlon, &
    maxlat,elon)


!     *  ROUTINE TO FIND EARTH LATITUDE/LONGITUDE FOR THE CORNER       *
!     *               POINTS OF AN ETA MODEL GRID                      *


!-----------------------------------------------------------------------
!                             D I M E N S I O N
!     & KHL0  (JM),KHH0  (JM), GLAT(IMJM),GLON(IMJM)

    REAL,ALLOCATABLE::GLAT(:),GLON(:)
    INTEGER,ALLOCATABLE::KHL0(:),KHH0(:)

    D A T A &
    PI/3.141592654/

    REAL :: DLMD,DPHD,WBD,SBD,TPH0D,TLM0D,minlat,wlon,maxlat,elon

!*******************************
    WBD=-(float(IM)-1.)*DLMD
    SBD=(-(float(JM)-1.)/2.)*DPHD

    DTR=PI/180.
    TPH0=TPH0D*DTR
    WB=WBD*DTR
    SB=SBD*DTR
    DLM=DLMD*DTR
    DPH=DPHD*DTR
    TDLM=DLM+DLM
    TDPH=DPH+DPH

    STPH0=SIN(TPH0)
    CTPH0=COS(TPH0)

    ALLOCATE(KHL0(JM),KHH0(JM))
    ALLOCATE(GLAT(IMJM),GLON(IMJM))

    DO 100 J=1,JM
        KHL0(J)=IM*(J-1)-(J-1)/2+1
        KHH0(J)=IM*J-J/2
    !     WRITE(6,9999) J, KHL0(J), KHH0(J)
    ! 999 FORMAT(2X,3(I10,1X))
    100 END DO
!--------------GEOGRAPHIC LAT AND LONG OF TLL GRID POINTS---------------
    TPH=SB-DPH
    maxlat=-999.
    minlat=99.
    DO J=1,JM
        KHL=KHL0(J)
        KHH=KHH0(J)
    
        TLM=WB-TDLM+MOD(J+1,2)*DLM
        TPH=TPH+DPH
        STPH=SIN(TPH)
        CTPH=COS(TPH)
        DO K=KHL,KHH
            TLM=TLM+TDLM
            SPH=CTPH0*STPH+STPH0*CTPH*COS(TLM)
            GLAT(K)=ASIN(SPH)
            CLM=CTPH*COS(TLM)/(COS(GLAT(K))*CTPH0)-TAN(GLAT(K))*TAN(TPH0)
            IF(CLM > 1.)      CLM=1.
            FACT=1.
            IF(TLM > 0.)      FACT=-1.
            GLON(K)=(-TLM0D*DTR+FACT*ACOS(CLM))/DTR

        ! p     at this point GLON is in DEGREES WEST
            if (GLON(K) < 0) GLON(K)=GLON(K)+360.
            if (GLON(K) > 360.) GLON(K)=GLON(K)-360.
            if (GLON(K) < 180) GLON(K)=-GLON(K)         ! make WH negative
            if (GLON(K) > 180) GLON(K)=360.-GLON(K)     ! make EH

            GLAT(K)=GLAT(K)/DTR

            if (glat(k) > maxlat) maxlat=glat(k)
            if (glat(k) < minlat) minlat=glat(k)

        enddo
    enddo

    DEALLOCATE(KHL0,KHH0)
    
   
    if (TPH0D >= 0) then
        wlon=glon(imjm-im+1)
        elon=glon(imjm)
    else
        wlon=glon(1)
        elon=glon(im)
    endif

!	write(6,*) 'raw lon values (w,e) ', wlon, elon
    if (tlm0d < 0 .AND. wlon > 0) wlon=wlon-360.
    if (tlm0d > 0 .AND. elon < 0) elon=elon+360.

    DEALLOCATE(GLAT,GLON)

    RETURN
    end subroutine CORNERS

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    SUBROUTINE SWAPVAL(A)

!      REVERSE ORDER OF BYTES IN INTEGER*2 WORD, or REAL*2

    INTEGER*2 ::   A

    CHARACTER(1) :: JTEMP(2)
    CHARACTER(1) :: KTEMP

    EQUIVALENCE (JTEMP(1),ITEMP)

    ITEMP    = A
    KTEMP    = JTEMP(2)
    JTEMP(2) = JTEMP(1)
    JTEMP(1) = KTEMP
    A     = ITEMP
    RETURN
    END SUBROUTINE SWAPVAL

! add subroutines

    subroutine conh12t(h1,imjm,lm,h2,imt,jmt)

! *** Routine for reordering thibue-height-point-1-dimensional
!        matrices for 2-dimensional indexing (imt,jmt).

    implicit none

    integer*4 :: imt,jmt,imjm,lm,i,j,k,l

    real*4 :: h1(imjm,lm),h2(imt,jmt,lm)
! ______________________________________________________________________________

    do l=1,lm
        k=0
        do j=1,jmt
            do i=1,imt,2
                k=k+1
                h2(i,j,l)=h1(k,l)
            enddo
        
            if (j < jmt) then
                do i=2,imt,2
                    k=k+1
                    h2(i,j,l)=h1(k,l)
                enddo
            endif
        enddo
    enddo

    return
    end subroutine conh12t

!===============================================================================

!===============================================================================

    subroutine conh12(h1,imjm,lm,h2,im,jm)

! *** Routine for reordering thibue-height-point-1-dimensional
!        matrices (im,jm) for 2-dimensional indexing.

    implicit none

    integer*4 :: im,jm,imjm,lm,i,j,k,l,imm1

    real*4 :: h1(imjm,lm),h2(im,jm,lm)
! ______________________________________________________________________________

    imm1=im-1
    do l=1,lm
        k=0
        do j=1,jm
            do i=1,imm1+mod(j,2)
                k=k+1
                h2(i,j,l)=h1(k,l)
            enddo
        enddo
        do j=2,jm-1,2
            h2(im,j,l)=h2(imm1,j,l)
        enddo
    enddo

    return
    end subroutine conh12

!===============================================================================

    subroutine ced_ij(RLAT,RLON,XPTS,YPTS,KGDS)

    integer :: kgds(200)


    IM=KGDS(2)
    JM=KGDS(3)
    RLAT1=KGDS(4)*1.E-3
    RLON1=KGDS(5)*1.E-3
    RLAT2=KGDS(7)*1.E-3
    RLON2=KGDS(8)*1.E-3
    ISCAN=MOD(KGDS(11)/128,2)
    JSCAN=MOD(KGDS(11)/64,2)
    NSCAN=MOD(KGDS(11)/32,2)
    HI=(-1.)**ISCAN
    HJ=(-1.)**(1-JSCAN)
    DLON=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(IM-1)
    DLAT=(RLAT2-RLAT1)/(JM-1)
    XMIN=0
    XMAX=IM+1
    IF(IM == NINT(360/ABS(DLON))) XMAX=IM+2
    YMIN=0
    YMAX=JM+1
    NRET=0
    LROT=0
    FILL=-999

    IF(ABS(RLON) <= 360 .AND. ABS(RLAT) <= 90) THEN
        XPTS=1+HI*MOD(HI*(RLON-RLON1)+3600,360.)/DLON
        YPTS=1+(RLAT-RLAT1)/DLAT

    ! p
        if (XMAX == 5401 .OR. XMAX == 2701) then

            if (XPTS > XMAX) then
                XPTS=XPTS-XMAX
                write(6,*) 'xpts now: ', xpts
            endif

            if (XPTS < XMIN) then
                XPTS=XPTS+XMAX
                write(6,*) 'xpts now: ', xpts
            endif

        endif
    ! p
        IF(XPTS >= XMIN .AND. XPTS <= XMAX .AND. &
        YPTS >= YMIN .AND. YPTS <= YMAX) THEN
            NRET=NRET+1
            IF(LROT == 1) THEN
                CROT=1
                SROT=0
            ENDIF
        ELSE
            XPTS=FILL
            YPTS=FILL
        ENDIF
    ELSE
        XPTS=FILL
        YPTS=FILL
    ENDIF

    return
    end subroutine ced_ij

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    subroutine get_sea30(isea)
    include 'ecommons.h'
    parameter(ID1=1,ID2=43204,IRECL=ID2)
    parameter(dlat=1./120.)
    parameter(r2d=57.29577951,dtr=0.017453293)
    integer :: ounit,iunit,nunit,iirecl
    parameter(OUNIT=50)
    character(12) :: fnameout
    INTEGER*1 :: LANDIN(ID1,ID2)
    INTEGER,ALLOCATABLE:: LAND2D(:,:)
    INTEGER :: GDS(200)
    character(2) :: uname
    CHARACTER(1) :: LANDIN_CHR(ID1,ID2)

    REAL, ALLOCATABLE:: GLATR(:,:),GLONR(:,:)
    real :: isea(imjm)
    real :: dum2d(im,jm)

    cshift=dlat/2.



!	details about target grid

!       compute lat lon on e-grid

    write(6,*) 'im,jm,wbd,sbd,tph0d,dtr ', im,jm,wbd,sbd,tph0d,dtr
    tph0=tph0d*dtr
    wb=wbd*dtr
    sb=sbd*dtr
    dlm=dlmd*dtr
    dph=dphd*dtr
    tdlm=dlm+dlm
    tdph=dph+dph
    tph=sb-dph
    ctph0=cos(tph0)
    stph0=sin(tph0)

    ALLOCATE(GLATR(IM,JM),GLONR(IM,JM))

    do j=1,jm
        tlm=wb-tdlm+mod(j+1,2)*dlm
        tph=tph+dph
        stph=sin(tph)
        ctph=cos(tph)
    
        do i=1,im
            tlm=tlm+tdlm
            sinphi=ctph0*stph+stph0*ctph*cos(tlm)
            glatr(i,j)=asin(sinphi)
            coslam=ctph*cos(tlm)/(cos(glatr(i,j))*ctph0) &
            -tan(glatr(i,j))*tan(tph0)
            coslam=min(coslam,1.)
            fact=1.
            if (tlm > 0.0) fact=-1.
            glonr(i,j)=-tlm0d*dtr+fact*acos(coslam)
        
        enddo
    enddo

    call corners(im,jm,imjm,tph0d,tlm0d,dlmd,dphd,apmn,almn,apmx,almx)

    rnlat=apmx+0.1
    slat=apmn-0.1
    west=almn-0.1
    east=almx+0.1

    write(6,*) glatr(im/2,jm),glatr(1,1),glonr(1,jm),glonr(im,jm),r2d

    write(6,*) 'west= ', west
    write(6,*) 'east= ', east
    write(6,*) 'north= ', rnlat
    write(6,*) 'south= ', slat

    IST=((west+180.-cshift)/dlat)+1
    IEND=((east+180.-cshift)/dlat)+1

    write(6,*) 'IST, IEND: ', IST,IEND

    wtrue=(ist-1)*dlat+cshift-180.
    etrue=(iend-1)*dlat+cshift-180.

    JST=-(RNLAT-90)/dlat+1
    JEND=-(SLAT-90)/dlat+1

    write(6,*) 'JST, JEND: ', JST, JEND

    rntrue=90.-(jst-1)*dlat-cshift
    strue=90.-(jend-1)*dlat-cshift

    ILIM=IEND-IST+1
    JLIM=JEND-JST+1

    write(6,*) 'ILIM, JLIM: ', ILIM, JLIM

!	STOP

    allocate(land2d(ILIM,JLIM))

    LAND2D=-99.

    JSUNIT=INT((JST-1)/600)+1
    JENDUNIT=INT((JEND-1)/600)+1

    write(6,*) 'true bounds: '
    write(6,*) 'N: ', rntrue
    write(6,*) 'W: ', wtrue
    write(6,*) 'S: ', strue
    write(6,*) 'E: ', etrue


!C	come up with an estimate of how many input 30 s points should
!C	be included.  Consider roughly enough to average over the size of
!C	the eta gridbox.

    ilimit=0.5*(dlmd/dlat)
    jlimit=0.5*(dphd/dlat)

    ilimit=max0(ilimit,1)
    jlimit=max0(jlimit,1)

!	ilimit=1
!	jlimit=1

    write(6,*) 'if possible, will go +/- ', ilimit, 'in I and +/-', &
    jlimit, 'in J'


    do J=JSUNIT,JENDUNIT
        write(uname,'(i2.2)') INT(J)
        fnameout='smask.30s.'//uname

        #ifdef DEC
        iirecl=IRECL/4
        #else
        iirecl=IRECL
        #endif
        IIUNIT=J+60
        open(unit=IIUNIT,file=fnameout, &
        access='direct',recl=IIRECL)
        write(6,*) 'opened unit,file ', IIUNIT,fnameout
    enddo

    GDS=0

    GDS(1)=0
    GDS(2)=ILIM
    GDS(3)=JLIM
    GDS(4)=INT(RNTRUE*1000)
    GDS(5)=INT(WTRUE*1000)
    GDS(6)=128
    GDS(7)=INT(STRUE*1000)
    GDS(8)=INT(ETRUE*1000)
    GDS(9)=INT(1./120.*1000)
    GDS(10)=INT(1./120.*1000)

    do J=1,10
        write(6,*) 'J,GDS(J): ', J,GDS(J)
    enddo

    do 80 N=JST,JEND

        IUNIT=INT((N-1)/600)+61
        IREC=N-(IUNIT-61)*600

        READ(IUNIT,REC=IREC,ERR=100) LANDIN_CHR

        DO IID1 = 1,ID1
            DO IID2 = 1,ID2
                LANDIN(IID1,IID2) = ICHAR(LANDIN_CHR(IID1,IID2))
            END DO
        END DO

        if (IST >= 1 .AND. IEND <= ID2 ) then

            do I=IST,IEND
                land2d(I-IST+1,N-JST+1)=LANDIN(1,I)
            enddo

        else

        !	write(6,*) 'special case!!!!'

            IF (IST < 1) then

            !	put the eastern part of LANDIN into LAND2D

                if (N == JST) then
                    write(6,*) 'LAND2D defined from : ', ID2+IST-(ID2+IST)+1, &
                    ID2-(ID2+IST)+1
                    write(6,*) '2nd part defines: ',ILIM-IEND+1, ILIM
                endif

                DO I=ID2+IST,ID2
                    LAND2D(I-(ID2+IST)+1,N-JST+1)=LANDIN(1,I)
                ENDDO

                DO I=1,IEND
                    LAND2D(ILIM-(IEND)+I,N-JST+1)=LANDIN(1,I)
                ENDDO

            ENDIF


        
            if (IEND > ID2) then

            !	Eastern part of target domain from western part of input data


                if (N == JST) then
                    write(6,*) 'LAND2D defined from : ', 1, ID2-IST+1
                    write(6,*) '2nd part defines: ',2+ID2-IST, IEND-ID2+ &
                    (ID2-IST+1)
                endif

                DO I=IST,ID2
                    LAND2D(I-IST+1,N-JST+1)=LANDIN(1,I)
                ENDDO

                DO I=1,IEND-ID2
                    LAND2D(I+(ID2-IST+1),N-JST+1)=LANDIN(1,I)
                ENDDO

            ENDIF
            
           
        !	STOP

        endif

    80 END DO

    ISEA=1

    KHL22=2*IM+1
    KHH22=IMJM-2*IM

    write(6,*) 'loop K from ', KHL22, KHH22
    do K=KHL22,KHH22
    !       convert K value to get I,J
        IMT=2*IM-1
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX <= IM)THEN
            I=IX
            J=2*JX-1
        ELSE
            I=IX-IM
            J=2*JX
        ENDIF

    ! DETERMINE LAT/LON OF TARGET (output) GRID POINT
    
        GLATD = GLATR(i,j)*r2d
        GLOND = GLONR(I,J)*r2d
        IF (GLOND < 0) GLOND=GLOND+360
        GLOND = 360. - GLOND

        call ced_ij(glatd,glond,x,y,gds)

        II=INT(X+0.5)
        JI=INT(Y+0.5)

        if (II >= ilimit .AND. II <= ILIM-ilimit .AND. &
        JI >= jlimit .AND. JI <= JLIM-jlimit) then

            ICOUNT=0
            sum=0.

            do JJJ=JI-jlimit,JI+jlimit
                do III=II-ilimit,II+ilimit
                    sum=sum+land2d(III,JJJ)
                    ICOUNT=ICOUNT+1
                enddo
            enddo
        

            frac=sum/float(icount)

            if (K == KHL22) then
                write(6,*) 'based this point on ', icount, '30 s data points'
                write(6,*) 'sea mask fraction (lower=fewer lakes) = ', fraclk
            endif

            if (frac > 0 .AND. frac < 1) then
            !!	write(6,*) 'K, fraction = ', K,frac
            endif

        !       If input data is 1 --> land --> isea=0
        !       If input data is 0 --> water --> isea=1

            if (frac > FRACLK) then
                isea(k)=0
            else
                isea(k)=1
            endif
            

        ELSE

            write(6,*) 'too close to the boundary...use n.n. ', K

            if (land2d(II,JI) == 1) then
                isea(k)=0
            elseif (land2d(II,JI) == 0) then
                isea(k)=1
            else
                write(6,*) 'bad land2d values!!!! ', II,JI,land2d(II,JI)
            endif

        ENDIF


    enddo

    icount=0

    do K=1,imjm
        if  (isea(K) == 1) then
            ICOUNT=ICOUNT+1
        endif
    enddo

    write(6,*) 'have ', ICOUNT , ' sea points out of ', &
    imjm, 'points total'

    call conh12(isea,imjm,1,dum2d,im,jm)
    write(6,*) 'sea within sea30 subroutine'
    do J=JM-1,1,-JM/40
        write(6,273)(dum2d(I,J),I=1,IM,IM/30)
    enddo
    273 format(71f2.0)

    RETURN

    STOP
    100	write(6,*) 'bad read!!!!'

    end subroutine get_sea30
