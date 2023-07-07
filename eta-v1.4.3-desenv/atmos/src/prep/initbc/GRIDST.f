C$$$  sub  PROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C sub  PROGRAM:  gridst   READ GLOBAL Grib SST
C                            
C
C ABSTRACT: SUBSTITUTE THE GRIBST SUB PROGRAM
C           READ THE SST FIELD IN A BINARY FORMAT 
C           INSTED A GRIB FORMAT
C           IN THIS CASE THE BINARY FILE IS GENERATE
C           FROM SST GRIB2 FILE WITH WGRIB2 PROGRAM: 
C           wgrib2 -d 1 sst2dvar_grb_0.5.grib2 -v \
C           -no_header big_endian -bin ssthiresgrd.bin
C
C PROGRAM HISTORY LOG:
C   08-11-03  JORGE LUIS GOMES
C
C
      subroutine gridst(insst,gsst)
c
      parameter(nx=4320,ny=2160,nxp1=nx+1)

      REAl out(nx,ny)
      dimension gsst(nxp1,ny)
c
C
c
C            OPEN UNIT FOR READING GRID FILE
       OPEN (INSST,form='unformatted',convert='big_endian')
       READ(INSST)((out(I,2161-J),I=1,4320),J=1,2160)
C
              do  jj= 1,ny
                do  kk= 1,nx
                 gsst(kk,jj) = out(kk,jj)
                end do
              end do
c
c...   add greenich to right side of grid
c
             do jj = 1,ny
               gsst(nxp1,jj) = gsst(1,jj)
             end do
         write(95)gsst
c
C
   30        CONTINUE
C
  999           CONTINUE
	write(6,*) 'leaving GRIDST'
             return
             END
