	program create

	parameter (ix=10800,jx=5400)
	parameter (iout=4501,jout=2101)

	real input(ix,jx)
	real output(iout,jout)

	open(unit=11,file='TOP2M_slm.ieee',form='unformatted',
     +		access='sequential')


	read(11) input


	do J=jx,1,-jx/50
C	write(6,64) (input(I,J),I=1,ix,ix/45)
	enddo

C
C	create output from J=3000 (9.983N) to J=5101 (80.0166N)
C	and from I=1 (179.983W) to I=4501 (29.96666 W)
c	

	write(6,*) 'output window'

	do J=5101,3000,-80
C	write(6,64) (input(I,J),I=5400,9901,4501/30)
	enddo

   64	format(50f2.0)


C
C       Write out the data subset
C

	do J=3000,5101
	JJ=J-2999
	do I=5400,9901
	II=I-5399
	output(II,JJ)=input(I,J)
	enddo
	enddo

	do J=2101,1,-80
	write(6,64) (output(I,J),I=1,4501,120)
	enddo

	open(unit=12,file='northfirst_2m.ieee',form='unformatted',
     +			access='sequential')

	do J=jout,1,-1
	write(12) (output(I,J),I=1,iout)
	enddo

	close(12)

	open(unit=13,file='southfirst_2m.ieee',form='unformatted',
     +			access='sequential')

	do J=1,jout
	write(13) (output(I,J),I=1,iout)
	enddo


	end
