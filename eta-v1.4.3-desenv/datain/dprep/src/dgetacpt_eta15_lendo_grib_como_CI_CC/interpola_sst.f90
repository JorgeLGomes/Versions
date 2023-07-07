PROGRAM interpola

      PARAMETER(nxin=4320,nyin=2160)   
      PARAMETER(nxout=1440,nyout=720)   

      REAL sst_in(nxin,nyin)
      REAL sst_out(nxout,nyout)
      
      CHARACTER date*8
      
      CALL getarg(1,date)
 
      print*, "Interpolating SST HIRES(0.1) to MRES(0.25) ..."
      
      OPEN (unit=10,file='gdas1.T00Z.sstgrd.'//date//'00',form='unformatted',convert='big_endian')
      
      OPEN (unit=20,file='sst_Eta_'//date//'.bin',form='unformatted',convert='big_endian')     
      
      READ(10)((sst_in(I,J),I=1,nxin),J=1,nyin)  
      
      Io=1
      
      DO I=1,nxin,3
     
      Jo=1

      DO J=1,nyin,3

!      print*, Io, Jo, I, J

      sst_out(Io,Jo)= ( sst_in(I,J) + sst_in(I,J+1) + sst_in(I,J+2)         &
                       +sst_in(I+1,J) + sst_in(I+1,J+1) + sst_in(I+1,J+2)   &
		       +sst_in(I+2,J) + sst_in(I+2,J+1) + sst_in(I+2,J+2) ) / 9

      Jo=Jo+1

      ENDDO

      Io=Io+1

      ENDDO
      

      WRITE(20)((sst_out(I,J),I=1,nxout),J=1,nyout)
      
      
END PROGRAM
