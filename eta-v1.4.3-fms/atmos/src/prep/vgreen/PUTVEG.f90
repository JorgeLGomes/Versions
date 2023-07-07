!   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       SUBROUTINE PUTVEG(GLAT,GLON,FVEG0,SM,SICE,FVEG1)
!   **************************************************************
!   * Interpolate NESDIS vegetation fraction product (five-year  *
!   * climatology with 0.144 degree resolution from 89.928S,180W * 
!   * to 89.928N, 180E)                                          *
!   * F. Chen 07/96                                              *
!   **************************************************************
 
    USE PARMCONF
    USE F77KINDS
!
    INTEGER(KIND=I4KIND), PARAMETER ::   L0  = 15680*10600
    REAL   (KIND=R4KIND), PARAMETER ::   RES = 0.00892857
!
    REAL   (KIND=R4KIND)   , DIMENSION(15680,10600)    :: FVEG0
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: FVEG1
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: GLAT
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: GLON
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: HLAT
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: HLON
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: SM
    REAL   (KIND=R4KIND)   , DIMENSION(IM,JM)          :: SICE   
        
!   This subroutine AVERAGES the value stored in FVEG0 TO the FVEG1 array.
!   Special effort is made to ensure that islands are taken care of when
!   not resolved on the 0.144 x 0.144 deg. grid.
!  *** Define the search limit for isolated island point 
        NLIM=10*(1+AINT(RES/DPHD))
!  *** Define the Eta averaging box size, NBOX=0 gives nearest neighbor
        NBOX1=NINT(MAX(DLMD,DPHD)/RES)
        NBOX=NBOX1
!        print*,' im,jm=',im,jm,'NBOX=',NBOX,'DPHD',DPHD,'NLIM',NLIM,
!     &         'DLMD',DLMD


        HLAT = GLAT
        DO J = 1,JM
        DO I = 1,IM

! 	added for workstation
        HLON(I,J) = 360.0 - GLON(I,J)
        IF(HLON(I,J) .GT. 360.) HLON(I,J) = GLON(I,J) - 360.
!

          IF((SM(I,J).GT.0.5).OR.(SICE(I,J).EQ.1.0)) THEN
            FVEG1(I,J) = 0.0
          ELSE
            DX=(HLON(I,J) - 240) / RES
            INDX = NINT(DX)
!  ***  Here, 179.928 is the starting longitude (180.072) + one grid cell 
!       width (0.144) so that the first index is one NOT 0.
            IF(INDX.LT.1) THEN 
              DX=(HLON(I,J) + 180.072) / RES
              INDX = NINT(DX)
            ENDIF
            DY=(HLAT(I,J) + 59.9911) / RES
            INDY = NINT(DY)
!  *** Get area-average value of FVEG1 from input grid FVEG0, for 
!      finer Eta grid, this becomes nearest-neighbour 
            ICON1=0
  100       CONTINUE
            IF(MOD(NBOX,2).EQ.0) THEN
              IF(DX.GT.REAL(INDX)) THEN
                NLO=MAX(1,INDX-NBOX/2+1)
                NHI=MIN(15680,INDX+NBOX/2)
              ELSE
                NLO=MAX(1,INDX-NBOX/2)
                NHI=MIN(15680,INDX+NBOX/2-1)
              ENDIF
              IF(DY.GT.REAL(INDY)) THEN
                MLO=MAX(1,INDY-NBOX/2+1)
                MHI=MIN(10600,INDY+NBOX/2) 
              ELSE
                MLO=MAX(1,INDY-NBOX/2)
                MHI=MIN(10600,INDY+NBOX/2-1)
              ENDIF
            ELSE
              NLO=MAX(1,INDX-NBOX/2)
              NHI=MIN(15680,INDX+NBOX/2)
              MLO=MAX(1,INDY-NBOX/2)
              MHI=MIN(10600,INDY+NBOX/2)
            ENDIF
!            NLO=MAX(1,INDX-NBOX)
!            NHI=MIN(15680,INDX+NBOX)
!            MLO=MAX(1,INDY-NBOX)
!            MHI=MIN(10600,INDY+NBOX)
            ICON=0
            VEGSUM=0.0
            DO N=NLO, NHI
              DO M=MLO, MHI
                IF(FVEG0(N,M).NE.0.0) THEN
                 ICON=ICON+1 
                 VEGSUM=VEGSUM+FVEG0(N,M)
                ENDIF
              ENDDO
            ENDDO
            IF(ICON.NE.0) THEN
              FVEG1(I,J)=VEGSUM/ICON
              NBOX=NBOX1
            ELSE
! *** Search for an isolated point 
              ICON1=ICON1+1
              NBOX=NBOX+1 
              IF(ICON1.LE.NLIM) THEN
                GO TO 100
              ELSE
                FVEG1(I,J)=0.55
                NBOX=NBOX1
                print *,"OH OH ",HLAT(I,J),HLON(I,J)
              ENDIF
            ENDIF
! ** End of land case
!          if(NBOX.NE.NBOX1) print*,' NBOX is NOT NBOX1,=',NBOX
          END IF  
          IF(FVEG1(I,J).LT.0.0.OR.FVEG1(I,J).GT.1.0) THEN 
           print*,' FVEG1 out of range, FVEG=',FVEG1(I,J)
          ENDIF
! ** End of IM,JM loops
        END DO  
        END DO  
        RETURN

        END
