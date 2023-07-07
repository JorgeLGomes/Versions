      subroutine es_ini
      use estab
      implicit none
      integer        :: it
      real           :: t,p1,p2,c1
!
! *** Create tables of the saturation vapour pressure with up to
!        two decimal figures of accuraccy:
!
       do it=15000,45000
         t=it*0.01
         p1 = 11.344-0.0303998*t
         p2 = 3.49149-1302.8844/t
         c1 = 23.832241-5.02808*alog10(t)
         esat(it) = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/t)
         es(it) = 610.78*exp(17.269*(t-273.16)/(t-35.86))
       enddo
      end subroutine es_ini
