!>--------------------------------------------------------------------------------------------------
!> CHECKIN.f90
!>
!> USE MODULES: -----
!>
!> DRIVER     : HZADV
!>              HZADV_LM1
!>
!> CALLS      : -----
!>-------------------------------------------------------------------------------------------------- 
    500 FORMAT(' ENTER HZADV')
!
    IEND = MYIE
!
    IF (MOD(MYPE+1,INPES) == 0) IEND = MYIE - 1
!
    DO 503 K=1,LM
        DO 502 J=MYJS,MYJE
            DO 502 I=MYIS,IEND
                IF (HTM(I,J,K) > 0.5) THEN
                    IF (T(I,J,K) < 150. .OR. T(I,J,K) > 335.) THEN
                        WRITE(6,500)
                        WRITE(6,504) I, J, K, T(I,J,K)
                        504 FORMAT(' I=',I3,' J=',I3,' K=',I2,' T=',E12.5)
                        STOP666
                    ELSE IF (Q(I,J,K) < -1.E-4 .OR. Q(I,J,K) > 30.E-3) THEN
                        WRITE(6,500)
                        WRITE(6,505) I, J, K, Q(I,J,K)
                        505 FORMAT(' I=',I3,' J=',I3,' K=',I2,' Q=',E12.5)
                        STOP666
                    ELSE IF (Q2(I,J,K) > 400.) THEN
                        WRITE(6,500)
                        WRITE(6,506) I, J, K, Q2(I,J,K)
                        506 FORMAT(' I=',I3,' J=',I3,' K=',I2,' Q2=',E12.5)
                        STOP666
                    ELSE IF (CWM(I,J,K) > 1.E-1) THEN
                        WRITE(6,500)
                        WRITE(6,507) I, J, K, CWM(I,J,K)
                        507 FORMAT(' I=',I3,' J=',I3,' K=',I2,' CWM=',E12.5)
                        STOP666
                    END IF
                END IF
    502 END DO
503 END DO
!
    DO 509 K=1,LM
        DO 508 J=MYJS,MYJE
            DO 508 I=MYIS,IEND
                IF (ABS(U(I,J,K)) > 125. .OR. ABS(V(I,J,K)) > 125.) THEN
                    WRITE(6,500)
                    WRITE(6,510) I, J, K, U(I,J,K), V(I,J,K)
                    510 FORMAT(' I=',I3,' J=',I3,' K=',I2,' U=',E12.5,' V=',E12.5)
                END IF
    508 END DO
509 END DO
