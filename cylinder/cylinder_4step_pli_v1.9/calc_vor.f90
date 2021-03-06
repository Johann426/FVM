    SUBROUTINE CALC_VOR
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!----------------------------------------------------------------------------
! Cell node vorticity
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Cell center vorticity
!----------------------------------------------------------------------------
    DO K=2,NKM1
    DO J=1,NJ-1
    DO I=1,NI-1
    VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
                +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
                -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))
                
    ENDDO
    ENDDO
    ENDDO

    DO K=2,NKM1
    DO J=1,NJ-1
    I=0
        VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-(3.d0*V(I,J,K)-3.d0*V(I+1,J,K)+V(I+2,J,K))) &
                    +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-(3.d0*U(I,J,K)-3.d0*U(I+1,J,K)+U(I+2,J,K))) &
                    -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))
    I=NI
        VOR(I,J,K) = XIX(I,J,K)*0.5D0*((3.d0*V(I,J,K)-3.d0*V(I-1,J,K)+V(I-2,J,K))-V(I-1,J,K)) &
                    +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                    -XIY(I,J,K)*0.5D0*((3.d0*U(I,J,K)-3.d0*U(I-1,J,K)+U(I-2,J,K))-U(I-1,J,K)) &
                    -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))
                    
    ENDDO
    ENDDO
    
    DO K=2,NKM1
    DO I=0,NI
    J=0 !branch cut region
        VOR(I,J,K) = VOR(I,NJ-1,K)
!        VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
!                    +ETX(I,J,K)*0.5D0*(V(I,NJ,K)-V(I,NJ-2,K)) &
!                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
!                    -ETY(I,J,K)*0.5D0*(U(I,NJ,K)-U(I,NJ-2,K))
    J=NJ !branch cut region
        VOR(I,J,K) = VOR(I,1,K)
!        VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
!                    +ETX(I,J,K)*0.5D0*(V(I,2,K)-V(I,0,K)) &
!                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
!                    -ETY(I,J,K)*0.5D0*(U(I,2,K)-U(I,0,K))
                    
    ENDDO
    ENDDO

    END SUBROUTINE CALC_VOR