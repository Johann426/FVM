    SUBROUTINE CALC_VOR
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!----------------------------------------------------------------------------
! Cell center vorticity
!----------------------------------------------------------------------------
    DO K=KS1,KE1
    DO J=JS1,JE1
    DO I=IS1,IE1
    VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
                +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
                -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))
                
    ENDDO
    ENDDO
    ENDDO

    DO K=KS1,KE1
    DO J=JS1,JE1
    I=IS
        VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-(3.d0*V(I,J,K)-3.d0*V(I+1,J,K)+V(I+2,J,K))) &
                    +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-(3.d0*U(I,J,K)-3.d0*U(I+1,J,K)+U(I+2,J,K))) &
                    -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))
    I=IE
		VOR(I,J,K) = XIX(I,J,K)*0.5D0*((3.d0*V(I,J,K)-3.d0*V(I-1,J,K)+V(I-2,J,K))-V(I-1,J,K)) &
					+ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
					-XIY(I,J,K)*0.5D0*((3.d0*U(I,J,K)-3.d0*U(I-1,J,K)+U(I-2,J,K))-U(I-1,J,K)) &
					-ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))
                    
    ENDDO
    ENDDO
    
    DO K=KS1,KE1
    DO I=IS1,IE1
    J=JS !branch cut region
        VOR(I,J,K) = VOR(I,NJ-1,K)
!        VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
!                    +ETX(I,J,K)*0.5D0*(V(I,NJ,K)-V(I,NJ-2,K)) &
!                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
!                    -ETY(I,J,K)*0.5D0*(U(I,NJ,K)-U(I,NJ-2,K))
    J=JE !branch cut region
        VOR(I,J,K) = VOR(I,1,K)
!        VOR(I,J,K) = XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
!                    +ETX(I,J,K)*0.5D0*(V(I,2,K)-V(I,0,K)) &
!                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
!                    -ETY(I,J,K)*0.5D0*(U(I,2,K)-U(I,0,K))
                    
    ENDDO
    ENDDO

    END SUBROUTINE CALC_VOR