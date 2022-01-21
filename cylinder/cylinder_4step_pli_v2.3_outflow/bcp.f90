    SUBROUTINE BCP

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	DUM_3D=P
!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	DO K=2,NKM1
	DO J=1,NJM1
	I=0
    P(I,J,K) = 1.D0/3.D0*((4.D0*P(I+1,J,k)-P(I+2,J,K))         &
             + Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))
    
    I=NI
    P(I,J,K) = 1.D0/3.D0*((4.D0*P(I-1,J,k)-P(I-2,J,K))         &
             - Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))
	ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=2,NKM1
	DO I=0,NI
    J=0
	    P(I,J,K) = P(I,NJM1,K)
	J=NJ
		P(I,J,K) = P(I,1,K)
	ENDDO
	ENDDO
	
    END SUBROUTINE BCP
