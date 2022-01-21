    SUBROUTINE BCP

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	INCLUDE 'mpif.h'
	
	DUM_3D=P

	CALL MPI_X(P)
	
!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	I=IS
    P(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*P(I+1,JS1:JE1,KS:KE)-P(I+2,JS1:JE1,KS:KE))         &
             + Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))
    
    I=IE
    P(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*P(I-1,JS1:JE1,KS:KE)-P(I-2,JS1:JE1,KS:KE))         &
             - Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))

    END SUBROUTINE BCP