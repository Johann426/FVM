    SUBROUTINE BCP

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
    DUM_3D=P
	
!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	I=IS
	! WALL BOUNDARY SURFACE
!!	AP(I+1,JS1:JE1,KS:KE) = AP(I+1,JS1:JE1,KS:KE)-AW(I+1,JS1:JE1,KS:KE)
!	AW(I+1,JS1:JE1,KS:KE) = 0.0D0
!!	P(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*P(I+1,JS1:JE1,KS:KE)-P(I+2,JS1:JE1,KS:KE))					&
!!			+Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))
	
	AP(I+1,JS1:JE1,KS:KE) = AP(I+1,JS1:JE1,KS:KE)-AW(I+1,JS1:JE1,KS:KE)
    AW(I+1,JS1:JE1,KS:KE) = 0.0D0
	P(I,JS1:JE1,KS:KE) = P(I+1,JS1:JE1,KS:KE)
	
    I=IE
	! OUTFLOW BOUNDARY SURFACE
!	AP(I-1,JS1:JE1,KS:KE) = AP(I-1,JS1:JE1,KS:KE)-AE(I-1,JS1:JE1,KS:KE)
    AE(I-1,JS1:JE1,KS:KE) = 0.0D0
!	P(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*P(I-1,JS1:JE1,KS:KE)-P(I-2,JS1:JE1,KS:KE))					&
!				-Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))
	P(I,JS1:JE1,KS:KE) = 3.D0*P(I-1,JS1:JE1,KS:KE)-3.D0*P(I-2,JS1:JE1,KS:KE)+P(I-3,JS1:JE1,KS:KE)
	
	CALL MPI_BARRIER(NMCW,ERR)
		
	CALL MPI_X(P)
	
    END SUBROUTINE BCP