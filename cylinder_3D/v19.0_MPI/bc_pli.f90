    SUBROUTINE BC_PLI

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
    DUM_3D=PLI
	
	CALL MPI_X(PLI)
	
!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	I=IS
	! WALL BOUNDARY SURFACE
    !AP(I+1,JS1:JE1,KS:KE) = AP(I+1,JS1:JE1,KS:KE)-AW(I+1,JS1:JE1,KS:KE)
    AW(I+1,JS1:JE1,KS:KE) = 0.0D0
    PLI(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*PLI(I+1,JS1:JE1,KS:KE)-PLI(I+2,JS1:JE1,KS:KE))					&
                +Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))

    I=IE
	! OUTFLOW BOUNDARY SURFACE
    !AP(I-1,JS1:JE1,KS:KE) = AP(I-1,JS1:JE1,KS:KE)-AE(I-1,JS1:JE1,KS:KE)
    AE(I-1,JS1:JE1,KS:KE) = 0.0D0
    PLI(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*PLI(I-1,JS1:JE1,KS:KE)-PLI(I-2,JS1:JE1,KS:KE))					&
 				-Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))

    END SUBROUTINE BC_PLI