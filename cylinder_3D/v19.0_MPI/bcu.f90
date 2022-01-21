	SUBROUTINE BCU(A,B)

	USE COMDAT_SHARED
    
	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	DOUBLE PRECISION, INTENT(INOUT) :: A(IS:IE,JS:JE,KS:KE)
	DOUBLE PRECISION, INTENT(IN)	:: B(IS:IE,JS:JE,KS:KE)
	    
    CALL MPI_X(A)
    
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	I=IS
		U_F = 0.D0	!VEL_ROT * -SIN(THETA(I,JS:JE,KS:KE))

    AP(I+1,JS:JE,KS:KE) = AP(I+1,JS:JE,KS:KE) +AW(I+1,JS:JE,KS:KE)
    SU(I+1,JS:JE,KS:KE) = SU(I+1,JS:JE,KS:KE) +2.D0*AW(I+1,JS:JE,KS:KE)*U_F
    AW(I+1,JS:JE,KS:KE) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	A(I,JS:JE,KS:KE) = -A(I+1,JS:JE,KS:KE) + 2.0D0*U_F

!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=KS,KE
	DO J=JS,JE
	I=IE
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    U_F = UIN
	    ELSE
	    U_F = 0.5D0*(B(I,J,K)+B(I-1,J,K)) - DTIME*1.D0	&
            *(0.5D0*(B(I,J,K)+B(I-1,J,K))-B(I-1,J,K))	&
            /(0.5D0*(  X(I,J,K)+  X(I-1,J,K)) - X(I-1,J,K))
	    ENDIF
    AP(IE1,J,K) = AP(IE1,J,K) +AE(IE1,J,K)
    SU(IE1,J,K) = SU(IE1,J,K) +2.D0*AE(IE1,J,K)*U_F
    AE(IE1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	A(I,J,K) = -A(I-1,J,K) + 2.0D0*U_F
    ENDDO
	ENDDO
	
	END SUBROUTINE BCU