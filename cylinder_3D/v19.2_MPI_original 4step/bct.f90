	SUBROUTINE BCT

	USE COMDAT_SHARED
    
	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	I=IS
		T_F = T_WALL

    AP(I+1,JS:JE,KS:KE) = AP(I+1,JS:JE,KS:KE) +AW(I+1,JS:JE,KS:KE)
    SU(I+1,JS:JE,KS:KE) = SU(I+1,JS:JE,KS:KE) +2.D0*AW(I+1,JS:JE,KS:KE)*T_F
    AW(I+1,JS:JE,KS:KE) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	T(I,JS:JE,KS:KE) = -T(I+1,JS:JE,KS:KE) + 2.0D0*T_F
	
!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IE
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    T_F = TIN
	    ELSE
	    T_F = 0.5D0*(T_N(I,J,K)+T_N(I-1,J,K)) - DTIME*1.D0										&
				*(XIX(I,J,K)*(T_N(I,J,K)-T_N(I-1,J,K))											&
				+ETX(I,J,K)*0.25D0*(T_N(I,J+1,K)+T_N(I-1,J+1,K)-T_N(I,J-1,K)-T_N(I-1,J-1,K))	&
				+ZTX(I,J,K)*0.25D0*(T_N(I,J,K+1)+T_N(I-1,J,K+1)-T_N(I,J,K-1)-T_N(I-1,J,K-1)))
	    ENDIF
    AP(IE1,J,K) = AP(IE1,J,K) +AE(IE1,J,K)
    SU(IE1,J,K) = SU(IE1,J,K) +2.D0*AE(IE1,J,K)*T_F
    AE(IE1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	T(I,J,K) = -T(I-1,J,K) + 2.0D0*T_F
    ENDDO
	ENDDO
	
	CALL MPI_BARRIER(NMCW,ERR)
	
	CALL MPI_X(T)
	
	END SUBROUTINE BCT