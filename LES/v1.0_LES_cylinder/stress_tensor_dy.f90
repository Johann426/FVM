	SUBROUTINE STRESS_TENSOR_DY
	
	USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE

	INCLUDE 'MPIF.H'
	
	DOUBLE PRECISION :: INC

	DOUBLE PRECISION :: C_AVG
	DOUBLE PRECISION :: LM_AVG
	DOUBLE PRECISION :: MM_AVG

	DOUBLE PRECISION :: T_LM_AVG
	DOUBLE PRECISION :: T_MM_AVG
	DOUBLE PRECISION :: T_INC
	
	FILT_WID(IS:IE,JS:JE,KS:KE) = DJAC(IS:IE,JS:JE,KS:KE)**(1.0D0/3.0D0)
		
	!	THE LARGE-SCALE STRAIN-RATE TENSOR
	LSS11(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIX(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - U(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETX(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - U(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTX(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - U(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS11(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*LSS11(IS1:IE1,JS1:JE1,KS1:KE1)
	
	LSS12(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIY(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - U(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETY(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - U(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTY(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - U(IS1:IE1,JS1:JE1,KS1-1:KE1-1))	&
		+ XIX(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - V(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETX(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - V(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTX(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - V(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS12(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*LSS12(IS1:IE1,JS1:JE1,KS1:KE1)
	
	LSS13(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIZ(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - U(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETZ(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - U(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)*(U(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - U(IS1:IE1,JS1:JE1,KS1-1:KE1-1))	&
		+ XIX(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - W(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETX(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - W(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTX(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - W(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS13(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*LSS13(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS21(IS1:IE1,JS1:JE1,KS1:KE1) = LSS12(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS22(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIY(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - V(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETY(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - V(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTY(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - V(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS22(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*LSS22(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS23(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIZ(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - V(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETZ(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - V(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)*(V(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - V(IS1:IE1,JS1:JE1,KS1-1:KE1-1))	&
		+ XIY(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - W(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETY(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - W(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTY(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - W(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS23(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*LSS23(IS1:IE1,JS1:JE1,KS1:KE1)
	
	LSS31(IS1:IE1,JS1:JE1,KS1:KE1) = LSS13(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS32(IS1:IE1,JS1:JE1,KS1:KE1) = LSS23(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS33(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIZ(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - W(IS1-1:IE1-1,JS1:JE1,KS1:KE1))	&
		+ ETZ(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - W(IS1:IE1,JS1-1:JE1-1,KS1:KE1))	&
		+ ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)*(W(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - W(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS33(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*LSS33(IS1:IE1,JS1:JE1,KS1:KE1)
	
	I=IS
	LSS11(I,JS:JE,KS:KE) = 3.D0*LSS11(I+1,JS:JE,KS:KE)-3.D0*LSS11(I+2,JS:JE,KS:KE)+LSS11(I+3,JS:JE,KS:KE)
	I=IE
	LSS11(I,JS:JE,KS:KE) = 3.D0*LSS11(I-1,JS:JE,KS:KE)-3.D0*LSS11(I-2,JS:JE,KS:KE)+LSS11(I-3,JS:JE,KS:KE)
	J=JS
	LSS11(IS:IE,J,KS:KE) = 3.D0*LSS11(IS:IE,J+1,KS:KE)-3.D0*LSS11(IS:IE,J+2,KS:KE)+LSS11(IS:IE,J+3,KS:KE)
	J=JE
	LSS11(IS:IE,J,KS:KE) = 3.D0*LSS11(IS:IE,J-1,KS:KE)-3.D0*LSS11(IS:IE,J-2,KS:KE)+LSS11(IS:IE,J-3,KS:KE)
	K=KS
	LSS11(IS:IE,JS:JE,K) = 3.D0*LSS11(IS:IE,JS:JE,K+1)-3.D0*LSS11(IS:IE,JS:JE,K+2)+LSS11(IS:IE,JS:JE,K+3)
	K=KE
	LSS11(IS:IE,JS:JE,K) = 3.D0*LSS11(IS:IE,JS:JE,K-1)-3.D0*LSS11(IS:IE,JS:JE,K-2)+LSS11(IS:IE,JS:JE,K-3)
	
!	DO K=KS1,KE1
!	DO J=JS1,JE1
!	I=0
!		LSS11(I,J,K) = 0.5D0*(2.0D0*XIX(I,J,K)*(U(I+1,J,K) - U(I,J,K))		&
!						           +ETX(I,J,K)*(U(I,J+1,K) - U(I,J-1,K))	&
!							       +ZTX(I,J,K)*(U(I,J,K+1) - U(I,J,K-1)) )
!	
!		LSS12(I,J,K) = 0.25D0*(2.0D0*XIY(I,J,K)*(U(I+1,J,K) - U(I,J,K))		&
!									+ETY(I,J,K)*(U(I,J+1,K) - U(I,J-1,K))	&
!									+ZTY(I,J,K)*(U(I,J,K+1) - U(I,J,K-1))	& 
!							  +2.0D0*XIX(I,J,K)*(V(I+1,J,K) - V(I,J,K))		&
!									+ETX(I,J,K)*(V(I,J+1,K) - V(I,J-1,K))	&
!									+ZTX(I,J,K)*(V(I,J,K+1) - V(I,J,K-1)) )
!
!		LSS13(I,J,K) = 0.25D0*(2.0D0*XIZ(I,J,K)*(U(I+1,J,K) - U(I,J,K))		&
!									+ETZ(I,J,K)*(U(I,J+1,K) - U(I,J-1,K))	&
!									+ZTZ(I,J,K)*(U(I,J,K+1) - U(I,J,K-1))	&
!							  +2.0D0*XIX(I,J,K)*(W(I+1,J,K) - W(I,J,K))		&
!									+ETX(I,J,K)*(W(I,J+1,K) - W(I,J-1,K))	&
!									+ZTX(I,J,K)*(W(I,J,K+1) - W(I,J,K-1)) )
!
!		LSS21(I,J,K) = LSS12(I,J,K)
!
!		LSS22(I,J,K) = 0.5D0*(2.0D0*XIY(I,J,K)*(V(I+1,J,K) - V(I,J,K))		&
!								   +ETY(I,J,K)*(V(I,J+1,K) - V(I,J-1,K))	&
!							       +ZTY(I,J,K)*(V(I,J,K+1) - V(I,J,K-1)) )
!
!		LSS23(I,J,K) = 0.25D0*(2.0D0*XIZ(I,J,K)*(V(I+1,J,K) - V(I,J,K))		&
!									+ETZ(I,J,K)*(V(I,J+1,K) - V(I,J-1,K))	&
!									+ZTZ(I,J,K)*(V(I,J,K+1) - V(I,J,K-1))	&
!							  +2.0D0*XIY(I,J,K)*(W(I+1,J,K) - W(I,J,K))		&
!									+ETY(I,J,K)*(W(I,J+1,K) - W(I,J-1,K))	&
!									+ZTY(I,J,K)*(W(I,J,K+1) - W(I,J,K-1)) )
!
!		LSS31(I,J,K) = LSS13(I,J,K)
!
!		LSS32(I,J,K) = LSS23(I,J,K)
!
!		LSS33(I,J,K) = 0.5D0*(2.0D0*XIZ(I,J,K)*(W(I+1,J,K) - W(I,J,K))		&
!								   +ETZ(I,J,K)*(W(I,J+1,K) - W(I,J-1,K))	&
!								   +ZTZ(I,J,K)*(W(I,J,K+1) - W(I,J,K-1)) )
!								   
!	I=NI
!		LSS11(I,J,K) = 0.5D0*(2.0D0*XIX(I,J,K)*(U(I,J,K)   - U(I-1,J,K))	&
!						           +ETX(I,J,K)*(U(I,J+1,K) - U(I,J-1,K))	&
!							       +ZTX(I,J,K)*(U(I,J,K+1) - U(I,J,K-1)) )
!	
!		LSS12(I,J,K) = 0.25D0*(2.0D0*XIY(I,J,K)*(U(I,J,K)   - U(I-1,J,K))	&
!									+ETY(I,J,K)*(U(I,J+1,K) - U(I,J-1,K))	&
!									+ZTY(I,J,K)*(U(I,J,K+1) - U(I,J,K-1))	&
!							  +2.0D0*XIX(I,J,K)*(V(I,J,K)   - V(I-1,J,K))	&
!									+ETX(I,J,K)*(V(I,J+1,K) - V(I,J-1,K))	&
!									+ZTX(I,J,K)*(V(I,J,K+1) - V(I,J,K-1)) )
!
!		LSS13(I,J,K) = 0.25D0*(2.0D0*XIZ(I,J,K)*(U(I,J,K)   - U(I-1,J,K))	&
!									+ETZ(I,J,K)*(U(I,J+1,K) - U(I,J-1,K))	&
!							        +ZTZ(I,J,K)*(U(I,J,K+1) - U(I,J,K-1))	&
!							  +2.0D0*XIX(I,J,K)*(W(I,J,K)   - W(I-1,J,K))	&
!									+ETX(I,J,K)*(W(I,J+1,K) - W(I,J-1,K))	&
!									+ZTX(I,J,K)*(W(I,J,K+1) - W(I,J,K-1)) )
!
!		LSS21(I,J,K) = LSS12(I,J,K)
!
!		LSS22(I,J,K) = 0.5D0*(2.0D0*XIY(I,J,K)*(V(I,J,K)   - V(I-1,J,K))	&
!								   +ETY(I,J,K)*(V(I,J+1,K) - V(I,J-1,K))	&
!							       +ZTY(I,J,K)*(V(I,J,K+1) - V(I,J,K-1)) )
!
!		LSS23(I,J,K) = 0.25D0*(2.0D0*XIZ(I,J,K)*(V(I,J,K)   - V(I-1,J,K))	&
!									+ETZ(I,J,K)*(V(I,J+1,K) - V(I,J-1,K))	&
!									+ZTZ(I,J,K)*(V(I,J,K+1) - V(I,J,K-1))	&
!							  +2.0D0*XIY(I,J,K)*(W(I,J,K)   - W(I-1,J,K))	&
!									+ETY(I,J,K)*(W(I,J+1,K) - W(I,J-1,K))	&
!									+ZTY(I,J,K)*(W(I,J,K+1) - W(I,J,K-1)) )
!
!		LSS31(I,J,K) = LSS13(I,J,K)
!
!		LSS32(I,J,K) = LSS23(I,J,K)
!
!		LSS33(I,J,K) = 0.5D0*(2.0D0*XIZ(I,J,K)*(W(I,J,K)   - W(I-1,J,K)) &
!								   +ETZ(I,J,K)*(W(I,J+1,K) - W(I,J-1,K)) &
!								   +ZTZ(I,J,K)*(W(I,J,K+1) - W(I,J,K-1)) )
!	ENDDO
!	ENDDO

	CALL MPI_X(LSS11)
	CALL MPI_X(LSS12)
	CALL MPI_X(LSS13)
	
	CALL MPI_X(LSS21)
	CALL MPI_X(LSS22)
	CALL MPI_X(LSS23)
	
	CALL MPI_X(LSS31)
	CALL MPI_X(LSS32)
	CALL MPI_X(LSS33)


	!   THE MAGNITUDE OF LARGE-SCALE STRAIN-RATE TENSOR

	MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) = DSQRT(2.D0*(													&
		+ LSS11(IS:IE,JS:JE,KS:KE)**2 + LSS12(IS:IE,JS:JE,KS:KE)**2	+ LSS13(IS:IE,JS:JE,KS:KE)**2	&
		+ LSS21(IS:IE,JS:JE,KS:KE)**2 + LSS22(IS:IE,JS:JE,KS:KE)**2	+ LSS23(IS:IE,JS:JE,KS:KE)**2	&
		+ LSS31(IS:IE,JS:JE,KS:KE)**2 + LSS32(IS:IE,JS:JE,KS:KE)**2	+ LSS33(IS:IE,JS:JE,KS:KE)**2 ) )
		
	CALL TEST_FILTER(U,U_HAT)
	CALL TEST_FILTER(V,V_HAT)
	CALL TEST_FILTER(W,W_HAT)
    
    UIN = 1.D0
    CALL BCU(U_HAT,U_N)
	UIN = 0.D0
    CALL BCU(U_HAT,V_N)
    CALL BCU(U_HAT,W_N)
	
	!	TEST-FILTERED RESOLVED STRAIN RATE TENSOR 
	LSS11_HAT(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIX(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - U_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETX(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - U_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTX(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - U_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS11_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*LSS11_HAT(IS1:IE1,JS1:JE1,KS1:KE1)
	
	LSS12_HAT(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIY(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - U_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETY(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - U_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTY(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - U_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))&
		+ XIX(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - V_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETX(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - V_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTX(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - V_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS12_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*LSS12_HAT(IS1:IE1,JS1:JE1,KS1:KE1)
	
	LSS13_HAT(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIZ(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - U_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETZ(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - U_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)*(U_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - U_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))&
		+ XIX(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - W_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETX(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - W_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTX(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - W_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS13_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*LSS13_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS21_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = LSS12_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS22_HAT(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIY(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - V_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETY(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - V_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTY(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - V_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS22_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*LSS22_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS23_HAT(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIZ(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - V_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETZ(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - V_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)*(V_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - V_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))&
		+ XIY(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - W_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETY(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - W_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTY(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - W_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS23_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*LSS23_HAT(IS1:IE1,JS1:JE1,KS1:KE1)
	
	LSS31_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = LSS13_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS32_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = LSS23_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	LSS33_HAT(IS1:IE1,JS1:JE1,KS1:KE1)																			&
		= XIZ(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1+1:IE1+1,JS1:JE1,KS1:KE1) - W_HAT(IS1-1:IE1-1,JS1:JE1,KS1:KE1))&
		+ ETZ(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1:IE1,JS1+1:JE1+1,KS1:KE1) - W_HAT(IS1:IE1,JS1-1:JE1-1,KS1:KE1))&
		+ ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)*(W_HAT(IS1:IE1,JS1:JE1,KS1+1:KE1+1) - W_HAT(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
	LSS33_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*LSS33_HAT(IS1:IE1,JS1:JE1,KS1:KE1)
	
	MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) = DSQRT(2.D0*(														&
		+ LSS11_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2 + LSS12_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2	+ LSS13_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
		+ LSS21_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2 + LSS22_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2	+ LSS23_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
		+ LSS31_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2 + LSS32_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2	+ LSS33_HAT(IS1:IE1,JS1:JE1,KS1:KE1)**2 ) )
	
!------------------------------------------------------------
!	COMPUTE L11
!------------------------------------------------------------

	DUM_IN(IS:IE,JS:JE,KS:KE) =	U(IS:IE,JS:JE,KS:KE) * U(IS:IE,JS:JE,KS:KE)
	
	CALL TEST_FILTER(DUM_IN,DUM_OUT)
	CALL TEST_FILTER(U,DUM_IN)

	DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) = DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)**2
	
	L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   = DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE M11
!------------------------------------------------------------
	
	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			4.0D0 * MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) * LSS11_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) * LSS11(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =		&
			FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) * FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1)*	&
			(M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1))

!------------------------------------------------------------
!  COMPUTE	L11*M11 
!	AND		M11*M11
!------------------------------------------------------------

	LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			L_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

	MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE L22
!------------------------------------------------------------

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			V(IS:IE,JS:JE,KS:KE) * V(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)
	CALL TEST_FILTER(V,DUM_IN)

	DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) = DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)**2
	
	L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   = DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE M22
!------------------------------------------------------------

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			4.0D0 * MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) * LSS22_HAT(IS1:IE1,JS1:JE1,KS1:KE1)
	
	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) * LSS22(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =		&
			FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) * FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1)*	&
			(M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1))

!------------------------------------------------------------
!  COMPUTE	L11*M11 + L22*M22 
!	AND		M11*M11 + M22*M22
!------------------------------------------------------------

	LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			L_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

	MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE L33
!------------------------------------------------------------

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			W(IS:IE,JS:JE,KS:KE) * W(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)
	CALL TEST_FILTER(W,DUM_IN)

	DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) = DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)**2
	
	L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   = DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)
			
!------------------------------------------------------------
!  COMPUTE M33
!------------------------------------------------------------

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			4.0D0 * MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) * LSS33_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) * LSS33(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) * FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) *		&
			(M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1))

!------------------------------------------------------------
!  COMPUTE	L11*M11 + L22*M22 + L33*M33 
!	AND		M11*M11 + M22*M22 + M33*M33
!------------------------------------------------------------

	LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1)  =	&
			LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			L_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)
		
	MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1)  =	&
			MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE L12
!------------------------------------------------------------

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			U(IS:IE,JS:JE,KS:KE) * V(IS:IE,JS:JE,KS:KE)
	
	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	CALL TEST_FILTER(U,DUM_IN)
	CALL TEST_FILTER(V,L_IJ)

	DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) * L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)
	
	L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
		DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE M12
!------------------------------------------------------------

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			4.0D0 * MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) * LSS12_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) * LSS12(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) * FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) *	&
			(M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1))

!------------------------------------------------------------
!  COMPUTE	L11*M11 + L22*M22 + L33*M33 + 2*L12*M12
!    AND	M11*M11 + M22*M22 + M33*M33 + 2*M12*M12
!------------------------------------------------------------

	LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			2.0D0 * L_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)
	
	MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			2.0D0 * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE L13
!------------------------------------------------------------

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			U(IS:IE,JS:JE,KS:KE) * W(IS:IE,JS:JE,KS:KE)
	
	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	CALL TEST_FILTER(U,DUM_IN)
	CALL TEST_FILTER(W,L_IJ)

	DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)  * L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

	L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE M13
!------------------------------------------------------------

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =		&
			4.0D0 * MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) * LSS13_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) * LSS13(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)
	
	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =		&
			FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) * FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1)* &
			(M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1))

!------------------------------------------------------------
!  COMPUTE	L11*M11 + L22*M22 + L33*M33 + 2*L12*M12 + 2*L13*M13
!	AND		M11*M11 + M22*M22 + M33*M33 + 2*M12*M12 + 2*M13*M13
!------------------------------------------------------------

	LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			2.0D0 * L_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

	MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			2.0D0 * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE L23
!------------------------------------------------------------

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			V(IS:IE,JS:JE,KS:KE) * W(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	CALL TEST_FILTER(V,DUM_IN)
	CALL TEST_FILTER(W,L_IJ)

	DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1) * L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

	L_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_IN(IS1:IE1,JS1:JE1,KS1:KE1)

!------------------------------------------------------------
!  COMPUTE M23
!------------------------------------------------------------

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			4.0D0 * MAGNITUDE_LSS_HAT(IS1:IE1,JS1:JE1,KS1:KE1) * LSS23_HAT(IS1:IE1,JS1:JE1,KS1:KE1)

	DUM_IN(IS:IE,JS:JE,KS:KE) =	&
			MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE) * LSS23(IS:IE,JS:JE,KS:KE)

	CALL TEST_FILTER(DUM_IN,DUM_OUT)

	M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)   =	&
			FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) * FILT_WID(IS1:IE1,JS1:JE1,KS1:KE1) *	&
			(M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) - DUM_OUT(IS1:IE1,JS1:JE1,KS1:KE1))


!------------------------------------------------------------
!  COMPUTE	L11*M11 + L22*M22 + L33*M33 + 2*L12*M12 + 2*L13*M13 + 2*L23*M23
!	AND		M11*M11 + M22*M22 + M33*M33 + 2*M12*M12 + 2*M13*M13 + 2*M23*M23
!------------------------------------------------------------

	LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			LM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			2.0D0 * L_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)

	MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) =	&
			MM_IJ(IS1:IE1,JS1:JE1,KS1:KE1) +	&
			2.0D0 * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1) * M_IJ(IS1:IE1,JS1:JE1,KS1:KE1)
	
	!	DYNAMIC (TIME-DEPENDENT,VOLUME-AVERAGED) EDDY VISCOSITY
	C_AVG  = 0.0D0

	LM_AVG = 0.0D0
	MM_AVG = 0.0D0
	INC    = 0.0D0

	DO K = KS1,KE1
	DO J = JS1,JE1
	DO I = IS1,IE1

		IF(LM_IJ(I,J,K) .LT. 0.0D0) THEN
			LM_AVG = LM_AVG + LM_IJ(I,J,K)
			MM_AVG = MM_AVG + MM_IJ(I,J,K)

			INC = INC + 1.0D0
		ENDIF

	ENDDO
	ENDDO
	ENDDO
	
	CALL MPI_ALLREDUCE(LM_AVG,T_LM_AVG,1,MPI_DOUBLE_PRECISION, MPI_SUM,NMCW,ERR)
	CALL MPI_ALLREDUCE(MM_AVG,T_MM_AVG,1,MPI_DOUBLE_PRECISION, MPI_SUM,NMCW,ERR)
	CALL MPI_ALLREDUCE(INC,T_INC,1,MPI_DOUBLE_PRECISION, MPI_SUM,NMCW,ERR)

	T_LM_AVG = T_LM_AVG / T_INC
	T_MM_AVG = T_MM_AVG / T_INC

	IF(T_MM_AVG.NE.0.D0) THEN
		C_AVG = -0.5D0 * T_LM_AVG/T_MM_AVG
	ENDIF
	
	!	THE EDDY VISCOSITY
	EDDY_VISCOS(IS:IE,JS:JE,KS:KE) = C_AVG*FILT_WID(IS:IE,JS:JE,KS:KE)**2*MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE)

	
!	SGS11 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS11(I,J,K)
!	SGS12 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS12(I,J,K)
!	SGS13 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS13(I,J,K)
!	
!	SGS21 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS21(I,J,K)
!	SGS22 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS22(I,J,K)
!	SGS23 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS23(I,J,K)
!
!	SGS31 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS31(I,J,K)
!	SGS32 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS32(I,J,K)
!	SGS33 =	-2.D0*EDDY_VISCOS(I,J,K)*LSS33(I,J,K)
	
	!	DIVERSENCE OF SUBGRID-SCALE STRESS TENSOR
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
	
	SGS_X(I,J,K) = 0.5D0*(	&
		+ XIX(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS11(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS11(I-1,J,K))	&
		+ ETX(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS11(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS11(I,J-1,K))	&
		+ ZTX(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS11(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS11(I,J,K-1))	&
	
		+ XIY(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS12(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS12(I-1,J,K))	&
		+ ETY(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS12(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS12(I,J-1,K))	&
		+ ZTY(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS12(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS12(I,J,K-1))	&
						
		+ XIZ(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS13(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS13(I-1,J,K))	&
		+ ETZ(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS13(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS13(I,J-1,K))	&
		+ ZTZ(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS13(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS13(I,J,K-1))	)

	SGS_Y(I,J,K) = 0.5D0*(	&
		+ XIX(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS21(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS21(I-1,J,K))	&
		+ ETX(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS21(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS21(I,J-1,K))	&
		+ ZTX(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS21(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS21(I,J,K-1))	&
	
		+ XIY(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS22(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS22(I-1,J,K))	&
		+ ETY(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS22(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS22(I,J-1,K))	&
		+ ZTY(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS22(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS22(I,J,K-1))	&
						
		+ XIZ(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS23(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS23(I-1,J,K))	&
		+ ETZ(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS23(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS23(I,J-1,K))	&
		+ ZTZ(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS23(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS23(I,J,K-1))	)
		          
	SGS_Z(I,J,K) = 0.5D0*(	&
		+ XIX(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS31(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS31(I-1,J,K))	&
		+ ETX(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS31(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS31(I,J-1,K))	&
		+ ZTX(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS31(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS31(I,J,K-1))	&

		+ XIY(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS32(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS32(I-1,J,K))	&
		+ ETY(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS32(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS32(I,J-1,K))	&
		+ ZTY(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS32(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS32(I,J,K-1))	&
					
		+ XIZ(I,J,K)*(-2.D0*EDDY_VISCOS(I+1,J,K)*LSS33(I+1,J,K) +2.D0*EDDY_VISCOS(I-1,J,K)*LSS33(I-1,J,K))	&
		+ ETZ(I,J,K)*(-2.D0*EDDY_VISCOS(I,J+1,K)*LSS33(I,J+1,K) +2.D0*EDDY_VISCOS(I,J-1,K)*LSS33(I,J-1,K))	&
		+ ZTZ(I,J,K)*(-2.D0*EDDY_VISCOS(I,J,K+1)*LSS33(I,J,K+1) +2.D0*EDDY_VISCOS(I,J,K-1)*LSS33(I,J,K-1))	)
	
	ENDDO
	ENDDO
	ENDDO
	
	END SUBROUTINE STRESS_TENSOR_DY