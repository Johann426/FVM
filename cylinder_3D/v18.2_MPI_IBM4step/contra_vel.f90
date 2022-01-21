	SUBROUTINE CONTRA_VEL

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
	
	INCLUDE 'mpif.h'
	

	!   CONTRA-VARIENT VELOCITY

    U1(IS:IE,JS:JE,KS:KE)	=XIX(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE)	&
							+XIY(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE)	&
							+XIZ(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
    
    U2(IS:IE,JS:JE,KS:KE)	=ETX(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE)	&
							+ETY(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE)	&
							+ETZ(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
    
    U3(IS:IE,JS:JE,KS:KE)	=ZTX(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE)	&
							+ZTY(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE)	&
							+ZTZ(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)


	!	MPI : COLUMN-DIRECTION BLOCK DIVISION

	CALL MPI_X(U1)
	
	CALL MPI_X(U2)
	
	CALL MPI_X(U3)
	
    !   VOLUME FLUX : CONTRA-VARIENT VELOCITY MULTIPLIED BY J^-1
    
	UE(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U1(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*U1(IS1+1:IE1+1,JS1:JE1,KS1:KE1))
    
    UW(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U1(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*U1(IS1-1:IE1-1,JS1:JE1,KS1:KE1))
    
    VN(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U2(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*U2(IS1:IE1,JS1+1:JE1+1,KS1:KE1))
    
    VS(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U2(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*U2(IS1:IE1,JS1-1:JE1-1,KS1:KE1))
    
    WT(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U3(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*U3(IS1:IE1,JS1:JE1,KS1+1:KE1+1))
    
    WB(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U3(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*U3(IS1:IE1,JS1:JE1,KS1-1:KE1-1))

	!	BOUNDARY CONDITION FOR VOLUME FLUX
	UW(IS1,JS1:JE1,KS1:KE1) = 0.D0
	
	END SUBROUTINE CONTRA_VEL