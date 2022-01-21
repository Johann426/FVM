	SUBROUTINE P_SOURCE(IN)

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
!------------------------------------------------------------
!   Define Argument variables
!------------------------------------------------------------
   	DOUBLE PRECISION, INTENT(IN)  :: IN(SX:EX,SY:EY,SZ:EZ)
	
!--------------------------------------------------------------------------------------------------------
!   PRESSURE GRADIENT
	S_U(SX1:EX1,SY1:EY1,SZ1:EZ1) = -0.5D0/DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)									&
	*(																									&
	+ DJAC(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*XIX(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*IN(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)&
	- DJAC(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*XIX(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*IN(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)&
	+ DJAC(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*ETX(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*IN(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)&
	- DJAC(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*ETX(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*IN(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)&
	+ DJAC(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*ZTX(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*IN(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)&
	- DJAC(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*ZTX(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*IN(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)&
	)

	S_V(SX1:EX1,SY1:EY1,SZ1:EZ1) = -0.5D0/DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)									&
	*(																									&
	+ DJAC(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*XIY(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*IN(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)&
	- DJAC(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*XIY(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*IN(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)&
	+ DJAC(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*ETY(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*IN(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)&
	- DJAC(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*ETY(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*IN(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)&
	+ DJAC(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*ZTY(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*IN(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)&
	- DJAC(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*ZTY(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*IN(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)&
    )
    
    S_W(SX1:EX1,SY1:EY1,SZ1:EZ1) = -0.5D0/DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)									&
	*(																									&
	+ DJAC(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*XIZ(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*IN(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)&
	- DJAC(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*XIZ(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*IN(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)&
	+ DJAC(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*ETZ(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*IN(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)&
	- DJAC(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*ETZ(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*IN(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)&
	+ DJAC(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*ZTZ(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*IN(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)&
	- DJAC(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*ZTZ(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*IN(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)&
    )
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE P_SOURCE