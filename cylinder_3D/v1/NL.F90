    SUBROUTINE NL(IN,OUT)
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!------------------------------------------------------------
!   Define Argument variables
!------------------------------------------------------------
   	DOUBLE PRECISION, INTENT(IN)  :: IN(SX:EX,SY:EY,SZ:EZ)
	DOUBLE PRECISION, INTENT(OUT) :: OUT(SX:EX,SY:EY,SZ:EZ)

    OUT(SX1:EX1,SY1:EY1,SZ1:EZ1)=0.5D0/DJAC(SX1:EX1,SY1:EY1,2:NKM1)*							&
	(                                                                                           &
    +DJAC(1+1:NIM1+1,1:NJM1,1:NKM1)*U1(1+1:NIM1+1,1:NJM1,1:NKM1)*IN(1+1:NIM1+1,1:NJM1,1:NKM1)   &
	-DJAC(1-1:NIM1-1,1:NJM1,1:NKM1)*U1(1-1:NIM1-1,1:NJM1,1:NKM1)*IN(1-1:NIM1-1,1:NJM1,1:NKM1)   &
	+DJAC(1:NIM1,1+1:NJM1+1,1:NKM1)*U2(1:NIM1,1+1:NJM1+1,1:NKM1)*IN(1:NIM1,1+1:NJM1+1,1:NKM1)   &
	-DJAC(1:NIM1,1-1:NJM1-1,1:NKM1)*U2(1:NIM1,1-1:NJM1-1,1:NKM1)*IN(1:NIM1,1-1:NJM1-1,1:NKM1)   &
	+DJAC(1:NIM1,1:NJM1,1+1:NKM1+1)*U3(1:NIM1,1:NJM1,1+1:NKM1+1)*IN(1:NIM1,1:NJM1,1+1:NKM1+1)   &
	-DJAC(1:NIM1,1:NJM1,1-1:NKM1-1)*U3(1:NIM1,1:NJM1,1-1:NKM1-1)*IN(1:NIM1,1:NJM1,1-1:NKM1-1)   &
	)
    END SUBROUTINE NL