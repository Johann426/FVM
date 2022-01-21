    SUBROUTINE NL(IN,OUT)
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!------------------------------------------------------------
!   Define Argument variables
!------------------------------------------------------------
   	DOUBLE PRECISION, INTENT(IN)  :: IN(SX:EX,SY:EY,2:NKM1)
	DOUBLE PRECISION, INTENT(OUT) :: OUT(SX:EX,SY:EY,2:NKM1)

    OUT(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0/DJAC(SX1:EX1,SY1:EY1,2:NKM1)*                               &
	(                                                                                               &
    +DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*U1(SX1+1:EX1+1,SY1:EY1,2:NKM1)*IN(SX1+1:EX1+1,SY1:EY1,2:NKM1) &
	-DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*U1(SX1-1:EX1-1,SY1:EY1,2:NKM1)*IN(SX1-1:EX1-1,SY1:EY1,2:NKM1) &
	+DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*U2(SX1:EX1,SY1+1:EY1+1,2:NKM1)*IN(SX1:EX1,SY1+1:EY1+1,2:NKM1) &
	-DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*U2(SX1:EX1,SY1-1:EY1-1,2:NKM1)*IN(SX1:EX1,SY1-1:EY1-1,2:NKM1) &
	)
    END SUBROUTINE NL