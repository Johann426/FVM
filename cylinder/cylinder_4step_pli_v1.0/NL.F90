    SUBROUTINE NL(IN,OUT)
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!------------------------------------------------------------
!   Define Argument variables
!------------------------------------------------------------
   	DOUBLE PRECISION, INTENT(IN)  :: IN(SX:EX,SY:EY,2:NKM1)
	DOUBLE PRECISION, INTENT(OUT) :: OUT(SX:EX,SY:EY,2:NKM1)

    OUT(sx1:ex1,sy1:ey1,2:NKM1)=0.5D0*                                      &
	(                                                                       &
    +DJAC(sx1+1:ex1+1,sy1:ey1,2:NKM1)*U1(sx1+1:ex1+1,sy1:ey1,2:NKM1)*IN(sx1+1:ex1+1,sy1:ey1,2:NKM1)    &
	-DJAC(sx1-1:ex1-1,sy1:ey1,2:NKM1)*U1(sx1-1:ex1-1,sy1:ey1,2:NKM1)*IN(sx1-1:ex1-1,sy1:ey1,2:NKM1)    &
	+DJAC(sx1:ex1,sy1+1:ey1+1,2:NKM1)*U2(sx1:ex1,sy1+1:ey1+1,2:NKM1)*IN(sx1:ex1,sy1+1:ey1+1,2:NKM1)    &
	-DJAC(sx1:ex1,sy1-1:ey1-1,2:NKM1)*U2(sx1:ex1,sy1-1:ey1-1,2:NKM1)*IN(sx1:ex1,sy1-1:ey1-1,2:NKM1)    &
	)
    END SUBROUTINE NL