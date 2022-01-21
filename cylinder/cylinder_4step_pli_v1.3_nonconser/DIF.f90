	SUBROUTINE DIF(IN,OUT)

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
!------------------------------------------------------------
!   Define Argument variables
!------------------------------------------------------------
	DOUBLE PRECISION, INTENT(IN)  :: IN(SX:EX,SY:EY,2:NKM1)
	DOUBLE PRECISION, INTENT(OUT) :: OUT(SX:EX,SY:EY,2:NKM1)
        
    CALL DIF_DIAGONAL(IN,DD)
    
	OUT(sx1:ex1,sy1:ey1,2:NKM1)=1.0D0/RE                            							    &
		*(DE(sx1:ex1,sy1:ey1,2:NKM1)*(IN(sx1+1:ex1+1,sy1:ey1,2:NKM1)-IN(sx1:ex1,sy1:ey1,2:NKM1))    &
		+ DW(sx1:ex1,sy1:ey1,2:NKM1)*(IN(sx1-1:ex1-1,sy1:ey1,2:NKM1)-IN(sx1:ex1,sy1:ey1,2:NKM1))    &
		+ DN(sx1:ex1,sy1:ey1,2:NKM1)*(IN(sx1:ex1,sy1+1:ey1+1,2:NKM1)-IN(sx1:ex1,sy1:ey1,2:NKM1))    &
		+ DS(sx1:ex1,sy1:ey1,2:NKM1)*(IN(sx1:ex1,sy1-1:ey1-1,2:NKM1)-IN(sx1:ex1,sy1:ey1,2:NKM1))    &
!		+ DT(sx1:ex1,sy1:ey1,2:NKM1)*(IN(sx1:ex1,sy1:ey1,2+1:NKM1+1)-IN(sx1:ex1,sy1:ey1,2:NKM1))    &
!		+ DB(sx1:ex1,sy1:ey1,2:NKM1)*(IN(sx1:ex1,sy1:ey1,2-1:NKM1-1)-IN(sx1:ex1,sy1:ey1,2:NKM1))    &
		+ DD                                                                                        &
		)
		
		
!------------------------------------------------------------	
	END SUBROUTINE DIF
!------------------------------------------------------------


