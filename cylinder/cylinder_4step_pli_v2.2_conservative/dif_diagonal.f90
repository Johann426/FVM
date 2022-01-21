    SUBROUTINE DIF_DIAGONAL(INN,OUT)
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!------------------------------------------------------------
!   Define Argument variables
!------------------------------------------------------------
	DOUBLE PRECISION, INTENT(IN)  :: INN(SX:EX,SY:EY,2:NKM1)
	DOUBLE PRECISION, INTENT(OUT) :: OUT(SX:EX,SY:EY,2:NKM1)
    
    OUT(SX1:EX1,SY1:EY1,2:NKM1) = 0.25D0*							                  &
    (                                                                                 &
        +DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*Q12(SX1+1:EX1+1,SY1:EY1,2:NKM1)             &
            *(INN(SX1+1:EX1+1,SY1+1:EY1+1,2:NKM1)-INN(SX1+1:EX1+1,SY1-1:EY1-1,2:NKM1))&
        -DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*Q12(SX1-1:EX1-1,SY1:EY1,2:NKM1)             &
            *(INN(SX1-1:EX1-1,SY1+1:EY1+1,2:NKM1)-INN(SX1-1:EX1-1,SY1-1:EY1-1,2:NKM1))&
        +DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*Q21(SX1:EX1,SY1+1:EY1+1,2:NKM1)             &
            *(INN(SX1+1:EX1+1,SY1+1:EY1+1,2:NKM1)-INN(SX1-1:EX1-1,SY1+1:EY1+1,2:NKM1))&
        -DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*Q21(SX1:EX1,SY1-1:EY1-1,2:NKM1)             &
            *(INN(SX1+1:EX1+1,SY1-1:EY1-1,2:NKM1)-INN(SX1-1:EX1-1,SY1-1:EY1-1,2:NKM1))&
    )
    
    END SUBROUTINE DIF_DIAGONAL