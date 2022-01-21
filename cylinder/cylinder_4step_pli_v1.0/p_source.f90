	SUBROUTINE P_SOURCE

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
!--------------------------------------------------------------------------------------------------------
!   PRESSURE GRADIENT
		S_U(SX1:EX1,SY1:EY1,2:NKM1) = -DJAC(SX1:EX1,SY1:EY1,2:NKM1)		&
			*( 0.5D0*XIX(SX1:EX1,SY1:EY1,2:NKM1)*(P(SX1+1:EX1+1,SY1  :EY1  ,2:NKM1)-P(SX1-1:EX1-1,SY1  :EY1  ,2:NKM1))	&
			  +0.5D0*ETX(SX1:EX1,SY1:EY1,2:NKM1)*(P(SX1  :EX1  ,SY1+1:EY1+1,2:NKM1)-P(SX1  :EX1  ,SY1-1:EY1-1,2:NKM1)) )

		S_V(SX1:EX1,SY1:EY1,2:NKM1) = -DJAC(SX1:EX1,SY1:EY1,2:NKM1)		&
			*( 0.5D0*XIY(SX1:EX1,SY1:EY1,2:NKM1)*(P(SX1+1:EX1+1,SY1  :EY1  ,2:NKM1)-P(SX1-1:EX1-1,SY1  :EY1  ,2:NKM1))	&
			  +0.5D0*ETY(SX1:EX1,SY1:EY1,2:NKM1)*(P(SX1  :EX1  ,SY1+1:EY1+1,2:NKM1)-P(SX1  :EX1  ,SY1-1:EY1-1,2:NKM1)) )
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE P_SOURCE