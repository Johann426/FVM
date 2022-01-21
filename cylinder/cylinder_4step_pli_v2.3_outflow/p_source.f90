	SUBROUTINE P_SOURCE

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
!--------------------------------------------------------------------------------------------------------
!   PRESSURE GRADIENT
    S_U(SX1:EX1,SY1:EY1,2:NKM1) = -0.5D0															&
    *(                                                                                              &
    + DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIX(SX1+1:EX1+1,SY1:EY1,2:NKM1)*P(SX1+1:EX1+1,SY1:EY1,2:NKM1)&
    - DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIX(SX1-1:EX1-1,SY1:EY1,2:NKM1)*P(SX1-1:EX1-1,SY1:EY1,2:NKM1)&
	+ DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETX(SX1:EX1,SY1+1:EY1+1,2:NKM1)*P(SX1:EX1,SY1+1:EY1+1,2:NKM1)&
    - DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETX(SX1:EX1,SY1-1:EY1-1,2:NKM1)*P(SX1:EX1,SY1-1:EY1-1,2:NKM1)&
!   + DJAC(SX1:EX1,SY1:EY1,2+1:NKM1+1)*ZTX(SX1:EX1,SY1:EY1,2+1:NKM1+1)*P(SX1:EX1,SY1:EY1,2+1:NKM1+1)&
!   - DJAC(SX1:EX1,SY1:EY1,2-1:NKM1-1)*ZTX(SX1:EX1,SY1:EY1,2-1:NKM1-1)*P(SX1:EX1,SY1:EY1,2-1:NKM1-1)&
    )

    S_V(SX1:EX1,SY1:EY1,2:NKM1) = -0.5D0															&
    *(                                                                                              &
    + DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIY(SX1+1:EX1+1,SY1:EY1,2:NKM1)*P(SX1+1:EX1+1,SY1:EY1,2:NKM1)&
    - DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIY(SX1-1:EX1-1,SY1:EY1,2:NKM1)*P(SX1-1:EX1-1,SY1:EY1,2:NKM1)&
    + DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETY(SX1:EX1,SY1+1:EY1+1,2:NKM1)*P(SX1:EX1,SY1+1:EY1+1,2:NKM1)&
    - DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETY(SX1:EX1,SY1-1:EY1-1,2:NKM1)*P(SX1:EX1,SY1-1:EY1-1,2:NKM1)&
!   + DJAC(SX1:EX1,SY1:EY1,2+1:NKM1+1)*ZTY(SX1:EX1,SY1:EY1,2+1:NKM1+1)*P(SX1:EX1,SY1:EY1,2+1:NKM1+1)&
!   - DJAC(SX1:EX1,SY1:EY1,2-1:NKM1-1)*ZTY(SX1:EX1,SY1:EY1,2-1:NKM1-1)*P(SX1:EX1,SY1:EY1,2-1:NKM1-1)&
    )
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE P_SOURCE