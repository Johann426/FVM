	SUBROUTINE UVP_CORRECTION

	USE ComDat_Shared	

	IMPLICIT NONE
	SAVE
!--------------------------------------------------------------------------------------------------------
!   VELOCITY CORRECTION
	CALL P_SOURCE
	U(SX1:EX1,SY1:EY1,2:NKM1) = U(SX1:EX1,SY1:EY1,2:NKM1)+DTIME*S_U(SX1:EX1,SY1:EY1,2:NKM1)
	V(SX1:EX1,SY1:EY1,2:NKM1) = V(SX1:EX1,SY1:EY1,2:NKM1)+DTIME*S_V(SX1:EX1,SY1:EY1,2:NKM1)
!--------------------------------------------------------------------------------------------------------
!   PRESSURE CORRECTION : CONTINUITY CONSTRAINT
    P(SX:EX,SY:EY,2:NKM1)=P(SX:EX,SY:EY,2:NKM1)+PLI(SX:EX,SY:EY,2:NKM1)
    DUM_A = P(NIM1,(NJ+1)/2,KMON)						! FAR FIELD	PRESSURE
	P(SX:EX,SY:EY,2:NKM1)=P(SX:EX,SY:EY,2:NKM1)-DUM_A
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE UVP_CORRECTION