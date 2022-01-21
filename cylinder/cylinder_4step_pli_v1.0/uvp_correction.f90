	SUBROUTINE UVP_CORRECTION

	USE ComDat_Shared	

	IMPLICIT NONE
	SAVE
!--------------------------------------------------------------------------------------------------------
!   U-VELOCITY CORRECTION : CONTINUITY CONSTRAINT
	DUM_1(1:NIM1,1:NJM1,2:NKM1)=-0.5D0*DTIME*(                                              &
        +XIX(1:NIM1,1:NJM1,2:NKM1)*(PLI(1+1:NI,1:NJM1,2:NKM1)-PLI(1-1:NIM2,1:NJM1,2:NKM1))	&
        +ETX(1:NIM1,1:NJM1,2:NKM1)*(PLI(1:NIM1,1+1:NJ,2:NKM1)-PLI(1:NIM1,1-1:NJM2,2:NKM1))	&
                                            )
	U(1:NIM1,1:NJM1,2:NKM1) = U(1:NIM1,1:NJM1,2:NKM1)+DUM_1(1:NIM1,1:NJM1,2:NKM1)

!   V-VELOCITY CORRECTION : CONTINUITY CONSTRAINT
	DUM_2(1:NIM1,1:NJM1,2:NKM1)=-0.5D0*DTIME*(                                              &
        +XIY(1:NIM1,1:NJM1,2:NKM1)*(PLI(1+1:NI,1:NJM1,2:NKM1)-PLI(1-1:NIM2,1:NJM1,2:NKM1))  &
        +ETY(1:NIM1,1:NJM1,2:NKM1)*(PLI(1:NIM1,1+1:NJ,2:NKM1)-PLI(1:NIM1,1-1:NJM2,2:NKM1))  &
								            )
	V(1:NIM1,1:NJM1,2:NKM1)=V(1:NIM1,1:NJM1,2:NKM1)+DUM_2(1:NIM1,1:NJM1,2:NKM1)
!--------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------
!   PRESSURE CORRECTION : CONTINUITY CONSTRAINT
	DUM_A = P(NIM1,(NJ+1)/2,KMON)						! FAR FIELD	PRESSURE
	P(0:NIM1,0:NJM1,2:NKM1)=P(0:NIM1,0:NJM1,2:NKM1)+PLI(0:NIM1,0:NJM1,2:NKM1)
	P(0:NIM1,0:NJM1,2:NKM1)=P(0:NIM1,0:NJM1,2:NKM1)-DUM_A
!--------------------------------------------------------------------------------------------------------

	END SUBROUTINE UVP_CORRECTION