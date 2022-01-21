	SUBROUTINE UVP_CORRECTION

	USE ComDat_Shared	

	IMPLICIT NONE
	SAVE
!--------------------------------------------------------------------------------------------------------
!   U-VELOCITY CORRECTION : CONTINUITY CONSTRAINT
	DUM_1(SX1:EX1,SY1:EY1,2:NKM1)=-0.5D0*DTIME															&
	*(DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIX(SX1+1:EX1+1,SY1:EY1,2:NKM1)*PLI(SX1+1:EX1+1,SY1:EY1,2:NKM1)  &
    - DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIX(SX1-1:EX1-1,SY1:EY1,2:NKM1)*PLI(SX1-1:EX1-1,SY1:EY1,2:NKM1)  &
	+ DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETX(SX1:EX1,SY1+1:EY1+1,2:NKM1)*PLI(SX1:EX1,SY1+1:EY1+1,2:NKM1)  &
	- DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETX(SX1:EX1,SY1-1:EY1-1,2:NKM1)*PLI(SX1:EX1,SY1-1:EY1-1,2:NKM1)  &
!   + DJAC(SX1:EX1,SY1:EY1,2+1:NKM1+1)*ZTX(SX1:EX1,SY1:EY1,2+1:NKM1+1)*PLI(SX1:EX1,SY1:EY1,2+1:NKM1+1)  &
!   - DJAC(SX1:EX1,SY1:EY1,2-1:NKM1-1)*ZTX(SX1:EX1,SY1:EY1,2-1:NKM1-1)*PLI(SX1:EX1,SY1:EY1,2-1:NKM1-1)  &
    )
	U(SX1:EX1,SY1:EY1,2:NKM1) = U(SX1:EX1,SY1:EY1,2:NKM1)+DUM_1(SX1:EX1,SY1:EY1,2:NKM1)/DJAC(SX1:EX1,SY1:EY1,2:NKM1)

!   V-VELOCITY CORRECTION : CONTINUITY CONSTRAINT
	DUM_2(SX1:EX1,SY1:EY1,2:NKM1)=-0.5D0*DTIME															&
	*(DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIY(SX1+1:EX1+1,SY1:EY1,2:NKM1)*PLI(SX1+1:EX1+1,SY1:EY1,2:NKM1)  &
    - DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIY(SX1-1:EX1-1,SY1:EY1,2:NKM1)*PLI(SX1-1:EX1-1,SY1:EY1,2:NKM1)  &
	+ DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETY(SX1:EX1,SY1+1:EY1+1,2:NKM1)*PLI(SX1:EX1,SY1+1:EY1+1,2:NKM1)  &
	- DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETY(SX1:EX1,SY1-1:EY1-1,2:NKM1)*PLI(SX1:EX1,SY1-1:EY1-1,2:NKM1)  &
!   + DJAC(SX1:EX1,SY1:EY1,2+1:NKM1+1)*ZTY(SX1:EX1,SY1:EY1,2+1:NKM1+1)*PLI(SX1:EX1,SY1:EY1,2+1:NKM1+1)  &
!   - DJAC(SX1:EX1,SY1:EY1,2-1:NKM1-1)*ZTY(SX1:EX1,SY1:EY1,2-1:NKM1-1)*PLI(SX1:EX1,SY1:EY1,2-1:NKM1-1)  &
    )
	V(SX1:EX1,SY1:EY1,2:NKM1)=V(SX1:EX1,SY1:EY1,2:NKM1)+DUM_2(SX1:EX1,SY1:EY1,2:NKM1)/DJAC(SX1:EX1,SY1:EY1,2:NKM1)
!--------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------
!   PRESSURE CORRECTION : CONTINUITY CONSTRAINT
    P(SX:EX,SY:EY,2:NKM1)=P(SX:EX,SY:EY,2:NKM1)+PLI(SX:EX,SY:EY,2:NKM1)
!--------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------
!   INTERPOLATION OF PRESSURE FIELD IN THE COLOCATED GRID
	IF(FLAG.EQ.0) THEN	
    K=KMON	
	DO J=1,NJ	! NODE VALUE
	DO I=1,NI
		PG(I,J,K) = 0.25D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K))
	ENDDO
	ENDDO
	
	DO I=1,NI	! BRANCH CUT
	J=1
		PG(I,J,K) = 0.5D0 * ( PG(I,J,K) + PG(I,EY,K) )
	J=NJ
		PG(I,J,K) = PG(I,SY1,K)
	ENDDO

	DO J=1,NJ-1	! CENTER VALUE
	DO I=1,NI-1
		P(I,J,K) = 0.25D0* ( PG(I+1,J+1,K) + PG(I+1,J,K) + PG(I,J+1,K) + PG(I,J,K) )
 	ENDDO
	ENDDO

	DO I=1,NI-1	! CENTER VALUE
	J=0
		P(I,J,K) = P(I,NJ-1,K)
	J=NJ
		P(I,J,K) = P(I,1,K)
	ENDDO
	DO J=0,NJ	! CENTER VALUE
	I=0
		P(I,J,K) = P(I+1,J,K)
	I=NI
		P(I,J,K) = P(I-1,J,K)
	ENDDO
	ENDIF
	
	DUM_A = P(NIM1,(NJ+1)/2,KMON)						! FAR FIELD	PRESSURE
	P(SX:EX,SY:EY,2:NKM1)=P(SX:EX,SY:EY,2:NKM1)-DUM_A
!--------------------------------------------------------------------------------------------------------

	END SUBROUTINE UVP_CORRECTION