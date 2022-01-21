    SUBROUTINE ITER_V
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE

!   MODIFIED : I-SWEEP

    AP = 1.D0+DE+DW
	AE = DE
	AW = DW
    SU = 2*RE*(-1.5D0*(NLX+DFEY) +0.5D0*(NLX_NM1+DFEY_NM1) +1.D0/RE*DFIY +S_V)
	
	CALL BCV
	
    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1,EX1
        A(I)=-AW(I,J,K)				! A(1) SHOULD BE 0(ZERO)
        D(I)=+AP(I,J,K)-A(I)*B(I-1)	! A(1) SHOULD BE 0(ZERO)
        B(I)=-AE(I,J,K)/D(I)		! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=+SU(I,J,K)
        U(I,J,K)=(C(I)-A(I)*U(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO I=EX1,SX1,-1
        V(I,J,K)=V(I,J,K)-B(I)*V(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO

!   MODIFIED : J-SWEEP
	
	AP = 1.D0+DN+DS
	AN = DN
	AS = DS
	SU = V

	CALL BCV
		
    DO I=SX1,EX1
    DO K=2,NKM1
    DO J=SY1,EY1
        A(I)=-AS(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(I)=AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO)
        B(I)=-AN(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=+SU(I,J,K)
        U(I,J,K)=(C(I)-A(I)*V(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO J=EY1,SY1,-1
        V(I,J,K)=V(I,J,K)-B(I)*U(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO

	V = V_N + V
	
    END SUBROUTINE ITER_V