    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE

!   MODIFIED : I-SWEEP

    AP = 1.D0+DE+DW
	AE = DE
	AW = DW
    SU = 2*PE*(-1.5D0*(ADV+DFET) +0.5D0*(ADV_NM1+DFET_NM1) +1.D0/PE*DFIT)
	
	CALL BCT
	
    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1,EX1
        A(I)=-AW(I,J,K)				! A(1) SHOULD BE 0(ZERO)
        D(I)=+AP(I,J,K)-A(I)*B(I-1)	! A(1) SHOULD BE 0(ZERO)
        B(I)=-AE(I,J,K)/D(I)		! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=+SU(I,J,K)
        T(I,J,K)=(C(I)-A(I)*T(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO I=EX1,SX1,-1
        T(I,J,K)=T(I,J,K)-B(I)*T(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO

!   MODIFIED : J-SWEEP
	
	AP = 1.D0+DN+DS
	AN = DN
	AS = DS
	SU = T

	CALL BCT
		
    DO I=SX1,EX1
    DO K=2,NKM1
    DO J=SY1,EY1
        A(I)=-AS(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(I)=+AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO)
        B(I)=-AN(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=+SU(I,J,K)
        T(I,J,K)=(C(I)-A(I)*T(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO J=EY1,SY1,-1
        T(I,J,K)=T(I,J,K)-B(I)*T(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO

	T = T_N + T
	
    END SUBROUTINE ITER_T