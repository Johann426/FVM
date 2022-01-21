	SUBROUTINE LISOLV(ISTART,JSTART,KSTART,NI2,NJ2,NK2,PHI)

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
		
	INTEGER, INTENT(IN) :: ISTART,JSTART,KSTART,NI2,NJ2,NK2
	DOUBLE PRECISION, INTENT(INOUT) :: PHI(IS:IE,JS:JE,KS:KE)

	DOUBLE PRECISION    :: A(0:MAX(NI,NJ)),B(0:MAX(NI,NJ)),C(0:MAX(NI,NJ)),D(0:MAX(NI,NJ))

!----------------------------------------------------------------------------
!	MODIFIED : I-SWEEP
    DO K=KSTART,NK2-1
    DO J=JSTART,NJ2-1
    DO I=ISTART,NI2-1
		A(I)=-AW(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
		D(I)=AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO)
		B(I)=-AE(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
		C(I)=AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)  &
			+AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)	&
			+SU(I,J,K)
		PHI(I,J,K)=(C(I)-A(I)*PHI(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO I=NI2-1,ISTART,-1
		PHI(I,J,K)=PHI(I,J,K)-B(I)*PHI(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!	MODIFIED : J-SWEEP
!----------------------------------------------------------------------------
	END SUBROUTINE LISOLV