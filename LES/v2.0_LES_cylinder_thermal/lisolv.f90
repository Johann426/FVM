	SUBROUTINE LISOLV(ISTART,JSTART,KSTART,NI2,NJ2,NK2,PHI)

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
		
	INTEGER, INTENT(IN) :: ISTART,JSTART,KSTART,NI2,NJ2,NK2
	DOUBLE PRECISION, INTENT(INOUT) :: PHI(IS:IE,JS:JE,KS:KE)

	DOUBLE PRECISION    :: A(0:MAX(NI,NJ)),B(0:MAX(NI,NJ)),C(0:MAX(NI,NJ)),D(0:MAX(NI,NJ))

!----------------------------------------------------------------------------
!	MODIFIED : I-SWEEP
    DO K=KSTART+1,NK2-1
    DO J=JSTART+1,NJ2-1
    DO I=ISTART+1,NI2-1
		A(I)=-AW(I,J,K)             ! A(1) SHOULD BE 0(ZERO) TO TREAT IMPLICITLY
		D(I)=AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO) TO TREAT IMPLICITLY
		B(I)=-AE(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO) TO TREAT IMPLICITLY
		C(I)=AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)  &
			+AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)	&
			+SU(I,J,K)
		PHI(I,J,K)=(C(I)-A(I)*PHI(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO) TO TREAT IMPLICITLY
    ENDDO
    DO I=NI2-1,ISTART+1,-1
		PHI(I,J,K)=PHI(I,J,K)-B(I)*PHI(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------
!	CALL MPI_BARRIER(NMCW,ERR)
!	CALL MPI_X(PHI)
!!----------------------------------------------------------------------------
!!	MODIFIED : J-SWEEP
!    DO K=KSTART+1,NK2-1
!    DO I=ISTART+1,NI2-1
!    DO J=JSTART+1,NJ2-1
!        A(J)=-AS(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
!        D(J)=AP(I,J,K)-A(J)*B(J-1)  ! A(1) SHOULD BE 0(ZERO)
!        B(J)=-AN(I,J,K)/D(J)        ! B(NI2-1) SHOULD BE 0(ZERO)
!        C(J)=AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)  &
!			+AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)	&
!            +SU(I,J,K)
!        PHI(I,J,K)=(C(J)-A(J)*PHI(I,J-1,K))/D(J) ! A(1) SHOULD BE 0(ZERO)
!    ENDDO
!    DO J=NJ2-1,JSTART+1,-1
!        PHI(I,J,K)=PHI(I,J,K)-B(J)*PHI(I,J+1,K)
!    ENDDO
!    ENDDO
!    ENDDO
!!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!	MODIFIED : K-SWEEP
!----------------------------------------------------------------------------
	END SUBROUTINE LISOLV