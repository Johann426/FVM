    SUBROUTINE LISOLV(ISTART,JSTART,KSTART,NI2,NJ2,NK2,PHI)

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
		
	INTEGER, INTENT(IN) :: ISTART,JSTART,KSTART,NI2,NJ2,NK2                   
	DOUBLE PRECISION, INTENT(INOUT) :: PHI(sx1-1:ex1+1,sy1-1:ey1+1,2:NKM1)

	INTEGER :: ISTM1,JSTM1,KSTM1,II,JJ,KK,NIMM1,NJMM1,NKMM1
	DOUBLE PRECISION    :: TERM
	DOUBLE PRECISION    :: A(0:MAX(NI,NJ)),B(0:MAX(NI,NJ)),C(0:MAX(NI,NJ)),D(0:MAX(NI,NJ))
	DOUBLE PRECISION    :: SWEEPDIR

!----------------------------------------------------------------------------
!   MODIFIED : I-SWEEP
    DO K=KSTART,NK2-1
    DO J=JSTART,NJ2-1
    DO I=ISTART,NI2-1
        A(I)=-AW(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(I)=AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO)
        B(I)=-AE(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)  &
!           +AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)	&
            +SU(I,J,K)
        PHI(I,J,K)=(C(I)-A(I)*PHI(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO I=NI2-1,ISTART,-1
        PHI(I,J,K)=PHI(I,J,K)-B(I)*PHI(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------

!-- (1)  E-W SWEEP  : ORIGINAL -----------
!	ISTM1=ISTART-1
!	A(ISTM1)=0.0d0
!	DO K=KSTART,NK2-1
!	DO J=JSTART,NJ2-1
!		C(ISTM1)=PPHI(ISTM1,J,K)
!	DO I=ISTART,NI2-1
!		A(I)=AE(I,J,K)
!		B(I)=AW(I,J,K)
!		C(I)=AN(I,J,K)*PPHI(I,J+1,K)+AS(I,J,K)*PPHI(I,J-1,K)	&
!!			+AT(I,J,K)*PPHI(I,J,K+1)+AB(I,J,K)*PPHI(I,J,K-1)	&
!			+SU(I,J,K) 
!		D(I)=AP(I,J,K)
!		TERM=1.0D0/(D(I)-B(I)*A(I-1))
!		A(I)=A(I)*TERM
!		C(I)=(C(I)+B(I)*C(I-1))*TERM
!	ENDDO
!
!	DO II=ISTART,NI2-1
!		I=NI2+ISTM1-II
!		PPHI(I,J,K)=A(I)*PPHI(I+1,J,K)+C(I)
!	ENDDO
!	ENDDO
!	ENDDO
!-------------------------------------------


!-- (2)  N-S SWEEP  : ORIGINAL -------------
!	JSTM1=JSTART-1
!	A(JSTM1)=0.0d0
!	DO K=KSTART,NK2-1
!	DO I=ISTART,NI2-1
!		C(JSTM1)=PPHI(I,JSTM1,K)
!	DO J=JSTART,NJ2-1
!		A(J)=AN(I,J,K)
!		B(J)=AS(I,J,K)
!		C(J)=AE(I,J,K)*PPHI(I+1,J,K)+AW(I,J,K)*PPHI(I-1,J,K)    &
!!			+AT(I,J,K)*PPHI(I,J,K+1)+AB(I,J,K)*PPHI(I,J,K-1)	&
!		    +SU(I,J,K) 
!		D(J)=AP(I,J,K)
!		TERM=1.0D0/(D(J)-B(J)*A(J-1))
!		A(J)=A(J)*TERM
!		C(J)=(C(J)+B(J)*C(J-1))*TERM
!	ENDDO
!
!	DO JJ=JSTART,NJ2-1
!		J=NJ2+JSTM1-JJ
!		PPHI(I,J,K)=A(J)*PPHI(I,J+1,K)+C(J)
!	ENDDO
!	ENDDO
!	ENDDO
!------------------------------------------

    END SUBROUTINE LISOLV

