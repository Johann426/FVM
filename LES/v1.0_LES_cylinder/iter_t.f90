    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
    
    INCLUDE 'mpif.h'

	DOUBLE PRECISION    :: A(0:MAX(NI,NJ,NK)),B(0:MAX(NI,NJ,NK)),C(0:MAX(NI,NJ,NK)),D(0:MAX(NI,NJ,NK))
	    
!---------------------------------------------------------------------------
!   APPROXIMATE FACTORIZATION SOLVER
	AE(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*XIX(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIX(IS1+1:IE1+1,JS1:JE1,KS1:KE1)**2)
	
	AW(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*XIX(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIX(IS1-1:IE1-1,JS1:JE1,KS1:KE1)**2)
	
	AN(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ETX(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETX(IS1:IE1,JS1+1:JE1+1,KS1:KE1)**2)
	
	AS(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ETX(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETX(IS1:IE1,JS1-1:JE1-1,KS1:KE1)**2)
								
	AT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ZTX(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTX(IS1:IE1,JS1:JE1,KS1+1:KE1+1)**2)
	
	AB(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ZTX(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTX(IS1:IE1,JS1:JE1,KS1-1:KE1-1)**2)
	
	AP(IS1:IE1,JS1:JE1,KS1:KE1) = 2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME	&
								+AE(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AW(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AN(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AS(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AT(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AB(IS1:IE1,JS1:JE1,KS1:KE1)

    DD(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*																	&
    (																										&
	+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIX(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*ETX(IS1+1:IE1+1,JS1:JE1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1+1:JE1+1,KS1:KE1)-T(IS1+1:IE1+1,JS1-1:JE1-1,KS1:KE1))	&
	-DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIX(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*ETX(IS1-1:IE1-1,JS1:JE1,KS1:KE1)	&
								*(T(IS1-1:IE1-1,JS1+1:JE1+1,KS1:KE1)-T(IS1-1:IE1-1,JS1-1:JE1-1,KS1:KE1))	&
	+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIX(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*ZTX(IS1+1:IE1+1,JS1:JE1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1+1:KE1+1)-T(IS1+1:IE1+1,JS1:JE1,KS1-1:KE1-1))	&
	-DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIX(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*ZTX(IS1-1:IE1-1,JS1:JE1,KS1:KE1)	&
								*(T(IS1-1:IE1-1,JS1:JE1,KS1+1:KE1+1)-T(IS1-1:IE1-1,JS1:JE1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETX(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*XIX(IS1:IE1,JS1+1:JE1+1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1+1:JE1+1,KS1:KE1)-T(IS1-1:IE1-1,JS1+1:JE1+1,KS1:KE1))	&
	-DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETX(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*XIX(IS1:IE1,JS1-1:JE1-1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1-1:JE1-1,KS1:KE1)-T(IS1-1:IE1-1,JS1-1:JE1-1,KS1:KE1))	&
	+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETX(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ZTX(IS1:IE1,JS1+1:JE1+1,KS1:KE1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1+1:KE1+1)-T(IS1:IE1,JS1+1:JE1+1,KS1-1:KE1-1))	&
	-DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETX(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ZTX(IS1:IE1,JS1-1:JE1-1,KS1:KE1)	&
								*(T(IS1:IE1,JS1-1:JE1-1,KS1+1:KE1+1)-T(IS1:IE1,JS1-1:JE1-1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTX(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*XIX(IS1:IE1,JS1:JE1,KS1+1:KE1+1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1+1:KE1+1)-T(IS1-1:IE1-1,JS1:JE1,KS1+1:KE1+1))	&
	-DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTX(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*XIX(IS1:IE1,JS1:JE1,KS1-1:KE1-1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1-1:KE1-1)-T(IS1-1:IE1-1,JS1:JE1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTX(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ETX(IS1:IE1,JS1:JE1,KS1+1:KE1+1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1+1:KE1+1)-T(IS1:IE1,JS1-1:JE1-1,KS1+1:KE1+1))	&
	-DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTX(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ETX(IS1:IE1,JS1:JE1,KS1-1:KE1-1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1-1:KE1-1)-T(IS1:IE1,JS1-1:JE1-1,KS1-1:KE1-1))	&
    )

	SU(IS1:IE1,JS1:JE1,KS1:KE1)	=																&
				+2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME*T(IS1:IE1,JS1:JE1,KS1:KE1)			&
				-3.D0*PE*ADV(IS1:IE1,JS1:JE1,KS1:KE1) +1.D0*PE*ADV_NM1(IS1:IE1,JS1:JE1,KS1:KE1)	&
				+DIFT(IS1:IE1,JS1:JE1,KS1:KE1)													&
				+DD(IS1:IE1,JS1:JE1,KS1:KE1)

    CALL BCT
!----------------------------------------------------------------------------
!   MODIFIED : I-SWEEP
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
		A(I)=-AW(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
		D(I)=AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO)
		B(I)=-AE(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
		C(I)=AN(I,J,K)*T(I,J+1,K)+AS(I,J,K)*T(I,J-1,K)	&
			+AT(I,J,K)*T(I,J,K+1)+AB(I,J,K)*T(I,J,K-1)	&
			+SU(I,J,K)
		T(I,J,K)=(C(I)-A(I)*T(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
	ENDDO
	DO I=NI-1,IS1,-1
		T(I,J,K)=T(I,J,K)-B(I)*T(I+1,J,K)
	ENDDO
	ENDDO
	ENDDO
!----------------------------------------------------------------------------
	AE(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*XIY(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIY(IS1+1:IE1+1,JS1:JE1,KS1:KE1)**2)
	
	AW(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*XIY(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIY(IS1-1:IE1-1,JS1:JE1,KS1:KE1)**2)
	
	AN(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ETY(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETY(IS1:IE1,JS1+1:JE1+1,KS1:KE1)**2)
	
	AS(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ETY(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETY(IS1:IE1,JS1-1:JE1-1,KS1:KE1)**2)
								
	AT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ZTY(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTY(IS1:IE1,JS1:JE1,KS1+1:KE1+1)**2)
	
	AB(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ZTY(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTY(IS1:IE1,JS1:JE1,KS1-1:KE1-1)**2)
	
	AP(IS1:IE1,JS1:JE1,KS1:KE1) = 2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME	&
								+AE(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AW(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AN(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AS(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AT(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AB(IS1:IE1,JS1:JE1,KS1:KE1)

    DD(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*																	&
    (																										&
	+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIY(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*ETY(IS1+1:IE1+1,JS1:JE1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1+1:JE1+1,KS1:KE1)-T(IS1+1:IE1+1,JS1-1:JE1-1,KS1:KE1))	&
	-DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIY(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*ETY(IS1-1:IE1-1,JS1:JE1,KS1:KE1)	&
								*(T(IS1-1:IE1-1,JS1+1:JE1+1,KS1:KE1)-T(IS1-1:IE1-1,JS1-1:JE1-1,KS1:KE1))	&
	+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIY(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*ZTY(IS1+1:IE1+1,JS1:JE1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1+1:KE1+1)-T(IS1+1:IE1+1,JS1:JE1,KS1-1:KE1-1))	&
	-DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIY(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*ZTY(IS1-1:IE1-1,JS1:JE1,KS1:KE1)	&
								*(T(IS1-1:IE1-1,JS1:JE1,KS1+1:KE1+1)-T(IS1-1:IE1-1,JS1:JE1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETY(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*XIY(IS1:IE1,JS1+1:JE1+1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1+1:JE1+1,KS1:KE1)-T(IS1-1:IE1-1,JS1+1:JE1+1,KS1:KE1))	&
	-DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETY(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*XIY(IS1:IE1,JS1-1:JE1-1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1-1:JE1-1,KS1:KE1)-T(IS1-1:IE1-1,JS1-1:JE1-1,KS1:KE1))	&
	+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETY(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ZTY(IS1:IE1,JS1+1:JE1+1,KS1:KE1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1+1:KE1+1)-T(IS1:IE1,JS1+1:JE1+1,KS1-1:KE1-1))	&
	-DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETY(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ZTY(IS1:IE1,JS1-1:JE1-1,KS1:KE1)	&
								*(T(IS1:IE1,JS1-1:JE1-1,KS1+1:KE1+1)-T(IS1:IE1,JS1-1:JE1-1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTY(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*XIY(IS1:IE1,JS1:JE1,KS1+1:KE1+1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1+1:KE1+1)-T(IS1-1:IE1-1,JS1:JE1,KS1+1:KE1+1))	&
	-DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTY(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*XIY(IS1:IE1,JS1:JE1,KS1-1:KE1-1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1-1:KE1-1)-T(IS1-1:IE1-1,JS1:JE1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTY(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ETY(IS1:IE1,JS1:JE1,KS1+1:KE1+1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1+1:KE1+1)-T(IS1:IE1,JS1-1:JE1-1,KS1+1:KE1+1))	&
	-DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTY(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ETY(IS1:IE1,JS1:JE1,KS1-1:KE1-1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1-1:KE1-1)-T(IS1:IE1,JS1-1:JE1-1,KS1-1:KE1-1))	&
    )

	SU(IS1:IE1,JS1:JE1,KS1:KE1)	=2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME*T(IS1:IE1,JS1:JE1,KS1:KE1)	&
								+DD(IS1:IE1,JS1:JE1,KS1:KE1)

    CALL BCT
!----------------------------------------------------------------------------
!   MODIFIED : J-SWEEP
    DO K=KS1,KE1
    DO I=IS1,IE1
    DO J=JS1,JE1
        A(J)=-AS(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(J)=AP(I,J,K)-A(J)*B(J-1)  ! A(1) SHOULD BE 0(ZERO)
        B(J)=-AN(I,J,K)/D(J)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(J)=AE(I,J,K)*T(I+1,J,K)+AW(I,J,K)*T(I-1,J,K)  &
			+AT(I,J,K)*T(I,J,K+1)+AB(I,J,K)*T(I,J,K-1)	&
            +SU(I,J,K)
        T(I,J,K)=(C(J)-A(J)*T(I,J-1,K))/D(J) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO J=JE1,JS1,-1
        T(I,J,K)=T(I,J,K)-B(J)*T(I,J+1,K)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------
	AE(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*XIZ(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIZ(IS1+1:IE1+1,JS1:JE1,KS1:KE1)**2)
	
	AW(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*XIZ(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIZ(IS1-1:IE1-1,JS1:JE1,KS1:KE1)**2)
	
	AN(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ETZ(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETZ(IS1:IE1,JS1+1:JE1+1,KS1:KE1)**2)
	
	AS(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ETZ(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETZ(IS1:IE1,JS1-1:JE1-1,KS1:KE1)**2)
								
	AT(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTZ(IS1:IE1,JS1:JE1,KS1+1:KE1+1)**2)
	
	AB(IS1:IE1,JS1:JE1,KS1:KE1) = 0.5D0*(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*ZTZ(IS1:IE1,JS1:JE1,KS1:KE1)**2	&
								+DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTZ(IS1:IE1,JS1:JE1,KS1-1:KE1-1)**2)
	
	AP(IS1:IE1,JS1:JE1,KS1:KE1) = 2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME	&
								+AE(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AW(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AN(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AS(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AT(IS1:IE1,JS1:JE1,KS1:KE1)					&
								+AB(IS1:IE1,JS1:JE1,KS1:KE1)

    DD(IS1:IE1,JS1:JE1,KS1:KE1) = 0.25D0*																	&
    (																										&
	+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIZ(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*ETZ(IS1+1:IE1+1,JS1:JE1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1+1:JE1+1,KS1:KE1)-T(IS1+1:IE1+1,JS1-1:JE1-1,KS1:KE1))	&
	-DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIZ(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*ETZ(IS1-1:IE1-1,JS1:JE1,KS1:KE1)	&
								*(T(IS1-1:IE1-1,JS1+1:JE1+1,KS1:KE1)-T(IS1-1:IE1-1,JS1-1:JE1-1,KS1:KE1))	&
	+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*XIZ(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*ZTZ(IS1+1:IE1+1,JS1:JE1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1+1:KE1+1)-T(IS1+1:IE1+1,JS1:JE1,KS1-1:KE1-1))	&
	-DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*XIZ(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*ZTZ(IS1-1:IE1-1,JS1:JE1,KS1:KE1)	&
								*(T(IS1-1:IE1-1,JS1:JE1,KS1+1:KE1+1)-T(IS1-1:IE1-1,JS1:JE1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETZ(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*XIZ(IS1:IE1,JS1+1:JE1+1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1+1:JE1+1,KS1:KE1)-T(IS1-1:IE1-1,JS1+1:JE1+1,KS1:KE1))	&
	-DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETZ(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*XIZ(IS1:IE1,JS1-1:JE1-1,KS1:KE1)	&
								*(T(IS1+1:IE1+1,JS1-1:JE1-1,KS1:KE1)-T(IS1-1:IE1-1,JS1-1:JE1-1,KS1:KE1))	&
	+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ETZ(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*ZTZ(IS1:IE1,JS1+1:JE1+1,KS1:KE1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1+1:KE1+1)-T(IS1:IE1,JS1+1:JE1+1,KS1-1:KE1-1))	&
	-DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ETZ(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*ZTZ(IS1:IE1,JS1-1:JE1-1,KS1:KE1)	&
								*(T(IS1:IE1,JS1-1:JE1-1,KS1+1:KE1+1)-T(IS1:IE1,JS1-1:JE1-1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTZ(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*XIZ(IS1:IE1,JS1:JE1,KS1+1:KE1+1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1+1:KE1+1)-T(IS1-1:IE1-1,JS1:JE1,KS1+1:KE1+1))	&
	-DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTZ(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*XIZ(IS1:IE1,JS1:JE1,KS1-1:KE1-1)	&
								*(T(IS1+1:IE1+1,JS1:JE1,KS1-1:KE1-1)-T(IS1-1:IE1-1,JS1:JE1,KS1-1:KE1-1))	&
	+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ZTZ(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*ETZ(IS1:IE1,JS1:JE1,KS1+1:KE1+1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1+1:KE1+1)-T(IS1:IE1,JS1-1:JE1-1,KS1+1:KE1+1))	&
	-DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ZTZ(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*ETZ(IS1:IE1,JS1:JE1,KS1-1:KE1-1)	&
								*(T(IS1:IE1,JS1+1:JE1+1,KS1-1:KE1-1)-T(IS1:IE1,JS1-1:JE1-1,KS1-1:KE1-1))	&
    )

	SU(IS1:IE1,JS1:JE1,KS1:KE1)	=2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME*T(IS1:IE1,JS1:JE1,KS1:KE1)	&
								+DD(IS1:IE1,JS1:JE1,KS1:KE1)

    CALL BCT
!----------------------------------------------------------------------------
!   MODIFIED : K-SWEEP
    DO J=JS1,JE1
    DO I=IS1,IE1
    DO K=KS1,KE1
        A(K)=-AB(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(K)=AP(I,J,K)-A(K)*B(K-1)  ! A(1) SHOULD BE 0(ZERO)
        B(K)=-AT(I,J,K)/D(K)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(K)=AE(I,J,K)*T(I+1,J,K)+AW(I,J,K)*T(I-1,J,K)  &
			+AN(I,J,K)*T(I,J+1,K)+AS(I,J,K)*T(I,J-1,K)	&
            +SU(I,J,K)
        T(I,J,K)=(C(K)-A(K)*T(I,J,K-1))/D(K) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO K=KE1,KS1,-1
        T(I,J,K)=T(I,J,K)-B(K)*T(I,J,K+1)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------
    END SUBROUTINE ITER_T