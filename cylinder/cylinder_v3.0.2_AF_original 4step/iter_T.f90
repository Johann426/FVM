    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE

    DOUBLE PRECISION    :: A(0:MAX(NI,NJ)),B(0:MAX(NI,NJ)),C(0:MAX(NI,NJ)),D(0:MAX(NI,NJ))

!---------------------------------------------------------------------------
!   APPROXIMATE FACTORIZATION SOLVER
	
	AE(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*XIX(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIX(SX1+1:EX1+1,SY1:EY1,2:NKM1)**2)
	
	AW(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*XIX(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIX(SX1-1:EX1-1,SY1:EY1,2:NKM1)**2)
	
	AN(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*ETX(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETX(SX1:EX1,SY1+1:EY1+1,2:NKM1)**2)
	
	AS(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*ETX(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETX(SX1:EX1,SY1-1:EY1-1,2:NKM1)**2)
	
	AP(SX1:EX1,SY1:EY1,2:NKM1) = 2.D0*PE*DJAC(SX1:EX1,SY1:EY1,2:NKM1)/DTIME	&
								+AE(SX1:EX1,SY1:EY1,2:NKM1)					&
								+AW(SX1:EX1,SY1:EY1,2:NKM1)					&
								+AN(SX1:EX1,SY1:EY1,2:NKM1)					&
								+AS(SX1:EX1,SY1:EY1,2:NKM1)
								
    DD(SX1:EX1,SY1:EY1,2:NKM1) = 0.25D0*																&
    (																									&
	+DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIX(SX1+1:EX1+1,SY1:EY1,2:NKM1)*ETX(SX1+1:EX1+1,SY1:EY1,2:NKM1)	&
								*(T(SX1+1:EX1+1,SY1+1:EY1+1,2:NKM1)-T(SX1+1:EX1+1,SY1-1:EY1-1,2:NKM1))	&
	-DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIX(SX1-1:EX1-1,SY1:EY1,2:NKM1)*ETX(SX1-1:EX1-1,SY1:EY1,2:NKM1)	&
								*(T(SX1-1:EX1-1,SY1+1:EY1+1,2:NKM1)-T(SX1-1:EX1-1,SY1-1:EY1-1,2:NKM1))	&
	+DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*XIX(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETX(SX1:EX1,SY1+1:EY1+1,2:NKM1)	&
								*(T(SX1+1:EX1+1,SY1+1:EY1+1,2:NKM1)-T(SX1-1:EX1-1,SY1+1:EY1+1,2:NKM1))	&
	-DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*XIX(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETX(SX1:EX1,SY1-1:EY1-1,2:NKM1)	&
								*(T(SX1+1:EX1+1,SY1-1:EY1-1,2:NKM1)-T(SX1-1:EX1-1,SY1-1:EY1-1,2:NKM1))	&
    )
    
    SU(sx1:ex1,sy1:ey1,2:NKM1)=2.D0*PE*DJAC(sx1:ex1,sy1:ey1,2:NKM1)/DTIME * T_N(sx1:ex1,sy1:ey1,2:NKM1)	&
                - 3.D0*PE*ADV    (sx1:ex1,sy1:ey1,2:NKM1) + 1.D0*PE*ADV_NM1(sx1:ex1,sy1:ey1,2:NKM1)		&
		        +		  DIFT   (sx1:ex1,sy1:ey1,2:NKM1)												&
		        + DD(sx1:ex1,sy1:ey1,2:NKM1)

    CALL BCT
!----------------------------------------------------------------------------
!   MODIFIED : I-SWEEP
    DO K=2,NK-1
    DO J=SY1,NJ-1
    DO I=SX1,NI-1
        A(I)=-AW(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(I)=AP(I,J,K)-A(I)*B(I-1)  ! A(1) SHOULD BE 0(ZERO)
        B(I)=-AE(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=AN(I,J,K)*T(I,J+1,K)+AS(I,J,K)*T(I,J-1,K)  &
!           +AT(I,J,K)*T(I,J,K+1)+AB(I,J,K)*T(I,J,K-1)	&
            +SU(I,J,K)
        T(I,J,K)=(C(I)-A(I)*T(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO I=NI-1,SX1,-1
        T(I,J,K)=T(I,J,K)-B(I)*T(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------
	AE(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*XIY(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIY(SX1+1:EX1+1,SY1:EY1,2:NKM1)**2)
	
	AW(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*XIY(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIY(SX1-1:EX1-1,SY1:EY1,2:NKM1)**2)
	
	AN(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*ETY(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETY(SX1:EX1,SY1+1:EY1+1,2:NKM1)**2)
	
	AS(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*ETY(SX1:EX1,SY1:EY1,2:NKM1)**2	&
								+DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETY(SX1:EX1,SY1-1:EY1-1,2:NKM1)**2)
	
	AP(SX1:EX1,SY1:EY1,2:NKM1) = 2.D0*PE*DJAC(SX1:EX1,SY1:EY1,2:NKM1)/DTIME	&
								+AE(SX1:EX1,SY1:EY1,2:NKM1)					&
								+AW(SX1:EX1,SY1:EY1,2:NKM1)					&
								+AN(SX1:EX1,SY1:EY1,2:NKM1)					&
								+AS(SX1:EX1,SY1:EY1,2:NKM1)
								
    DD(SX1:EX1,SY1:EY1,2:NKM1) = 0.25D0*																&
    (																									&
	+DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*XIY(SX1+1:EX1+1,SY1:EY1,2:NKM1)*ETY(SX1+1:EX1+1,SY1:EY1,2:NKM1)	&
								*(V(SX1+1:EX1+1,SY1+1:EY1+1,2:NKM1)-V(SX1+1:EX1+1,SY1-1:EY1-1,2:NKM1))	&
	-DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*XIY(SX1-1:EX1-1,SY1:EY1,2:NKM1)*ETY(SX1-1:EX1-1,SY1:EY1,2:NKM1)	&
								*(V(SX1-1:EX1-1,SY1+1:EY1+1,2:NKM1)-V(SX1-1:EX1-1,SY1-1:EY1-1,2:NKM1))	&
	+DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*XIY(SX1:EX1,SY1+1:EY1+1,2:NKM1)*ETY(SX1:EX1,SY1+1:EY1+1,2:NKM1)	&
								*(V(SX1+1:EX1+1,SY1+1:EY1+1,2:NKM1)-V(SX1-1:EX1-1,SY1+1:EY1+1,2:NKM1))	&
	-DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*XIY(SX1:EX1,SY1-1:EY1-1,2:NKM1)*ETY(SX1:EX1,SY1-1:EY1-1,2:NKM1)	&
								*(V(SX1+1:EX1+1,SY1-1:EY1-1,2:NKM1)-V(SX1-1:EX1-1,SY1-1:EY1-1,2:NKM1))	&
    )
    
    SU(sx1:ex1,sy1:ey1,2:NKM1)=2.D0*PE*DJAC(sx1:ex1,sy1:ey1,2:NKM1)/DTIME*T(sx1:ex1,sy1:ey1,2:NKM1)   &
                              +DD(sx1:ex1,sy1:ey1,2:NKM1)
                              
    CALL BCT
!----------------------------------------------------------------------------
!   MODIFIED : J-SWEEP
    DO K=2,NK-1
    DO I=SX1,NI-1
    DO J=SY1,NJ-1
        A(J)=-AS(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(J)=AP(I,J,K)-A(J)*B(J-1)  ! A(1) SHOULD BE 0(ZERO)
        B(J)=-AN(I,J,K)/D(J)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(J)=AE(I,J,K)*T(I+1,J,K)+AW(I,J,K)*T(I-1,J,K)  &
!           +AT(I,J,K)*T(I,J,K+1)+AB(I,J,K)*T(I,J,K-1)	&
            +SU(I,J,K)
        T(I,J,K)=(C(J)-A(J)*T(I,J-1,K))/D(J) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO J=NJ-1,SY1,-1
        T(I,J,K)=T(I,J,K)-B(J)*T(I,J+1,K)
    ENDDO
    ENDDO
    ENDDO
!----------------------------------------------------------------------------
    END SUBROUTINE ITER_T