    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : P
!   THE NEIGHBOR CELLS : E,W,N,S
!   THE EXPLICIT TERMS : Su
    AP = fp+DE+DW+DN+DS
    AE = DE
	AW = DW
	AN = DN
	AS = DS
	SU = 0
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   TDMA solver iteration
    NITER = 0
801 NITER = NITER + 1
    CALL DIF_DIAGONAL(T,DD) ! DIAGONAL DIF OF U^* : TREATED AS EXPLICIT

!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(sx1:ex1,sy1:ey1,2:NKM1)=fp(sx1:ex1,sy1:ey1,2:NKM1) * T_N(sx1:ex1,sy1:ey1,2:NKM1)             &
                - 3.D0*PE*ADV    (sx1:ex1,sy1:ey1,2:NKM1) + 1.D0*PE*ADV_NM1(sx1:ex1,sy1:ey1,2:NKM1) &
		        + 1.D0*PE*DIFT   (sx1:ex1,sy1:ey1,2:NKM1)                                           &
		        + DD(sx1:ex1,sy1:ey1,2:NKM1)

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BCT
	CALL LISOLV(1,1,2,NI,NJ,NK,T) ! WITHOUT OVER-RELAXATION
        
!---- U-VEL CONVERGE TEST : RESIDUAL OF TEMPERATURE ------------------------------
	RESIT = 0.D0
	DUM = 0.D0
	DO K=2,NKM1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*T(I+1,J  ,K) + AW(I,J,K)*T(I-1,J  ,K) +		&
				AN(I,J,K)*T(I  ,J+1,K) + AS(I,J,K)*T(I  ,J-1,K) +		&
				SU(I,J,K) - AP(I,J,K)*T(I,J,K) 
		RESIT = RESIT + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*T(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIT = RESIT / DUM ! SCALED RESIDUAL
    !WRITE(81,*) NITER,RESIT
	IF(RESIT.GT.RESIMAX .AND. NITER.LE.NSWPT) GOTO 801
!---------------------------------------------------------------------------
    END SUBROUTINE ITER_T