    SUBROUTINE ITER_V
    
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
802 NITER = NITER + 1
    CALL DIF_DIAGONAL(V,DD) ! DIAGONAL DIF OF V^* : TREATED AS EXPLICIT

!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(sx1:ex1,sy1:ey1,2:NKM1)=fp(sx1:ex1,sy1:ey1,2:NKM1) * V_N(sx1:ex1,sy1:ey1,2:NKM1)             &
                - 3.D0*RE*NLY    (sx1:ex1,sy1:ey1,2:NKM1) + 1.D0*RE*NLY_NM1(sx1:ex1,sy1:ey1,2:NKM1) &
		        + 1.D0*RE*DIFY   (sx1:ex1,sy1:ey1,2:NKM1) + 2.D0*RE*S_V(sx1:ex1,sy1:ey1,2:NKM1)     &
		        + DD(sx1:ex1,sy1:ey1,2:NKM1)

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BCV
	CALL LISOLV(1,1,2,NI,NJ,NK,V) ! WITHOUT OVER-RELAXATION
        
!---- V-VEL CONVERGE TEST : RESIDUAL OF V-VEL ------------------------------
	RESIV = 0.0D0
	DUM = 0.D0
	DO K=2,NKM1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*V(I+1,J  ,K) + AW(I,J,K)*V(I-1,J  ,K) +		&
				AN(I,J,K)*V(I  ,J+1,K) + AS(I,J,K)*V(I  ,J-1,K) +		&
				SU(I,J,K) - AP(I,J,K)*V(I,J,K) 
		RESIV = RESIV + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*V(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIV = RESIV / DUM ! SCALED RESIDUAL
	!WRITE(82,*) NITER,RESIV
	IF(RESIV.GT.RESIMAX .AND. NITER.LE.NSWPV) GOTO 802
!---------------------------------------------------------------------------
    END SUBROUTINE ITER_V