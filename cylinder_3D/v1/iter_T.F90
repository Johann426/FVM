    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : T
!   THE NEIGHBOR CELLS : E,W,N,S,T,B
!   THE EXPLICIT TERMS : Su
    AP = DE+DW+DN+DS+DT+DB
    AE = DE
	AW = DW
	AN = DN
	AS = DS
	AT = DT
	AB = DB
	SU = 0
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   TDMA solver iteration
    NITER = 0
801 NITER = NITER + 1
    CALL DIF_DIAGONAL(T,DD) ! DIAGONAL DIF OF U^* : TREATED AS EXPLICIT

!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(SX1:EX1,SY1:EY1,SZ1:EZ1)=FP(SX1:EX1,SY1:EY1,SZ1:EZ1) * T_N(SX1:EX1,SY1:EY1,SZ1:EZ1)             &
                - 3.D0*PE*ADV    (SX1:EX1,SY1:EY1,SZ1:EZ1) + 1.D0*PE*ADV_NM1(SX1:EX1,SY1:EY1,SZ1:EZ1) &
		        + 1.D0*PE*DIFT   (SX1:EX1,SY1:EY1,SZ1:EZ1)                                           &
		        + DD(SX1:EX1,SY1:EY1,SZ1:EZ1)

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BCT
	CALL LISOLV(1,1,1,NI,NJ,NK,T) ! WITHOUT OVER-RELAXATION
        
!---- U-VEL CONVERGE TEST : RESIDUAL OF TEMPERATURE ------------------------------
	RESIT = 0.D0
	DUM = 0.D0
	DO K=SZ1,EZ1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*T(I+1,J,K) + AW(I,J,K)*T(I-1,J,K) +		&
				AN(I,J,K)*T(I,J+1,K) + AS(I,J,K)*T(I,J-1,K) +		&
				AT(I,J,K)*T(I,J,K+1) + AB(I,J,K)*T(I,J,K-1) +		&
				SU(I,J,K) - AP(I,J,K)*T(I,J,K) 
		RESIT = RESIT + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*T(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIT = RESIT / DUM	! SCALED RESIDUAL
	IF(RESIT.GT.RESIMAX .AND. NITER.LT.NSWPT) GOTO 801
!---------------------------------------------------------------------------
    END SUBROUTINE ITER_T