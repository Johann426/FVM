    SUBROUTINE ITER_U
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : P
!   THE NEIGHBOR CELLS : E,W,N,S,T,B
!   THE EXPLICIT TERMS : Su
    AP = fp+DE+DW+DN+DS+DT+DB
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
    CALL DIF_DIAGONAL(U,DD) ! DIAGONAL DIF OF U^* : TREATED AS EXPLICIT

!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(SX1:EX1,SY1:EY1,SZ1:EZ1)=FP(SX1:EX1,SY1:EY1,SZ1:EZ1) * U_N(SX1:EX1,SY1:EY1,SZ1:EZ1)             &
                - 3.D0*RE*NLX    (SX1:EX1,SY1:EY1,SZ1:EZ1) + 1.D0*RE*NLX_NM1(SX1:EX1,SY1:EY1,SZ1:EZ1) &
		        + 1.D0*RE*DIFX   (SX1:EX1,SY1:EY1,SZ1:EZ1) + 2.D0*RE*S_U(SX1:EX1,SY1:EY1,SZ1:EZ1)     &
		        + DD(SX1:EX1,SY1:EY1,SZ1:EZ1)

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BCU
	CALL LISOLV(1,1,1,NI,NJ,NK,U) ! WITHOUT OVER-RELAXATION

!---- U-VEL CONVERGE TEST : RESIDUAL OF U-VEL ------------------------------
	RESIU = 0.D0
	DUM = 0.D0
	DO K=SZ1,EZ1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*U(I+1,J,K) + AW(I,J,K)*U(I-1,J,K) +		&
				AN(I,J,K)*U(I,J+1,K) + AS(I,J,K)*U(I,J-1,K) +		&
				AT(I,J,K)*U(I,J,K+1) + AB(I,J,K)*U(I,J,K-1) +		&
				SU(I,J,K) - AP(I,J,K)*U(I,J,K) 
		RESIU = RESIU + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*U(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIU = RESIU / DUM ! SCALED RESIDUAL
	IF(RESIU.GT.RESIMAX .AND. NITER.LT.NSWPU) GOTO 801
!---------------------------------------------------------------------------
    END SUBROUTINE ITER_U