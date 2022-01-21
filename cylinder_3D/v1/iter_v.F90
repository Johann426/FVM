    SUBROUTINE ITER_V
    
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
802 NITER = NITER + 1
    CALL DIF_DIAGONAL(V,DD) ! DIAGONAL DIF OF V^* : TREATED AS EXPLICIT

!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(SX1:EX1,SY1:EY1,SZ1:EZ1)=FP(SX1:EX1,SY1:EY1,SZ1:EZ1) * V_N(SX1:EX1,SY1:EY1,SZ1:EZ1)             &
                - 3.D0*RE*NLY    (SX1:EX1,SY1:EY1,SZ1:EZ1) + 1.D0*RE*NLY_NM1(SX1:EX1,SY1:EY1,SZ1:EZ1) &
		        + 1.D0*RE*DIFY   (SX1:EX1,SY1:EY1,SZ1:EZ1) + 2.D0*RE*S_V(SX1:EX1,SY1:EY1,SZ1:EZ1)     &
		        + DD(SX1:EX1,SY1:EY1,SZ1:EZ1)

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BCV
	CALL LISOLV(1,1,1,NI,NJ,NK,V) ! WITHOUT OVER-RELAXATION
        
!---- V-VEL CONVERGE TEST : RESIDUAL OF V-VEL ------------------------------
	RESIV = 0.0D0
	DUM = 0.D0
	DO K=SZ1,EZ1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*V(I+1,J,K) + AW(I,J,K)*V(I-1,J,K) +		&
				AN(I,J,K)*V(I,J+1,K) + AS(I,J,K)*V(I,J-1,K) +		&
				AT(I,J,K)*V(I,J,K+1) + AB(I,J,K)*V(I,J,K-1) +		&
				SU(I,J,K) - AP(I,J,K)*V(I,J,K) 
		RESIV = RESIV + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*V(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIV = RESIV / DUM ! SCALED RESIDUAL
	IF(RESIV.GT.RESIMAX .AND. NITER.LT.NSWPV) GOTO 802
!---------------------------------------------------------------------------
    END SUBROUTINE ITER_V