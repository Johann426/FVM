    SUBROUTINE ITER_PLI
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : PLI
!   THE NEIGHBOR CELLS : E,W,N,S
!   THE EXPLICIT TERMS : Su
    AP = DE+DW+DN+DS
    AE = DE
	AW = DW
	AN = DN
	AS = DS
	SU = 0
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   TDMA solver iteration
    NITER = 0 
804	NITER = NITER+1
    CALL DIF_DIAGONAL(PLI,DD) ! DIAGONAL DIF OF PLI : TREATED AS EXPLICIT
    
!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(sx1:ex1,sy1:ey1,2:NKM1) = -1.0D0/DTIME/DJAC(sx1:ex1,sy1:ey1,2:NKM1) *(   &
                        +UE(sx1:ex1,sy1:ey1,2:NKM1) -UW(sx1:ex1,sy1:ey1,2:NKM1) &
                        +VN(sx1:ex1,sy1:ey1,2:NKM1) -VS(sx1:ex1,sy1:ey1,2:NKM1))&
	                    -DD

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BC_PLI
	CALL LISOLV(1,1,2,NI,NJ,NK,PLI)
    CALL BC_PLI
        
!---- PLI CONVERGE TEST : RESIDUAL OF PLI ----------------------------------
	RESIP = 0.0D0
    DUM = 0.D0
	DO K=2,NKM1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*PLI(I+1,J  ,K) + AW(I,J,K)*PLI(I-1,J  ,K) +   &
				AN(I,J,K)*PLI(I  ,J+1,K) + AS(I,J,K)*PLI(I  ,J-1,K) +	&
				SU(I,J,K) - AP(I,J,K)*PLI(I,J,K) 
		RESIP = RESIP + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*PLI(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIP = RESIP / DUM ! SCALED RESIDUAL
	IF(RESIP.GT.RESIMAX .AND. NITER.LT.NSWPP) GOTO 804
!---------------------------------------------------------------------------

    END SUBROUTINE ITER_PLI