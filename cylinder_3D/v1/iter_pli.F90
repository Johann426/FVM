    SUBROUTINE ITER_PLI
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : PLI
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
804	NITER = NITER+1
    CALL DIF_DIAGONAL(PLI,DD) ! DIAGONAL DIF OF PLI : TREATED AS EXPLICIT
    
!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(SX1:EX1,SY1:EY1,SZ1:EZ1) = -1.0D0/DTIME/DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1) *(		&
                        +UE(SX1:EX1,SY1:EY1,SZ1:EZ1) -UW(SX1:EX1,SY1:EY1,SZ1:EZ1)	&
                        +VN(SX1:EX1,SY1:EY1,SZ1:EZ1) -VS(SX1:EX1,SY1:EY1,SZ1:EZ1)	&
                        +WT(SX1:EX1,SY1:EY1,SZ1:EZ1) -WB(SX1:EX1,SY1:EY1,SZ1:EZ1))	&
	                    +DD

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BC_PLI
	CALL LISOLV(1,1,1,NI,NJ,NK,PLI)
    CALL BC_PLI
        
!---- PLI CONVERGE TEST : RESIDUAL OF PLI ----------------------------------
	RESIP = 0.0D0
    DUM = 0.D0
	DO K=SZ1,EZ1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*PLI(I+1,J,K) + AW(I,J,K)*PLI(I-1,J,K) +   &
				AN(I,J,K)*PLI(I,J+1,K) + AS(I,J,K)*PLI(I,J-1,K) +	&
				AT(I,J,K)*PLI(I,J,K+1) + AB(I,J,K)*PLI(I,J,K-1) +	&
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