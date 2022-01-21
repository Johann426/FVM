    SUBROUTINE ITER_PLI
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE

	INCLUDE 'mpif.h'
	
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : P
!   THE NEIGHBOR CELLS : E,W,N,S,T,B
!   THE EXPLICIT TERMS : Su
    AP(IS:IE,JS:JE,KS:KE) = +DE(IS:IE,JS:JE,KS:KE)+DW(IS:IE,JS:JE,KS:KE)	&
							+DN(IS:IE,JS:JE,KS:KE)+DS(IS:IE,JS:JE,KS:KE)	&
							+DT(IS:IE,JS:JE,KS:KE)+DB(IS:IE,JS:JE,KS:KE)
							
    AE(IS:IE,JS:JE,KS:KE) = DE(IS:IE,JS:JE,KS:KE)
	AW(IS:IE,JS:JE,KS:KE) = DW(IS:IE,JS:JE,KS:KE)
	AN(IS:IE,JS:JE,KS:KE) = DN(IS:IE,JS:JE,KS:KE)
	AS(IS:IE,JS:JE,KS:KE) = DS(IS:IE,JS:JE,KS:KE)
	AT(IS:IE,JS:JE,KS:KE) = DT(IS:IE,JS:JE,KS:KE)
	AB(IS:IE,JS:JE,KS:KE) = DB(IS:IE,JS:JE,KS:KE)
	SU(IS:IE,JS:JE,KS:KE) = 0
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   TDMA solver iteration
    NITER = 0 
804	NITER = NITER+1
    CALL DIF_DIAGONAL(PLI,DD) ! DIAGONAL DIF OF PLI : TREATED AS EXPLICIT
    
!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(IS1:IE1,JS1:JE1,KS1:KE1) = -1.0D0/DTIME/DJAC(IS1:IE1,JS1:JE1,KS1:KE1) *(		&
                        +UE(IS1:IE1,JS1:JE1,KS1:KE1) -UW(IS1:IE1,JS1:JE1,KS1:KE1)	&
                        +VN(IS1:IE1,JS1:JE1,KS1:KE1) -VS(IS1:IE1,JS1:JE1,KS1:KE1)	&
                        +WT(IS1:IE1,JS1:JE1,KS1:KE1) -WB(IS1:IE1,JS1:JE1,KS1:KE1))	&
	                    +DD

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BC_PLI
	CALL LISOLV(IS1,JS1,KS1,IE,JE,KE,PLI)
    CALL BC_PLI
        
!---- PLI CONVERGE TEST : RESIDUAL OF PLI ----------------------------------
	RESIP = 0.0D0
    DUM = 0.D0
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
	    RESOR = AE(I,J,K)*PLI(I+1,J,K) + AW(I,J,K)*PLI(I-1,J,K) +   &
				AN(I,J,K)*PLI(I,J+1,K) + AS(I,J,K)*PLI(I,J-1,K) +	&
				AT(I,J,K)*PLI(I,J,K+1) + AB(I,J,K)*PLI(I,J,K-1) +	&
				SU(I,J,K) - AP(I,J,K)*PLI(I,J,K)
		RESIP = RESIP + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*PLI(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESOR = RESIP / DUM	!	RESIP = RESIP / DUM ! SCALED RESIDUAL
	
	CALL MPI_ALLREDUCE(RESOR,RESIP,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ERR)
	
	IF(RESIP.GT.RESIMAX .AND. NITER.LT.NSWPP) GOTO 804
!---------------------------------------------------------------------------
	
	CALL MPI_X(PLI)
	
    END SUBROUTINE ITER_PLI