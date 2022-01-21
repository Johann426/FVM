    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
    
    INCLUDE 'mpif.h'
    
	!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION

    AP(IS1:IE1,JS1:JE1,KS1:KE1) = 2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME				&
    							+ DE(IS1:IE1,JS1:JE1,KS1:KE1)+DW(IS1:IE1,JS1:JE1,KS1:KE1)	&
								+ DN(IS1:IE1,JS1:JE1,KS1:KE1)+DS(IS1:IE1,JS1:JE1,KS1:KE1)	&
								+ DT(IS1:IE1,JS1:JE1,KS1:KE1)+DB(IS1:IE1,JS1:JE1,KS1:KE1)
	
	AE(IS1:IE1,JS1:JE1,KS1:KE1) = DE(IS1:IE1,JS1:JE1,KS1:KE1)
	AW(IS1:IE1,JS1:JE1,KS1:KE1) = DW(IS1:IE1,JS1:JE1,KS1:KE1)
	AN(IS1:IE1,JS1:JE1,KS1:KE1) = DN(IS1:IE1,JS1:JE1,KS1:KE1)
	AS(IS1:IE1,JS1:JE1,KS1:KE1) = DS(IS1:IE1,JS1:JE1,KS1:KE1)
	AT(IS1:IE1,JS1:JE1,KS1:KE1) = DT(IS1:IE1,JS1:JE1,KS1:KE1)
	AB(IS1:IE1,JS1:JE1,KS1:KE1) = DB(IS1:IE1,JS1:JE1,KS1:KE1)
	
	!   TDMA solver
    
    NITER = 0 
805	NITER = NITER+1

    CALL DIF_DIAGONAL(T,DD) ! DIAGONAL DIF OF T : TREATED AS EXPLICIT

	!	SUMMATION OF EXPLICIT COMPONENTS
	
    SU(IS1:IE1,JS1:JE1,KS1:KE1)	=															&
				+2.D0*PE*DJAC(IS1:IE1,JS1:JE1,KS1:KE1)/DTIME*T_N(IS1:IE1,JS1:JE1,KS1:KE1)	&
                -3.D0*PE*ADV (IS1:IE1,JS1:JE1,KS1:KE1) +1.D0*PE*ADV_NM1(IS1:IE1,JS1:JE1,KS1:KE1)	&
		        +DIFT(IS1:IE1,JS1:JE1,KS1:KE1)												&
		        +DD(IS1:IE1,JS1:JE1,KS1:KE1)

    CALL BCT
	
	CALL LISOLV(IS,JS,KS,IE,JE,KE,T) ! WITHOUT OVER-RELAXATION
	
	!	T CONVERGE TEST : RESIDUAL OF T
	
	RESIT = 0.0D0
    DUM = 0.D0
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
	    RESOR = AE(I,J,K)*T(I+1,J,K) + AW(I,J,K)*T(I-1,J,K) +   &
				AN(I,J,K)*T(I,J+1,K) + AS(I,J,K)*T(I,J-1,K) +	&
				AT(I,J,K)*T(I,J,K+1) + AB(I,J,K)*T(I,J,K-1) +	&
				SU(I,J,K) - AP(I,J,K)*T(I,J,K)
		RESIT = RESIT + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*T(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESOR = RESIT / DUM	!	RESIT = RESIT / DUM ! SCALED RESIDUAL
	
	CALL MPI_ALLREDUCE(RESOR,RESIT,1,MPI_DOUBLE_PRECISION,MPI_MAX,NMCW,ERR)
	
	IF(RESIT.GT.RESIMAX .AND. NITER.LT.NSWPT) GOTO 805
!---------------------------------------------------------------------------
    END SUBROUTINE ITER_T