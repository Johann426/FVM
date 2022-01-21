    SUBROUTINE ITER_T
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
    
    INCLUDE 'mpif.h'
    
	!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION

    AP(IS1:IE1,JS1:JE1,KS1:KE1) = FP	+DE(IS1:IE1,JS1:JE1,KS1:KE1)+DW(IS1:IE1,JS1:JE1,KS1:KE1)	&
										+DN(IS1:IE1,JS1:JE1,KS1:KE1)+DS(IS1:IE1,JS1:JE1,KS1:KE1)	&
										+DT(IS1:IE1,JS1:JE1,KS1:KE1)+DB(IS1:IE1,JS1:JE1,KS1:KE1)
	
	AE(IS1:IE1,JS1:JE1,KS1:KE1) = DE(IS1:IE1,JS1:JE1,KS1:KE1)
	AW(IS1:IE1,JS1:JE1,KS1:KE1) = DW(IS1:IE1,JS1:JE1,KS1:KE1)
	AN(IS1:IE1,JS1:JE1,KS1:KE1) = DN(IS1:IE1,JS1:JE1,KS1:KE1)
	AS(IS1:IE1,JS1:JE1,KS1:KE1) = DS(IS1:IE1,JS1:JE1,KS1:KE1)
	AT(IS1:IE1,JS1:JE1,KS1:KE1) = DT(IS1:IE1,JS1:JE1,KS1:KE1)
	AB(IS1:IE1,JS1:JE1,KS1:KE1) = DB(IS1:IE1,JS1:JE1,KS1:KE1)
	
	!   TDMA solver
    
    CALL DIF_DIAGONAL(T,DD) ! DIAGONAL DIF OF T : TREATED AS EXPLICIT

	!	SUMMATION OF EXPLICIT COMPONENTS
	
    SU(IS1:IE1,JS1:JE1,KS1:KE1)	= FP*T_N(IS1:IE1,JS1:JE1,KS1:KE1)									&
                -3.D0*PE*ADV (IS1:IE1,JS1:JE1,KS1:KE1) +1.D0*PE*ADV_NM1(IS1:IE1,JS1:JE1,KS1:KE1)	&
		        +1.D0*PE*DIFT(IS1:IE1,JS1:JE1,KS1:KE1)												&
		        +DD(IS1:IE1,JS1:JE1,KS1:KE1)

    CALL BCT
	
	CALL LISOLV(IS,JS,KS,IE,JE,KE,T) ! WITHOUT OVER-RELAXATION
	
    END SUBROUTINE ITER_T