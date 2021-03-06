!***************************************************************************!
!	 < Incompressible Navier-Stokes Equation in Curvilinear Coordinates >	!
!	  Three dimensional viscous flow past a circular cylinder				!
!																			!
!	- FVM : Cell center scheme												!
!	- Fractional 4-Step method												!
!	- Time-Dependent, incompressible, laminar flow							!
!	- Constant materials													!
!	- Based on colocated grid arrangement(a non-staggerd grid)				!
!																			!
!	  3D CYLINDER FORCED EXTERNAL FLOW  v19.2, Jin-Wook Lee					!
!																			!
!							Last Modified : Jul.05, 2010					!
!***************************************************************************!

	PROGRAM MAIN

	USE ComDat_Shared

	IMPLICIT NONE
	SAVE
    
	INCLUDE 'mpif.h'

	CALL DEFINE

	CALL GRID

	CALL METRIC

	CALL INIT   ! WITH B.C.
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ITERATION FOR THE TIME ADVANCEMENT
777	NITER_TIME = NITER_TIME+1
	RTIME      = RTIME + DTIME 
	
	CALL CPU_TIME(TIME_BEGIN)
	
	U_N(IS:IE,JS:JE,KS:KE) = U(IS:IE,JS:JE,KS:KE)
	V_N(IS:IE,JS:JE,KS:KE) = V(IS:IE,JS:JE,KS:KE)
	W_N(IS:IE,JS:JE,KS:KE) = W(IS:IE,JS:JE,KS:KE)
	T_N(IS:IE,JS:JE,KS:KE) = T(IS:IE,JS:JE,KS:KE)

!   CONTRA-VARIENT VELOCITY : CONSERVATIVE FORM
    CALL CONTRA_VEL

!   ADVECTIVE EXPLICIT TERM : CONSERVATIVE FORM
    NLX_NM1(IS1:IE1,JS1:JE1,KS1:KE1) = NLX(IS1:IE1,JS1:JE1,KS1:KE1)
    NLY_NM1(IS1:IE1,JS1:JE1,KS1:KE1) = NLY(IS1:IE1,JS1:JE1,KS1:KE1)
    NLZ_NM1(IS1:IE1,JS1:JE1,KS1:KE1) = NLZ(IS1:IE1,JS1:JE1,KS1:KE1)
    ADV_NM1(IS1:IE1,JS1:JE1,KS1:KE1) = ADV(IS1:IE1,JS1:JE1,KS1:KE1)
    CALL NL(U_N,NLX)
    CALL NL(V_N,NLY)
    CALL NL(W_N,NLZ)
    CALL NL(T_N,ADV)	!   THERMAL ADVECTIVE EXPLICIT TERM : CONSERVATIVE FORM

!   DIFFUSIVE EXPLICIT TERM : CONSERVATIVE FORM
    CALL DIF(U_N,DIFX)
    CALL DIF(V_N,DIFY)
    CALL DIF(W_N,DIFZ)
	CALL DIF(T_N,DIFT)	!   THERMAL DIFFUSIVE EXPLICIT TERM : CONSERVATIVE FORM
	
!   PRESSURE GRADIENT EXPLICIT TERM : NON-CONSERVATIVE FORM
    CALL P_SOURCE(P)

!   INTERMEDIATE VELOCITY : LINE GAUSS-GEIDEL ITERATION
	
	UIN = 1.D0
    CALL ITER_U(U,U_N,NLX,NLX_NM1,DIFX,S_U)
	
    UIN = 0.D0
    CALL ITER_U(V,V_N,NLY,NLY_NM1,DIFY,S_V)
    	
    CALL ITER_U(W,W_N,NLZ,NLZ_NM1,DIFZ,S_W)
	
	U(IS1:IE1,JS1:JE1,KS1:KE1) = U(IS1:IE1,JS1:JE1,KS1:KE1)-DTIME*S_U(IS1:IE1,JS1:JE1,KS1:KE1)/DJAC(IS1:IE1,JS1:JE1,KS1:KE1)
	V(IS1:IE1,JS1:JE1,KS1:KE1) = V(IS1:IE1,JS1:JE1,KS1:KE1)-DTIME*S_V(IS1:IE1,JS1:JE1,KS1:KE1)/DJAC(IS1:IE1,JS1:JE1,KS1:KE1)
	W(IS1:IE1,JS1:JE1,KS1:KE1) = W(IS1:IE1,JS1:JE1,KS1:KE1)-DTIME*S_W(IS1:IE1,JS1:JE1,KS1:KE1)/DJAC(IS1:IE1,JS1:JE1,KS1:KE1)
    
    UIN = 1.D0
    CALL BCU(U,U_N)
	UIN = 0.D0
    CALL BCU(V,V_N)
    CALL BCU(W,W_N)
        
!   PRESSURE POISSON EQUATION : LINE GAUSS-GEIDEL ITERATION
    CALL CONTRA_VEL

    CALL ITER_P   ! WITH BCP
	
	CALL BCP
	
!   CORRECTION OF VELOCITY AND PRESSURE : CONTINUITY CONSTRAINT	
	CALL P_SOURCE(P)
	U(IS1:IE1,JS1:JE1,KS1:KE1) = U(IS1:IE1,JS1:JE1,KS1:KE1)+DTIME*S_U(IS1:IE1,JS1:JE1,KS1:KE1)/DJAC(IS1:IE1,JS1:JE1,KS1:KE1)
	V(IS1:IE1,JS1:JE1,KS1:KE1) = V(IS1:IE1,JS1:JE1,KS1:KE1)+DTIME*S_V(IS1:IE1,JS1:JE1,KS1:KE1)/DJAC(IS1:IE1,JS1:JE1,KS1:KE1)
	W(IS1:IE1,JS1:JE1,KS1:KE1) = W(IS1:IE1,JS1:JE1,KS1:KE1)+DTIME*S_W(IS1:IE1,JS1:JE1,KS1:KE1)/DJAC(IS1:IE1,JS1:JE1,KS1:KE1)
	
	DUM_B = 0.D0
	IF(JS.LT.NJ/2 .AND. JE.GT.NJ/2) THEN
		DO K=KS1,KE1
			DUM_B = DUM_B + 0.25D0*(P(IE,NJ/2,K)+P(IE1,NJ/2,K)+P(IE1,NJ/2+1,K)+P(IE,NJ/2+1,K))
		ENDDO
			DUM_B = DUM_B/DBLE(KE1-KS1+1)
	ENDIF
	
	CALL MPI_ALLREDUCE(DUM_B,DUM_A,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)	! REF. PRESSURE ON THE INLET BNDR
	
	P(IS:IE,JS:JE,KS:KE)=P(IS:IE,JS:JE,KS:KE) - DUM_A / NDIMK

!   CORRECTION FOR BOUNDARY VALUES AND BRANCH CUT
    UIN = 1.D0
    CALL BCU(U,U_N)
	UIN = 0.D0
    CALL BCU(V,V_N)
    CALL BCU(W,W_N)
	CALL BCP
	
!   INTERMIDIATE VELOCITY : LINE GAUSS-GEIDEL ITERATION
    CALL ITER_T ! WITH BCT
    
    CALL BCT
	
!   CALCULATE COEFFICIENTS
	CALL CALCCOEF
	
	CALL MPI_BARRIER(NMCW,ERR)
	
	CALL WRITEDATA
	
	CALL CPU_TIME(TIME_END)

	IF(RANK.EQ.0) WRITE(*,109) TIME_END-TIME_BEGIN
109 FORMAT(TR1'Time of operation = ',F6.3,TR4'seconds')

	IF(RTIME.LT.RTIME_MAX) GOTO 777
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
!   POST-PROCESS SUBROUTINE
    CALL POST_PROCESS
    IF(SIGN) GOTO 777
    
	CALL MPI_FINALIZE(ERR)
	
	END PROGRAM MAIN