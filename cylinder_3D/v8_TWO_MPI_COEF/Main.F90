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
!	  3D CYLINDER FORCED EXTERNAL FLOW  v5.0, Jin-Wook Lee					!
!																			!
!							Last Modified : Feb.23, 2010					!
!***************************************************************************!

	PROGRAM MAIN

	USE ComDat_Shared

	IMPLICIT NONE
	SAVE
    
	INCLUDE 'mpif.h'

	CALL DEFINE
	
    CALL CPU_TIME(TIME_BEGIN)
    
	CALL GRID

	CALL METRIC

	CALL INIT   ! WITH B.C.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ITERATION FOR THE TIME ADVANCEMENT
777	NITER_TIME = NITER_TIME+1
	RTIME      = RTIME + DTIME 

	U_N(IS:IE,JS:JE,KS:KE) = U(IS:IE,JS:JE,KS:KE)
	V_N(IS:IE,JS:JE,KS:KE) = V(IS:IE,JS:JE,KS:KE)
	W_N(IS:IE,JS:JE,KS:KE) = W(IS:IE,JS:JE,KS:KE)
	T_N(IS:IE,JS:JE,KS:KE) = T(IS:IE,JS:JE,KS:KE)
	PLI(IS:IE,JS:JE,KS:KE) = 0.0D0

!   fp = 2RE/DT
    fp = 2.D0*RE/DTIME

!   CONTRA-VARIENT VELOCITY : CONSERVATIVE FORM
    FLAG = 0
    CALL CONTRA_VEL

!   ADVECTIVE EXPLICIT TERM : CONSERVATIVE FORM
    NLX_NM1 = NLX
    NLY_NM1 = NLY
    NLZ_NM1 = NLZ
    CALL NL(U_N,NLX)
    CALL NL(V_N,NLY)
    CALL NL(W_N,NLZ)

!   DIFFUSIVE EXPLICIT TERM : CONSERVATIVE FORM
    CALL DIF(U_N,DIFX)
    CALL DIF(V_N,DIFY)
    CALL DIF(W_N,DIFZ)

!   PRESSURE GRADIENT EXPLICIT TERM : NON-CONSERVATIVE FORM
    CALL P_SOURCE(P)

!   INTERMEDIATE VELOCITY : LINE GAUSS-GEIDEL ITERATION
	UIN = 1.D0
    CALL ITER_U(U,U_N,NLX,NLX_NM1,DIFX,S_U)
	
    UIN = 0.D0
    CALL ITER_U(V,V_N,NLY,NLY_NM1,DIFY,S_V)
    	
    CALL ITER_U(W,W_N,NLZ,NLZ_NM1,DIFZ,S_W)

    UIN = 1.D0
    CALL BCU(U,U_N)
	UIN = 0.D0
    CALL BCU(V,V_N)
    CALL BCU(W,W_N)

!   PRESSURE POISSON EQUATION : LINE GAUSS-GEIDEL ITERATION
    FLAG = 1        ! 0 : LINEAR INTERPOLATION / 1 : RHIE & CHOW INTERPOLATION
    CALL CONTRA_VEL
    CALL ITER_PLI   ! WITH BC_PLI
	
!   CORRECTION OF VELOCITY AND PRESSURE : CONTINUITY CONSTRAINT	
	CALL P_SOURCE(PLI)
	U(IS1:IE1,JS1:JE1,KS1:KE1) = U(IS1:IE1,JS1:JE1,KS1:KE1)+DTIME*S_U(IS1:IE1,JS1:JE1,KS1:KE1)
	V(IS1:IE1,JS1:JE1,KS1:KE1) = V(IS1:IE1,JS1:JE1,KS1:KE1)+DTIME*S_V(IS1:IE1,JS1:JE1,KS1:KE1)
	W(IS1:IE1,JS1:JE1,KS1:KE1) = W(IS1:IE1,JS1:JE1,KS1:KE1)+DTIME*S_W(IS1:IE1,JS1:JE1,KS1:KE1)
    P(IS1:IE1,JS1:JE1,KS1:KE1) = P(IS1:IE1,JS1:JE1,KS1:KE1)+PLI(IS1:IE1,JS1:JE1,KS1:KE1)
	
!	CALL MPI_X(U)
!	CALL MPI_X(V)
!	CALL MPI_X(W)

	CALL MPI_X(P)
	
    DUM_A = P(IE1,(JE+1)/2,(KS+KE)/2)
	P(IS:IE,JS:JE,KS:KE)=P(IS:IE,JS:JE,KS:KE) - DUM_A

!   CORRECTION FOR BOUNDARY VALUES AND BRANCH CUT
    UIN = 1.D0
    CALL BCU(U,U_N)
	UIN = 0.D0
    CALL BCU(V,V_N)
    CALL BCU(W,W_N)
	CALL BCP

!   fp = 2PE/DT
    fp = 2.D0*PE/DTIME
	
!   CONTRA-VARIENT VELOCITY : CONSERVATIVE FORM
    FLAG = 0
    CALL CONTRA_VEL
    
!   THERMAL ADVECTIVE EXPLICIT TERM : CONSERVATIVE FORM
    ADV_NM1 = ADV
    CALL NL(T,ADV)
	
!   THERMAL DIFFUSIVE EXPLICIT TERM : CONSERVATIVE FORM
    CALL DIF(T,DIFT)
	
!   INTERMIDIATE VELOCITY : LINE GAUSS-GEIDEL ITERATION
    CALL ITER_T ! WITH BCT
    CALL BCT
	
!   CALCULATE COEFFICIENTS
    CALL CALCCOEF
	
    CALL WRITEDATA
    
    IF(RTIME.LT.RTIME_MAX) GOTO 777
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
!   POST-PROCESS SUBROUTINE
    CALL POST_PROCESS
    IF(SIGN) GOTO 777
    
	CALL CPU_TIME(TIME_END)
    
    PRINT *, ""
    
	WRITE(*,109) TIME_END-TIME_BEGIN
	
109 FORMAT(TR1'Time of operation = ',F10.0,TR4'seconds')

	CALL MPI_FINALIZE(ERR)

	END PROGRAM MAIN