!                                                                           !
!    < Two dimensional viscous flow past a circular cylinder >              !
!                                                                           !
!   - FVM : Cell center scheme                                              !
!   - Fractional 4-Step method                                              !
!	- Time-Dependent, incompressible, laminar flow                          !
!   - Constant materials                                                    !
!	- Based on colocated grid arrangement(a non-staggerd grid)              !
!                                                                           !
!     2D CYLINDER FORCED EXTERNAL FLOW  v2.3, Jin-Wook Lee                  !
!                                                                           !
!                           Last Modified : May.20, 2010                    !
!                                                                           !
	PROGRAM MAIN

	USE ComDat_Shared

	IMPLICIT NONE
	SAVE
	
	PRINT *, ""
	
    PRINT *, "Input restart parameter. (0=start : 1=restart)"
    
    READ (*,*)  IRESTART
    
    IF(IRESTART.EQ.0) THEN
	    WRITE(*,*) 'Input Reynolds Number.'
        READ (*,*)  RE
        WRITE(*,*) 'Input Prandtl  Number.'
        READ (*,*)  PR
        WRITE(*,*) 'Input rotating velocity.'
        READ (*,*)  VEL_ROT
        WRITE(*,*) 'Input time step size.'
        READ (*,*)  DTIME
	ENDIF
	
    CALL CPU_TIME(TIME_BEGIN)
    
	CALL GRID
    
	CALL METRIC
    
	CALL INIT   ! WITH B.C.
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ITERATION FOR THE TIME ADVANCEMENT
777	NITER_TIME = NITER_TIME+1

	RTIME      = RTIME + DTIME 

	U_N(:,:,:) = U(:,:,:)
	V_N(:,:,:) = V(:,:,:)
	T_N(:,:,:) = T(:,:,:)
	
	PLI(:,:,:) = 0.0D0
 	
!   CONTRA-VARIENT VELOCITY : CONSERVATIVE FORM
    FLAG = 0
    CALL CONTRA_VEL

!   ADVECTIVE EXPLICIT TERM : CONSERVATIVE FORM
    NLX_NM1 = NLX
    NLY_NM1 = NLY
    ADV_NM1 = ADV
    CALL NL(U_N,NLX)
    CALL NL(V_N,NLY)
    CALL NL(T_N,ADV)	!   THERMAL ADVECTIVE EXPLICIT TERM : CONSERVATIVE FORM

!   DIFFUSIVE EXPLICIT TERM : CONSERVATIVE FORM
    CALL DIF(U_N,DIFX)
    CALL DIF(V_N,DIFY)
    CALL DIF(T_N,DIFT)	!   THERMAL DIFFUSIVE EXPLICIT TERM : CONSERVATIVE FORM
    
!   PRESSURE GRADIENT EXPLICIT TERM : NON-CONSERVATIVE FORM
    CALL P_SOURCE

!   INTERMEDIATE VELOCITY : LINE GAUSS-GEIDEL ITERATION
    CALL ITER_U ! WITH BCU
    CALL ITER_V ! WITH BCV
	CALL BCU
	CALL BCV
	CALL BCP
	
!   PRESSURE POISSON EQUATION : LINE GAUSS-GEIDEL ITERATION
    FLAG = 1        ! 0 : LINEAR INTERPOLATION / 1 : RHIE & CHOW INTERPOLATION
    CALL CONTRA_VEL
    CALL ITER_PLI   ! WITH BC_PLI
	CALL BC_PLI
	
!   CORRECTION OF VELOCITY AND PRESSURE
	CALL UVP_CORRECTION
!   CORRECTION FOR BOUNDARY VALUES AND BRANCH CUT
	CALL BCU
	CALL BCV
	CALL BCP

!   INTERMIDIATE VELOCITY : LINE GAUSS-GEIDEL ITERATION
    CALL ITER_T ! WITH BCT
    CALL BCT

!   CALCULATE COEFFICIENTS
    CALL CALCCOEF
    
    CALL WRITEDATA
    
    IF(RTIME.LT.RTIME_MAX) GOTO 777
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!   POST-PROCESS SUBROUTINE
    CALL POST_PROCESS
    IF(SIGN) GOTO 777
    
    CALL CPU_TIME(TIME_END)
    
	PRINT *, ""
	
	WRITE(*,109) TIME_END-TIME_BEGIN
	
109 FORMAT(TR1'Time of operation = ',F10.0,TR4'seconds')

    PAUSE
    
	END PROGRAM MAIN