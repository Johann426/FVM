!***************************************************************************!
!    < Incompressible Navier-Stokes Equation in Curvilinear Coordinates >   !
!     3 dimensional viscous flow past a circular cylinder					!
!                                                                           !
!   - FVM : Cell center scheme                                              !
!   - Fractional 4-Step method                                              !
!	- Time-Dependent, incompressible, laminar flow                          !
!   - Constant materials                                                    !
!	- Based on colocated grid arrangement(a non-staggerd grid)              !
!                                                                           !
!     3D CYLINDER FORCED EXTERNAL FLOW  v1.0, Jin-Wook Lee                  !
!                                                                           !
!                           Last Modified : Feb.23, 2010                    !
!***************************************************************************!

	PROGRAM MAIN

	USE ComDat_Shared

	IMPLICIT NONE
	SAVE

!----------------------------------------------------------------------------
!	PARALLEL PROGRAM USING MPI(MASSAGE PASSING INTERFACE)
!    INCLUDE 'mpif.h'
!    
!    INTEGER :: err, rank, nprc, count, status(MPI_STATUS_SIZE)
!    
!	CALL MPI_INIT(err)
!	
!    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
!    
!    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,err)
!    
!    PRINT*, "PROCESSOR",RANK,"OF",NPRC,":CAN I TAKE YOUR ORDER?"
!----------------------------------------------------------------------------
 
    WRITE(*,*) 'Input restart parameter. (0=start : 1=restart)'
    READ (*,*)  IRESTART
    IF(IRESTART.EQ.0) THEN
	    WRITE(*,*) 'Input Reynolds Number.'
        RE =100		!READ (*,*)  RE
        WRITE(*,*) 'Input Prandtl  Number.'
        PR =0.7		!READ (*,*)  PR
        WRITE(*,*) 'Input rotating velocity.'
		VEL_ROT =0	!READ (*,*)  VEL_ROT
        WRITE(*,*) 'Input time step size.'
        DTIME =1.0d-3	!READ (*,*)  DTIME
	ELSE
        RTIME_MAX = RTIME_MAX * 10.D0
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
	W_N(:,:,:) = W(:,:,:)
	T_N(:,:,:) = T(:,:,:)
	PLI(:,:,:) = 0.0D0
 	
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
    CALL ITER_U ! WITH BCU
    CALL ITER_V ! WITH BCV
    CALL ITER_W ! WITH BCW
	CALL BCU
	CALL BCV
	CALL BCP

!   PRESSURE POISSON EQUATION : LINE GAUSS-GEIDEL ITERATION
    FLAG = 1        ! 0 : LINEAR INTERPOLATION / 1 : RHIE & CHOW INTERPOLATION
    CALL CONTRA_VEL
    CALL ITER_PLI   ! WITH BC_PLI

!   CORRECTION OF VELOCITY AND PRESSURE : CONTINUITY CONSTRAINT	
	CALL P_SOURCE(PLI)
	U = U+DTIME*S_U
	V = V+DTIME*S_V
	W = W+DTIME*S_W
    P=P+PLI

    DUM_A = P(NIM1,(NJ+1)/2,NK/2)
	P=P-DUM_A

!   CORRECTION FOR BOUNDARY VALUES AND BRANCH CUT
	CALL BCU
	CALL BCV
	CALL BCW
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!   POST-PROCESS SUBROUTINE
    CALL POST_PROCESS
    IF(SIGN) GOTO 777
    
!	CALL MPI_finalize(err)

    CALL CPU_TIME(TIME_END)
    
	WRITE(*,109) TIME_END-TIME_BEGIN
109 FORMAT(TR1'Time of operation = ',F10.0,TR4'seconds')
    PAUSE
	    
	END PROGRAM MAIN