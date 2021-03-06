	MODULE ComDat_Shared

	USE ComDat_Param
    
	IMPLICIT NONE
	SAVE

	INTEGER :: I,J,K
	INTEGER :: FLAG
	INTEGER :: NUM_GRID
    INTEGER :: NUM_WRITE
    INTEGER :: ALLOC_ERR
	
	INTEGER :: NITER,N
	INTEGER :: NITER_TIME
	
	INTEGER :: sx,sx1,ex,ex1,sy,sy1,ey,ey1
	INTEGER :: NIM1,NJM1,NKM1
	INTEGER :: NIM2,NJM2,NKM2
	
	LOGICAL :: SIGN
		
	CHARACTER(11) :: FILENAME
	CHARACTER(13) :: CONTOUR
		
	DOUBLE PRECISION :: PI
	
	DOUBLE PRECISION :: RESOR   !,RESORU,RESORV,RESORP
	DOUBLE PRECISION :: RESIU,RESIV,RESIP,RESIT
	
	DOUBLE PRECISION :: CFL_MAX,CFL_LOCAL
	
	DOUBLE PRECISION :: RTIME
	DOUBLE PRECISION :: TIME_BEGIN,TIME_END
	
	DOUBLE PRECISION ::	DUM,DUM_A,DUM_B,DUM_K

	DOUBLE PRECISION :: U(0:NI,0:NJ,2:2),V(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: P(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: T(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: PLI(0:NI,0:NJ,2:2)
    
    DOUBLE PRECISION :: U_F,V_F
	DOUBLE PRECISION :: T_F
    
! N-TIME STEP VALUES ARE NEEDED UNTIL N+1 TIME STEP VALUES ARE CONVERGED IN JACOBI ITERATION METHOD
    DOUBLE PRECISION :: U_N(0:NI,0:NJ,2:2),V_N(0:NI,0:NJ,2:2)   
    DOUBLE PRECISION :: T_N(0:NI,0:NJ,2:2)
!--------------------------------------------------------------------
!   CONTRA-VARIENT VELOCITY
	DOUBLE PRECISION :: U1(0:NI,0:NJ,2:2),U2(0:NI,0:NJ,2:2)     
	DOUBLE PRECISION :: UE(1:NI,1:NJ,2:2),UW(1:NI,1:NJ,2:2)
	DOUBLE PRECISION :: VN(1:NI,1:NJ,2:2),VS(1:NI,1:NJ,2:2)
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!   VORTICITY
    DOUBLE PRECISION :: VOR(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: VOR_G(0:NI,0:NJ,2:2) 
!--------------------------------------------------------------------
	DOUBLE PRECISION :: UG(0:NI,0:NJ,2:2),VG(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: TG(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: PG(0:NI,0:NJ,2:2)
	
	DOUBLE PRECISION :: NLX(0:NI,0:NJ,2:2),NLY(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: NLX_NM1(0:NI,0:NJ,2:2),NLY_NM1(0:NI,0:NJ,2:2)

	DOUBLE PRECISION :: DIFX(0:NI,0:NJ,2:2),DIFY(0:NI,0:NJ,2:2)

	DOUBLE PRECISION :: ADV(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: ADV_NM1(0:NI,0:NJ,2:2)
	
    DOUBLE PRECISION :: DIFT(0:NI,0:NJ,2:2)
    
	DOUBLE PRECISION :: DUM_3D(0:NI,0:NJ,2:2) 
	DOUBLE PRECISION ::	DUM_0 (0:NI,0:NJ,2:2) 
	DOUBLE PRECISION ::	DUM_1 (0:NI,0:NJ,2:2) 
	DOUBLE PRECISION ::	DUM_2 (0:NI,0:NJ,2:2) 

	DOUBLE PRECISION :: AP(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: AE(0:NI,0:NJ,2:2),AW(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: AN(0:NI,0:NJ,2:2),AS(0:NI,0:NJ,2:2)

	DOUBLE PRECISION :: SU(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: FP(0:NI,0:NJ,2:2)

    DOUBLE PRECISION :: S_U(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: S_V(0:NI,0:NJ,2:2)
    
    DOUBLE PRECISION :: DE(0:NI,0:NJ,2:2),DW(0:NI,0:NJ,2:2),DN(0:NI,0:NJ,2:2),DS(0:NI,0:NJ,2:2) 
    DOUBLE PRECISION :: DD(0:NI,0:NJ,2:2)

!---------------------------------------------------------------------------
!   GEOMETRY
!---------------------------------------------------------------------------
	DOUBLE PRECISION :: X(0:NI,0:NJ,2:2) ,Y(0:NI,0:NJ,2:2),Z(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: XG(1:NI,1:NJ,2:2),YG(1:NI,1:NJ,2:2),ZG(1:NI,1:NJ,2:2)
    
    DOUBLE PRECISION :: RADIUS(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: THETA(0:NI,0:NJ,2:2)
	!DOUBLE PRECISION :: DTHETA(0:NI,0:NJ,2:2)
	
!---------------------------------------------------------------------------
!   METRIC_CENTER & METRIC_NODE
!---------------------------------------------------------------------------
	DOUBLE PRECISION, ALLOCATABLE :: XXI(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: YXI(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: XET(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: YET(:,:,:)
	
	DOUBLE PRECISION:: DJAC(0:NI,0:NJ,2:2)
	
	DOUBLE PRECISION :: XIX(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: XIY(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: ETX(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: ETY(0:NI,0:NJ,2:2)

	DOUBLE PRECISION :: Q11(0:NI,0:NJ,2:2),Q12(0:NI,0:NJ,2:2)
	DOUBLE PRECISION :: Q21(0:NI,0:NJ,2:2),Q22(0:NI,0:NJ,2:2)

!---------------------------------------------------------------------------
!   COMMON DATA FOR CIRCULAR CYLINDER
!   SURFACE_AREA
	DOUBLE PRECISION :: AREA(1,1:NJ-1,2:2)
    DOUBLE PRECISION :: SUM_AREA

!   FORM_DRAG_LIFT
	DOUBLE PRECISION :: P_DRAG
	DOUBLE PRECISION :: P_LIFT

!   FRICTION_DRAG_LIFT
	DOUBLE PRECISION :: DUDN(1,1:NJ,2:2)
	DOUBLE PRECISION :: DTDN(1,1:NJ,2:2)
	DOUBLE PRECISION :: T_DRAG
	DOUBLE PRECISION :: T_LIFT

!   CALCCOEF
	DOUBLE PRECISION :: ST
	DOUBLE PRECISION :: NU,CD,CL
	DOUBLE PRECISION :: CD_P,CD_V,CL_P,CL_V,CPB
	DOUBLE PRECISION :: SUM_NU
	

    DOUBLE PRECISION :: L1,L2,L3
    DOUBLE PRECISION :: TH1,TH2
    DOUBLE PRECISION :: DN1,DN2
    DOUBLE PRECISION :: LA,LB,LA_,LB_
    DOUBLE PRECISION :: X1,Y1

!   POST-PROCESSING
    INTEGER :: IPOST
    INTEGER :: TOGGLE
    INTEGER :: NITER_POS
    DOUBLE PRECISION :: U_TA(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: V_TA(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: P_TA(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: T_TA(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: VOR_TA(0:NI,0:NJ,2:2)
    DOUBLE PRECISION :: NU_TA,CD_TA,CL_TA,CPB_TA,ST_TA,CL_AMP_TA
    DOUBLE PRECISION :: RMS_CD,RMS_CL
    DOUBLE PRECISION :: CL1,CL_AMP,CL_NM1
    DOUBLE PRECISION :: TIME1,TIME2
    
!---------------------------------------------------------------------------
    ENDMODULE ComDat_Shared