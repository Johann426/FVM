	MODULE ComDat_Shared

	USE ComDat_Param
    
	IMPLICIT NONE
	SAVE

	INTEGER :: I,J,K
	INTEGER :: FLAG,SIGN
	INTEGER :: NUM_GRID
    INTEGER :: ALLOC_ERR
	
	INTEGER :: NITER,N,MT 
	INTEGER :: NITER_TIME
	
	INTEGER :: SX,SX1,EX,EX1,SY,SY1,EY,EY1,SZ,SZ1,EZ,EZ1
	INTEGER :: NIM1,NJM1,NKM1
	INTEGER :: NIM2,NJM2,NKM2
	
	DOUBLE PRECISION :: pi
	
	DOUBLE PRECISION :: RESOR   !,RESORU,RESORV,RESORP
	DOUBLE PRECISION :: RESIU,RESIV,RESIP,RESIT
	
	DOUBLE PRECISION :: CFL_MAX,CFL_LOCAL
	
	DOUBLE PRECISION :: RTIME
	DOUBLE PRECISION :: TIME_BEGIN,TIME_END
	
	DOUBLE PRECISION ::	DUM,DUM_A,DUM_B
	DOUBLE PRECISION :: DUM_3D(0:NI,0:NJ,0:NK)
	
	DOUBLE PRECISION :: U(0:NI,0:NJ,0:NK),V(0:NI,0:NJ,0:NK),W(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: P(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: T(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: PLI(0:NI,0:NJ,0:NK)
    
    DOUBLE PRECISION :: U_F(0:NI,0:NJ,0:NK),V_F(0:NI,0:NJ,0:NK),W_F(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: T_F(0:NI,0:NJ,0:NK)
    
! N-TIME STEP VALUES ARE NEEDED UNTIL N+1 TIME STEP VALUES ARE CONVERGED IN JACOBI ITERATION METHOD
    DOUBLE PRECISION :: U_N(0:NI,0:NJ,0:NK),V_N(0:NI,0:NJ,0:NK),W_N(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: T_N(0:NI,0:NJ,0:NK)
!--------------------------------------------------------------------
!   CONTRA-VARIENT VELOCITY
	DOUBLE PRECISION :: U1(0:NI,0:NJ,0:NK),U2(0:NI,0:NJ,0:NK),U3(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: UE(0:NI,0:NJ,0:NK),UW(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: VN(0:NI,0:NJ,0:NK),VS(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: WT(0:NI,0:NJ,0:NK),WB(0:NI,0:NJ,0:NK)
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!   VORTICITY
    DOUBLE PRECISION :: VOR(0:NI,0:NJ,0:NK)
!--------------------------------------------------------------------
	DOUBLE PRECISION :: UG(0:NI,0:NJ,0:NK),VG(0:NI,0:NJ,0:NK),WG(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: TG(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: PG(0:NI,0:NJ,0:NK)
	
	DOUBLE PRECISION :: NLX(0:NI,0:NJ,0:NK),NLY(0:NI,0:NJ,0:NK),NLZ(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: NLX_NM1(0:NI,0:NJ,0:NK),NLY_NM1(0:NI,0:NJ,0:NK),NLZ_NM1(0:NI,0:NJ,0:NK)
	
	DOUBLE PRECISION :: ADV(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: ADV_NM1(0:NI,0:NJ,0:NK)
	
	DOUBLE PRECISION :: DIFX(0:NI,0:NJ,0:NK),DIFY(0:NI,0:NJ,0:NK),DIFZ(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: DIFT(0:NI,0:NJ,0:NK)
		
	DOUBLE PRECISION :: AP(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: AE(0:NI,0:NJ,0:NK),AW(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: AN(0:NI,0:NJ,0:NK),AS(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: AT(0:NI,0:NJ,0:NK),AB(0:NI,0:NJ,0:NK)
	
	DOUBLE PRECISION :: SU(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: FP(0:NI,0:NJ,0:NK)

    DOUBLE PRECISION :: S_U(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: S_V(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: S_W(0:NI,0:NJ,0:NK)
    
    DOUBLE PRECISION :: DE(0:NI,0:NJ,0:NK),DW(0:NI,0:NJ,0:NK),DN(0:NI,0:NJ,0:NK),DS(0:NI,0:NJ,0:NK) 
    DOUBLE PRECISION :: DT(0:NI,0:NJ,0:NK),DB(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: DD(0:NI,0:NJ,0:NK)

!---------------------------------------------------------------------------
!   GEOMETRY
!---------------------------------------------------------------------------
	DOUBLE PRECISION :: X(0:NI,0:NJ,0:NK) ,Y(0:NI,0:NJ,0:NK),Z(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: XG(0:NI,0:NJ,0:NK),YG(0:NI,0:NJ,0:NK),ZG(0:NI,0:NJ,0:NK)
    
    DOUBLE PRECISION :: RADIUS(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: THETA(0:NI,0:NJ,0:NK)
	
!---------------------------------------------------------------------------
!   METRIC_CENTER & METRIC_NODE
!---------------------------------------------------------------------------
	DOUBLE PRECISION :: DJAC(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: XXI(0:NI,0:NJ,0:NK),YXI(0:NI,0:NJ,0:NK),ZXI(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: XET(0:NI,0:NJ,0:NK),YET(0:NI,0:NJ,0:NK),ZET(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: XZT(0:NI,0:NJ,0:NK),YZT(0:NI,0:NJ,0:NK),ZZT(0:NI,0:NJ,0:NK)
    
	DOUBLE PRECISION :: XIX(0:NI,0:NJ,0:NK),XIY(0:NI,0:NJ,0:NK),XIZ(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: ETX(0:NI,0:NJ,0:NK),ETY(0:NI,0:NJ,0:NK),ETZ(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: ZTX(0:NI,0:NJ,0:NK),ZTY(0:NI,0:NJ,0:NK),ZTZ(0:NI,0:NJ,0:NK)

	DOUBLE PRECISION :: Q11(0:NI,0:NJ,0:NK),Q12(0:NI,0:NJ,0:NK),Q13(0:NI,0:NJ,0:NK)
	DOUBLE PRECISION :: Q21(0:NI,0:NJ,0:NK),Q22(0:NI,0:NJ,0:NK),Q23(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: Q31(0:NI,0:NJ,0:NK),Q32(0:NI,0:NJ,0:NK),Q33(0:NI,0:NJ,0:NK)
	
!---------------------------------------------------------------------------
!   COMMON DATA FOR CIRCULAR CYLINDER
!   SURFACE_AREA
	DOUBLE PRECISION :: AREA(1,1:NJ-1,1:NK-1)
    DOUBLE PRECISION :: SUM_AREA

!   FORM_DRAG_LIFT
	DOUBLE PRECISION :: P_DRAG
	DOUBLE PRECISION :: P_LIFT

!   FRICTION_DRAG_LIFT
	DOUBLE PRECISION :: DUDN(1,1:NJ-1,1:NK-1)
	DOUBLE PRECISION :: DTDN(1,1:NJ-1,1:NK-1)
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
    DOUBLE PRECISION :: U_TA(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: V_TA(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: P_TA(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: T_TA(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: VOR_TA(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: NU_TA,CD_TA,CL_TA,CPB_TA,ST_TA,CL_AMP_TA
    DOUBLE PRECISION :: RMS_CD,RMS_CL
    DOUBLE PRECISION :: CL1,CL_AMP,CL_NM1
    DOUBLE PRECISION :: TIME1,TIME2
    
    CHARACTER(20) :: FILENAME,FILENAME2
    
!---------------------------------------------------------------------------
    ENDMODULE ComDat_Shared