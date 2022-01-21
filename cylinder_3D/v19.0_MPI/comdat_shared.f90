	MODULE ComDat_Shared

	USE ComDat_Param
    
	IMPLICIT NONE
	SAVE

	!   MPI(Massage Passing Interface) VARIABLES
    INTEGER :: RANK, NPRC, ERR, NMCW 
    INTEGER :: AROW, CLMN, JPLN, KPLN
    INTEGER :: IPRV, INXT
    INTEGER :: SIZ
    INTEGER :: DIMS(0:NDIMS-1),COORDS(0:NDIMS-1)
    INTEGER :: NBR_S,NBR_N,NBR_B,NBR_T
    INTEGER :: ISEND1, ISEND2, ISEND3, ISEND4
    INTEGER :: IRECV1, IRECV2, IRECV3, IRECV4
    INTEGER, ALLOCATABLE :: DSIZ(:), DISP(:)

    LOGICAL :: REORDER(0:NDIMS-1)
    LOGICAL :: PERIODS(0:NDIMS-1)
    
	INTEGER :: I,J,K
	INTEGER :: FLAG
	INTEGER :: NUM_GRID
	INTEGER :: NUM_WRITE
	INTEGER :: ALLOC_ERR
	INTEGER :: INTA,INTB
	
	INTEGER :: N,NITER
	INTEGER :: NITER_TIME
	
	INTEGER :: NIM1,NJM1,NKM1
	INTEGER :: NIM2,NJM2,NKM2
	INTEGER :: IS,JS,KS,IE,JE,KE,IS1,JS1,KS1,IE1,JE1,KE1
		
	LOGICAL :: SIGN
	
    CHARACTER(15) :: FILENAME
    CHARACTER(14) :: FILENAME2
    CHARACTER(11) :: CONTOUR
    
	DOUBLE PRECISION :: PI
	
	DOUBLE PRECISION :: RESOR
	DOUBLE PRECISION :: RESIP,RESIU,RESIT
	
	DOUBLE PRECISION :: CFL_MAX,CFL_LOCAL
	
	DOUBLE PRECISION :: RTIME
	DOUBLE PRECISION :: TIME_BEGIN,TIME_END
	
	DOUBLE PRECISION ::	DUM,DUM_A,DUM_B,DUM_C,DUM_D
	
	DOUBLE PRECISION :: U_F,T_F
	
	DOUBLE PRECISION, ALLOCATABLE :: U(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: V(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: W(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: P(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: T(:,:,:)
    
    DOUBLE PRECISION, ALLOCATABLE :: UG(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: VG(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: WG(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: PG(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: TG(:,:,:)
	
	DOUBLE PRECISION, ALLOCATABLE :: VOR_X(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: VOR_Y(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: VOR_Z(:,:,:)
	
	! N-TIME STEP VALUES ARE NEEDED FOR PROPER BOUNDARY VALUE
    DOUBLE PRECISION, ALLOCATABLE :: U_N(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: V_N(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: W_N(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: T_N(:,:,:)
	
	DOUBLE PRECISION, ALLOCATABLE :: PLI(:,:,:)

	DOUBLE PRECISION, ALLOCATABLE :: DUM_3D(:,:,:)
!--------------------------------------------------------------------
!   CONTRA-VARIENT VELOCITY
	DOUBLE PRECISION, ALLOCATABLE :: U1(:,:,:),U2(:,:,:),U3(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: UE(:,:,:),UW(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: VN(:,:,:),VS(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: WT(:,:,:),WB(:,:,:)
!--------------------------------------------------------------------

	DOUBLE PRECISION, ALLOCATABLE :: NLX(:,:,:),NLY(:,:,:),NLZ(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: NLX_NM1(:,:,:),NLY_NM1(:,:,:),NLZ_NM1(:,:,:)
	
	DOUBLE PRECISION, ALLOCATABLE :: ADV(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: ADV_NM1(:,:,:)
	
	DOUBLE PRECISION, ALLOCATABLE :: DIFX(:,:,:),DIFY(:,:,:),DIFZ(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: DIFT(:,:,:)
		
    DOUBLE PRECISION, ALLOCATABLE :: S_U(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: S_V(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: S_W(:,:,:)
    
	DOUBLE PRECISION, ALLOCATABLE :: AP(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: AE(:,:,:),AW(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: AN(:,:,:),AS(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: AT(:,:,:),AB(:,:,:)
	
    DOUBLE PRECISION, ALLOCATABLE :: DE(:,:,:),DW(:,:,:),DN(:,:,:),DS(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: DT(:,:,:),DB(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: DD(:,:,:)

	DOUBLE PRECISION, ALLOCATABLE :: SU(:,:,:)
!---------------------------------------------------------------------------
!   GEOMETRY
!---------------------------------------------------------------------------
	DOUBLE PRECISION :: X_BUF
	DOUBLE PRECISION :: Y_BUF
	DOUBLE PRECISION :: Z_BUF
	
	DOUBLE PRECISION, ALLOCATABLE :: X(:,:,:) ,Y(:,:,:),Z(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: XG(:,:,:),YG(:,:,:),ZG(:,:,:)
    
    DOUBLE PRECISION, ALLOCATABLE :: RADIUS(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: THETA(:,:,:)
	
!---------------------------------------------------------------------------
!   METRIC_CENTER & METRIC_NODE
!---------------------------------------------------------------------------
	DOUBLE PRECISION, ALLOCATABLE :: DJAC(:,:,:)
	
	DOUBLE PRECISION, ALLOCATABLE :: XXI(:,:,:),YXI(:,:,:),ZXI(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: XET(:,:,:),YET(:,:,:),ZET(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: XZT(:,:,:),YZT(:,:,:),ZZT(:,:,:)
    
	DOUBLE PRECISION, ALLOCATABLE :: XIX(:,:,:),XIY(:,:,:),XIZ(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: ETX(:,:,:),ETY(:,:,:),ETZ(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: ZTX(:,:,:),ZTY(:,:,:),ZTZ(:,:,:)

	DOUBLE PRECISION, ALLOCATABLE :: Q11(:,:,:),Q12(:,:,:),Q13(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE :: Q21(:,:,:),Q22(:,:,:),Q23(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Q31(:,:,:),Q32(:,:,:),Q33(:,:,:)
	
!---------------------------------------------------------------------------
!   COMMON DATA FOR CIRCULAR CYLINDER
!   SURFACE_AREA,DUDN
	DOUBLE PRECISION, ALLOCATABLE :: DUDN(:,:,:)	!(1,1:NJ-1,1:NK-1)
	DOUBLE PRECISION, ALLOCATABLE :: DTDN(:,:,:)	!(1,1:NJ-1,1:NK-1)
	DOUBLE PRECISION, ALLOCATABLE :: AREA1(:,:,:)	!(1,1:NJ-1,1:NK-1)
	DOUBLE PRECISION, ALLOCATABLE :: AREA2(:,:,:)	!(1,1:NJ-1,1:NK-1)
    DOUBLE PRECISION :: SUM_AREA
!   FRICTION_DRAG_LIFT
	DOUBLE PRECISION :: T_DRAG
	DOUBLE PRECISION :: T_LIFT_Y
	DOUBLE PRECISION :: T_LIFT_Z
!   FORM_DRAG_LIFT
	DOUBLE PRECISION :: P_DRAG
	DOUBLE PRECISION :: P_LIFT_Y
	DOUBLE PRECISION :: P_LIFT_Z
!   CALCCOEF
	DOUBLE PRECISION :: ST
	DOUBLE PRECISION :: NU,CD,CL
	DOUBLE PRECISION :: CD_P,CD_V,CL_P,CL_V,CPB
	DOUBLE PRECISION :: SUM_NU
!	LINE,ANGLE,DN
    DOUBLE PRECISION :: D_N
    DOUBLE PRECISION :: TH1
	DOUBLE PRECISION :: L1,L2,L3
	
!   POST-PROCESSING
    INTEGER :: TOGGLE
    INTEGER :: NITER_POS
!    DOUBLE PRECISION :: U_TA(0:NI,0:NJ,0:NK)
!    DOUBLE PRECISION :: V_TA(0:NI,0:NJ,0:NK)
!    DOUBLE PRECISION :: P_TA(0:NI,0:NJ,0:NK)
!    DOUBLE PRECISION :: T_TA(0:NI,0:NJ,0:NK)
    DOUBLE PRECISION :: NU_TA,CD_TA,CL_TA,CPB_TA,ST_TA,CL_AMP_TA
    DOUBLE PRECISION :: RMS_CD,RMS_CL
    DOUBLE PRECISION :: CL1,CL_AMP,CL_NM1
    DOUBLE PRECISION :: TIME1,TIME2
    
!---------------------------------------------------------------------------
    ENDMODULE ComDat_Shared