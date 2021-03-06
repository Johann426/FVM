	SUBROUTINE INIT

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	!	ALLOCATE THE VARIABLES IN THE PARALLEL PRECESS

	ALLOCATE(RADIUS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(U(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(V(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(W(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(P(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(UG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(VG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(WG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(PG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(VOR_X(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VOR_Y(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VOR_Z(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VOR_M(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(U_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(V_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(W_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(DUM_3D(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(U1(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(U2(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(U3(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(UE(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(UW(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VN(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(WT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(WB(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(NLX(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(NLY(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(NLZ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(NLX_NM1(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(NLY_NM1(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(NLZ_NM1(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(DIFX(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(DIFY(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(DIFZ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(S_U(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S_V(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S_W(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(SU(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(THETA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
    
	ALLOCATE(AREA1(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	ALLOCATE(AREA2(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	
	ALLOCATE(LAMDA2(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LAMDAG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	
	!	LES VARIABLES
	IF(USE_LES.EQ.1) THEN
	
	ALLOCATE(WID1(IS1:IE1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	ALLOCATE(WID2(IS1:IE1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	ALLOCATE(WID3(IS1:IE1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	
	ALLOCATE(S11(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S12(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S13(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S22(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S23(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S33(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(A11(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A12(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A13(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(A21(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A22(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A23(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(A31(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A32(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A33(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(S11A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S12A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S13A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S22A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S23A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S33A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(A11A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A12A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A13A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(A21A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A22A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A23A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(A31A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A32A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(A33A(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(PHIG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(PHIA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(UA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(WA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(A_IJ_A_IJ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(S_IJ_S_IJ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(CV(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(EDDY_VISCOS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(SGS_X(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(SGS_Y(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(SGS_Z(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ENDIF
	
	
	NITER = 0
	
	NITER_TIME = 0
	
	NUM_WRITE = 0
	
	RTIME = 0.0D0
	
	!   INITIALIZING
	U(IS:IE,JS:JE,KS:KE)   = 1.0D0
	V(IS:IE,JS:JE,KS:KE)   = 0.0D0
    W(IS:IE,JS:JE,KS:KE)   = 0.0D0
	
	P(IS:IE,JS:JE,KS:KE)   = 0.0D0
    
    NLX(IS:IE,JS:JE,KS:KE) = 0.0D0
	NLY(IS:IE,JS:JE,KS:KE) = 0.0D0
	NLZ(IS:IE,JS:JE,KS:KE) = 0.0D0
	
	!   RESTART ROUTINE
	IF(IRESTART.EQ.0.AND.RANK.EQ.0)	THEN
		PRINT *,'Program start...'
		
		OPEN(UNIT=11,FILE='CD_CL.PLT',STATUS='UNKNOWN')
	    WRITE(11,*) 'variables=TIME,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,cfl_max'
	    CLOSE(11)
		
		OPEN(UNIT=12,FILE='CV.PLT',STATUS='UNKNOWN')
	    WRITE(12,*) 'variables=TIME,Cv'
	    CLOSE(12)
	    
	ENDIF

	!	READ RESTART PARAMETER AND NONDIMENSIONAL NUMBER
	
	IF(IRESTART.EQ.1) THEN
		
		PRINT *,'Reading restart files...'
    	
    	FILENAME2 = 'RESTART_OO.DAT'
    	
    	WRITE(FILENAME2(9:10),'(I2.2)') INT(RANK)
		
		OPEN(UNIT=16+RANK,FILE=FILENAME2,STATUS='OLD')
		READ(16+RANK,122) RE,PE,DTIME,RTIME,NITER_TIME
122		FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,I16)
		DO K=KS,KE
		DO J=JS,JE
		DO I=IS,IE
			READ(16+RANK,128) NLX(I,J,K),NLY(I,J,K),NLZ(I,J,K)
			READ(16+RANK,129) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K)
128		FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9)
129		FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
		ENDDO
		ENDDO
		ENDDO
		CLOSE(16+RANK)
	ENDIF

	!   ANGULAR VECTOR (RADIAN)
	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
		RADIUS(I,J,K) = DSQRT(X(I,J,K)**2+Y(I,J,K)**2)
	ENDDO
    ENDDO
	ENDDO
	
	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
      THETA(I,J,K) = DATAN(Y(I,J,K)/X(I,J,K))
 	  IF     (THETA(I,J,K).GE.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).GT.0.0D0) THEN
		THETA(I,J,K) = THETA(I,J,K)
 	  
	  ELSEIF (THETA(I,J,K).LT.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).GT.0.0D0) THEN
		THETA(I,J,K) = PI+THETA(I,J,K)

 	  ELSEIF (THETA(I,J,K).GE.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).LE.0.0D0) THEN
		THETA(I,J,K) = PI+THETA(I,J,K)

	  ELSEIF (THETA(I,J,K).LT.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).LT.0.0D0) THEN
		THETA(I,J,K) = 2.0D0*PI+THETA(I,J,K)
	  ENDIF
	ENDDO
	ENDDO
    ENDDO

	!   SURFACE AREA ON THE CELL'S CENTER
    DUM = 0.D0	!SUM_AREA = 0.D0
    
    DO K=KS1,KE1
    DO J=JS1,JE1
    I=IS1
    L1 =DSQRT((XG(I,J+1,K)-XG(I,J,K))**2+(YG(I,J+1,K)-YG(I,J,K))**2+(ZG(I,J+1,K)-ZG(I,J,K))**2)
	L2 =DSQRT((XG(I,J,K+1)-XG(I,J,K))**2+(YG(I,J,K+1)-YG(I,J,K))**2+(ZG(I,J,K+1)-ZG(I,J,K))**2)
	L3 =DSQRT((XG(I,J+1,K)-XG(I,J,K+1))**2+(YG(I,J+1,K)-YG(I,J,K+1))**2+(ZG(I,J+1,K)-ZG(I,J,K+1))**2)
	TH1=DACOS((L3**2-L2**2-L1**2)/(-2.D0*L1*L2))
	AREA1(I,J,K) = 0.5D0*L1*L2*DSIN(TH1)
    
    L1 =DSQRT((XG(I,J+1,K+1)-XG(I,J+1,K))**2+(YG(I,J+1,K+1)-YG(I,J+1,K))**2+(ZG(I,J+1,K+1)-ZG(I,J+1,K))**2)
	L2 =DSQRT((XG(I,J+1,K+1)-XG(I,J,K+1))**2+(YG(I,J+1,K+1)-YG(I,J,K+1))**2+(ZG(I,J+1,K+1)-ZG(I,J,K+1))**2)
!	L3 =DSQRT((XG(I,J+1,K)-XG(I,J,K+1))**2+(YG(I,J+1,K)-YG(I,J,K+1))**2+(ZG(I,J+1,K)-ZG(I,J,K+1))**2)
	TH1=DACOS((L3**2-L2**2-L1**2)/(-2.D0*L1*L2))
	AREA2(I,J,K) = 0.5D0*L1*L2*DSIN(TH1)
	
		DUM = DUM + AREA1(I,J,K)+AREA2(I,J,K)	!SUM_AREA = SUM_AREA + AREA(I,J,K)
    ENDDO
    ENDDO

	CALL MPI_ALLREDUCE(DUM,SUM_AREA,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR) ! SUM_AREA = DUM
	
   				IF(RANK.EQ.0) THEN
    WRITE (*,2) NI,NJ,NK,NUM_GRID,SUM_AREA
    WRITE (*,3) RE,PR,PE,DTIME
2   FORMAT(TR1,'NI =',I4,TR3,'NJ =',I4,TR3,'NK =',I4,TR4,'Cells =',I7,TR3,'Surface Area =',F10.6)
3   FORMAT(TR1,'Re =',F10.2,TR4,'Pr =',F7.2,TR4,'Pe =',F7.2,TR4,'dtime =',F10.4,TR12)
				ENDIF	! PROCESS ZERO ONLY
	
	IF(USE_LES.EQ.1) THEN
	
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
	LX= MAX(XG(I,J,K),XG(I+1,J,K),XG(I,J+1,K),XG(I+1,J+1,K),XG(I,J,K+1),XG(I+1,J,K+1),XG(I,J+1,K+1),XG(I+1,J+1,K+1))
	LY= MAX(YG(I,J,K),YG(I+1,J,K),YG(I,J+1,K),YG(I+1,J+1,K),YG(I,J,K+1),YG(I+1,J,K+1),YG(I,J+1,K+1),YG(I+1,J+1,K+1))
	LZ= MAX(ZG(I,J,K),ZG(I+1,J,K),ZG(I,J+1,K),ZG(I+1,J+1,K),ZG(I,J,K+1),ZG(I+1,J,K+1),ZG(I,J+1,K+1),ZG(I+1,J+1,K+1))
	
	DUM = DJAC(I,J,K)**(1.0D0/3.0D0)
	
	WID1(I,J,K)=DJAC(I,J,K)**(1.0D0/3.0D0)	!DUM*LX/(LX*LY*LZ)	! ACTUALLY WID1=DX
	WID2(I,J,K)=DJAC(I,J,K)**(1.0D0/3.0D0)	!DUM*LY/(LX*LY*LZ)	! ACTUALLY WID2=DY
	WID3(I,J,K)=DJAC(I,J,K)**(1.0D0/3.0D0)	!DUM*LZ/(LX*LY*LZ)	! ACTUALLY WID3=DZ
	IF(WID1(I,J,K)*WID2(I,J,K)*WID3(I,J,K).GT.100) WRITE(*,*) 'OH!!!'
	ENDDO
	ENDDO
	ENDDO

	ENDIF
!--------------------------------------------------------------------------------------
	END SUBROUTINE INIT