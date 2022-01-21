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
	ALLOCATE(T(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(UG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(VG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(WG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(PG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	ALLOCATE(TG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(VOR_X(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VOR_Y(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VOR_Z(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(VOR_M(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(U_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(V_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(W_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(T_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
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
	
	ALLOCATE(ADV(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(ADV_NM1(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(DIFX(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(DIFY(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(DIFZ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(DIFT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
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
	
	ALLOCATE(FILT_WID(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(LSS11(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS12(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS13(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)        

	ALLOCATE(LSS21(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS22(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS23(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)    
        
	ALLOCATE(LSS31(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS32(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS33(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)            

	ALLOCATE(MAGNITUDE_LSS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(U_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(V_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(W_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(LSS11_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS12_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS13_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)        

	ALLOCATE(LSS21_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS22_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS23_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)    
        
	ALLOCATE(LSS31_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS32_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LSS33_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)            

	ALLOCATE(MAGNITUDE_LSS_HAT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)

	ALLOCATE(DUM_IN(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(DUM_OUT(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(L_IJ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(M_IJ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(LM_IJ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(MM_IJ(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
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
	
	T(IS:IE,JS:JE,KS:KE)   = 0.0D0
    
	P(IS:IE,JS:JE,KS:KE)   = 0.0D0
    
    NLX(IS:IE,JS:JE,KS:KE) = 0.0D0
	NLY(IS:IE,JS:JE,KS:KE) = 0.0D0
	NLZ(IS:IE,JS:JE,KS:KE) = 0.0D0
	
	ADV(IS:IE,JS:JE,KS:KE) = 0.0D0

	!   RESTART ROUTINE
	IF(IRESTART.EQ.0.AND.RANK.EQ.0)	THEN
		PRINT *,'Program start...'
		
		OPEN(UNIT=11,FILE='NU_CD_CL.PLT',STATUS='UNKNOWN')
	    WRITE(11,*) 'variables=TIME,NU,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,cfl_max'
	    CLOSE(11)
	ENDIF

	!	READ RESTART PARAMETER AND NONDIMENSIONAL NUMBER
	
	IF(IRESTART.EQ.1) THEN
		
		PRINT *,'Reading restart files...'
    	
    	FILENAME2 = 'RESTART_OO.DAT'
    	
    	WRITE(FILENAME2(9:10),'(I2.2)') INT(RANK)
		
		OPEN(UNIT=16+RANK,FILE=FILENAME2,STATUS='OLD')
		READ(16+RANK,122) RE,PE,DTIME,RTIME,NITER_TIME
122		FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
		DO K=KS,KE
		DO J=JS,JE
		DO I=IS,IE
			READ(16+RANK,128) NLX(I,J,K),NLY(I,J,K),NLZ(I,J,K),ADV(I,J,K)
			READ(16+RANK,129) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K)
128		FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
129		FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
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
!--------------------------------------------------------------------------------------
	END SUBROUTINE INIT