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
	
	ALLOCATE(U_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(V_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(W_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	ALLOCATE(T_N(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(PLI(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
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
    
	ALLOCATE(DUDN(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	ALLOCATE(DTDN(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	ALLOCATE(AREA1(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	ALLOCATE(AREA2(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
	
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
	PLI(IS:IE,JS:JE,KS:KE) = 0.0D0
    
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
		RTIME_MAX = RTIME_MAX * 2.D0
		
		PRINT *,'Reading restart files...'
    	
    	FILENAME2 = 'RESTART_OO.DAT'
    	
    	WRITE(FILENAME2(9:10),'(I2.2)') INT4(RANK)
		
		OPEN(UNIT=16+RANK,FILE=FILENAME2,STATUS='OLD')
		READ(16+RANK,*) RE,PR,DTIME
		READ(16+RANK,*) RTIME,NITER_TIME
		K=KS
		DO J=JS,JE
		DO I=IS,IE
			READ(16+RANK,*) NLX(I,J,K),NLY(I,J,K),ADV(I,J,K)				! RESTART FROM 2-DIMENSIONAL DATA
			READ(16+RANK,*) U(I,J,K),V(I,J,K),P(I,J,K),T(I,J,K)				! RESTART FROM 2-DIMENSIONAL DATA
		ENDDO
		ENDDO
		DO K=KS,KE
		DO J=JS,JE
		DO I=IS,IE
!			READ(16+RANK,*) NLX(I,J,K),NLY(I,J,K),NLZ(I,J,K),ADV(I,J,K)
!			READ(16+RANK,*) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K)
			NLX(I,J,K)=NLX(I,J,KS)	! RESTART FROM 2-DIMENSIONAL DATA
			NLY(I,J,K)=NLY(I,J,KS)	! RESTART FROM 2-DIMENSIONAL DATA
			ADV(I,J,K)=ADV(I,J,KS)	! RESTART FROM 2-DIMENSIONAL DATA
			U(I,J,K)=U(I,J,KS)		! RESTART FROM 2-DIMENSIONAL DATA
			V(I,J,K)=V(I,J,KS)		! RESTART FROM 2-DIMENSIONAL DATA
			P(I,J,K)=P(I,J,KS)		! RESTART FROM 2-DIMENSIONAL DATA
			T(I,J,K)=T(I,J,KS)		! RESTART FROM 2-DIMENSIONAL DATA
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

!	CONTOUR = 'THET_OO.PLT'
!	
!	WRITE(CONTOUR(6:7),'(I2.2)') INT4(RANK)
!	
!	OPEN(UNIT=34+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
!	WRITE(34+RANK,*) 'title=""'
!	WRITE(34+RANK,*) 'variables=x,y,z,THETA,J,Q11,Q12'
!	WRITE(34+RANK,*) 'zone t="',RTIME,'",i=',IE-IS+1,',j=',JE-JS+1,',k=',KE-KS+1
!	DO K=KS,KE
!	DO J=JS,JE
!	DO I=IS,IE
!		WRITE(34+RANK,*) X(I,J,K),Y(I,J,K),Z(I,J,K),THETA(I,J,K),DJAC(I,J,K),Q11(I,J,K),Q12(I,J,K)
!	ENDDO
!	ENDDO
!	ENDDO
!	CLOSE(34+RANK)


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

	CALL MPI_ALLREDUCE(DUM,SUM_AREA,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR) ! SUM_AREA = DUM
	
	REF_A = 2.D0*PI	! REFERANCE AREA FOR CIRCULAR CYLINDER
	
   				IF(RANK.EQ.0) THEN
    WRITE (*,2) NI,NJ,NK,NUM_GRID,SUM_AREA
    WRITE (*,3) RE,PR,DTIME
2   FORMAT(TR1'NI =',I4,TR3'NJ =',I4,TR3'NK =',I4,TR4'Cells =',I7,TR3'Surface Area =',F10.6)
3   FORMAT(TR1'Re =',F10.2,TR4'Pr =',F7.2,TR4'dtime =',F10.4,TR12)
				ENDIF	! PROCESS ZERO ONLY
!--------------------------------------------------------------------------------------
	END SUBROUTINE INIT