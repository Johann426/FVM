	SUBROUTINE POST_PROCESS
	
	USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE

	IF(IPOST.EQ.0 .AND. POS_LES.EQ.0) THEN
		ALLOCATE(U_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(V_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(W_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(P_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(T_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		U_AVG = 0.D0
		V_AVG = 0.D0
		W_AVG = 0.D0
		P_AVG = 0.D0
		T_AVG = 0.D0
		ALLOCATE(UU_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(UV_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(UW_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(VV_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(VW_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(WW_AVG(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(U_RMS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(V_RMS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(W_RMS(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		UU_AVG = 0.D0
        UV_AVG = 0.D0
        UW_AVG = 0.D0
        VV_AVG = 0.D0
        VW_AVG = 0.D0
        WW_AVG = 0.D0
        U_RMS = 0.D0
		V_RMS = 0.D0
		W_RMS = 0.D0
        ALLOCATE(SGS11(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(SGS12(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(SGS13(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(SGS22(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(SGS23(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        ALLOCATE(SGS33(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
        SGS11 = 0.D0
        SGS12 = 0.D0
        SGS13 = 0.D0
        SGS22 = 0.D0
        SGS23 = 0.D0
        SGS33 = 0.D0
		ALLOCATE(SH_S (IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
		ALLOCATE(NU_S (IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
		ALLOCATE(CP_S (IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
		SH_S (:,:,:) = 0.D0
		NU_S (:,:,:) = 0.D0
		CP_S (:,:,:) = 0.D0
		RTIME_MAX  = RTIME_END
		NITER_POS = 0
		IPOST = 1
		SIGN = .TRUE.
		DTIME =	DTIME_POST
		
	ELSE IF(IPOST.EQ.1 .AND. POS_LES.EQ.0) THEN
		U_AVG(:,:,:) = U_AVG(:,:,:) / DBLE(NITER_POS)
		V_AVG(:,:,:) = V_AVG(:,:,:) / DBLE(NITER_POS)
		W_AVG(:,:,:) = W_AVG(:,:,:) / DBLE(NITER_POS)
		P_AVG(:,:,:) = P_AVG(:,:,:) / DBLE(NITER_POS)
		T_AVG(:,:,:) = T_AVG(:,:,:) / DBLE(NITER_POS)
		
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS1
	! (A,B,C) : VECTOR, NORMAL VECTOR OF THE DISCRETIZED CYLINDER SURFACE
	DUM_A = YG(I,J,K)*(ZG(I,J+1,K)-ZG(I,J,K+1))+YG(I,J+1,K)*(ZG(I,J,K+1)-ZG(I,J,K))+YG(I,J,K+1)*(ZG(I,J,K)-ZG(I,J+1,K))
	DUM_B = ZG(I,J,K)*(XG(I,J+1,K)-XG(I,J,K+1))+ZG(I,J+1,K)*(XG(I,J,K+1)-XG(I,J,K))+ZG(I,J,K+1)*(XG(I,J,K)-XG(I,J+1,K))
	DUM_C = XG(I,J,K)*(YG(I,J+1,K)-YG(I,J,K+1))+XG(I,J+1,K)*(YG(I,J,K+1)-YG(I,J,K))+XG(I,J,K+1)*(YG(I,J,K)-YG(I,J+1,K))
	! -D : SCALAR, LENGTH OF THE NORMAL VECTOR
	DUM_D = XG(I,J,K)*(YG(I,J+1,K)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K))+XG(I,J+1,K)*(YG(I,J,K+1)*ZG(I,J,K)-YG(I,J,K)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J,K)*ZG(I,J+1,K)-YG(I,J+1,K)*ZG(I,J,K))
	DUM_D = -1.D0*DUM_D
	
	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
	
	D_N = DABS(DUM_A*X(I,J,K)+DUM_B*Y(I,J,K)+DUM_C*Z(I,J,K)+DUM_D) / DUM

	DUDN = DSQRT((V_AVG(I,J,K)*DUM_C/DUM-W_AVG(I,J,K)*DUM_B/DUM)**2	&
				+(U_AVG(I,J,K)*DUM_C/DUM-W_AVG(I,J,K)*DUM_A/DUM)**2	&
				+(U_AVG(I,J,K)*DUM_B/DUM-V_AVG(I,J,K)*DUM_A/DUM)**2)
	
	DUDN = DUDN/D_N
	
	DTDN = (T_AVG(I,J,K)-T_WALL) / D_N
	
	SH_S(I,J,K) = DUDN*AREA1(I,J,K)/(AREA1(I,J,K)+AREA2(I,J,K))
	NU_S(I,J,K)	= DTDN*AREA1(I,J,K)/(AREA1(I,J,K)+AREA2(I,J,K))
	CP_S(I,J,K)	= 0.125D0*(P_AVG(I,J,K)+P_AVG(I,J-1,K)+P_AVG(I-1,J,K)+P_AVG(I-1,J-1,K)			&
						+P_AVG(I,J,K+1)+P_AVG(I,J-1,K+1)+P_AVG(I-1,J,K+1)+P_AVG(I-1,J-1,K+1))
	
	DUM_A = YG(I,J+1,K+1)*(ZG(I,J,K+1)-ZG(I,J,K))+YG(I,J,K+1)*(ZG(I,J,K)-ZG(I,J+1,K+1))+YG(I,J,K)*(ZG(I,J+1,K+1)-ZG(I,J,K+1))
	DUM_B = ZG(I,J+1,K+1)*(XG(I,J,K+1)-XG(I,J,K))+ZG(I,J,K+1)*(XG(I,J,K)-XG(I,J+1,K+1))+ZG(I,J,K)*(XG(I,J+1,K+1)-XG(I,J,K+1))
	DUM_C = XG(I,J+1,K+1)*(YG(I,J,K+1)-YG(I,J,K))+XG(I,J,K+1)*(YG(I,J,K)-YG(I,J+1,K+1))+XG(I,J,K)*(YG(I,J+1,K+1)-YG(I,J,K+1))
	! D : LENGTH OF THE NORMAL VECTOR
	DUM_D = XG(I,J+1,K+1)*(YG(I,J,K+1)*ZG(I,J,K)-YG(I,J,K)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J,K)*ZG(I,J+1,K+1)-YG(I,J+1,K+1)*ZG(I,J,K))+XG(I,J,K)*(YG(I,J+1,K+1)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K+1))
	DUM_D = -1.D0*DUM_D
	
	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
	
	D_N = DABS(DUM_A*X(I,J,K)+DUM_B*Y(I,J,K)+DUM_C*Z(I,J,K)+DUM_D) / DUM
	
	DUDN = DSQRT((V_AVG(I,J,K)*DUM_C/DUM-W_AVG(I,J,K)*DUM_B/DUM)**2	&
				+(U_AVG(I,J,K)*DUM_C/DUM-W_AVG(I,J,K)*DUM_A/DUM)**2	&
				+(U_AVG(I,J,K)*DUM_B/DUM-V_AVG(I,J,K)*DUM_A/DUM)**2)
	
	DUDN = DUDN/D_N
	
	DTDN = (T_AVG(I,J,K)-T_WALL) / D_N
	
	SH_S(I,J,K) = SH_S(I,J,K) + DUDN*AREA2(I,J,K)/(AREA1(I,J,K)+AREA2(I,J,K))
	NU_S(I,J,K)	= NU_S(I,J,K) + DTDN*AREA2(I,J,K)/(AREA1(I,J,K)+AREA2(I,J,K))
	!CP_S(I,J,K)	= CP_S(I,J,K)
	ENDDO
	ENDDO
		
		DO K=KS1,KE
		DO J=JS1,JE
		DO I=IS1,IE
		PG(I,J,K) = 0.125D0*(P_AVG(I,J,K)+P_AVG(I,J-1,K)+P_AVG(I-1,J,K)+P_AVG(I-1,J-1,K)		&
						+P_AVG(I,J,K-1)+P_AVG(I,J-1,K-1)+P_AVG(I-1,J,K-1)+P_AVG(I-1,J-1,K-1))
		UG(I,J,K) = 0.125D0*(U_AVG(I,J,K)+U_AVG(I,J-1,K)+U_AVG(I-1,J,K)+U_AVG(I-1,J-1,K)		&
						+U_AVG(I,J,K-1)+U_AVG(I,J-1,K-1)+U_AVG(I-1,J,K-1)+U_AVG(I-1,J-1,K-1))
		VG(I,J,K) = 0.125D0*(V_AVG(I,J,K)+V_AVG(I,J-1,K)+V_AVG(I-1,J,K)+V_AVG(I-1,J-1,K)		&
						+V_AVG(I,J,K-1)+V_AVG(I,J-1,K-1)+V_AVG(I-1,J,K-1)+V_AVG(I-1,J-1,K-1))
		WG(I,J,K) = 0.125D0*(W_AVG(I,J,K)+W_AVG(I,J-1,K)+W_AVG(I-1,J,K)+W_AVG(I-1,J-1,K)		&
						+W_AVG(I,J,K-1)+W_AVG(I,J-1,K-1)+W_AVG(I-1,J,K-1)+W_AVG(I-1,J-1,K-1))
		TG(I,J,K) = 0.125D0*(T_AVG(I,J,K)+T_AVG(I,J-1,K)+T_AVG(I-1,J,K)+T_AVG(I-1,J-1,K)		&
						+T_AVG(I,J,K-1)+T_AVG(I,J-1,K-1)+T_AVG(I-1,J,K-1)+T_AVG(I-1,J-1,K-1))
		ENDDO
		ENDDO
		ENDDO
	
		CONTOUR = 'Avrg000.PLT'
		
		WRITE(CONTOUR(5:7),'(I3.3)') INT(RANK)
		
		OPEN(UNIT=1000+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
		WRITE(1000+RANK,*) 'TITLE="','RE',RE,'PR',PR,'"'
		WRITE(1000+RANK,*) 'variables=x,y,z,p,u,v,w,T'
	    WRITE(1000+RANK,68) IE-IS,JE-JS,KE-KS,RTIME
68		FORMAT(TR1,'ZONE T="WHOLE", I=',I6,', J=',I6,', K=',I6,', SOLUTIONTIME=',F10.4)
		DO K=KS1,KE
		DO J=JS1,JE
		DO I=IS1,IE
			WRITE(1000+RANK,73) XG(I,J,K),YG(I,J,K),ZG(I,J,K),PG(I,J,K),UG(I,J,K),VG(I,J,K),WG(I,J,K),TG(I,J,K)
73			FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
		ENDDO
		ENDDO
		ENDDO
		CLOSE(1000+RANK)
		
        !-----------------------------------------------------------------------------
   		UU_AVG(:,:,:) = UU_AVG(:,:,:) / DBLE(NITER_POS)
		UV_AVG(:,:,:) = UV_AVG(:,:,:) / DBLE(NITER_POS)
        UW_AVG(:,:,:) = UW_AVG(:,:,:) / DBLE(NITER_POS)
        VV_AVG(:,:,:) = VV_AVG(:,:,:) / DBLE(NITER_POS)
        VW_AVG(:,:,:) = VW_AVG(:,:,:) / DBLE(NITER_POS)
        WW_AVG(:,:,:) = WW_AVG(:,:,:) / DBLE(NITER_POS)
        
        SGS11(:,:,:) = SGS11(:,:,:) / DBLE(NITER_POS)
        SGS12(:,:,:) = SGS12(:,:,:) / DBLE(NITER_POS)
        SGS13(:,:,:) = SGS13(:,:,:) / DBLE(NITER_POS)
        SGS22(:,:,:) = SGS22(:,:,:) / DBLE(NITER_POS)
        SGS23(:,:,:) = SGS23(:,:,:) / DBLE(NITER_POS)
        SGS33(:,:,:) = SGS33(:,:,:) / DBLE(NITER_POS)
        
		U_RMS(:,:,:) = DSQRT(UU_AVG(:,:,:))
		V_RMS(:,:,:) = DSQRT(VV_AVG(:,:,:))
		W_RMS(:,:,:) = DSQRT(WW_AVG(:,:,:))
		
		DO K=KS1,KE
		DO J=JS1,JE
		DO I=IS1,IE
		UG(I,J,K) = 0.125D0*(U_RMS(I,J,K)+U_RMS(I,J-1,K)+U_RMS(I-1,J,K)+U_RMS(I-1,J-1,K)		&
						+U_RMS(I,J,K-1)+U_RMS(I,J-1,K-1)+U_RMS(I-1,J,K-1)+U_RMS(I-1,J-1,K-1))
		VG(I,J,K) = 0.125D0*(V_RMS(I,J,K)+V_RMS(I,J-1,K)+V_RMS(I-1,J,K)+V_RMS(I-1,J-1,K)		&
						+V_RMS(I,J,K-1)+V_RMS(I,J-1,K-1)+V_RMS(I-1,J,K-1)+V_RMS(I-1,J-1,K-1))
		WG(I,J,K) = 0.125D0*(W_RMS(I,J,K)+W_RMS(I,J-1,K)+W_RMS(I-1,J,K)+W_RMS(I-1,J-1,K)		&
						+W_RMS(I,J,K-1)+W_RMS(I,J-1,K-1)+W_RMS(I-1,J,K-1)+W_RMS(I-1,J-1,K-1))
		ENDDO
		ENDDO
		ENDDO
		
		CONTOUR = 'Urms000.PLT'
		
		WRITE(CONTOUR(5:7),'(I3.3)') INT(RANK)
		
		OPEN(UNIT=1000+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
		WRITE(1000+RANK,*) 'TITLE="','RE',RE,'PR',PR,'"'
		WRITE(1000+RANK,*) 'variables=x,y,z,u_rms,v_rms,w_rms'
	    WRITE(1000+RANK,104) IE-IS,JE-JS,KE-KS,RTIME
104		FORMAT(TR1,'ZONE T="WHOLE", I=',I6,', J=',I6,', K=',I6,', SOLUTIONTIME=',F10.4)
		DO K=KS1,KE
		DO J=JS1,JE
		DO I=IS1,IE
			WRITE(1000+RANK,109) XG(I,J,K),YG(I,J,K),ZG(I,J,K),UG(I,J,K),VG(I,J,K),WG(I,J,K)
109			FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
		ENDDO
		ENDDO
		ENDDO
		CLOSE(1000+RANK)

		UU_AVG(:,:,:) = UU_AVG(:,:,:) - U_AVG(:,:,:)*U_AVG(:,:,:) + SGS11(:,:,:)
        UV_AVG(:,:,:) = UV_AVG(:,:,:) - U_AVG(:,:,:)*V_AVG(:,:,:) + SGS12(:,:,:)
        UW_AVG(:,:,:) = UW_AVG(:,:,:) - U_AVG(:,:,:)*W_AVG(:,:,:) + SGS13(:,:,:)
        VV_AVG(:,:,:) = VV_AVG(:,:,:) - V_AVG(:,:,:)*V_AVG(:,:,:) + SGS22(:,:,:)
        VW_AVG(:,:,:) = VW_AVG(:,:,:) - V_AVG(:,:,:)*W_AVG(:,:,:) + SGS23(:,:,:)
        WW_AVG(:,:,:) = WW_AVG(:,:,:) - W_AVG(:,:,:)*W_AVG(:,:,:) + SGS33(:,:,:)
		
		DO K=KS1,KE
		DO J=JS1,JE
		DO I=IS1,IE
		UG(I,J,K) = 0.125D0*(UU_AVG(I,J,K)+UU_AVG(I,J-1,K)+UU_AVG(I-1,J,K)+UU_AVG(I-1,J-1,K)		&
						+UU_AVG(I,J,K-1)+UU_AVG(I,J-1,K-1)+UU_AVG(I-1,J,K-1)+UU_AVG(I-1,J-1,K-1))
		VG(I,J,K) = 0.125D0*(VV_AVG(I,J,K)+VV_AVG(I,J-1,K)+VV_AVG(I-1,J,K)+VV_AVG(I-1,J-1,K)		&
						+VV_AVG(I,J,K-1)+VV_AVG(I,J-1,K-1)+VV_AVG(I-1,J,K-1)+VV_AVG(I-1,J-1,K-1))
		WG(I,J,K) = 0.125D0*(WW_AVG(I,J,K)+WW_AVG(I,J-1,K)+WW_AVG(I-1,J,K)+WW_AVG(I-1,J-1,K)		&
						+WW_AVG(I,J,K-1)+WW_AVG(I,J-1,K-1)+WW_AVG(I-1,J,K-1)+WW_AVG(I-1,J-1,K-1))
        U_RMS(I,J,K) = 0.125D0*(UV_AVG(I,J,K)+UV_AVG(I,J-1,K)+UV_AVG(I-1,J,K)+UV_AVG(I-1,J-1,K)		&
						+UV_AVG(I,J,K-1)+UV_AVG(I,J-1,K-1)+UV_AVG(I-1,J,K-1)+UV_AVG(I-1,J-1,K-1))
		V_RMS(I,J,K) = 0.125D0*(UW_AVG(I,J,K)+UW_AVG(I,J-1,K)+UW_AVG(I-1,J,K)+UW_AVG(I-1,J-1,K)		&
						+UW_AVG(I,J,K-1)+UW_AVG(I,J-1,K-1)+UW_AVG(I-1,J,K-1)+UW_AVG(I-1,J-1,K-1))
		W_RMS(I,J,K) = 0.125D0*(VW_AVG(I,J,K)+VW_AVG(I,J-1,K)+VW_AVG(I-1,J,K)+VW_AVG(I-1,J-1,K)		&
						+VW_AVG(I,J,K-1)+VW_AVG(I,J-1,K-1)+VW_AVG(I-1,J,K-1)+VW_AVG(I-1,J-1,K-1))
		ENDDO
		ENDDO
		ENDDO
		
		CONTOUR = 'Strs000.PLT'
		
		WRITE(CONTOUR(5:7),'(I3.3)') INT(RANK)
		
		OPEN(UNIT=1000+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
		WRITE(1000+RANK,*) 'TITLE="','RE',RE,'PR',PR,'"'
		WRITE(1000+RANK,*) 'variables=x,y,z,uu,uv,uw,vv,vw,ww'
	    WRITE(1000+RANK,244) IE-IS,JE-JS,KE-KS,RTIME
244		FORMAT(TR1,'ZONE T="WHOLE", I=',I6,', J=',I6,', K=',I6,', SOLUTIONTIME=',F10.4)
		DO K=KS1,KE
		DO J=JS1,JE
		DO I=IS1,IE
			WRITE(1000+RANK,249) XG(I,J,K),YG(I,J,K),ZG(I,J,K),UG(I,J,K),U_RMS(I,J,K),V_RMS(I,J,K),VG(I,J,K),W_RMS(I,J,K),WG(I,J,K)
249			FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
		ENDDO
		ENDDO
		ENDDO
		CLOSE(1000+RANK)
        
        
        
		THETA(I,J,K)= 180.D0/PI		*THETA(I,J,K)
		SH_S(:,:,:)	= 2.D0/DSQRT(RE)*SH_S(:,:,:)	! SURFACE FRICTION COEF.
		CP_S(:,:,:) = 2.D0			*CP_S(:,:,:)	! SURFACE PRESSURE COEF.
		NU_S(:,:,:) =				-NU_S(:,:,:)	! SURFACE NUSSELT NUMBER
		
		CONTOUR = 'Surf000.PLT'
		
		WRITE(CONTOUR(5:7),'(I3.3)') INT(RANK)
		
		OPEN(UNIT=100+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
		WRITE(100+RANK,*) 'TITLE="SURFACE COEF."'
		WRITE(100+RANK,*) 'variables="x","y","z","Theta","Cf","Cp","Nu"'
		WRITE(100+RANK,132) 1,JE1-JS1+1,KE1-KS1+1
132		FORMAT(TR1,'ZONE T="POST", I=',I6,', j=',I6,', k=',I6)
		DO K=KS1,KE1
		DO J=JS1,JE1
		I=IS1
			WRITE(100+RANK,137) X(I,J,K),Y(I,J,K),Z(I,J,K),THETA(I,J,K),SH_S(I,J,K),CP_S(I,J,K),NU_S(I,J,K)
137			FORMAT(TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9,TR1,D16.9)
		ENDDO
		ENDDO
		CLOSE(100+RANK)
		
		SIGN = .FALSE.
	
	ELSE IF(POS_LES.EQ.1) THEN
		
		SIGN = .FALSE.

	ENDIF
!---------------------------------------------------------------------------------------------
	END SUBROUTINE POST_PROCESS