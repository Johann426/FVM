    SUBROUTINE CALCCOEF

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	INCLUDE 'mpif.h'
		
	!   CFL(Courant-Friedrichs-Lewy) CALCULATION
	CFL_MAX = 0.0D0
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
		CFL_LOCAL = DTIME *(ABS(XIX(I,J,K)*U(I,J,K) +XIY(I,J,K)*V(I,J,K) +XIZ(I,J,K)*W(I,J,K))  &
     		               +ABS(ETX(I,J,K)*U(I,J,K) +ETY(I,J,K)*V(I,J,K) +ETZ(I,J,K)*W(I,J,K))	&
						   +ABS(ZTX(I,J,K)*U(I,J,K) +ZTY(I,J,K)*V(I,J,K) +ZTZ(I,J,K)*W(I,J,K))  )
		CFL_MAX = DMAX1(CFL_MAX,CFL_LOCAL)
	ENDDO
	ENDDO
	ENDDO
	
	CALL MPI_ALLREDUCE(CFL_MAX,DUM,1,MPI_DOUBLE_PRECISION,MPI_MAX,NMCW,ERR)
	
	CFL_MAX = DUM
	
!	!	USER-DEFINED CFL(Courant-Friedrichs-Lewy) CONDITION
!	IF(CFL_MAX.GT.CFL_USE) THEN
!        DTIME = DTIME + 0.05D0*DTIME    &
!                *DBLE(IPOST-1)  ! CANCELED BY POST-PROCESSING PARAMETER
!    ELSE
!        DTIME = DTIME + 0.05D0*DTIME    &
!                *DBLE(1-IPOST)  ! CANCELED BY POST-PROCESSING PAREMETER
!	ENDIF
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   CPB CALCULATION
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS1
		PG(I,J,K) = 0.125D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K)+P(I,J,K+1)+P(I,J-1,K+1)+P(I-1,J,K+1)+P(I-1,J-1,K+1))
		UG(I,J,K) = 0.125D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K)+U(I,J,K+1)+U(I,J-1,K+1)+U(I-1,J,K+1)+U(I-1,J-1,K+1))
		VG(I,J,K) = 0.125D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K)+V(I,J,K+1)+V(I,J-1,K+1)+V(I-1,J,K+1)+V(I-1,J-1,K+1))
	I=IE
		PG(I,J,K) = 0.125D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K)+P(I,J,K+1)+P(I,J-1,K+1)+P(I-1,J,K+1)+P(I-1,J-1,K+1))
		UG(I,J,K) = 0.125D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K)+U(I,J,K+1)+U(I,J-1,K+1)+U(I-1,J,K+1)+U(I-1,J-1,K+1))
		VG(I,J,K) = 0.125D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K)+V(I,J,K+1)+V(I,J-1,K+1)+V(I-1,J,K+1)+V(I-1,J-1,K+1))
	ENDDO
    ENDDO
    
    CPB = 0.D0
    DUM_A = 0.D0
    DUM_B = 0.D0
    
    IF(JS1.EQ.1) THEN
		DO K=KS1,KE1
		I=1
		!DUM_B = DUM_B + PG(I,JS1,(KS+KE)/2)			! BASE LINE SUCTION PRESSURE
		DUM_B = DUM_B + PG(I,JS1,K)
		ENDDO
	ENDIF

	DUM_B = DUM_B/DBLE(KE1-KS1+1)

	IF(JS1.EQ.1) THEN
		CPB = CPB + 2.0D0 * ( DUM_A - DUM_B )
	ELSE
		CPB = 0.D0
	ENDIF
	
	CALL MPI_ALLREDUCE(CPB,DUM,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
	
	CPB = DUM / DBLE(NDIMK)

!---------------------------------------------------------------------------
!	SURFACE COEF.
	P_DRAG = 0.D0
	P_LIFT_Y = 0.D0
	P_LIFT_Z = 0.D0
    
    T_DRAG = 0.D0
	T_LIFT_Y = 0.D0
	T_LIFT_Z = 0.D0
	
	SUM_NU = 0.D0
	
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

	DUDN = DSQRT((V(I,J,K)*DUM_C/DUM-W(I,J,K)*DUM_B/DUM)**2	&
				+(U(I,J,K)*DUM_C/DUM-W(I,J,K)*DUM_A/DUM)**2	&
				+(U(I,J,K)*DUM_B/DUM-V(I,J,K)*DUM_A/DUM)**2)
	
	DUDN = DUDN/D_N
		
	P_DRAG = P_DRAG - 2.0D0*P(I,J,K)*AREA1(I,J,K)*DUM_A/DUM
	
	P_LIFT_Y = P_LIFT_Y - 2.0D0*P(I,J,K)*AREA1(I,J,K)*DUM_B/DUM
	
	T_DRAG = T_DRAG + (2.0D0/RE)*DUDN*AREA1(I,J,K) *U(I,J,K)/DSQRT(U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)
					
	T_LIFT_Y = T_LIFT_Y + (2.0D0/RE)*DUDN*AREA1(I,J,K) *V(I,J,K)/DSQRT(U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)
							
	DTDN = (T(I,J,K)-T_WALL) / D_N
	
	SUM_NU  = SUM_NU + DTDN*AREA1(I,J,K)
	
	YPLUS(I,J,K)= DSQRT(RE*DABS(DUDN))*D_N*AREA1(I,J,K)/(AREA1(I,J,K)+AREA2(I,J,K))
	
    ENDDO
	ENDDO
    
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS1
	! (A,B,C) : NORMAL VECTOR OF THE DISCRETIZED CYLINDER SURFACE
	DUM_A = YG(I,J+1,K+1)*(ZG(I,J,K+1)-ZG(I,J,K))+YG(I,J,K+1)*(ZG(I,J,K)-ZG(I,J+1,K+1))+YG(I,J,K)*(ZG(I,J+1,K+1)-ZG(I,J,K+1))
	DUM_B = ZG(I,J+1,K+1)*(XG(I,J,K+1)-XG(I,J,K))+ZG(I,J,K+1)*(XG(I,J,K)-XG(I,J+1,K+1))+ZG(I,J,K)*(XG(I,J+1,K+1)-XG(I,J,K+1))
	DUM_C = XG(I,J+1,K+1)*(YG(I,J,K+1)-YG(I,J,K))+XG(I,J,K+1)*(YG(I,J,K)-YG(I,J+1,K+1))+XG(I,J,K)*(YG(I,J+1,K+1)-YG(I,J,K+1))
	! D : LENGTH OF THE NORMAL VECTOR
	DUM_D = XG(I,J+1,K+1)*(YG(I,J,K+1)*ZG(I,J,K)-YG(I,J,K)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J,K)*ZG(I,J+1,K+1)-YG(I,J+1,K+1)*ZG(I,J,K))+XG(I,J,K)*(YG(I,J+1,K+1)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K+1))
	DUM_D = -1.D0*DUM_D
	
	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
	
	D_N = DABS(DUM_A*X(I,J,K)+DUM_B*Y(I,J,K)+DUM_C*Z(I,J,K)+DUM_D) / DUM
	
	DUDN = DSQRT((V(I,J,K)*DUM_C/DUM-W(I,J,K)*DUM_B/DUM)**2	&
				+(U(I,J,K)*DUM_C/DUM-W(I,J,K)*DUM_A/DUM)**2	&
				+(U(I,J,K)*DUM_B/DUM-V(I,J,K)*DUM_A/DUM)**2)

	DUDN = DUDN/D_N
	
	P_DRAG = P_DRAG - 2.0D0*P(I,J,K)*AREA2(I,J,K)*DUM_A/DUM
	
	P_LIFT_Y = P_LIFT_Y - 2.0D0*P(I,J,K)*AREA2(I,J,K)*DUM_B/DUM
	
	T_DRAG = T_DRAG + (2.0D0/RE)*DUDN*AREA2(I,J,K) *U(I,J,K)/DSQRT(U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)
					
	T_LIFT_Y = T_LIFT_Y + (2.0D0/RE)*DUDN*AREA2(I,J,K) *V(I,J,K)/DSQRT(U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)
	
	DTDN = (T(I,J,K)-T_WALL) / D_N
	
	SUM_NU  = SUM_NU + DTDN*AREA2(I,J,K)
	
	YPLUS(I,J,K)= YPLUS(I,J,K)+ DSQRT(RE*DABS(DUDN))*D_N*AREA2(I,J,K)/(AREA1(I,J,K)+AREA2(I,J,K))
	
    ENDDO
	ENDDO
	
    CALL MPI_ALLREDUCE(P_DRAG,CD_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
	CALL MPI_ALLREDUCE(P_LIFT_Y,CL_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
	
	CALL MPI_ALLREDUCE(T_DRAG,CD_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
	CALL MPI_ALLREDUCE(T_LIFT_Y,CL_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
	
	CALL MPI_ALLREDUCE(SUM_NU,NU,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)


	IF(RANK.EQ.0) THEN

	CD_P = CD_P / REF_A	!P_DRAG
	CL_P = CL_P / REF_A	!P_LIFT

	CD_V = CD_V / REF_A	!T_DRAG
	CL_V = CL_V / REF_A	!T_LIFT

	CD  = CD_P + CD_V
	CL  = CL_P + CL_V
    
    NU  = -NU / SUM_AREA	!-SUM_NU / SUM_AREA
	
	OPEN(UNIT=11,FILE='NU_CD_CL.PLT',STATUS='OLD',ACCESS='APPEND')
	WRITE (11,*) RTIME,NU,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,CFL_MAX ! PLOTTING
	CLOSE(11)

!   PRINTING
    PRINT *, ''
    WRITE (*,101) RTIME,NITER_TIME,CFL_MAX
	WRITE (*,102) NU
	WRITE (*,103) CD,CD_P,CD_V
	WRITE (*,104) CL,CL_P,CL_V,CPB
101 FORMAT(TR6,'Time',F10.3,TR12,'ITER',I10,TR12,'CFL',F10.3)
102 FORMAT(TR1,'NU =',F10.6)
103 FORMAT(TR1,'CD =',F10.6,TR4,'CD_P =',F10.6,TR4,'CD_V =',F10.6)
104 FORMAT(TR1,'CL =',F10.6,TR4,'CL_P =',F10.6,TR4,'CL_V =',F10.6,TR4,'CPB =',F10.6)

	ENDIF	! PROCESS ZERO ONLY

!!----------------------------------------------------------------------------------
!!	MASS BALANCE CHECK
!	MASS_INLET	= 0.D0
!	MASS_OUTLET = 0.D0
!	
!	DO K=KS1,KE1
!	DO J=JS1,JE1
!	I=IE
!	! (A,B,C) : VECTOR, NORMAL VECTOR OF THE DISCRETIZED CYLINDER SURFACE
!	DUM_A = YG(I,J,K)*(ZG(I,J+1,K)-ZG(I,J,K+1))+YG(I,J+1,K)*(ZG(I,J,K+1)-ZG(I,J,K))+YG(I,J,K+1)*(ZG(I,J,K)-ZG(I,J+1,K))
!	DUM_B = ZG(I,J,K)*(XG(I,J+1,K)-XG(I,J,K+1))+ZG(I,J+1,K)*(XG(I,J,K+1)-XG(I,J,K))+ZG(I,J,K+1)*(XG(I,J,K)-XG(I,J+1,K))
!	DUM_C = XG(I,J,K)*(YG(I,J+1,K)-YG(I,J,K+1))+XG(I,J+1,K)*(YG(I,J,K+1)-YG(I,J,K))+XG(I,J,K+1)*(YG(I,J,K)-YG(I,J+1,K))
!	! -D : SCALAR, LENGTH OF THE NORMAL VECTOR
!	DUM_D = XG(I,J,K)*(YG(I,J+1,K)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K))+XG(I,J+1,K)*(YG(I,J,K+1)*ZG(I,J,K)-YG(I,J,K)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J,K)*ZG(I,J+1,K)-YG(I,J+1,K)*ZG(I,J,K))
!	DUM_D = -1.D0*DUM_D
!	
!	IF(DUM_D.LT.0.D0) THEN
!		DUM_A = -1.D0*DUM_A
!		DUM_B = -1.D0*DUM_B
!		DUM_C = -1.D0*DUM_C
!		DUM_D = -1.D0*DUM_D
!	ENDIF
!	
!	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
!	
!    L1 =DSQRT((XG(I,J+1,K)-XG(I,J,K))**2+(YG(I,J+1,K)-YG(I,J,K))**2+(ZG(I,J+1,K)-ZG(I,J,K))**2)
!	L2 =DSQRT((XG(I,J,K+1)-XG(I,J,K))**2+(YG(I,J,K+1)-YG(I,J,K))**2+(ZG(I,J,K+1)-ZG(I,J,K))**2)
!	L3 =DSQRT((XG(I,J+1,K)-XG(I,J,K+1))**2+(YG(I,J+1,K)-YG(I,J,K+1))**2+(ZG(I,J+1,K)-ZG(I,J,K+1))**2)
!	TH1=DACOS((L3**2-L2**2-L1**2)/(-2.D0*L1*L2))
!	DAREA = 0.5D0*L1*L2*DSIN(TH1)
!	
!	IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
!		MASS_INLET	= MASS_INLET +(U(I,J,K)*DUM_A+V(I,J,K)*DUM_B+W(I,J,K)*DUM_C)/DUM*DAREA	!(u,0,0,) .dot. unit(a,b,c)
!	ELSE
!		MASS_OUTLET	= MASS_OUTLET +(U(I,J,K)*DUM_A+V(I,J,K)*DUM_B+W(I,J,K)*DUM_C)/DUM*DAREA	!(u,0,0,) .dot. unit(a,b,c)
!	ENDIF
!	
!	ENDDO
!	ENDDO
!	
!	DO K=KS1,KE1
!	DO J=JS1,JE1
!	I=IE
!	! (A,B,C) : NORMAL VECTOR OF THE DISCRETIZED CYLINDER SURFACE
!	DUM_A = YG(I,J+1,K+1)*(ZG(I,J+1,K)-ZG(I,J,K+1))+YG(I,J+1,K)*(ZG(I,J,K+1)-ZG(I,J+1,K+1))+YG(I,J,K+1)*(ZG(I,J+1,K+1)-ZG(I,J+1,K))
!	DUM_B = ZG(I,J+1,K+1)*(XG(I,J+1,K)-XG(I,J,K+1))+ZG(I,J+1,K)*(XG(I,J,K+1)-XG(I,J+1,K+1))+ZG(I,J,K+1)*(XG(I,J+1,K+1)-XG(I,J+1,K))
!	DUM_C = XG(I,J+1,K+1)*(YG(I,J+1,K)-YG(I,J,K+1))+XG(I,J+1,K)*(YG(I,J,K+1)-YG(I,J+1,K+1))+XG(I,J,K+1)*(YG(I,J+1,K+1)-YG(I,J+1,K))
!	! D : LENGTH OF THE NORMAL VECTOR
!	DUM_D = XG(I,J+1,K+1)*(YG(I,J+1,K)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K))+XG(I,J+1,K)*(YG(I,J,K+1)*ZG(I,J+1,K+1)-YG(I,J+1,K+1)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J+1,K+1)*ZG(I,J+1,K)-YG(I,J+1,K)*ZG(I,J+1,K+1))
!	DUM_D = -1.D0*DUM_D
!	
!	IF(DUM_D.LT.0.D0) THEN
!		DUM_A = -1.D0*DUM_A
!		DUM_B = -1.D0*DUM_B
!		DUM_C = -1.D0*DUM_C
!		DUM_D = -1.D0*DUM_D
!	ENDIF
!	
!	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
!	
!	L1 =DSQRT((XG(I,J+1,K+1)-XG(I,J+1,K))**2+(YG(I,J+1,K+1)-YG(I,J+1,K))**2+(ZG(I,J+1,K+1)-ZG(I,J+1,K))**2)
!	L2 =DSQRT((XG(I,J+1,K+1)-XG(I,J,K+1))**2+(YG(I,J+1,K+1)-YG(I,J,K+1))**2+(ZG(I,J+1,K+1)-ZG(I,J,K+1))**2)
!	L3 =DSQRT((XG(I,J+1,K)-XG(I,J,K+1))**2+(YG(I,J+1,K)-YG(I,J,K+1))**2+(ZG(I,J+1,K)-ZG(I,J,K+1))**2)
!	TH1=DACOS((L3**2-L2**2-L1**2)/(-2.D0*L1*L2))
!	DAREA = 0.5D0*L1*L2*DSIN(TH1)
!		
!	IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
!		MASS_INLET	= MASS_INLET +(U(I,J,K)*DUM_A+V(I,J,K)*DUM_B+W(I,J,K)*DUM_C)/DUM*DAREA	!(u,0,0,) .dot. unit(a,b,c)
!	ELSE
!		MASS_OUTLET	= MASS_OUTLET +(U(I,J,K)*DUM_A+V(I,J,K)*DUM_B+W(I,J,K)*DUM_C)/DUM*DAREA	!(u,0,0,) .dot. unit(a,b,c)
!	ENDIF
!	
!	ENDDO
!	ENDDO
!	
!	CALL MPI_ALLREDUCE(MASS_INLET,DUM_A,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
!	CALL MPI_ALLREDUCE(MASS_OUTLET,DUM_B,1,MPI_DOUBLE_PRECISION,MPI_SUM,NMCW,ERR)
!	IF(RANK.EQ.0) WRITE (*,100) DUM_A,-DUM_B,100.D0-(-DUM_B-DUM_A)/DUM_A*100.D0
!100 FORMAT(TR1,'MASS INLET =',F10.6,TR4,'MASS OUTLET =',F10.6,TR4,'BALANCE =',F10.6,'%')
!!----------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!   SUBROUTINE FOR POST-PROCESSING
    IF(IPOST.EQ.1) THEN
		NITER_POS = NITER_POS + 1
		U_AVG(IS:IE,JS:JE,KS:KE) = U_AVG(IS:IE,JS:JE,KS:KE) + U(IS:IE,JS:JE,KS:KE)
		V_AVG(IS:IE,JS:JE,KS:KE) = V_AVG(IS:IE,JS:JE,KS:KE) + V(IS:IE,JS:JE,KS:KE)
		W_AVG(IS:IE,JS:JE,KS:KE) = W_AVG(IS:IE,JS:JE,KS:KE) + W(IS:IE,JS:JE,KS:KE)
		P_AVG(IS:IE,JS:JE,KS:KE) = P_AVG(IS:IE,JS:JE,KS:KE) + P(IS:IE,JS:JE,KS:KE)
		T_AVG(IS:IE,JS:JE,KS:KE) = T_AVG(IS:IE,JS:JE,KS:KE) + T(IS:IE,JS:JE,KS:KE)
        
		UU_AVG(IS:IE,JS:JE,KS:KE) = UU_AVG(IS:IE,JS:JE,KS:KE) + U(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE)
		UV_AVG(IS:IE,JS:JE,KS:KE) = UV_AVG(IS:IE,JS:JE,KS:KE) + U(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE)
        UW_AVG(IS:IE,JS:JE,KS:KE) = UW_AVG(IS:IE,JS:JE,KS:KE) + U(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
        VV_AVG(IS:IE,JS:JE,KS:KE) = VV_AVG(IS:IE,JS:JE,KS:KE) + V(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE)
		VW_AVG(IS:IE,JS:JE,KS:KE) = VW_AVG(IS:IE,JS:JE,KS:KE) + V(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
        WW_AVG(IS:IE,JS:JE,KS:KE) = WW_AVG(IS:IE,JS:JE,KS:KE) + W(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
        
    	SGS11(IS:IE,JS:JE,KS:KE) =	SGS11(IS:IE,JS:JE,KS:KE) - 2.D0*EDDY_VISCOS(IS:IE,JS:JE,KS:KE)*LSS11(IS:IE,JS:JE,KS:KE)
    	SGS12(IS:IE,JS:JE,KS:KE) =	SGS12(IS:IE,JS:JE,KS:KE) - 2.D0*EDDY_VISCOS(IS:IE,JS:JE,KS:KE)*LSS12(IS:IE,JS:JE,KS:KE)
    	SGS13(IS:IE,JS:JE,KS:KE) =	SGS13(IS:IE,JS:JE,KS:KE) - 2.D0*EDDY_VISCOS(IS:IE,JS:JE,KS:KE)*LSS13(IS:IE,JS:JE,KS:KE)
    	SGS22(IS:IE,JS:JE,KS:KE) =	SGS22(IS:IE,JS:JE,KS:KE) - 2.D0*EDDY_VISCOS(IS:IE,JS:JE,KS:KE)*LSS22(IS:IE,JS:JE,KS:KE)
    	SGS23(IS:IE,JS:JE,KS:KE) =	SGS23(IS:IE,JS:JE,KS:KE) - 2.D0*EDDY_VISCOS(IS:IE,JS:JE,KS:KE)*LSS23(IS:IE,JS:JE,KS:KE)
    	SGS33(IS:IE,JS:JE,KS:KE) =	SGS33(IS:IE,JS:JE,KS:KE) - 2.D0*EDDY_VISCOS(IS:IE,JS:JE,KS:KE)*LSS33(IS:IE,JS:JE,KS:KE)
	ENDIF
!---------------------------------------------------------------------------------------------


	END SUBROUTINE CALCCOEF