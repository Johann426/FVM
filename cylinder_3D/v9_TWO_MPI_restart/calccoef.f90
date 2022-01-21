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
	
	CALL MPI_ALLREDUCE(CFL_MAX,DUM,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_NEW,ERR)
	
	CFL_MAX = DUM
	
!	!	USER-DEFINED CFL(Courant-Friedrichs-Lewy) CONDITION
!	IF(CFL_MAX.GT.CFL_USE) THEN
!        DTIME = DTIME + 0.05D0*DTIME    &
!                *DBLE(IPOST-1)  ! CANCELED BY POST-PROCESSING PARAMETER
!    ELSE
!        DTIME = DTIME + 0.05D0*DTIME    &
!                *DBLE(1-IPOST)  ! CANCELED BY POST-PROCESSING PAREMETER
!	ENDIF
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   CPB CALCULATION
	DO K=KS1,KE
	DO J=JS1,JE
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
    
    DO K=KS1,KE1
	I=NI
	    DUM_A = PG(I,(NJ+1)/2,(KS+KE)/2)		! Far field pressure
	I=1
	    DUM_B = PG(I,1,(KS+KE)/2)			! Base suction pressure
	    
	CPB  = CPB + 2.0D0 * ( DUM_A - DUM_B )
	ENDDO
	
	CPB = CPB / DBLE(KE1-KS1+1)
	
	CALL MPI_ALLREDUCE(CPB,DUM,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	
	CPB = DUM / DBLE(NPRC)

!---------------------------------------------------------------------------
!   CALCULATION OF DUDN
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
	! (A,B,C) : NORMAL VECTOR OF THE DISCRETIZED CYLINDER SURFACE
	DUM_A = YG(I,J,K)*(ZG(I,J+1,K)-ZG(I,J,K+1))+YG(I,J+1,K)*(ZG(I,J,K+1)-ZG(I,J,K))+YG(I,J,K+1)*(ZG(I,J,K)-ZG(I,J+1,K))
	DUM_B = ZG(I,J,K)*(XG(I,J+1,K)-XG(I,J,K+1))+ZG(I,J+1,K)*(XG(I,J,K+1)-XG(I,J,K))+ZG(I,J,K+1)*(XG(I,J,K)-XG(I,J+1,K))
	DUM_C = XG(I,J,K)*(YG(I,J+1,K)-YG(I,J,K+1))+XG(I,J+1,K)*(YG(I,J,K+1)-YG(I,J,K))+XG(I,J,K+1)*(YG(I,J,K)-YG(I,J+1,K))
	! D : LENGTH OF THE NORMAL VECTOR
	DUM_D = XG(I,J,K)*(YG(I,J+1,K)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K))+XG(I,J+1,K)*(YG(I,J,K+1)*ZG(I,J,K)-YG(I,J,K)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J,K)*ZG(I,J+1,K)-YG(I,J+1,K)*ZG(I,J,K))
	DUM_D = -1.D0*DUM_D
	
	IF(-1.D0*DUM_D.GT.0.D0) THEN
		DUM_A = -1.D0*DUM_A
		DUM_B = -1.D0*DUM_B
		DUM_C = -1.D0*DUM_C
		DUM_D = -1.D0*DUM_D
	ENDIF
	
	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
	
	D_N = DABS(DUM_A*X(I,J,K)+DUM_B*Y(I,J,K)+DUM_C*Z(I,J,K)+DUM_D) / DUM
	
	P_DRAG = P_DRAG - 2.0D0*P(I,J,K)*AREA1(I,J,K)*DUM_A/DUM
	
	P_LIFT_Y = P_LIFT_Y - 2.0D0*P(I,J,K)*AREA1(I,J,K)*DUM_B/DUM
	
	P_LIFT_Z = P_LIFT_Z - 2.0D0*P(I,J,K)*AREA1(I,J,K)*DUM_C/DUM
	
	T_DRAG = T_DRAG + (2.0D0/RE)*U(I,J,K)*DABS(DUM_B)/DUM/ D_N*AREA1(I,J,K)		! (u,0,0) x (A,B,C)
	
	T_LIFT_Y = T_LIFT_Y + (2.0D0/RE)*V(I,J,K)*DABS(DUM_A)/DUM/ D_N*AREA1(I,J,K)	! (0,v,0) x (A,B,C)
	
	T_LIFT_Z = T_LIFT_Z + (2.0D0/RE)*W(I,J,K)*DABS(DUM_B)/DUM/ D_N*AREA1(I,J,K)	! (0,0,w) x (A,B,C)
	
	DTDN(I,J,K) = (T(I,J,K)-T_WALL) / D_N
	
	SUM_NU  = SUM_NU + DTDN(I,J,K)*AREA1(I,J,K)
	
    ENDDO
	ENDDO
    
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS1
	! (A,B,C) : NORMAL VECTOR OF THE DISCRETIZED CYLINDER SURFACE
	DUM_A = YG(I,J+1,K+1)*(ZG(I,J+1,K)-ZG(I,J,K+1))+YG(I,J+1,K)*(ZG(I,J,K+1)-ZG(I,J+1,K+1))+YG(I,J,K+1)*(ZG(I,J+1,K+1)-ZG(I,J+1,K))
	DUM_B = ZG(I,J+1,K+1)*(XG(I,J+1,K)-XG(I,J,K+1))+ZG(I,J+1,K)*(XG(I,J,K+1)-XG(I,J+1,K+1))+ZG(I,J,K+1)*(XG(I,J+1,K+1)-XG(I,J+1,K))
	DUM_C = XG(I,J+1,K+1)*(YG(I,J+1,K)-YG(I,J,K+1))+XG(I,J+1,K)*(YG(I,J,K+1)-YG(I,J+1,K+1))+XG(I,J,K+1)*(YG(I,J+1,K+1)-YG(I,J+1,K))
	! D : LENGTH OF THE NORMAL VECTOR
	DUM_D = XG(I,J+1,K+1)*(YG(I,J+1,K)*ZG(I,J,K+1)-YG(I,J,K+1)*ZG(I,J+1,K))+XG(I,J+1,K)*(YG(I,J,K+1)*ZG(I,J+1,K+1)-YG(I,J+1,K+1)*ZG(I,J,K+1))+XG(I,J,K+1)*(YG(I,J+1,K+1)*ZG(I,J+1,K)-YG(I,J+1,K)*ZG(I,J+1,K+1))
	DUM_D = -1.D0*DUM_D
	
	IF(-1.D0*DUM_D.GT.0.D0) THEN
		DUM_A = -1.D0*DUM_A
		DUM_B = -1.D0*DUM_B
		DUM_C = -1.D0*DUM_C
		DUM_D = -1.D0*DUM_D
	ENDIF
	
	DUM = DSQRT(DUM_A**2+DUM_B**2+DUM_C**2)
	
	D_N = DABS(DUM_A*X(I,J,K)+DUM_B*Y(I,J,K)+DUM_C*Z(I,J,K)+DUM_D) / DUM
	
	P_DRAG = P_DRAG - 2.0D0*P(I,J,K)*AREA2(I,J,K)*DUM_A/DUM
	
	P_LIFT_Y = P_LIFT_Y - 2.0D0*P(I,J,K)*AREA2(I,J,K)*DUM_B/DUM
	
	P_LIFT_Z = P_LIFT_Z - 2.0D0*P(I,J,K)*AREA2(I,J,K)*DUM_C/DUM
	
	T_DRAG = T_DRAG + (2.0D0/RE)*U(I,J,K)*DABS(DUM_B)/DUM/ D_N*AREA2(I,J,K)		! (u,0,0) x (A,B,C)
	
	T_LIFT_Y = T_LIFT_Y + (2.0D0/RE)*V(I,J,K)*DABS(DUM_A)/DUM/ D_N*AREA2(I,J,K)	! (0,v,0) x (A,B,C)
	
	T_LIFT_Z = T_LIFT_Z + (2.0D0/RE)*W(I,J,K)*DABS(DUM_B)/DUM/ D_N*AREA2(I,J,K)	! (0,0,w) x (A,B,C)
	
	DTDN(I,J,K) = (T(I,J,K)-T_WALL) / D_N
	
	SUM_NU  = SUM_NU + DTDN(I,J,K)*AREA2(I,J,K)
	
    ENDDO
	ENDDO
	
    CALL MPI_ALLREDUCE(P_DRAG,CD_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	CALL MPI_ALLREDUCE(P_LIFT_Y,CL_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
!	CALL MPI_ALLREDUCE(P_LIFT_Z,CL_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	
	CALL MPI_ALLREDUCE(T_DRAG,CD_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	CALL MPI_ALLREDUCE(T_LIFT_Y,CL_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
!	CALL MPI_ALLREDUCE(T_LIFT_Z,CL_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	
	CALL MPI_ALLREDUCE(SUM_NU,NU,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	
				IF(RANK.EQ.0) THEN
	
	CD_P = CD_P / REF_A	!P_DRAG
	CL_P = CL_P / REF_A	!P_LIFT

	CD_V = CD_V / REF_A	!T_DRAG
	CL_V = CL_V / REF_A	!T_LIFT

	CL_NM1 = CL
	
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
101 FORMAT(TR6'Time',F10.3,TR12'ITER',I10,TR12'CFL',F10.3)
102 FORMAT(TR1'NU =',F10.6)
103 FORMAT(TR1'CD =',F10.6,TR4'CD_P =',F10.6,TR4'CD_V =',F10.6)
104 FORMAT(TR1'CL =',F10.6,TR4'CL_P =',F10.6,TR4'CL_V =',F10.6,TR4'CPB =',F10.6)
	
				ENDIF	! PROCESS ZERO ONLY
	
!------------------------------------------------------------------------------------------
!   SUBROUTINE FOR POST-PROCESSING
    IF(IPOST.AND.RANK.EQ.0) THEN
        NITER_POS = NITER_POS + 1
        DUM = DBLE(TOGGLE)*(CL-CL_NM1)
        ! STROUHAL NUMBER AND CL_AMPLITUDE
        IF(DUM.LT.0.D0) THEN
        TOGGLE  = -1 * TOGGLE
            ST  = 1.D0/(RTIME-TIME2)
            TIME2 = TIME1
            TIME1 = RTIME
            CL_AMP  = ABS(CL-CL1)/2.D0
            CL1 = CL
        ENDIF
        
    ! SUMMATION FOR TIME-AVERAGED VALUE
!        U_TA(:,:,:)    = U_TA(:,:,:) + U(:,:,:)
!        V_TA(:,:,:)    = V_TA(:,:,:) + V(:,:,:)
!        P_TA(:,:,:)    = P_TA(:,:,:) + P(:,:,:)
!        T_TA(:,:,:)    = T_TA(:,:,:) + T(:,:,:)
        NU_TA   = NU_TA + NU
        CD_TA   = CD_TA + CD
        CL_TA   = CL_TA + CL
        ST_TA   = ST_TA + ST
        CPB_TA  = CPB_TA + CPB
        CL_AMP_TA = CL_AMP_TA + CL_AMP
    ! SUMMATION FOR TIME-AVERAGED RMS VALUE
        RMS_CD  = RMS_CD + CD**2
        RMS_CL  = RMS_CL + CL**2
    ! TIME-AVERAGED VALUE
    
		OPEN(UNIT=12,FILE='POST_NU_CD_CL.PLT',STATUS='OLD',ACCESS='APPEND')
        WRITE(12,*) RTIME,NU_TA/DBLE(NITER_POS),CD_TA/DBLE(NITER_POS),CL_TA/DBLE(NITER_POS)     &
            ,CL_AMP_TA/DBLE(NITER_POS),CPB_TA/DBLE(NITER_POS),ST_TA/DBLE(NITER_POS)             &
            ,SQRT(RMS_CD/DBLE(NITER_POS)),SQRT(RMS_CL/DBLE(NITER_POS))
		CLOSE(12)
    
        WRITE (*,105) ST,CL_AMP
105     FORMAT(TR1'ST =',F10.6,TR4'CL_AMP =',F10.6)
    ENDIF !179
    
    CALL MPI_BCAST(ST,1,MPI_DOUBLE_PRECISION,0,COMM_NEW,ERR)
    
!---------------------------------------------------------------------------------------------
	END SUBROUTINE CALCCOEF