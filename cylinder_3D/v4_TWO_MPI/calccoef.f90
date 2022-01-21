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
	
	!	USER-DEFINED CFL(Courant-Friedrichs-Lewy) CONDITION
	IF(CFL_MAX.GT.CFL_USE) THEN
        DTIME = DTIME + 0.05D0*DTIME    &
                *DBLE(IPOST-1)  ! CANCELED BY POST-PROCESSING PARAMETER
    ELSE
        DTIME = DTIME + 0.05D0*DTIME    &
                *DBLE(1-IPOST)  ! CANCELED BY POST-PROCESSING PAREMETER
	ENDIF
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   CPB CALCULATION
	DO K=KS1,KE
	DO J=JS1,JE
	I=IS1
		PG(I,J,K) = 0.25D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K))
		UG(I,J,K) = 0.25D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K))
		VG(I,J,K) = 0.25D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K))
	I=IE
		PG(I,J,K) = 0.25D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K))
		UG(I,J,K) = 0.25D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K))
		VG(I,J,K) = 0.25D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K))
	ENDDO
    ENDDO
    
    CPB = 0.D0
    
    DO K=KS1,KE1
	I=IE
	    DUM_A = PG(I,(NJ+1)/2,NK/2)		! Far field pressure
	I=IS1
	    DUM_B = PG(I,1,NK/2)			! Base suction pressure
	    
	CPB  = CPB + 2.0D0 * ( DUM_A - DUM_B )
	ENDDO
	
	CPB = CPB / DBLE(KE1-KS1+1)
	
	CALL MPI_ALLREDUCE(CPB,DUM,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	
	CPB = DUM / DBLE(NPRC)

	!   CALCULATION OF DUDN
	P_DRAG = 0.D0
	P_LIFT = 0.D0
    
    T_DRAG = 0.D0
	T_LIFT = 0.D0
	
	SUM_NU = 0.D0
	
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS1
!	L1 = ((X(I,J,K)   -XG(I,J  ,K))**2.D0+(Y(I,J,K)   -YG(I,J  ,K))**2.D0)**0.5D0
!	L2 = ((X(I,J,K)   -XG(I,J+1,K))**2.D0+(Y(I,J,K)   -YG(I,J+1,K))**2.D0)**0.5D0
!	L3 = ((XG(I,J+1,K)-XG(I,J  ,K))**2.D0+(YG(I,J+1,K)-YG(I,J  ,K))**2.D0)**0.5D0
	DUM_A = 0.5D0*(XG(I,J,K)+XG(I,J,K+1))
	DUM_B = 0.5D0*(YG(I,J,K)+YG(I,J,K+1))
	DUM_C = 0.5D0*(XG(I,J+1,K)+XG(I,J+1,K+1))
	DUM_D = 0.5D0*(YG(I,J+1,K)+YG(I,J+1,K))
	L1 = ((X(I,J,K)-DUM_A)**2.D0	+(Y(I,J,K)-DUM_B)**2.D0)**0.5D0
	L2 = ((X(I,J,K)-DUM_C)**2.D0+(Y(I,J,K)-DUM_D)**2.D0)**0.5D0
	L3 = ((DUM_C-DUM_A)**2.D0+(DUM_D-DUM_B)**2.D0)**0.5D0
	TH1= DACOS( (L1**2+L2**2-L3**2)/(2.D0*L1*L2) )
	DN1= L1*L2/L3*DSIN(TH1)
!	LA  = (YG(I,J+1,K)-YG(I,J,K))/(XG(I,J+1,K)-XG(I,J,K))
!	LB  = YG(I,J,K)-LA*XG(I,J,K)
!	LA_ = -1.D0/LA
!	LB_ = Y(I,J,K)-LA_*X(I,J,K)
    LA  = (DUM_D-DUM_B)/(DUM_C-DUM_A)
    LB  = DUM_B-LA*DUM_A
    LA_ = -1.D0/LA
    LB_ = Y(I,J,K)-LA_*X(I,J,K)
    X1 = (LB_-LB)/(LA-LA_)
    Y1 = LA*(LB_-LB)/(LA-LA_)+LB
    DN2=((Y(I,J,K)-Y1)**2.D0+(X(I,J,K)-X1)**2.D0)**0.5D0
    IF(ABS(DN2-DN1).GT.1.0D-15) THEN
        WRITE(*,*) 'POST-PROCESSING ERROR IS HIGH'
    ENDIF
    TH2= DATAN((Y(I,J,K)-Y1)/(X(I,J,K)-X1))
      IF     (TH2.GE.0.0D0.AND.DASIN((Y(I,J,K)-Y1)/DN2).GT.0.0D0) THEN
		TH2 = TH2
 	  
	  ELSEIF (TH2.LT.0.0D0.AND.DASIN((Y(I,J,K)-Y1)/DN2).GT.0.0D0) THEN
		TH2 = PI+TH2

 	  ELSEIF (TH2.GE.0.0D0.AND.DASIN((Y(I,J,K)-Y1)/DN2).LE.0.0D0) THEN
		TH2 = PI+TH2

	  ELSEIF (TH2.LT.0.0D0.AND.DASIN((Y(I,J,K)-Y1)/DN2).LT.0.0D0) THEN
		TH2 = 2.0D0*PI+TH2
	  ENDIF
    DUDN(I,J,K) = &!<--- DU/DN = D(Usin0-Vcos0)/DN
	            (U(I,J,K)*DSIN(TH2)-V(I,J,K)*DCOS(TH2))	& !DU
	            /DN1
    DTDN(I,J,K) = &!<--- DT/DN = DT/DN
	            (T(I,J,K)-T_WALL)						& !DT
	            /DN1
!   FORM DRAG FORCE
		P_DRAG = P_DRAG - 2.0D0*P(I,J,K)*AREA(I,J,K)*DCOS(THETA(I,J,K))
!   FORM LIFT FORCE
		P_LIFT = P_LIFT - 2.0D0*P(I,J,K)*AREA(I,J,K)*DSIN(THETA(I,J,K))
!   FRICTION DRAG FORCE
		T_DRAG = T_DRAG + (2.0D0/RE)*DUDN(I,J,K)*AREA(I,J,K)*DSIN(THETA(I,J,K))
!   FRICTION LIFT FORCE
		T_LIFT = T_LIFT - (2.0D0/RE)*DUDN(I,J,K)*AREA(I,J,K)*DCOS(THETA(I,J,K))
		
		SUM_NU  = SUM_NU + DTDN(I,J,K)*AREA(I,J,K)
    ENDDO
	ENDDO
    
    CALL MPI_ALLREDUCE(P_DRAG,CD_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	CALL MPI_ALLREDUCE(P_LIFT,CL_P,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	CALL MPI_ALLREDUCE(T_DRAG,CD_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
	CALL MPI_ALLREDUCE(T_LIFT,CL_V,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_NEW,ERR)
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
	
	WRITE (11,*) RTIME,NU,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,CFL_MAX ! PLOTTING

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
    IF(IPOST) THEN
        CALL CALC_VOR
        NITER_POS = NITER_POS + 1
        
        ! STROUHAL NUMBER AND CL_AMPLITUDE
        IF(DBLE(TOGGLE)*(CL-CL_NM1).LT.0.0D0) THEN
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
!        VOR_TA(:,:,:)  = VOR_TA(:,:,:) + VOR(:,:,:)
!        NU_TA   = NU_TA + NU
!        CD_TA   = CD_TA + CD
!        CL_TA   = CL_TA + CL
!        ST_TA   = ST_TA + ST
!        CPB_TA  = CPB_TA + CPB
        CL_AMP_TA = CL_AMP_TA + CL_AMP
    ! SUMMATION FOR TIME-AVERAGED RMS VALUE
!        RMS_CD  = RMS_CD + CD**2.D0
!        RMS_CL  = RMS_CL + CL**2.D0
    ! TIME-AVERAGED VALUE
    
!        WRITE(12,*) RTIME,NU_TA/DBLE(NITER_POS),CD_TA/DBLE(NITER_POS),CL_TA/DBLE(NITER_POS)     &
!            ,CL_AMP_TA/DBLE(NITER_POS),CPB_TA/DBLE(NITER_POS),ST_TA/DBLE(NITER_POS)             &
!            ,SQRT(RMS_CD/DBLE(NITER_POS)),SQRT(RMS_CL/DBLE(NITER_POS))
    
        WRITE (*,105) ST,CL_AMP
105     FORMAT(TR1'ST =',F10.6,TR4'CL_AMP =',F10.6)
    ENDIF
!---------------------------------------------------------------------------------------------
	END SUBROUTINE CALCCOEF