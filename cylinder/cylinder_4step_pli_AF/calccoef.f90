    SUBROUTINE CALCCOEF

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
!---------------------------------------------------------------------------------
!   CFL(Courant-Friedrichs-Lewy) CALCULATION
	CFL_MAX = 0.0D0
	DO I=sx1,ex1
	DO J=sy1,ey1
	DO K=2,NK-1        
		CFL_LOCAL = DTIME *(ABS(XIX(I,J,K)*U(I,J,K) + XIY(I,J,K)*V(I,J,K))  & 
     		              + ABS(ETX(I,J,K)*U(I,J,K) + ETY(I,J,K)*V(I,J,K))  )
		CFL_MAX = DMAX1(CFL_MAX,CFL_LOCAL)
	ENDDO
	ENDDO
	ENDDO
!   USER-DEFINED CFL(Courant-Friedrichs-Lewy) CONDITION
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
	DO K=2,NKM1
	DO I=1,3
	DO J=1,nj
		PG(I,J,K) = 0.25D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K))
		UG(I,J,K) = 0.25D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K))
		VG(I,J,K) = 0.25D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K))
	ENDDO
	ENDDO
    ENDDO
    
	I=NI
	    DUM_A = PG(I,NJ/2,2)        ! Far field pressure
	I=1
	    DUM_B = PG(I,1,2)           ! Base suction pressure
	CPB  = 2.0D0 * ( DUM_A - DUM_B  )
!---------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!   CALCULATION OF DUDN on orthogonal grid
	P_DRAG = 0.D0
	P_LIFT = 0.D0
    
    T_DRAG = 0.D0
	T_LIFT = 0.D0
	
	SUM_NU = 0.D0
	
	DO K=2,NKM1
	DO J=1,NJ-1
	I=1
	L1 = ((X(I,J,K)   -XG(I,J  ,K))**2.D0+(Y(I,J,K)   -YG(I,J  ,K))**2.D0)**0.5D0
    L2 = ((X(I,J,K)   -XG(I,J+1,K))**2.D0+(Y(I,J,K)   -YG(I,J+1,K))**2.D0)**0.5D0
    L3 = ((XG(I,J+1,K)-XG(I,J  ,K))**2.D0+(YG(I,J+1,K)-YG(I,J  ,K))**2.D0)**0.5D0
    TH1= DACOS( (L1**2+L2**2-L3**2)/(2.D0*L1*L2) )
    DN1= L1*L2/L3*DSIN(TH1)
    LA  = (YG(I,J+1,K)-YG(I,J,K))/(XG(I,J+1,K)-XG(I,J,K))
    LB  = YG(I,J,K)-LA*XG(I,J,K)
    LA_ = -1.D0/LA
    LB_ = Y(I,J,K)-LA_*X(I,J,K)
    X1 = (LB_-LB)/(LA-LA_)
    Y1 = LA*(LB_-LB)/(LA-LA_)+LB
    DN2=((Y(I,J,K)-Y1)**2.D0+(X(I,J,K)-X1)**2.D0)**0.5D0
    IF(ABS(DN2-DN1).GT.1.0D-15) THEN
        WRITE(*,*) 'POST-PROCESSING ERROR IS HIGH'
        PAUSE
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
	            (U(I,J,K)*DSIN(TH2)-V(I,J,K)*DCOS(TH2)-VEL_ROT) & !DU
	            /DN1
    DTDN(I,J,K) = &!<--- DT/DN = DT/DN
	            (T(I,J,K)-T_WALL)                               & !DT
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
!-------------------------------------------------------------------------------
    CL_NM1 = CL
    
	CD_P = P_DRAG
	CL_P = P_LIFT

	CD_V = T_DRAG
	CL_V = T_LIFT

	CD  = CD_P + CD_V
	CL  = CL_P + CL_V
    
    NU  = -SUM_NU / SUM_AREA
    
	WRITE (11,*) RTIME,NU,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,CFL_MAX ! PLOTTING
!------------------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------------------

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
        U_TA(:,:,:)    = U_TA(:,:,:) + U(:,:,:)
        V_TA(:,:,:)    = V_TA(:,:,:) + V(:,:,:)
        P_TA(:,:,:)    = P_TA(:,:,:) + P(:,:,:)
        T_TA(:,:,:)    = T_TA(:,:,:) + T(:,:,:)
        VOR_TA(:,:,:)  = VOR_TA(:,:,:) + VOR(:,:,:)
        NU_TA   = NU_TA + NU
        CD_TA   = CD_TA + CD
        CL_TA   = CL_TA + CL
        ST_TA   = ST_TA + ST
        CPB_TA  = CPB_TA + CPB
        CL_AMP_TA = CL_AMP_TA + CL_AMP
    ! SUMMATION FOR TIME-AVERAGED RMS VALUE
        RMS_CD  = RMS_CD + CD**2.D0
        RMS_CL  = RMS_CL + CL**2.D0
    ! TIME-AVERAGED VALUE
        WRITE(12,*) RTIME,NU_TA/DBLE(NITER_POS),CD_TA/DBLE(NITER_POS),CL_TA/DBLE(NITER_POS)     &
            ,CL_AMP_TA/DBLE(NITER_POS),CPB_TA/DBLE(NITER_POS),ST_TA/DBLE(NITER_POS)             &
            ,SQRT(RMS_CD/DBLE(NITER_POS)),SQRT(RMS_CL/DBLE(NITER_POS))
    
        WRITE (*,105) ST,CL_AMP
105     FORMAT(TR1'ST =',F10.6,TR4'CL_AMP =',F10.6)
    ENDIF
!---------------------------------------------------------------------------------------------
	END SUBROUTINE CALCCOEF