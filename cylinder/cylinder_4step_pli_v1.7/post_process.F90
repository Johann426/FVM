    SUBROUTINE POST_PROCESS
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------------------------
!   POST-PROCESSING FOR NU,CD,CL
    IF(IPOST.EQ.0) THEN
        OPEN(UNIT=12,FILE='POST_NU_CD_CL.PLT',STATUS='UNKNOWN')
        WRITE(12,*) 'variables=TIME,NU,CD,CL,CL_AMP,CPB,ST,,RMS_CD,RMS_CL'
        CLOSE(12)
		RTIME_MAX  = RTIME_MAX * 2.D0
		NITER_POS = 0
        IRESTART = 1
        TOGGLE = 1
        IPOST = 1
        SIGN = .TRUE.
    ELSE IF(IPOST.EQ.1.AND.IRESTART.EQ.1) THEN
        RTIME_MAX  = RTIME_MAX + (1.D0/ST)*3.D0 ! 3 CYCLE
        NITER_W = (1.D0/ST)/DTIME/24 ! 24 FRAME
        IRESTART = 0
        SIGN = .TRUE.
    ELSE IF(IPOST.EQ.1.AND.IRESTART.EQ.0) THEN
        SIGN = .FALSE.
    ! TIME-AVERAGED VALUE
        U(:,:,:)    = U_TA(:,:,:) / DBLE(NITER_POS)
        V(:,:,:)    = V_TA(:,:,:) / DBLE(NITER_POS)
        P(:,:,:)    = P_TA(:,:,:) / DBLE(NITER_POS)
        T(:,:,:)    = T_TA(:,:,:) / DBLE(NITER_POS)
        VOR(:,:,:)  = VOR_TA(:,:,:) / DBLE(NITER_POS)
        NU  = NU_TA / DBLE(NITER_POS)
        CD  = CD_TA / DBLE(NITER_POS)
        CL  = CL_TA / DBLE(NITER_POS)
        ST  = ST_TA / DBLE(NITER_POS)
        CPB = CPB_TA / DBLE(NITER_POS)
        CL_AMP = CL_AMP_TA/DBLE(NITER_POS)
    ! TIME-AVERAGED RMS VALUE
        RMS_CD  = SQRT(RMS_CD / DBLE(NITER_POS))
        RMS_CL  = SQRT(RMS_CL / DBLE(NITER_POS))
    ! TIME-AVERAGED FIELD RESULTS
        OPEN(UNIT=41,FILE='TIME-AVERAGED FIELD.PLT',STATUS='UNKNOWN')
        WRITE(41,*) 'title=""'
        WRITE(41,*) 'variables=x,y,u,v,P,Vorticity,T'
        WRITE(41,*) 'zone t="',RTIME,'",i=',NI+1,',j=',NJ+1
        K=KMON
        DO J=0,NJ
        DO I=0,NI
        WRITE(41,*) X(I,J,K),Y(I,J,K),U(I,J,K),V(I,J,K),P(I,J,K),VOR(I,J,K),T(I,J,K)
   	    ENDDO
   	    ENDDO
   	    CLOSE(41)
    ! TIME-AVERAGED NON-DIMENSIONAL COEFFICIENTS RESULTS
        PRINT *, ''
        PRINT *, 'TIME-AVERAGED NON-DIMENSIONAL COEFFICIENTS RESULTS'
	    WRITE (*,106) NU,CD,CL
    	WRITE (*,107) CL_AMP,CPB,ST
	    WRITE (*,108) RMS_CD,RMS_CL
106     FORMAT(TR1'NU =',F10.6,TR4'CD =',F10.6,TR4'CL =',F10.6)
107     FORMAT(TR1'CL_AMP =',F10.6,TR4'CPB =',F10.6,TR4'ST =',F10.6)
108     FORMAT(TR1'RMS_CD =',F10.6,TR4'RMS_CL =',F10.6)
        
        OPEN(UNIT=13,FILE='TIME_AVERAGED_RESULTS.DAT',STATUS='UNKNOWN')
        WRITE(13,*) 'variables=TIME,NU,CD,CL,CL_AMP,CPB,ST,RMS_CD,RMS_CL'
        WRITE(13,*) RTIME,NU,CD,CL,CL_AMP,CPB,ST,RMS_CD,RMS_CL
        CLOSE(13)
    ENDIF
!---------------------------------------------------------------------------------------------
    END SUBROUTINE POST_PROCESS