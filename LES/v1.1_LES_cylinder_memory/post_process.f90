	SUBROUTINE POST_PROCESS
	
	USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE
!---------------------------------------------------------------------------------------------
!   POST-PROCESSING FOR NU,CD,CL
	IF(IPOST.EQ.0) THEN
!		ALLOCATE(U_TA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
!		ALLOCATE(V_TA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
!		ALLOCATE(W_TA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
!		ALLOCATE(P_TA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
!		ALLOCATE(T_TA(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
		ALLOCATE(SH_S(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
		ALLOCATE(NU_S(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
		ALLOCATE(CP_S(IS1,JS1:JE1,KS1:KE1),STAT=ALLOC_ERR)
		SH_S(:,:,:) = 0.D0
		NU_S(:,:,:) = 0.D0
		CP_S(:,:,:) = 0.D0
		RTIME_MAX  = RTIME_END
		NITER_POS = 0
		IPOST = 1
		SIGN = .TRUE.
	ELSE IF(IPOST.EQ.1) THEN
!		U(:,:,:)    = U_TA(:,:,:) / DBLE(NITER_POS)
!		V(:,:,:)    = V_TA(:,:,:) / DBLE(NITER_POS)
!		W(:,:,:)    = W_TA(:,:,:) / DBLE(NITER_POS)
!		P(:,:,:)    = P_TA(:,:,:) / DBLE(NITER_POS)
!		T(:,:,:)    = T_TA(:,:,:) / DBLE(NITER_POS)
		SH_S(:,:,:) = SH_S(:,:,:) / DBLE(NITER_POS)
		NU_S(:,:,:) = NU_S(:,:,:) / DBLE(NITER_POS)
		CP_S(:,:,:) = CP_S(:,:,:) / DBLE(NITER_POS)
		
		CONTOUR = 'POST_OO.PLT'
		WRITE(CONTOUR(6:7),'(I2.2)') INT(RANK)
		OPEN(UNIT=100+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
		WRITE(100+RANK,*) 'TITLE="SURFACE COEF."'
		WRITE(100+RANK,*) 'variables="x","y","z","Theta","Shear Strain","Nu","Cp"'
		WRITE(100+RANK,24) 1,JE1-JS1+1,KE1-KS1+1
24		FORMAT(TR1,'ZONE T="POST", I=',I6,', j=',I6,', k=',I6)
		DO K=KS1,KE1
		DO J=JS1,JE1
		I=IS1
			WRITE(100+RANK,*) X(I,J,K),Y(I,J,K),Z(I,J,K),THETA(I,J,K)*180.D0/PI,SH_S(I,J,K),-NU_S(I,J,K),CP_S(I,J,K)
		ENDDO
		ENDDO
		CLOSE(100+RANK)
		SIGN = .FALSE.
	ENDIF
!---------------------------------------------------------------------------------------------
	END SUBROUTINE POST_PROCESS