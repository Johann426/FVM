	SUBROUTINE BCW

	USE COMDAT_SHARED
    
	IMPLICIT NONE
	SAVE
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	DO K=0,NK
	DO J=0,NJ
	I=0
		W_F(SX1,J,K) = 0.0D0

    AP(SX1,J,K) = AP(SX1,J,K) +AW(SX1,J,K)
    SU(SX1,J,K) = SU(SX1,J,K) +2.D0*AW(SX1,J,K)*W_F(SX1,J,K)
    AW(SX1,J,K) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	W(I,J,K) = -W(I+1,J,K) + 2.0D0*W_F(SX1,J,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=0,NK
	DO J=0,NJ
	I=NI
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    W_F(NIM1,J,K) = WIN
	    ELSE
	    W_F(NIM1,J,K) = 0.5D0*(W_N(I,J,K)+W_N(I-1,J,K)) - DTIME*Uin       &
                    *(0.5D0*(W_N(I,J,K)+W_N(I-1,J,K))-W_N(I-1,J,K)) &
                    /(0.5D0*(  X(I,J,K)+  X(I-1,J,K)) - X(I-1,J,K))	
	    ENDIF
    AP(NIM1,J,K) = AP(NIM1,J,K) +AE(NIM1,J,K)
    SU(NIM1,J,K) = SU(NIM1,J,K) +2.D0*AE(NIM1,J,K)*W_F(NIM1,J,K)
    AE(NIM1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	W(I,J,K) = -W(I-1,J,K) + 2.0D0*W_F(NIM1,J,K)
    ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=0,NK
	DO I=0,NI
    J=0
	    W(I,J,K) = W(I,NJ-1,K)
	J=NJ
		W(I,J,K) = W(I,1,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
	DO J=0,NJ
	DO I=0,NI
	K=0
!		AP(I,J,SZ1) = AP(I,J,SZ1)-AB(I,J,SZ1)
!		AB(I,J,SZ1) = 0.D0
!		W(I,J,K) = W(I,J,K+1)
		W(I,J,K) = W(I,J,NKM1)
	K=NK
!		AP(I,J,EZ1) = AP(I,J,EZ1)-AT(I,J,EZ1)
!		AT(I,J,EZ1) = 0.D0
!		W(I,J,K) = W(I,J,K-1)
		W(I,J,K) = W(I,J,1)
	ENDDO
	ENDDO
	
	END SUBROUTINE BCW