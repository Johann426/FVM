	SUBROUTINE BCT

	USE COMDAT_SHARED
    
	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)
	
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	DO K=KS,KE
	DO J=JS,JE
	I=IS
		T_F = T_WALL

    AP(IS1,J,K) = AP(IS1,J,K) +AW(IS1,J,K)
    SU(IS1,J,K) = SU(IS1,J,K) +2.D0*AW(IS1,J,K)*T_F
    AW(IS1,J,K) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	T(I,J,K) = -T(I+1,J,K) + 2.0D0*T_F
	ENDDO
	ENDDO
!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=KS,KE
	DO J=JS,JE
	I=IE
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    T_F = TIN
	    ELSE
	    T_F = 0.5D0*(T_N(I,J,K)+T_N(I-1,J,K)) - DTIME*1.D0	&
			*(0.5D0*(T_N(I,J,K)+T_N(I-1,J,K))-T_N(I-1,J,K))	&
			/(0.5D0*(  X(I,J,K)+  X(I-1,J,K)) - X(I-1,J,K))	
	    ENDIF
    AP(IE1,J,K) = AP(IE1,J,K) +AE(IE1,J,K)
    SU(IE1,J,K) = SU(IE1,J,K) +2.D0*AE(IE1,J,K)*U_F
    AE(IE1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	T(I,J,K) = -T(I-1,J,K) + 2.0D0*T_F
    ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=KS,KE
	DO I=IS,IE
    J=JS
	    T(I,J,K) = T(I,JE1,K)
	J=JE
		T(I,J,K) = T(I,JS1,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
!	DO J=0,NJ
!	DO I=0,NI
!	K=0
!!		AP(I,J,SZ1) = AP(I,J,SZ1)-AB(I,J,SZ1)
!!		AB(I,J,SZ1) = 0.D0
!!		T(I,J,K) = T(I,J,K+1)	! JSMMETRIC BOUNDARY
!		T(I,J,K) = T(I,J,NKM1)	! PERIODIC BOUNDARY
!	K=NK
!!		AP(I,J,EZ1) = AP(I,J,EZ1)-AT(I,J,EZ1)
!!		AT(I,J,EZ1) = 0.D0
!!		T(I,J,K) = T(I,J,K-1)	! JSMMETRIC BOUNDARY
!		T(I,J,K) = T(I,J,1)		! PERIODIC BOUNDARY
!	ENDDO
!	ENDDO
	
	DO J=JS,JE
	
	DO I=IS,IE
	
	IF(RANK.EQ.0)		CALL MPI_ISEND(T(I,J,KS1),1,MPI_DOUBLE_PRECISION,NPRC-1,1,MPI_COMM_WORLD,ISEND2,ERR)
	
	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(T(I,J,KE1),1,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ISEND1,ERR)

	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(T(I,J,KE) ,1,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,IRECV2,ERR)
		
	IF(RANK.EQ.0)		CALL MPI_IRECV(T(I,J,KS) ,1,MPI_DOUBLE_PRECISION,NPRC-1,1,MPI_COMM_WORLD,IRECV1,ERR)
	
	CALL MPI_WAIT(ISEND1,STATUS,ERR)
	
	CALL MPI_WAIT(ISEND2,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV1,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV2,STATUS,ERR)
	
	ENDDO
	
	ENDDO
	
	END SUBROUTINE BCT