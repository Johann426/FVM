	SUBROUTINE BCU(A,B)

	USE COMDAT_SHARED
    
	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)
	
	DOUBLE PRECISION, INTENT(INOUT) :: A(IS:IE,JS:JE,KS:KE)
	DOUBLE PRECISION, INTENT(IN)	:: B(IS:IE,JS:JE,KS:KE)
	
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	DO K=KS,KE
	DO J=JS,JE
	I=IS
		U_F = 0.D0	!VEL_ROT * -SIN(THETA(I,J,K))

    AP(IS1,J,K) = AP(IS1,J,K) +AW(IS1,J,K)
    SU(IS1,J,K) = SU(IS1,J,K) +2.D0*AW(IS1,J,K)*U_F
    AW(IS1,J,K) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	A(I,J,K) = -A(I+1,J,K) + 2.0D0*U_F
	ENDDO
	ENDDO
!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=KS,KE
	DO J=JS,JE
	I=IE
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    U_F = UIN
	    ELSE
	    U_F = 0.5D0*(B(I,J,K)+B(I-1,J,K)) - DTIME*1.D0	&
            *(0.5D0*(B(I,J,K)+B(I-1,J,K))-B(I-1,J,K))	&
            /(0.5D0*(  X(I,J,K)+  X(I-1,J,K)) - X(I-1,J,K))
	    ENDIF
    AP(IE1,J,K) = AP(IE1,J,K) +AE(IE1,J,K)
    SU(IE1,J,K) = SU(IE1,J,K) +2.D0*AE(IE1,J,K)*U_F
    AE(IE1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	A(I,J,K) = -A(I-1,J,K) + 2.0D0*U_F
    ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=KS,KE
	DO I=IS,IE
    J=JS
	    A(I,J,K) = A(I,JE1,K)
	J=JE
		A(I,J,K) = A(I,JS1,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
!	DO J=JS,JE
!	DO I=IS,IE
!	K=0
!!		AP(I,J,SZ1) = AP(I,J,SZ1)-AB(I,J,SZ1)
!!		AB(I,J,SZ1) = 0.D0
!!		A(I,J,K) = A(I,J,K+1)	! SYMMETRIC BOUNDARY
!		A(I,J,K) = A(I,J,NKM1)	! PERIODIC BOUNDARY
!	K=NK
!!		AP(I,J,EZ1) = AP(I,J,EZ1)-AT(I,J,EZ1)
!!		AT(I,J,EZ1) = 0.D0
!!		A(I,J,K) = A(I,J,K-1)	! SYMMETRIC BOUNDARY
!		A(I,J,K) = A(I,J,1)		! PERIODIC BOUNDARY
!	ENDDO
!	ENDDO
	
	DO J=JS,JE
	
	IF(RANK.EQ.0)		CALL MPI_ISEND(A(IS,J,KS1),IE-IS+1,MPI_DOUBLE_PRECISION,NPRC-1,1,COMM_NEW,ISEND2,ERR)
	
	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(A(IS,J,KE1),IE-IS+1,MPI_DOUBLE_PRECISION,0,1,COMM_NEW,ISEND1,ERR)

	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(A(IS,J,KE) ,IE-IS+1,MPI_DOUBLE_PRECISION,0,1,COMM_NEW,IRECV2,ERR)
		
	IF(RANK.EQ.0)		CALL MPI_IRECV(A(IS,J,KS) ,IE-IS+1,MPI_DOUBLE_PRECISION,NPRC-1,1,COMM_NEW,IRECV1,ERR)
	
	CALL MPI_WAIT(ISEND1,STATUS,ERR)
	
	CALL MPI_WAIT(ISEND2,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV1,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV2,STATUS,ERR)
	
	ENDDO
	
	END SUBROUTINE BCU