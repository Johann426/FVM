	SUBROUTINE BCT

	USE COMDAT_SHARED
    
	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)

!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
    T(IS:IE,JS,KS:KE) = T(IS:IE,JE1,KS:KE)

	T(IS:IE,JE,KS:KE) = T(IS:IE,JS1,KS:KE)

!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
	CALL MPI_X(T)

	IF(RANK.EQ.0)		CALL MPI_ISEND(T(IS,JS,KS1),1,APLN,NPRC-1,1,COMM_NEW,ISEND1,ERR)
	
	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(T(IS,JS,KE1),1,APLN,0,1,COMM_NEW,ISEND2,ERR)

	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(T(IS,JS,KE) ,1,APLN,0,1,COMM_NEW,IRECV1,ERR)
		
	IF(RANK.EQ.0)		CALL MPI_IRECV(T(IS,JS,KS) ,1,APLN,NPRC-1,1,COMM_NEW,IRECV2,ERR)
	
	CALL MPI_WAIT(ISEND1,STATUS,ERR)
	
	CALL MPI_WAIT(ISEND2,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV1,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV2,STATUS,ERR)
		
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	I=IS
		T_F = T_WALL

    AP(I+1,JS:JE,KS:KE) = AP(I+1,JS:JE,KS:KE) +AW(I+1,JS:JE,KS:KE)
    SU(I+1,JS:JE,KS:KE) = SU(I+1,JS:JE,KS:KE) +2.D0*AW(I+1,JS:JE,KS:KE)*T_F
    AW(I+1,JS:JE,KS:KE) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	T(I,JS:JE,KS:KE) = -T(I+1,JS:JE,KS:KE) + 2.0D0*T_F
	
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
    SU(IE1,J,K) = SU(IE1,J,K) +2.D0*AE(IE1,J,K)*T_F
    AE(IE1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	T(I,J,K) = -T(I-1,J,K) + 2.0D0*T_F
    ENDDO
	ENDDO

	END SUBROUTINE BCT