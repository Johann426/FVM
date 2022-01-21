    SUBROUTINE BCP

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)
	
	DUM_3D=P

!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
    J=JS
	    P(IS:IE,J,KS:KE) = P(IS:IE,JE1,KS:KE)
	J=JE
		P(IS:IE,J,KS:KE) = P(IS:IE,JS1,KS:KE)
	
!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
	CALL MPI_X(P)

	IF(RANK.EQ.0)		CALL MPI_ISEND(P(IS,JS,KS1),1,APLN,NPRC-1,1,COMM_NEW,ISEND2,ERR)
	
	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(P(IS,JS,KE1),1,APLN,0,1,COMM_NEW,ISEND1,ERR)

	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(P(IS,JS,KE) ,1,APLN,0,1,COMM_NEW,IRECV2,ERR)
		
	IF(RANK.EQ.0)		CALL MPI_IRECV(P(IS,JS,KS) ,1,APLN,NPRC-1,1,COMM_NEW,IRECV1,ERR)
	
	CALL MPI_WAIT(ISEND1,STATUS,ERR)
	
	CALL MPI_WAIT(ISEND2,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV1,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV2,STATUS,ERR)

!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	I=IS
    P(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*P(I+1,JS1:JE1,KS:KE)-P(I+2,JS1:JE1,KS:KE))         &
             + Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))
    
    I=IE
    P(I,JS1:JE1,KS:KE) = 1.D0/3.D0*((4.D0*P(I-1,JS1:JE1,KS:KE)-P(I-2,JS1:JE1,KS:KE))         &
             - Q12(I,JS1:JE1,KS:KE)/Q11(I,JS1:JE1,KS:KE) *0.5D0*(DUM_3D(I,JS1+1:JE1+1,KS:KE)-DUM_3D(I,JS1-1:JE1-1,KS:KE)))
	
    END SUBROUTINE BCP