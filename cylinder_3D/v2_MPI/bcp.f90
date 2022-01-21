    SUBROUTINE BCP

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)
	
	DUM_3D=P
!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	DO K=KS,KE
	DO J=JS1,JE1
	I=IS
    P(I,J,K) = 1.D0/3.D0*((4.D0*P(I+1,J,k)-P(I+2,J,K))         &
             + Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))
    
    I=IE
    P(I,J,K) = 1.D0/3.D0*((4.D0*P(I-1,J,k)-P(I-2,J,K))         &
             - Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))
	ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=KS,KE
	DO I=IS,IE
    J=JS
	    P(I,J,K) = P(I,JE1,K)
	J=JE
		P(I,J,K) = P(I,JS1,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
!	DO J=0,NJ
!	DO I=0,NI
!	K=0
!		P(I,J,K) = P(I,J,K+1)	!U(I,J,K) = U(I,J,NKM1)
!	K=NK
!		P(I,J,K) = P(I,J,K-1)	!U(I,J,K) = U(I,J,1)
!	ENDDO
!	ENDDO
	
	DO J=JS,JE
	
	DO I=IS,IE
	
	IF(RANK.EQ.0)		CALL MPI_ISEND(P(I,J,KS1),1,MPI_DOUBLE_PRECISION,NPRC-1,TAG1,MPI_COMM_WORLD,ISEND2,ERR)
	
	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(P(I,J,KE1),1,MPI_DOUBLE_PRECISION,0,TAG2,MPI_COMM_WORLD,ISEND1,ERR)

	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(P(I,J,KE) ,1,MPI_DOUBLE_PRECISION,0,TAG1,MPI_COMM_WORLD,IRECV2,ERR)
		
	IF(RANK.EQ.0)		CALL MPI_IRECV(P(I,J,KS) ,1,MPI_DOUBLE_PRECISION,NPRC-1,TAG2,MPI_COMM_WORLD,IRECV1,ERR)
	
	CALL MPI_WAIT(ISEND1,STATUS,ERR)
	
	CALL MPI_WAIT(ISEND2,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV1,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV2,STATUS,ERR)
	
	ENDDO
	
	ENDDO
	
    END SUBROUTINE BCP