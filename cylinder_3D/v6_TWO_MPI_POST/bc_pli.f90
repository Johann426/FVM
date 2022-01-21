    SUBROUTINE BC_PLI

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)
	
    DUM_3D=PLI

!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=KS,KE
	DO I=IS,IE
    J=JS
	    PLI(I,J,K) = PLI(I,JE1,K)
	J=JE
		PLI(I,J,K) = PLI(I,JS1,K)
	ENDDO
	ENDDO

!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
	CALL MPI_X(PLI)
	
	IF(RANK.EQ.0)		CALL MPI_ISEND(PLI(IS,JS,KS1),1,APLN,NPRC-1,1,COMM_NEW,ISEND2,ERR)
	
	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(PLI(IS,JS,KE1),1,APLN,0,1,COMM_NEW,ISEND1,ERR)

	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(PLI(IS,JS,KE) ,1,APLN,0,1,COMM_NEW,IRECV2,ERR)
		
	IF(RANK.EQ.0)		CALL MPI_IRECV(PLI(IS,JS,KS) ,1,APLN,NPRC-1,1,COMM_NEW,IRECV1,ERR)
	
	CALL MPI_WAIT(ISEND1,STATUS,ERR)
	
	CALL MPI_WAIT(ISEND2,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV1,STATUS,ERR)
	
	CALL MPI_WAIT(IRECV2,STATUS,ERR)

!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	DO K=KS,KE
	DO J=JS1,JE1
	I=IS
	! WALL BOUNDARY SURFACE
    !AP(I+1,J,K) = AP(I+1,J,K)-AW(I+1,J,K)
    AW(I+1,J,K) = 0.0D0
    PLI(I,J,K) = 1.D0/3.D0*((4.D0*PLI(I+1,J,k)-PLI(I+2,J,K))					&
                +Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))

    I=IE
	! OUTFLOW BOUNDARY SURFACE
    !AP(I-1,J,K) = AP(I-1,J,K)-AE(I-1,J,K)
    AE(I-1,J,K) = 0.0D0
    PLI(I,J,K) = 1.D0/3.D0*((4.D0*PLI(I-1,J,k)-PLI(I-2,J,K))					&
 				-Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))
	ENDDO
	ENDDO
		
!	DO J=JS,JE
!	
!	IF(RANK.EQ.0)		CALL MPI_ISEND(PLI(IS,J,KS1),IE-IS+1,MPI_DOUBLE_PRECISION,NPRC-1,1,COMM_NEW,ISEND2,ERR)
!	
!	IF(RANK.EQ.NPRC-1)	CALL MPI_ISEND(PLI(IS,J,KE1),IE-IS+1,MPI_DOUBLE_PRECISION,0,1,COMM_NEW,ISEND1,ERR)
!
!	IF(RANK.EQ.NPRC-1)	CALL MPI_IRECV(PLI(IS,J,KE) ,IE-IS+1,MPI_DOUBLE_PRECISION,0,1,COMM_NEW,IRECV2,ERR)
!		
!	IF(RANK.EQ.0)		CALL MPI_IRECV(PLI(IS,J,KS) ,IE-IS+1,MPI_DOUBLE_PRECISION,NPRC-1,1,COMM_NEW,IRECV1,ERR)
!	
!	CALL MPI_WAIT(ISEND1,STATUS,ERR)
!	
!	CALL MPI_WAIT(ISEND2,STATUS,ERR)
!	
!	CALL MPI_WAIT(IRECV1,STATUS,ERR)
!	
!	CALL MPI_WAIT(IRECV2,STATUS,ERR)
!	
!	ENDDO
	
    END SUBROUTINE BC_PLI