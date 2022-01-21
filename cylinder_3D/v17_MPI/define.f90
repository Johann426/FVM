	SUBROUTINE DEFINE
	
	USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE

    INCLUDE 'mpif.h'
	
	!	LAUNCHING PARALLEL PROCESS USING MPI(MASSAGE PASSING INTERFACE)
    
	CALL MPI_INIT(ERR)
	
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,ERR)
    
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPRC,ERR)
    
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,ERR)
	
	CALL MPI_TYPE_CONTIGUOUS(IE-IS+1,MPI_DOUBLE_PRECISION,AROW,ERR)
	
	CALL MPI_TYPE_COMMIT(AROW,ERR)	! A ROW DATA IN THE CELL-CENTER BASED COMPUTATIONAL DOMAIN
	
	CALL MPI_TYPE_VECTOR(JE-JS+1,1,IE-IS+1,MPI_DOUBLE_PRECISION,CLMN,ERR)
	
	CALL MPI_TYPE_COMMIT(CLMN,ERR)	! A COLUMN DATA IN THE CELL-CENTER BASED COMPUTATIONAL DOMAIN
	
	CALL MPI_TYPE_CONTIGUOUS((NI+1)*(NJ+1),MPI_DOUBLE_PRECISION,APLN,ERR)
	
	CALL MPI_TYPE_COMMIT(APLN,ERR)	! A PLANE DATA IN THE CELL-CENTER BASED COMPUTATIONAL DOMAIN
	
	ALLOCATE(DISP(0:NPRC-1),STAT=ALLOC_ERR)
	
	ALLOCATE(DSIZ(0:NPRC-1),STAT=ALLOC_ERR)

	PI=4.0D0*DATAN(1.0D0)
	
	NIM1 = NI-1
	NJM1 = NJ-1
	NKM1 = NK-1
	
	NUM_GRID=(NI+1)*(NJ+1)*(NK-2)

	! DEFINE THE INTERVAL PARAMETER
	IS =0
	IS1=IS +1
	IE =NI
	IE1=IE -1
	
	JS =0
	JS1=JS +1
	JE =NJ
	JE1=JE -1
	
	KS =0
	KS1=KS +1
	KE =NK
	KE1=KE -1

	INTA = (NK+1)/NPRC
	INTB = MOD(NK+1,NPRC)

	KS =RANK*INTA+MIN(RANK,INTB)
	KE =KS+INTA-1
		IF(INTB.GT.RANK) KE=KE+1
	KS = KS-1	! BOUNDARY BETWEEN PROCESSORS
	KE = KE+1
		IF(RANK.EQ.0)		KS=KS+1
		IF(RANK.EQ.NPRC-1)	KE=KE-1
	KS1=KS +1
	KE1=KE -1

	INXT = RANK +1
	IPRV = RANK -1
	IF(RANK.EQ.0)		IPRV = MPI_PROC_NULL
	IF(RANK.EQ.NPRC-1)	INXT = MPI_PROC_NULL

	SIZ = (IE-IS+1)*(JE-JS+1)*(KE-KS+1)
	CALL MPI_ALLGATHER(SIZ,1,MPI_INTEGER,DSIZ,1,MPI_INTEGER,MPI_COMM_WORLD,ERR)
	
	DISP(0) = 0
	DO N=1,NPRC-1
		DISP(N) = DISP(N-1) + DSIZ(N-1)
	ENDDO
	
    WRITE(*,92) RANK,NPRC,KS,KE,DSIZ(RANK),DISP(RANK)
	WRITE(*,93) RANK,NPRC,IPRV,INXT
92	FORMAT(TR1,'PROCESSOR',I2,TR1,'OF',I2,TR4,' : KMIN =',I4,', KMAX =',I4,', ELMS =',I7,', DIPS =',I7)
93	FORMAT(TR1,'PROCESSOR',I2,TR1,'OF',I2,TR4,' : PREV =',I4,', NEXT =',I4)
	
	END SUBROUTINE DEFINE