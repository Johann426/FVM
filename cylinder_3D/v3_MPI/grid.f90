	SUBROUTINE GRID

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	IF(RANK.EQ.0) THEN
!-------------------------------------------------
!   READING GRID DATA
	OPEN(10,FILE='80X81.dat')
!	DO K=1,NK
	K=1
	DO J=1,NJ
	DO I=1,NI
		!READ(10,*) XG(I,J,K),YG(I,J,K),ZG(I,J,K)
		READ(10,*) XG(I,J,K),YG(I,J,K)
	ENDDO
	ENDDO
!	ENDDO
	CLOSE(10)
	
!-------------------------------------------
! ASSIGN THE CELL NODE LOCATION.
	DO K=1,NK
	DO J=1,NJ
    DO I=1,NI
	    XG(I,J,K)=XG(I,J,1)
	    YG(I,J,K)=YG(I,J,1)
		ZG(I,J,K)=DBLE(K-1)/DBLE(NK-1)*2.D0*PI
    ENDDO
	ENDDO
	ENDDO
	
	REF_A = 1.D0*(ZG(IS1,JS1,KE)-ZG(IS1,JS1,KS1))
	
!-------------------------------------------
! ASSIGN THE CELL CENTER LOCATION.
	DO K=1,NKM1
	DO J=1,NJM1
	DO I=1,NIM1
		X(I,J,K)=0.125D0* (XG(I+1,J+1,K)+XG(I+1,J,K)+XG(I,J+1,K)+XG(I,J,K)+XG(I+1,J+1,K+1)+XG(I+1,J,K+1)+XG(I,J+1,K+1)+XG(I,J,K+1))
		Y(I,J,K)=0.125D0* (YG(I+1,J+1,K)+YG(I+1,J,K)+YG(I,J+1,K)+YG(I,J,K)+YG(I+1,J+1,K+1)+YG(I+1,J,K+1)+YG(I,J+1,K+1)+YG(I,J,K+1))
		Z(I,J,K)=0.125D0* (ZG(I+1,J+1,K)+ZG(I+1,J,K)+ZG(I,J+1,K)+ZG(I,J,K)+ZG(I+1,J+1,K+1)+ZG(I+1,J,K+1)+ZG(I,J+1,K+1)+ZG(I,J,K+1))
 	ENDDO
	ENDDO
	ENDDO
!-----------------------------------------------------------------------------
!   ALIGNING CELL CENTER LOCATION ON THE BOUNDARY FACES
    DO K=1,NK-1
	DO J=1,NJ-1
	I=0
		X(I,J,K) = 0.25D0* ( XG(I+1,J+1,K)+XG(I+1,J,K)+XG(I+1,J+1,K+1)+XG(I+1,J,K+1) )
		Y(I,J,K) = 0.25D0* ( YG(I+1,J+1,K)+YG(I+1,J,K)+YG(I+1,J+1,K+1)+YG(I+1,J,K+1) )
		Z(I,J,K) = 0.25D0* ( ZG(I+1,J+1,K)+ZG(I+1,J,K)+ZG(I+1,J+1,K+1)+ZG(I+1,J,K+1) )
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K)-Z(I,J,K))

	I=NI
		X(I,J,K) = 0.25D0* ( XG(I,J+1,K)+XG(I,J,K)+XG(I,J+1,K+1)+XG(I,J,K+1) )
		Y(I,J,K) = 0.25D0* ( YG(I,J+1,K)+YG(I,J,K)+YG(I,J+1,K+1)+YG(I,J,K+1) )
		Z(I,J,K) = 0.25D0* ( ZG(I,J+1,K)+ZG(I,J,K)+ZG(I,J+1,K+1)+ZG(I,J,K+1) )
!Re-positioning
		X(I,J,K) = X(I,J,K) + (X(I,J,K)-X(I-1,J,K))
		Y(I,J,K) = Y(I,J,K) + (Y(I,J,K)-Y(I-1,J,K))
		Z(I,J,K) = Z(I,J,K) + (Z(I,J,K)-Z(I-1,J,K))
	ENDDO
    ENDDO
    
    DO K=1,NK-1
	DO I=1,NI-1
	J=0
		X(I,J,K) = X(I,NJM1,K)
		Y(I,J,K) = Y(I,NJM1,K)
		Z(I,J,K) = Z(I,NJM1,K)
	J=NJ
		X(I,J,K) = X(I,1,K)
		Y(I,J,K) = Y(I,1,K)
	    Z(I,J,K) = Z(I,1,K)
	ENDDO
    ENDDO
    
	DO J=1,NJ-1
    DO I=1,NI-1
	K=0
		X(I,J,K) = 0.25D0* ( XG(I,J,K+1)+XG(I+1,J,K+1)+XG(I+1,J+1,K+1)+XG(I,J+1,K+1) )
		Y(I,J,K) = 0.25D0* ( YG(I,J,K+1)+YG(I+1,J,K+1)+YG(I+1,J+1,K+1)+YG(I,J+1,K+1) )
		Z(I,J,K) = 0.25D0* ( ZG(I,J,K+1)+ZG(I+1,J,K+1)+ZG(I+1,J+1,K+1)+ZG(I,J+1,K+1) )
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I,J,K+1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I,J,K+1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I,J,K+1)-Z(I,J,K))

	K=NK
		X(I,J,K) = 0.25D0* ( XG(I,J,K)+XG(I+1,J,K)+XG(I+1,J+1,K)+XG(I,J+1,K) )
		Y(I,J,K) = 0.25D0* ( YG(I,J,K)+YG(I+1,J,K)+YG(I+1,J+1,K)+YG(I,J+1,K) )
		Z(I,J,K) = 0.25D0* ( ZG(I,J,K)+ZG(I+1,J,K)+ZG(I+1,J+1,K)+ZG(I,J+1,K) )
!Re-positioning
		X(I,J,K) = X(I,J,K) + (X(I,J,K)-X(I,J,K-1))
		Y(I,J,K) = Y(I,J,K) + (Y(I,J,K)-Y(I,J,K-1))
		Z(I,J,K) = Z(I,J,K) + (Z(I,J,K)-Z(I,J,K-1))
	ENDDO
    ENDDO
!-----------------------------------------------------------------------------
!   ALIGNING CELL CENTER LOCATION ON THE BOUNDARY LINES
	DO I=1,NI-1
	J=0
	K=0
        X(I,J,K) = X(I,NJ-1,K)
		Y(I,J,K) = Y(I,NJ-1,K)
        Z(I,J,K) = Z(I,NJ-1,K)
	K=NK
        X(I,J,K) = X(I,NJ-1,K)
		Y(I,J,K) = Y(I,NJ-1,K)
        Z(I,J,K) = Z(I,NJ-1,K)
    J=NJ
        X(I,J,K) = X(I,1,K)
		Y(I,J,K) = Y(I,1,K)
		Z(I,J,K) = Z(I,1,K)
    K=0
        X(I,J,K) = X(I,1,K)
		Y(I,J,K) = Y(I,1,K)
		Z(I,J,K) = Z(I,1,K)
	ENDDO
	
    DO J=1,NJ-1
    I=0
    K=0
		X(I,J,K) = 0.5D0* ( XG(I+1,J+1,K+1) + XG(I+1,J,K+1) )
		Y(I,J,K) = 0.5D0* ( YG(I+1,J+1,K+1) + YG(I+1,J,K+1) )
		Z(I,J,K) = 0.5D0* ( ZG(I+1,J+1,K+1) + ZG(I+1,J,K+1) )
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K+1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K+1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K+1)-Z(I,J,K))
    K=NK
        X(I,J,K) = 0.5D0* ( XG(I+1,J+1,K) + XG(I+1,J,K) )
		Y(I,J,K) = 0.5D0* ( YG(I+1,J+1,K) + YG(I+1,J,K) )
		Z(I,J,K) = 0.5D0* ( ZG(I+1,J+1,K) + ZG(I+1,J,K) )
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K-1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K-1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K-1)-Z(I,J,K))
    I=NI
        X(I,J,K) = 0.5D0* ( XG(I,J+1,K) + XG(I,J,K) )
		Y(I,J,K) = 0.5D0* ( YG(I,J+1,K) + YG(I,J,K) )
		Z(I,J,K) = 0.5D0* ( ZG(I,J+1,K) + ZG(I,J,K) )
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I-1,J,K-1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I-1,J,K-1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I-1,J,K-1)-Z(I,J,K))
    K=0
		X(I,J,K) = 0.5D0* ( XG(I,J+1,K+1) + XG(I,J,K+1) )
		Y(I,J,K) = 0.5D0* ( YG(I,J+1,K+1) + YG(I,J,K+1) )
		Z(I,J,K) = 0.5D0* ( ZG(I,J+1,K+1) + ZG(I,J,K+1) )
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I-1,J,K+1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I-1,J,K+1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I-1,J,K+1)-Z(I,J,K))
    ENDDO
    
	DO K=1,NK-1
	I=0
	J=0
        X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
        Z(I,J,K) =  Z(I,NJ-1,K)
	J=NJ
		X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
	I=NI
		X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
	J=0
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    ENDDO
!-----------------------------------------------------------------------------
!   ALIGNING CELL CENTER LOCATION AT THE BOUNDARY POINTS
    I=0
    J=0
    K=0
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    I=NI
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    J=NJ
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
    I=0
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
    K=NK
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
    J=0
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    I=NI
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    J=NJ
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)

	ENDIF

!----------------------------------------------------------------------------	
	DO K=1,NK
	
	DO J=1,NJ
	
	CALL MPI_BCAST(XG(1,J,K),NI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
	
	CALL MPI_BCAST(YG(1,J,K),NI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
	
	CALL MPI_BCAST(ZG(1,J,K),NI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
	
	ENDDO
	
	ENDDO
	
	DO K=0,NK
	
	DO J=0,NJ
	
	CALL MPI_BCAST(X(0,J,K),NI+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
	
	CALL MPI_BCAST(Y(0,J,K),NI+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
	
	CALL MPI_BCAST(Z(0,J,K),NI+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
	
	ENDDO
	
	ENDDO
	

!----------------------------------------------------------------------------	

!----------------------------------------------------------------------------
!   WRITE THE CELL CENTER GRID COORDINATE
	IF(RANK.EQ.NPRC-1) THEN
	OPEN(unit=10,file='GRID-CCENT.PLT')
	WRITE(10,*) 'TITLE="CELL CENTER GRID"'
	WRITE(10,*) 'variables=X,Y,Z'
	WRITE(10,*) 'zone t="',RTIME,'",,i=',NI+1,',j=',NJ+1,',K=',NK+1
	DO K=0,NK
	DO J=0,NJ
	DO I=0,NI
		WRITE(10,*) X(I,J,K),Y(I,J,K),Z(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(10)
	ENDIF
!----------------------------------------------------------------------------
    END SUBROUTINE GRID