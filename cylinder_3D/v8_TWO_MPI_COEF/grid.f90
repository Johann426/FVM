	SUBROUTINE GRID

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	ALLOCATE(XG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(YG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(ZG(IS1:IE,JS1:JE,KS1:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(X(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(Y(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	ALLOCATE(Z(IS:IE,JS:JE,KS:KE),STAT=ALLOC_ERR)
	
	FILENAME = '000x000x000_OO.dat'
	
	WRITE(FILENAME(1:3),'(I3.3)') INT4(NI)
	
	WRITE(FILENAME(5:7),'(I3.3)') INT4(NJ)
	
	WRITE(FILENAME(9:11),'(I3.3)') INT4(NK)
	
	WRITE(FILENAME(13:14),'(I2.2)') INT4(RANK)
	
	!   READING GRID DATA
	OPEN(10+RANK,FILE=FILENAME)
	DO K=KS1,KE
	DO J=JS1,JE
	DO I=IS1,IE
		READ(10+RANK,*) XG(I,J,K),YG(I,J,K),ZG(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(10+RANK)
	
	REF_A = 1.D0*(ZG(IS1,JS1,KE)-ZG(IS1,JS1,KS1))	! FOR CIRCULAR CYLINDER
			
	!	ASSIGN THE CELL CENTER LOCATION.
	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
		X(I,J,K)=0.125D0* (XG(I+1,J+1,K)+XG(I+1,J,K)+XG(I,J+1,K)+XG(I,J,K)+XG(I+1,J+1,K+1)+XG(I+1,J,K+1)+XG(I,J+1,K+1)+XG(I,J,K+1))
		Y(I,J,K)=0.125D0* (YG(I+1,J+1,K)+YG(I+1,J,K)+YG(I,J+1,K)+YG(I,J,K)+YG(I+1,J+1,K+1)+YG(I+1,J,K+1)+YG(I,J+1,K+1)+YG(I,J,K+1))
		Z(I,J,K)=0.125D0* (ZG(I+1,J+1,K)+ZG(I+1,J,K)+ZG(I,J+1,K)+ZG(I,J,K)+ZG(I+1,J+1,K+1)+ZG(I+1,J,K+1)+ZG(I,J+1,K+1)+ZG(I,J,K+1))
 	ENDDO
	ENDDO
	ENDDO

	!   ALIGNING CELL CENTER LOCATION ON THE BOUNDARY FACES
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS
		X(I,J,K) = 0.25D0* ( XG(I+1,J+1,K)+XG(I+1,J,K)+XG(I+1,J+1,K+1)+XG(I+1,J,K+1) )
		Y(I,J,K) = 0.25D0* ( YG(I+1,J+1,K)+YG(I+1,J,K)+YG(I+1,J+1,K+1)+YG(I+1,J,K+1) )
		Z(I,J,K) = 0.25D0* ( ZG(I+1,J+1,K)+ZG(I+1,J,K)+ZG(I+1,J+1,K+1)+ZG(I+1,J,K+1) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K)-Z(I,J,K))

	I=IE
		X(I,J,K) = 0.25D0* ( XG(I,J+1,K)+XG(I,J,K)+XG(I,J+1,K+1)+XG(I,J,K+1) )
		Y(I,J,K) = 0.25D0* ( YG(I,J+1,K)+YG(I,J,K)+YG(I,J+1,K+1)+YG(I,J,K+1) )
		Z(I,J,K) = 0.25D0* ( ZG(I,J+1,K)+ZG(I,J,K)+ZG(I,J+1,K+1)+ZG(I,J,K+1) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) + (X(I,J,K)-X(I-1,J,K))
		Y(I,J,K) = Y(I,J,K) + (Y(I,J,K)-Y(I-1,J,K))
		Z(I,J,K) = Z(I,J,K) + (Z(I,J,K)-Z(I-1,J,K))
	ENDDO
    ENDDO
    
	DO K=KS1,KE1
	DO I=IS1,IE1
	J=JS
		X(I,J,K) = X(I,NJM1,K)
		Y(I,J,K) = Y(I,NJM1,K)
		Z(I,J,K) = Z(I,NJM1,K)
	J=JE
		X(I,J,K) = X(I,1,K)
		Y(I,J,K) = Y(I,1,K)
	    Z(I,J,K) = Z(I,1,K)
	ENDDO
    ENDDO
    
	DO J=JS1,JE1
	DO I=IS1,IE1
	K=KS
		X(I,J,K) = 0.25D0* ( XG(I,J,K+1)+XG(I+1,J,K+1)+XG(I+1,J+1,K+1)+XG(I,J+1,K+1) )
		Y(I,J,K) = 0.25D0* ( YG(I,J,K+1)+YG(I+1,J,K+1)+YG(I+1,J+1,K+1)+YG(I,J+1,K+1) )
		Z(I,J,K) = 0.25D0* ( ZG(I,J,K+1)+ZG(I+1,J,K+1)+ZG(I+1,J+1,K+1)+ZG(I,J+1,K+1) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I,J,K+1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I,J,K+1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I,J,K+1)-Z(I,J,K))

	K=KE
		X(I,J,K) = 0.25D0* ( XG(I,J,K)+XG(I+1,J,K)+XG(I+1,J+1,K)+XG(I,J+1,K) )
		Y(I,J,K) = 0.25D0* ( YG(I,J,K)+YG(I+1,J,K)+YG(I+1,J+1,K)+YG(I,J+1,K) )
		Z(I,J,K) = 0.25D0* ( ZG(I,J,K)+ZG(I+1,J,K)+ZG(I+1,J+1,K)+ZG(I,J+1,K) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) + (X(I,J,K)-X(I,J,K-1))
		Y(I,J,K) = Y(I,J,K) + (Y(I,J,K)-Y(I,J,K-1))
		Z(I,J,K) = Z(I,J,K) + (Z(I,J,K)-Z(I,J,K-1))
	ENDDO
    ENDDO

	!   ALIGNING CELL CENTER LOCATION ON THE BOUNDARY LINES
	DO I=IS1,IE1
	J=JS
	K=KS
        X(I,J,K) = X(I,NJ-1,K)
		Y(I,J,K) = Y(I,NJ-1,K)
        Z(I,J,K) = Z(I,NJ-1,K)
	K=KE
        X(I,J,K) = X(I,NJ-1,K)
		Y(I,J,K) = Y(I,NJ-1,K)
        Z(I,J,K) = Z(I,NJ-1,K)
    J=JE
        X(I,J,K) = X(I,1,K)
		Y(I,J,K) = Y(I,1,K)
		Z(I,J,K) = Z(I,1,K)
    K=KS
        X(I,J,K) = X(I,1,K)
		Y(I,J,K) = Y(I,1,K)
		Z(I,J,K) = Z(I,1,K)
	ENDDO
	
    DO J=JS1,JE1
    I=IS
    K=KS
		X(I,J,K) = 0.5D0* ( XG(I+1,J+1,K+1) + XG(I+1,J,K+1) )
		Y(I,J,K) = 0.5D0* ( YG(I+1,J+1,K+1) + YG(I+1,J,K+1) )
		Z(I,J,K) = 0.5D0* ( ZG(I+1,J+1,K+1) + ZG(I+1,J,K+1) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K+1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K+1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K+1)-Z(I,J,K))
    K=KE
        X(I,J,K) = 0.5D0* ( XG(I+1,J+1,K) + XG(I+1,J,K) )
		Y(I,J,K) = 0.5D0* ( YG(I+1,J+1,K) + YG(I+1,J,K) )
		Z(I,J,K) = 0.5D0* ( ZG(I+1,J+1,K) + ZG(I+1,J,K) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K-1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K-1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K-1)-Z(I,J,K))
    I=IE
        X(I,J,K) = 0.5D0* ( XG(I,J+1,K) + XG(I,J,K) )
		Y(I,J,K) = 0.5D0* ( YG(I,J+1,K) + YG(I,J,K) )
		Z(I,J,K) = 0.5D0* ( ZG(I,J+1,K) + ZG(I,J,K) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I-1,J,K-1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I-1,J,K-1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I-1,J,K-1)-Z(I,J,K))
    K=KS
		X(I,J,K) = 0.5D0* ( XG(I,J+1,K+1) + XG(I,J,K+1) )
		Y(I,J,K) = 0.5D0* ( YG(I,J+1,K+1) + YG(I,J,K+1) )
		Z(I,J,K) = 0.5D0* ( ZG(I,J+1,K+1) + ZG(I,J,K+1) )
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I-1,J,K+1)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I-1,J,K+1)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I-1,J,K+1)-Z(I,J,K))
    ENDDO
    
	DO K=KS1,KE1
	I=IS
	J=JS
        X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
        Z(I,J,K) =  Z(I,NJ-1,K)
	J=JE
		X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
	I=IE
		X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
	J=JS
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    ENDDO

	!   ALIGNING CELL CENTER LOCATION AT THE BOUNDARY POINTS
    I=IS
    J=JS
    K=KS
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    I=IE
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    J=JE
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
    I=IS
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
    K=KE
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)
    J=JS
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    I=IE
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)
		Z(I,J,K) =  Z(I,NJ-1,K)
    J=JE
        X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
		Z(I,J,K) =  Z(I,1,K)

	CALL MPI_X(X)
	
	CALL MPI_X(Y)
	
	CALL MPI_X(Z)
	
    END SUBROUTINE GRID