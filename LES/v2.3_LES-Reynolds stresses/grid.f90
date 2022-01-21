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

	!!----------------------------------------------------------------
	!FILENAME = '00000x00000.dat'
	!
	!WRITE(FILENAME(1:5),'(I5.5)') INT(NI)
	!
	!WRITE(FILENAME(7:11),'(I5.5)') INT(NJ)
	!
	!DO K=1,NK
	!OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD')
	!DO J=1,NJ
	!DO I=1,NI
	!	READ(10,*) X_BUF,Y_BUF
	!	Z_BUF = 6.283185307 * DBLE(k-1)/DBLE(NK-1)
	!	IF(J.GE.JS1.AND.J.LE.JE .AND. K.GE.KS1.AND.K.LE.KE) THEN
	!		XG(I,J,K) = X_BUF
	!		YG(I,J,K) = Y_BUF
	!		ZG(I,J,K) = Z_BUF
	!	ENDIF
	!ENDDO
	!ENDDO
	!	CLOSE(10)
	!ENDDO
	!!----------------------------------------------------------------
	
	FILENAME = '000x000x000.dat'
	
	WRITE(FILENAME(1:3),'(I3.3)') INT(NI)
	
	WRITE(FILENAME(5:7),'(I3.3)') INT(NJ)
	
	WRITE(FILENAME(9:11),'(I3.3)') INT(NK)

	OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD')
	DO K=1,NK
	DO J=1,NJ
	DO I=1,NI
		READ(10,*) X_BUF,Y_BUF,Z_BUF
		IF(J.GE.JS1.AND.J.LE.JE .AND. K.GE.KS1.AND.K.LE.KE) THEN
			XG(I,J,K) = X_BUF
			YG(I,J,K) = Y_BUF
			ZG(I,J,K) = Z_BUF
		ENDIF
	ENDDO
	ENDDO
	ENDDO
	CLOSE(10)
   
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
		X(I,J,K) = 0.25D0* ( XG(I+1,J+1,K+1) + XG(I+1,J+1,K) + XG(I,J+1,K+1) + XG(I,J+1,K)	)
		Y(I,J,K) = 0.25D0* ( YG(I+1,J+1,K+1) + YG(I+1,J+1,K) + YG(I,J+1,K+1) + YG(I,J+1,K)	)
		Z(I,J,K) = 0.25D0* ( ZG(I+1,J+1,K+1) + ZG(I+1,J+1,K) + ZG(I,J+1,K+1) + ZG(I,J+1,K)	)
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I,J+1,K)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I,J+1,K)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I,J+1,K)-Z(I,J,K))
!		X(I,J,K) = X(I,NJM1,K)
!		Y(I,J,K) = Y(I,NJM1,K)
!		Z(I,J,K) = Z(I,NJM1,K)

	J=JE
		X(I,J,K) = 0.25D0* ( XG(I+1,J,K+1) + XG(I+1,J,K) + XG(I,J,K+1) + XG(I,J,K)	)
		Y(I,J,K) = 0.25D0* ( YG(I+1,J,K+1) + YG(I+1,J,K) + YG(I,J,K+1) + YG(I,J,K)	)
		Z(I,J,K) = 0.25D0* ( ZG(I+1,J,K+1) + ZG(I+1,J,K) + ZG(I,J,K+1) + ZG(I,J,K)	)
	!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I,J-1,K)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I,J-1,K)-Y(I,J,K))
		Z(I,J,K) = Z(I,J,K) - (Z(I,J-1,K)-Z(I,J,K))
!		X(I,J,K) = X(I,1,K)
!		Y(I,J,K) = Y(I,1,K)
!	    Z(I,J,K) = Z(I,1,K)
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

!	!   ALIGNING CELL CENTER LOCATION ON THE BOUNDARY LINES
!	DO I=IS1,IE1
!	J=JS
!	K=KS
!        X(I,J,K) = X(I,NJ-1,K)
!		Y(I,J,K) = Y(I,NJ-1,K)
!        Z(I,J,K) = Z(I,NJ-1,K)
!	K=KE
!        X(I,J,K) = X(I,NJ-1,K)
!		Y(I,J,K) = Y(I,NJ-1,K)
!        Z(I,J,K) = Z(I,NJ-1,K)
!    J=JE
!        X(I,J,K) = X(I,1,K)
!		Y(I,J,K) = Y(I,1,K)
!		Z(I,J,K) = Z(I,1,K)
!    K=KS
!        X(I,J,K) = X(I,1,K)
!		Y(I,J,K) = Y(I,1,K)
!		Z(I,J,K) = Z(I,1,K)
!	ENDDO
!	
!    DO J=JS1,JE1
!    I=IS
!    K=KS
!		X(I,J,K) = 0.5D0* ( XG(I+1,J+1,K+1) + XG(I+1,J,K+1) )
!		Y(I,J,K) = 0.5D0* ( YG(I+1,J+1,K+1) + YG(I+1,J,K+1) )
!		Z(I,J,K) = 0.5D0* ( ZG(I+1,J+1,K+1) + ZG(I+1,J,K+1) )
!	!Re-positioning
!		X(I,J,K) = X(I,J,K) - (X(I+1,J,K+1)-X(I,J,K))
!		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K+1)-Y(I,J,K))
!		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K+1)-Z(I,J,K))
!    K=KE
!        X(I,J,K) = 0.5D0* ( XG(I+1,J+1,K) + XG(I+1,J,K) )
!		Y(I,J,K) = 0.5D0* ( YG(I+1,J+1,K) + YG(I+1,J,K) )
!		Z(I,J,K) = 0.5D0* ( ZG(I+1,J+1,K) + ZG(I+1,J,K) )
!	!Re-positioning
!		X(I,J,K) = X(I,J,K) - (X(I+1,J,K-1)-X(I,J,K))
!		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K-1)-Y(I,J,K))
!		Z(I,J,K) = Z(I,J,K) - (Z(I+1,J,K-1)-Z(I,J,K))
!    I=IE
!        X(I,J,K) = 0.5D0* ( XG(I,J+1,K) + XG(I,J,K) )
!		Y(I,J,K) = 0.5D0* ( YG(I,J+1,K) + YG(I,J,K) )
!		Z(I,J,K) = 0.5D0* ( ZG(I,J+1,K) + ZG(I,J,K) )
!	!Re-positioning
!		X(I,J,K) = X(I,J,K) - (X(I-1,J,K-1)-X(I,J,K))
!		Y(I,J,K) = Y(I,J,K) - (Y(I-1,J,K-1)-Y(I,J,K))
!		Z(I,J,K) = Z(I,J,K) - (Z(I-1,J,K-1)-Z(I,J,K))
!    K=KS
!		X(I,J,K) = 0.5D0* ( XG(I,J+1,K+1) + XG(I,J,K+1) )
!		Y(I,J,K) = 0.5D0* ( YG(I,J+1,K+1) + YG(I,J,K+1) )
!		Z(I,J,K) = 0.5D0* ( ZG(I,J+1,K+1) + ZG(I,J,K+1) )
!	!Re-positioning
!		X(I,J,K) = X(I,J,K) - (X(I-1,J,K+1)-X(I,J,K))
!		Y(I,J,K) = Y(I,J,K) - (Y(I-1,J,K+1)-Y(I,J,K))
!		Z(I,J,K) = Z(I,J,K) - (Z(I-1,J,K+1)-Z(I,J,K))
!    ENDDO
!    
!	DO K=KS1,KE1
!	I=IS
!	J=JS
!        X(I,J,K) =  X(I,NJ-1,K)
!		Y(I,J,K) =  Y(I,NJ-1,K)
!        Z(I,J,K) =  Z(I,NJ-1,K)
!	J=JE
!		X(I,J,K) =  X(I,1,K)
!		Y(I,J,K) =  Y(I,1,K)
!		Z(I,J,K) =  Z(I,1,K)
!	I=IE
!		X(I,J,K) =  X(I,1,K)
!		Y(I,J,K) =  Y(I,1,K)
!		Z(I,J,K) =  Z(I,1,K)
!	J=JS
!		X(I,J,K) =  X(I,NJ-1,K)
!		Y(I,J,K) =  Y(I,NJ-1,K)
!		Z(I,J,K) =  Z(I,NJ-1,K)
!    ENDDO
!
!	!   ALIGNING CELL CENTER LOCATION AT THE BOUNDARY POINTS
!    I=IS
!    J=JS
!    K=KS
!		X(I,J,K) =  X(I,NJ-1,K)
!		Y(I,J,K) =  Y(I,NJ-1,K)
!		Z(I,J,K) =  Z(I,NJ-1,K)
!    I=IE
!		X(I,J,K) =  X(I,NJ-1,K)
!		Y(I,J,K) =  Y(I,NJ-1,K)
!		Z(I,J,K) =  Z(I,NJ-1,K)
!    J=JE
!        X(I,J,K) =  X(I,1,K)
!		Y(I,J,K) =  Y(I,1,K)
!		Z(I,J,K) =  Z(I,1,K)
!    I=IS
!        X(I,J,K) =  X(I,1,K)
!		Y(I,J,K) =  Y(I,1,K)
!		Z(I,J,K) =  Z(I,1,K)
!    K=KE
!        X(I,J,K) =  X(I,1,K)
!		Y(I,J,K) =  Y(I,1,K)
!		Z(I,J,K) =  Z(I,1,K)
!    J=JS
!		X(I,J,K) =  X(I,NJ-1,K)
!		Y(I,J,K) =  Y(I,NJ-1,K)
!		Z(I,J,K) =  Z(I,NJ-1,K)
!    I=IE
!		X(I,J,K) =  X(I,NJ-1,K)
!		Y(I,J,K) =  Y(I,NJ-1,K)
!		Z(I,J,K) =  Z(I,NJ-1,K)
!    J=JE
!        X(I,J,K) =  X(I,1,K)
!		Y(I,J,K) =  Y(I,1,K)
!		Z(I,J,K) =  Z(I,1,K)
	
	CALL MPI_X(X)
!	CALL MPI_X(X)
	
	CALL MPI_X(Y)
!	CALL MPI_X(Y)
	
	CALL MPI_X(Z)
!	CALL MPI_X(Z)

	!	Z-DIRECTIONAL PERIODIC CELL-CENTER POSITION
	IF(KS.EQ.0) THEN
	K=KS
	
	X(IS:IE,JS:JE,K) = X(IS:IE,JS:JE,K+1)
	Y(IS:IE,JS:JE,K) = Y(IS:IE,JS:JE,K+1)
	DO J=JS1,JE1
	DO I=IS1,IE1
		Z(I,J,K) = 0.25D0* ( ZG(I,J,K+1)+ZG(I+1,J,K+1)+ZG(I+1,J+1,K+1)+ZG(I,J+1,K+1) )
	!Re-positioning
		Z(I,J,K) = Z(I,J,K) - (Z(I,J,K+1)-Z(I,J,K))
    ENDDO
    ENDDO
	DO J=JS1,JE1
	I=IS
		Z(I,J,K) = Z(I+1,J,K)
	I=IE
		Z(I,J,K) = Z(I-1,J,K)
	ENDDO
	DO I=IS1,IE1
	J=JS
		Z(I,J,K) = Z(I,J+1,K)
	J=JE
		Z(I,J,K) = Z(I,J-1,K)
	ENDDO
	I=IS
	J=JS
		Z(I,J,K) = Z(I+1,J+1,K)
	I=IE
		Z(I,J,K) = Z(I-1,J+1,K)
	J=JE
		Z(I,J,K) = Z(I-1,J-1,K)
	I=IS
		Z(I,J,K) = Z(I+1,J-1,K)
		
    ENDIF
	
	IF(KE.EQ.NK) THEN
	K=KE
	
	X(IS:IE,JS:JE,K) = X(IS:IE,JS:JE,K-1)
	Y(IS:IE,JS:JE,K) = Y(IS:IE,JS:JE,K-1)
	DO J=JS1,JE1
	DO I=IS1,IE1
		Z(I,J,K) = 0.25D0* ( ZG(I,J,K)+ZG(I+1,J,K)+ZG(I+1,J+1,K)+ZG(I,J+1,K) )
	!Re-positioning
		Z(I,J,K) = Z(I,J,K) + (Z(I,J,K)-Z(I,J,K-1))
	ENDDO
    ENDDO
	DO J=JS1,JE1
	I=IS
		Z(I,J,K) = Z(I+1,J,K)
	I=IE
		Z(I,J,K) = Z(I-1,J,K)
	ENDDO
	DO I=IS1,IE1
	J=JS
		Z(I,J,K) = Z(I,J+1,K)
	J=JE
		Z(I,J,K) = Z(I,J-1,K)
	ENDDO
	I=IS
	J=JS
		Z(I,J,K) = Z(I+1,J+1,K)
	I=IE
		Z(I,J,K) = Z(I-1,J+1,K)
	J=JE
		Z(I,J,K) = Z(I-1,J-1,K)
	I=IS
		Z(I,J,K) = Z(I+1,J-1,K)
    
    ENDIF
    
    	
!	CONTOUR = 'GRID_OO.PLT'
!	
!	WRITE(CONTOUR(6:7),'(I2.2)') INT4(RANK)
!	
!	OPEN(UNIT=100+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
!	WRITE(100+RANK,*) 'TITLE="','RE',RE,'PR',PR,'"'
!	WRITE(100+RANK,*) 'variables=x,y,z'
!	WRITE(100+RANK,33) IE-IS+1,JE-JS+1,KE-KS+1
!33	FORMAT(TR1'ZONE T="WHOLE", I=',I6,', J=',I6,', K=',I6)
!	DO K=KS,KE
!	DO J=JS,JE
!	DO I=IS,IE
!		WRITE(100+RANK,*) X(I,J,K),Y(I,J,K),Z(I,J,K)
!	ENDDO
!	ENDDO
!	ENDDO
!	CLOSE(100+RANK)
	
    END SUBROUTINE GRID