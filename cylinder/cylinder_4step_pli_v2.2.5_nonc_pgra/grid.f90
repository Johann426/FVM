	SUBROUTINE GRID

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	SX =0
	SX1=1
	EX =NI
	EX1=NI-1

	SY =0
	SY1=1
	EY =NJ
	EY1=NJ-1

    NIM1=NI-1
    NJM1=NJ-1
    NKM1=NK-1

    NIM2=NI-2
    NJM2=NJ-2
    NKM2=NK-2   
    
    NUM_GRID=(NI+1)*(NJ+1)*(NK-2)

	FILENAME = '000x000.dat'
	
	WRITE(FILENAME(1:3),'(I3.3)') INT4(NI)
	
	WRITE(FILENAME(5:7),'(I3.3)') INT4(NJ)
!-------------------------------------------
!   READING GRID DATA
    OPEN(10,FILE=FILENAME)
	DO K=2,NKM1
	DO J=1,NJ 
	DO I=1,NI
         READ(10,*) XG(I,J,K),YG(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(10)
	
!-------------------------------------------
! ASSIGN THE CELL NODE LOCATION.
	DO K=2,NKM1
	DO J=1,NJ
    DO I=1,NI
	    XG(I,J,K)=XG(I,J,K) 
	    YG(I,J,K)=YG(I,J,K) 
    ENDDO
	ENDDO
	ENDDO

!-------------------------------------------
! ASSIGN THE CELL CENTER LOCATION.
	DO K=2, NKM1
	DO J=1,NJ-1
	DO I=1,NI-1
		X(I,J,K) = 0.25D0* ( XG(I+1,J+1,K) + XG(I+1,J,K) + XG(I,J+1,K) + XG(I,J,K) )
		Y(I,J,K) = 0.25D0* ( YG(I+1,J+1,K) + YG(I+1,J,K) + YG(I,J+1,K) + YG(I,J,K) )
 	ENDDO
	ENDDO
	
	DO I=1,NI-1
	J=0
		X(I,J,K) = x(i,nj-1,k)
		Y(I,J,K) = y(i,nj-1,k)
	J=NJ
		X(I,J,K) = x(i,1,k)
		Y(I,J,K) = y(i,1,k)
	
	ENDDO

	DO J=1,NJ-1
	I=0
		X(I,J,K) = ( XG(I+1,J+1,K) + XG(I+1,J,K) ) * 0.5D0
		Y(I,J,K) = ( YG(I+1,J+1,K) + YG(I+1,J,K) ) * 0.5D0
!Re-positioning
		X(I,J,K) = X(I,J,K) - (X(I+1,J,K)-X(I,J,K))
		Y(I,J,K) = Y(I,J,K) - (Y(I+1,J,K)-Y(I,J,K))

	I=NI
		X(I,J,K) = ( XG(I,J+1,K)   + XG(I,J,K)	) * 0.5D0
		Y(I,J,K) = ( YG(I,J+1,K)   + YG(I,J,K)	) * 0.5D0
!Re-positioning
		X(I,J,K) = X(I,J,K) + (X(I,J,K)-X(I-1,J,K))
		Y(I,J,K) = Y(I,J,K) + (Y(I,J,K)-Y(I-1,J,K))

	ENDDO
    
!--------------------------------------	
	I=0
	J=0
        X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)

!--------------------------------------
	J=NJ
		X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)

!--------------------------------------
	I=NI
	J=0
		X(I,J,K) =  X(I,NJ-1,K)
		Y(I,J,K) =  Y(I,NJ-1,K)

!--------------------------------------
	J=NJ
		X(I,J,K) =  X(I,1,K)
		Y(I,J,K) =  Y(I,1,K)
    ENDDO

!----------------------------------------------------------------------------
!!   WRITE THE CELL CENTER GRID COORDINATE
!	OPEN(unit=10,file='GRID-CCENT.PLT')
!	WRITE(10,*) 'TITLE="CELL CENTER GRID"'
!	WRITE(10,*) 'variables=X,Y'
!	WRITE(10,*) 'zone t="',RTIME,'",,i=',NI+1,',j=',NJ+1
!	K=KMON
!	DO J=0,NJ
!	DO I=0,NI
!		WRITE(10,*) X(I,J,K),Y(I,J,K)
!	ENDDO
!	ENDDO
!	CLOSE(10)
!----------------------------------------------------------------------------
!!----------------------------------------------------------------------------
!!   WRITE THE CELL NODE GRID COORDINATE
!	OPEN(unit=10,file='GRID-NODE.PLT')
!	WRITE(10,*) 'TITLE="CELL NODE GRID"'
!	WRITE(10,*) 'variables=X,Y'
!	WRITE(10,*) 'zone t="',RTIME,'",,i=',NI,',j=',NJ
!	K=KMON
!	DO J=1,NJ
!	DO I=1,NI
!		WRITE(10,*) XG(I,J,K),YG(I,J,K)
!	ENDDO
!	ENDDO
!   CLOSE(10)
!!----------------------------------------------------------------------------
    END SUBROUTINE GRID