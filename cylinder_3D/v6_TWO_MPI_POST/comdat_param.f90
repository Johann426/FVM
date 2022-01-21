	MODULE ComDat_Param                                  

	IMPLICIT NONE
	SAVE

	INTEGER, PARAMETER :: NI	= 80
	INTEGER, PARAMETER :: NJ	= 81
	INTEGER, PARAMETER :: NK	= 9

!	INTEGER, PARAMETER :: IMON	= 100
!	INTEGER, PARAMETER :: JMON	= 100
	INTEGER, PARAMETER :: KMON	=   2
	
	INTEGER, PARAMETER :: NSWPU=1
	INTEGER, PARAMETER :: NSWPP=20
	INTEGER, PARAMETER :: NSWPT=1

    INTEGER :: NITER_W   = 1E6
    INTEGER :: RTIME_MAX = 1E1
    
    DOUBLE PRECISION, PARAMETER :: RESIMAX = 1.0D-6
    DOUBLE PRECISION, PARAMETER :: CFL_USE = 0.29D0
    
!	DOUBLE PRECISION, PARAMETER :: UIN     = 1.0D0
	DOUBLE PRECISION, PARAMETER :: TIN     = 0.0D0
    DOUBLE PRECISION, PARAMETER :: T_WALL  = 1.0D0
    
    INTEGER :: IRESTART       
    DOUBLE PRECISION :: RE
    DOUBLE PRECISION :: PR
    DOUBLE PRECISION :: PE
    DOUBLE PRECISION :: UIN
    DOUBLE PRECISION :: DTIME
!	DOUBLE PRECISION :: VEL_ROT
        
	END MODULE ComDat_Param