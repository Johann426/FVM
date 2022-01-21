	MODULE ComDat_Param                                  

	IMPLICIT NONE
	SAVE

	INTEGER, PARAMETER :: NI	= 301
	INTEGER, PARAMETER :: NJ	= 301
	INTEGER, PARAMETER :: NK	= 64

!	INTEGER, PARAMETER :: IMON	= 100
!	INTEGER, PARAMETER :: JMON	= 100
	INTEGER, PARAMETER :: KMON	=   2
	
	INTEGER, PARAMETER :: NSWPU=1
	INTEGER, PARAMETER :: NSWPP=20
	INTEGER, PARAMETER :: NSWPT=1

    INTEGER :: NITER_W   = 1E4
    INTEGER :: NITER_R   = 50
    INTEGER :: RTIME_MAX = 300
    
    INTEGER :: IRESTART
    
    DOUBLE PRECISION, PARAMETER :: RESIMAX = 1.0D-6
    DOUBLE PRECISION, PARAMETER :: CFL_USE = 0.29D0
    
    DOUBLE PRECISION :: UIN 
!	DOUBLE PRECISION, PARAMETER :: UIN     = 1.D0
	DOUBLE PRECISION, PARAMETER :: TIN     = 0.D0
    DOUBLE PRECISION, PARAMETER :: T_WALL  = 1.D0
    
    DOUBLE PRECISION :: RE = 100.D0
    DOUBLE PRECISION :: PR = 0.7D0
    DOUBLE PRECISION :: PE = 70.D0
    DOUBLE PRECISION :: DTIME = 1.0D-3
!	DOUBLE PRECISION :: VEL_ROT
        
	END MODULE ComDat_Param