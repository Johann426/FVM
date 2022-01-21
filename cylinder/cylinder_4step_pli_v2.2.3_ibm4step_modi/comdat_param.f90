	MODULE ComDat_Param                                  

	IMPLICIT NONE                                         
	SAVE                                                  

	INTEGER, PARAMETER :: NI	= 80
	INTEGER, PARAMETER :: NJ	= 81
	INTEGER, PARAMETER :: NK	= 3

	INTEGER, PARAMETER :: KMON	=   2

	INTEGER, PARAMETER :: NSWPU=2
	INTEGER, PARAMETER :: NSWPV=2
	INTEGER, PARAMETER :: NSWPP=10
	INTEGER, PARAMETER :: NSWPT=1

	DOUBLE PRECISION, PARAMETER :: RESIMAX = 1.0D-6
    DOUBLE PRECISION, PARAMETER :: CFL_USE = 0.29D0
    
	DOUBLE PRECISION, PARAMETER :: UIN     = 1.0D0
	DOUBLE PRECISION, PARAMETER :: VIN     = 0.0D0
	DOUBLE PRECISION, PARAMETER :: TIN     = 0.0D0
    DOUBLE PRECISION, PARAMETER :: T_WALL  = 1.0D0
    
!    DOUBLE PRECISION, PARAMETER :: RE      = 100.D0
!    DOUBLE PRECISION, PARAMETER :: PR      = 1.D0
!    DOUBLE PRECISION, PARAMETER :: PE      = 100.D0
!    DOUBLE PRECISION, PARAMETER :: DTIME   = 1.0D-3
!    DOUBLE PRECISION, PARAMETER :: VEL_ROT = 0.0D0

    INTEGER :: NITER_W   = 1E6
    INTEGER :: NITER_R   = 2E4
    INTEGER :: IRESTART
    
    DOUBLE PRECISION :: RTIME_MAX = 200
    
    DOUBLE PRECISION :: RE
    DOUBLE PRECISION :: PR
    DOUBLE PRECISION :: PE
    DOUBLE PRECISION :: DTIME
    DOUBLE PRECISION :: VEL_ROT
        
	END MODULE ComDat_Param