!----------------------------------------------------------------------------------
!:sdoc+:
!
!
! GOME2LER_Module
!
! Module containing the load/destruction routines to handel
! the surface reflectivity model data in NetCDF format. 
!
! PUBLIC DATA:
!  GOME2LER:  Data structure containing the GOME2 surface
!              reflectivity data.
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure GOME2_LER_Type.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CSEM initialisation.
!
! CREATION HISTORY:
!       Written by:     Ming  Chen, 28-March-2023
!                       ming.chen@noaa.gov
!       Modified by:    Quanhua Liu, 08-April-2023
!
MODULE GOME2LER_Module
  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE CRTM_Module
  USE netcdf

  ! Disable implicit typing
  IMPLICIT NONE
 
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! The shared data
  PUBLIC :: GOME2LER_ATLAS
  ! Procedures
  PUBLIC :: GOME2LER_Atlas_Load
  PUBLIC :: GOME2LER_Mode
  PUBLIC :: GOME2LER_Mode_pc
  PUBLIC :: GOME2LER_Atlas_CleanUp
  PUBLIC :: GOME2LER_Atlas_IsLoaded
!  PUBLIC :: GOME2LER_Atlas_Setup
  PUBLIC :: GOME2LER_Atlas_Channels
  PUBLIC :: GOME2LER_Atlas_Close
  ! ---------------------------------
  !GOME2LER Atalas data type definition
  ! ---------------------------------
  !:tdoc+:
  TYPE :: GOME2_LER_Type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = 1
    INTEGER(Long) :: Version = 1

    ! dimention parameters
    INTEGER(Long) ::  n_Longitude  = -1
    INTEGER(Long) ::  n_Latitude   = -1
    INTEGER(Long) ::  n_Month      = -1
    INTEGER(Long) ::  n_Wavelength = -1

    ! Derived type components
    REAL(fp), ALLOCATABLE :: longitude(:)      
    REAL(fp), ALLOCATABLE :: latitude(:)       
    REAL(fp), ALLOCATABLE :: wavelength(:)      
    REAL(fp), ALLOCATABLE :: mode_LER(:,:,:)   
    REAL(fp), ALLOCATABLE :: mode_LER_pc(:,:,:,:)  
  END TYPE GOME2_LER_Type
  !:tdoc-:

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 512

  ! Group name of the target set 
  INTEGER,      PARAMETER :: MAX_VAR_DIMS = 5

  ! Dimension names
  CHARACTER(*), PARAMETER :: Longitude_DimName      = 'longitude'
  CHARACTER(*), PARAMETER :: Latitude_DimName       = 'latitude'
  CHARACTER(*), PARAMETER :: Wavelength_DimName     = 'wavelength'
  CHARACTER(*), PARAMETER :: Time_DimName           = 'month'

  ! Variable names. Case sensitive.
  CHARACTER(*), PARAMETER :: Longitude_VarName      = 'longitude'
  CHARACTER(*), PARAMETER :: Latitude_VarName       = 'latitude'
  CHARACTER(*), PARAMETER :: Wavelength_VarName     = 'wavelength'
  CHARACTER(*), PARAMETER :: Time_VarName           = 'month'

 
  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FillValue_ATTNAME = '_FillValue'

  INTEGER,      PARAMETER :: IP_FILLVALUE = -1
  REAL(fp),     PARAMETER :: FP_FILLVALUE = -999.0_fp
  
 
 
  ! --------------------------------------------------
  ! Status Control variable 
  ! -------------------------------------------------
  LOGICAL,SAVE :: Atlas_IsLoaded  = .FALSE.
  ! --------------------------------------------------
  ! The shared microwave water surface emissivity data
  ! --------------------------------------------------
  TYPE(GOME2_LER_Type),  SAVE ::   GOME2LER_ATLAS


  LOGICAL, SAVE :: GOME2LER_atlas_init
  INTEGER, SAVE :: GOME2LER_atlas_month


CONTAINS


  FUNCTION GOME2LER_Atlas_Initialized(imonth) RESULT( Atlas_Status )
    INTEGER :: imonth 
    LOGICAL :: Atlas_Status
    Atlas_Status=(GOME2LER_atlas_init .AND. (imonth == GOME2LER_atlas_month)) 
  END FUNCTION GOME2LER_Atlas_Initialized

  SUBROUTINE GOME2LER_Atlas_Close( )
    INTEGER :: Error_Status
    IF (GOME2LER_atlas_init) THEN
      PRINT*,'Clean GOME2LER Atlas ...'
      Error_Status = GOME2LER_Atlas_CleanUP()
      GOME2LER_atlas_init = .FALSE.
    ENDIF
   
  END SUBROUTINE GOME2LER_Atlas_Close

  SUBROUTINE GOME2LER_Atlas_Channels( Wavelength )                 
    REAL(fp), ALLOCATABLE, INTENT(OUT)  :: Wavelength(:)
    IF(.NOT. GOME2LER_atlas_init) THEN
      PRINT*,"GOME2LER Atlas has not been initilizaed, Call GOME2LER_Atlas_Setup ...."
      RETURN
    ENDIF
    Wavelength = GOME2LER_ATLAS%wavelength
  END SUBROUTINE GOME2LER_Atlas_Channels

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      GOME2LER_Mode
!
! PURPOSE:
!       Function to provide the surface reflectivity values of multiple 
!       channels 
!
!       This function is dedicated to using GOME2LER monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status =  GOME2LER_Mode(           &
!                         Latitude,              & ! input
!                         Longitude,             & ! input
!                         Wavelength,            & ! input
!                         Reflectivity)            ! output
!
!
! INPUTS:
!
!       Latitude:        User's latitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Longitude:       User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       imonth:          User's month index (1-12)
!                        UNITS:      N/A
!                        TYPE:       integer
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Wavelength  :    Wavelength of sensor channels
!                        UNITS:      nm
!                        TYPE:       float
!                        DIMENSION:  Array of n_Channels
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Reflectivity :   Reflectivityvalues of sensor channels
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Array of n_Channels
!                        ATTRIBUTES: INTENT(OUT)
!
!
!:sdoc-:
!----------------------------------------------------------------------------------



  FUNCTION GOME2LER_Mode(            &
    & Latitude,                           & ! input
    & Longitude,                          & ! input
    & imonth,                             & ! input
    & reflectivity)                       & ! output
    RESULT ( Error_Status )

    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    INTEGER,  INTENT(IN)      :: imonth
!    REAL(fp), ALLOCATABLE, INTENT(OUT)     :: reflectivity(:)
    REAL(fp)  :: reflectivity(:)
    ! Function result
    INTEGER :: Error_Status
    ! local
    INTEGER :: iLat, iLon, nchan

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'GOME2LER_Mode'
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
!    IF(.NOT. GOME2LER_Atlas_Initialized(imonth)) THEN
!      PRINT*,'GOME2LER Atlas has not been initilizaed for the month ', imonth
!      Error_Status = GOME2LER_Atlas_Setup(imonth)
!    ENDIF
    nchan = GOME2LER_ATLAS%n_Wavelength
!    ALLOCATE(reflectivity(nchan), STAT = Error_Status)
    
    CALL MAP_GEOLOC_INDEX(Latitude, Longitude, iLat, iLon)
    reflectivity = GOME2LER_ATLAS%mode_LER(iLat,iLon,:)
  
  END FUNCTION GOME2LER_Mode
!
  FUNCTION GOME2LER_Mode_pc(            &
    & Latitude,                           & ! input
    & Longitude,                          & ! input
    & imonth,                             & ! input
    & Zenith_Angle,                       & ! input
    & reflectivity)                       & ! output
    RESULT ( Error_Status )

    REAL(fp), INTENT(IN)      :: Latitude, Longitude, Zenith_Angle
    INTEGER,  INTENT(IN)      :: imonth
!    REAL(fp), ALLOCATABLE, INTENT(OUT)     :: reflectivity(:)
    REAL(fp)  :: reflectivity(:)
    ! Function result
    INTEGER :: Error_Status
    ! local
    INTEGER :: iLat, iLon, nchan

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'GOME2LER_Mode'
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
!    IF(.NOT. GOME2LER_Atlas_Initialized(imonth)) THEN
!      PRINT*,'GOME2LER Atlas has not been initilizaed for the month ', imonth
!      Error_Status = GOME2LER_Atlas_Setup(imonth)
!    ENDIF
    nchan = GOME2LER_ATLAS%n_Wavelength
!    ALLOCATE(reflectivity(nchan), STAT = Error_Status)
    
    CALL MAP_GEOLOC_INDEX(Latitude, Longitude, iLat, iLon)
    reflectivity = GOME2LER_ATLAS%mode_LER(iLat,iLon,:) + GOME2LER_ATLAS%mode_LER_pc(1,iLat,iLon,:) &
       + GOME2LER_ATLAS%mode_LER_pc(2,iLat,iLon,:)*abs(Zenith_Angle) &
       + GOME2LER_ATLAS%mode_LER_pc(3,iLat,iLon,:)*abs(Zenith_Angle)**2
  END FUNCTION GOME2LER_Mode_pc

!
  SUBROUTINE MAP_GEOLOC_INDEX(Lat, Lon, iLat, iLon)

    REAL(fp),        INTENT(IN)  :: Lat,  Lon
    INTEGER(long) ,  INTENT(OUT) :: iLat, iLon
 
    REAL(fp) :: aLat, aLon, rLat, rLon
   
    INTEGER(long)  :: jLat, jLon
    INTEGER(long)  :: nLat, nLon
  
    !nLat = GOME2LER_ATLAS%n_Latitude
    !nLon = GOME2LER_ATLAS%n_Longitude
    nLat = 720 ;  nLon = 1440

    rLat = 180.0/REAL(nLat) ; rLon = 360.0/REAL(nLon)
    
    aLat = lat + 90.0 
    
    aLon = Lon
    IF ( Lon > 180.0 ) aLon = Lon - 360.0 
    alon = aLon + 180.0
    
    jLat = nint(aLat/rLat+0.5) 
    jLon = nint(aLon/rLon+0.5)

    iLat = max(min(nLat, jLat),1)
    iLon = max(min(jLon, nLon),1)

  END SUBROUTINE MAP_GEOLOC_INDEX
!
  FUNCTION GOME2LER_ATLAS_Load( &
    Filename         , &  ! Input
    mm               , &  ! Input
    wavelength       , &  ! Output
    File_Path        , &  ! Optional input
    Version)           &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    INTEGER     , OPTIONAL, INTENT(IN) :: Version             
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'GOME2LER_Atlas_Load'
    ! Local variables
    CHARACTER(ML) :: Message
    CHARACTER(LEN=256) :: GOME2LER_ATLAS_File
    
    INTEGER :: FileID
    INTEGER :: varid  
    INTEGER :: n_Longitude
    INTEGER :: n_Latitude
    INTEGER :: n_Wavelength, mm
    REAL(fp), Allocatable :: wavelength(:)
    REAL, Allocatable :: xx(:,:,:,:), xx_pc(:,:,:,:,:)
    ! Setup 
    Error_Status = SUCCESS
    ! ...Assign the filename to local variable
    GOME2LER_Atlas_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) GOME2LER_Atlas_File = TRIM(ADJUSTL(File_Path))//TRIM(GOME2LER_Atlas_File)

    IF(Atlas_IsLoaded)  THEN
       PRINT*, 'GOME2LER_Atlas already loaded, reloading ...'
       Error_Status = GOME2LER_Atlas_CleanUP()
    ENDIF
    Error_Status = INQ_EmisCoeff_File( &
          GOME2LER_Atlas_File,   &
          n_Longitude,           &
          n_Latitude,            &
          n_Wavelength,          &
          FileID)  
    
    ! Allocate the output structure
    Error_Status  = GOME2LER_Atlas_Create(n_Longitude, n_Latitude, n_Wavelength )
 
    IF ( Error_Status /= SUCCESS ) THEN
      PRINT*,'Error in creating GOME2LER_ATLAS_Type structure.'
      Error_Status = check( nf90_close(FileID) )
      RETURN
    END IF

    ! Read the GOME2LER_Atlas data file
    IF(PRESENT(Version)) GOME2LER_ATLAS%Version = Version

    Error_Status = check(nf90_inq_varid(FileID, Longitude_VarName, varid))      
    Error_Status = check(nf90_get_var(FileID,   varid,  &
               GOME2LER_ATLAS%longitude,  count=(/n_Longitude/)))

    Error_Status = check(nf90_inq_varid(FileID, Latitude_VarName, varid))      
    Error_Status = check(nf90_get_var(FileID,   varid,  &
               GOME2LER_ATLAS%latitude,  count=(/n_Latitude/)))

    Error_Status = check(nf90_inq_varid(FileID, Wavelength_VarName, varid))      
    Error_Status = check(nf90_get_var(FileID,   varid,  &
               GOME2LER_ATLAS%wavelength,  count=(/n_Wavelength/)))
    Allocate( wavelength( n_Wavelength ) )
    wavelength = GOME2LER_ATLAS%wavelength
   
    Error_Status = check(nf90_inq_varid(FileID,'mode_LER', varid))      
!    Error_Status = check(nf90_get_var(FileID,   varid,  &
!               GOME2LER_ATLAS%mode_LER,  count=(/ n_Latitude, n_Longitude, n_Wavelength/)))
    Allocate( xx(n_Latitude, n_Longitude, n_Wavelength, 12) )
    Error_Status = check(nf90_get_var(FileID,   varid,  &
        xx(:,:,:,mm), start=(/1,1,1,mm/), count=(/ n_Latitude, n_Longitude, n_Wavelength, 1/)))
    GOME2LER_ATLAS%mode_LER(:,:,:) = xx(:,:,:,mm)
    Deallocate( xx )

    Error_Status = check(nf90_inq_varid(FileID,'polynomial_coefficients_mode_LER', varid))      
    Allocate( xx_pc(3,n_Latitude, n_Longitude, n_Wavelength, 12) )
    Error_Status = check(nf90_get_var(FileID,   varid,  &
        xx_pc(:,:,:,:,mm), start=(/1,1,1,1,mm/), count=(/ 3, n_Latitude, n_Longitude, n_Wavelength, 1/)))
    GOME2LER_ATLAS%mode_LER_pc(:,:,:,:) = xx_pc(:,:,:,:,mm)
    Deallocate( xx_pc )
   
!integer, parameter :: numLons = 10, numLats = 5, numTimes = 3 
!real, dimension(numLons, numLats, numTimes) :: rhValues
!        !Read the values at the last time by passing an array section 
!        status = nf90_get_var(ncid, rhVarId, rhValues(:, :, 3),  &
!        start = (/ 1, 1, numTimes /),  count = (/ numLats, numLons, 1 /))

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in loading the GOME2LER structure.'
      PRINT*,Message
      Error_Status = GOME2LER_Atlas_CleanUp()    
      Error_Status = check(nf90_close(FileID) )
      RETURN
    END IF
   
    Atlas_IsLoaded = .TRUE. 
 
 
  END FUNCTION GOME2LER_ATLAS_Load
  
  FUNCTION GOME2LER_Atlas_Create( &
       n_Longitude,           &
       n_Latitude,            &
       n_Wavelength)          &
  RESULT (alloc_stat)
    INTEGER,      INTENT(IN) ::  n_Longitude
    INTEGER,      INTENT(IN) ::  n_Latitude
    INTEGER,      INTENT(IN) ::  n_Wavelength
    ! Local variables
    INTEGER :: alloc_stat

    ! Perform the allocation
    print *,' create dim '
    ALLOCATE( GOME2LER_ATLAS%longitude(n_Longitude), STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN
    ALLOCATE( GOME2LER_ATLAS%latitude(n_Latitude), STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN
    ALLOCATE( GOME2LER_ATLAS%wavelength(n_Wavelength), STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN
    ALLOCATE( GOME2LER_ATLAS%mode_LER(n_Latitude, n_Longitude, n_Wavelength), STAT = alloc_stat )
    ALLOCATE( GOME2LER_ATLAS%mode_LER_pc(3,n_Latitude, n_Longitude, n_Wavelength), STAT = alloc_stat )    
    IF ( alloc_stat /= 0 ) THEN
       PRINT*, "Allocation failed ..."
       RETURN
    ENDIF
    ! Initialise
   
    GOME2LER_ATLAS%n_longitude  = n_Longitude    
    GOME2LER_ATLAS%n_latitude   = n_Latitude    
    GOME2LER_ATLAS%n_Wavelength = n_Wavelength    
   
    GOME2LER_ATLAS%mode_LER = 0.0   

    ! Set allocation indicator
    GOME2LER_ATLAS%Is_Allocated = .TRUE.
END FUNCTION GOME2LER_Atlas_Create


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       GOME2LER_Atlas_CleanUp
!
! PURPOSE:
!       Function to deallocate the public data structure GOME2LER_ATLAS containing
!       the GOME2 surface reflectivity  data.
!
! CALLING SEQUENCE:
!       Error_Status = GOME2LER_Atlas_CleanUp( Process_ID = Process_ID )
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION GOME2LER_Atlas_CleanUp( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'GOME2LER_Atlas_CleanUp'
    ! Local variables
    CHARACTER(ML) :: pid_msg

    ! Setup
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF
    
 
    ! Destroy the structure
    IF ( GOME2LER_ATLAS%Is_Allocated ) THEN
       DEALLOCATE(GOME2LER_ATLAS%longitude, &
                  GOME2LER_ATLAS%latitude, &
                  GOME2LER_ATLAS%wavelength, &
                  GOME2LER_ATLAS%mode_LER, GOME2LER_ATLAS%mode_LER_pc, STAT=err_stat)

       GOME2LER_ATLAS%n_Longitude   = -1
       GOME2LER_ATLAS%n_Latitude   = -1
       GOME2LER_ATLAS%n_Wavelength = -1
    ELSE
      RETURN
    END IF
    
    Atlas_IsLoaded = .FALSE.

  END FUNCTION GOME2LER_Atlas_CleanUp

  
!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       GOME2LER_Atlas_IsLoaded
!
! PURPOSE:
!       Function to test if GOME2 surface LER data has
!       been loaded into the public data structure GOME2LER_ATLAS.
!
! CALLING SEQUENCE:
!       status =  GOME2LER_Atlas_IsLoaded()
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION GOME2LER_Atlas_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = Atlas_IsLoaded
  END FUNCTION GOME2LER_Atlas_IsLoaded


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################


  FUNCTION INQ_EmisCoeff_File( &
       FILE_NAME,             &
       n_Longitude,           &
       n_Latitude,            &
       n_Wavelength,          &
       FileID) & 
    RESULT ( Error_status)
  
    CHARACTER(*), INTENT(IN)  :: FILE_NAME
    INTEGER,      INTENT(OUT),   OPTIONAL  :: FileID
    INTEGER,      INTENT(OUT) ::  n_Longitude
    INTEGER,      INTENT(OUT) ::  n_Latitude
    INTEGER,      INTENT(OUT) ::  n_Wavelength
    INTEGER   :: Error_status

    LOGICAL :: Existence
    INTEGER :: ncid,  dimid
 
    Error_Status = SUCCESS
    INQUIRE( FILE = TRIM( FILE_NAME), EXIST = Existence )
    IF ( .NOT. Existence ) THEN
      PRINT*,'File '//TRIM( FILE_NAME )//' not found.'
      Error_Status = FAILURE
      RETURN
    END IF

    !Open the NetCDF file:
    Error_Status = check(nf90_open( TRIM( FILE_NAME ), nf90_nowrite, ncid))
    IF(Error_Status /= SUCCESS) THEN
       Error_Status = check( nf90_close(ncid) ) 
       PRINT*, 'open NC file fail .....' 
       RETURN    
    END IF

    Error_Status = check(nf90_inq_dimid(ncid, Longitude_DimName,  dimid))
    Error_Status = check(nf90_inquire_dimension(ncid, dimid, len=n_Longitude))

    Error_Status = check(nf90_inq_dimid(ncid, Latitude_DimName,  dimid))
    Error_Status = check(nf90_inquire_dimension(ncid, dimid, len=n_Latitude))

    Error_Status = check(nf90_inq_dimid(ncid, Wavelength_DimName,  dimid))
    Error_Status = check(nf90_inquire_dimension(ncid, dimid, len=n_Wavelength))

    IF(PRESENT(FileID)) THEN
      FileID = ncid
      RETURN
    ENDIF
    
    Error_Status = check( nf90_close(ncid) )

  END FUNCTION  INQ_EmisCoeff_File


  
  FUNCTION check(status) RESULT ( Err_Status )

    INTEGER, INTENT ( IN) :: status
    INTEGER :: err_status
      err_status = SUCCESS
      IF (status /= nf90_noerr) THEN
         PRINT *, trim(nf90_strerror(status))
         err_status = FAILURE
      END IF
  END FUNCTION check


END MODULE GOME2LER_Module
