!
! CRTM_ActiveSensor_Module
!
! This module calculates attenuated reflectivity for active sensors within the CRTM. 
! The module is based on the CRTM_AOD_Module. This module designed for radar reflectivity,
! lidar reflectivity and scatterometer for surface. The first application is for radar 
! attenuated reflectivity in dbZ here.
!
!
! CREATION HISTORY:  September 22, 2023
!       Contributed by:  Quanhua Liu
!                        Benjamin Johnson
!                        Isaac Moradi
!                        Yingtao Ma
!
! Inputs are the same as the CRTM model: sensor ID, geometryInfo (e.g. zenith and azimuth angles),
!   Atmosphere and Surface data types.
!  The output will be 
!     RTSolution(channel_idx, profile_idx)%Reflectivity_Attenuated(1:n_Layers).
!     Missing value or zero value is set to -9999.0 for the dbZ values.
!
MODULE CRTM_ActiveSensor_Module


  ! ------------
  ! Module usage
  ! ------------
  USE Type_Kinds,                 ONLY: fp, LLong
  USE ODPS_CoordinateMapping,     ONLY: Geopotential_Height
  USE Message_Handler,            ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,            ONLY: SET,NOT_SET,ZERO,ONE, RT_VMOM, &
                                        MAX_N_LAYERS        , &
                                        MAX_N_PHASE_ELEMENTS, &
                                        MAX_N_LEGENDRE_TERMS, &
                                        MAX_N_STOKES        , &
                                        MAX_N_ANGLES        , &
                                        MAX_N_AZIMUTH_FOURIER, &
                                        MAX_SOURCE_ZENITH_ANGLE, &
                                        MAX_N_STREAMS, &
                                        AIRCRAFT_PRESSURE_THRESHOLD, &
                                        MIN_COVERAGE_THRESHOLD, &
                                        SCATTERING_ALBEDO_THRESHOLD
  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      CRTM_Atmosphere_IsValid
  USE CRTM_Surface_Define,        ONLY: CRTM_Surface_type, &
                                        CRTM_Surface_IsValid
  USE CRTM_ChannelInfo_Define,  ONLY: CRTM_ChannelInfo_type, &
                                      CRTM_ChannelInfo_n_Channels
  USE CRTM_Options_Define,      ONLY: CRTM_Options_type
  USE CRTM_AtmOptics_Define,    ONLY: CRTM_AtmOptics_type      , &
                                      CRTM_AtmOptics_Associated, &
                                      CRTM_AtmOptics_Create    , &
                                      CRTM_AtmOptics_Destroy   , &
                                      CRTM_AtmOptics_Zero
  USE CRTM_GeometryInfo_Define,   ONLY: CRTM_GeometryInfo_type, &
                                        CRTM_GeometryInfo_SetValue, &
                                        CRTM_GeometryInfo_GetValue
  USE CRTM_GeometryInfo,          ONLY: CRTM_GeometryInfo_Compute
  USE CRTM_Geometry_Define,       ONLY: CRTM_Geometry_type, &
                                        CRTM_Geometry_IsValid
  USE CRTM_Predictor_Define,      ONLY: CRTM_Predictor_type      , &
                                        CRTM_Predictor_Associated, &
                                        CRTM_Predictor_Destroy   , &
                                        CRTM_Predictor_Create
  USE CRTM_Predictor,             ONLY: CRTM_PVar_type => iVar_type, &
                                        CRTM_Compute_Predictors
  USE CRTM_AtmAbsorption,         ONLY: CRTM_AAvar_type => iVar_type, &
                                        CRTM_Compute_AtmAbsorption
  USE CRTM_AerosolScatter,      ONLY: CRTM_Compute_AerosolScatter, &
                                      CRTM_Compute_AerosolScatter_TL, &
                                      CRTM_Compute_AerosolScatter_AD
  USE CRTM_SfcOptics_Define,      ONLY: CRTM_SfcOptics_type      , &
                                        CRTM_SfcOptics_Associated, &
                                        CRTM_SfcOptics_Create    , &
                                        CRTM_SfcOptics_Destroy
  USE CRTM_SfcOptics,             ONLY: CRTM_Compute_SurfaceT
  USE CRTM_RTSolution_Define,     ONLY: CRTM_RTSolution_type   , &
                                        CRTM_RTSolution_Destroy, &
                                        CRTM_RTSolution_Zero,    &
                                        CRTM_RTSolution_Inspect
  USE CRTM_AncillaryInput_Define, ONLY: CRTM_AncillaryInput_type
  USE CRTM_AerosolCoeff,        ONLY: CRTM_AerosolCoeff_IsLoaded
  USE CRTM_CloudCoeff,            ONLY: CRTM_CloudCoeff_IsLoaded
  USE CRTM_TauCoeff,              ONLY: TC
  USE CRTM_SpcCoeff,              ONLY: SC, &
                                        SpcCoeff_IsVisibleSensor, &
                                        SpcCoeff_IsMicrowaveSensor, &
                                        SpcCoeff_IsInfraredSensor, &
                                        SpcCoeff_IsUltravioletSensor
  USE CRTM_CloudScatter,          ONLY: CRTM_Compute_CloudScatter

  ! ...CloudScatter
  USE CSvar_Define, ONLY: CSvar_type, &
                          CSvar_Associated, &
                          CSvar_Destroy   , &
                          CSvar_Create  
  ! Internal variable definition modules
  ! ...AerosolScatter
  USE ASvar_Define, ONLY: ASvar_type, &
                          ASvar_Associated, &
                          ASvar_Destroy   , &
                          ASvar_Create

  ! ...AtmOptics
  USE AOvar_Define, ONLY: AOvar_type, &
                          AOvar_Associated, &
                          AOvar_Destroy   , &
                          AOvar_Create
  ! ...Radiative transfer
  USE RTV_Define,   ONLY: RTV_type, &
                          RTV_Associated, &
                          RTV_Destroy   , &
                          RTV_Create
  ! Active Sensors
  USE ActiveSensor_Model, ONLY: Radar_Solution

  
  ! -----------------------
  ! Disable implicit typing
  ! -----------------------
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Public procedures
  PUBLIC :: CRTM_ActiveSensor

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_ActiveSensor
!
! PURPOSE:
!       Function that calculates layer total optical depth profile at nadir.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_ActiveSensor( Atmosphere       , &
!                                         ChannelInfo      , &
!                                         RTSolution       , &
!                                         Options = Options  )
!
! INPUTS:
!       Atmosphere:     Structure containing the Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Rank-1 (n_Profiles)
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelInfo:    Structure returned from the CRTM_Init() function
!                       that contains the satellite/sensor channel index
!                       information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_ChannelInfo_type
!                       DIMENSION:  Rank-1 (n_Sensors)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:     Structure containing the layer aerosol optical
!                       profile for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Options:        Options structure containing the optional arguments
!                       for the CRTM.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Options_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       - Many of the components of the Options optional input structure
!         are not used in this function. Consult the CRTM User Guide for
!         which Options components are usable for AOD calculations.
!
!:sdoc-:
!--------------------------------------------------------------------------------
  FUNCTION CRTM_ActiveSensor( &
    Atmosphere , &  ! Input, M
    Surface    , &  ! Input, M
    Geometry   , &  ! Input, M
    ChannelInfo, &  ! Input, n_Sensors
    RTSolution , &  ! Output, L x M
    Options    ) &  ! Optional input, M
  RESULT( Error_Status )

    ! Arguments
    TYPE(CRTM_Surface_type),           INTENT(IN)     :: Surface(:)        ! M
    TYPE(CRTM_Geometry_type),          INTENT(IN)     :: Geometry(:)       ! M
    TYPE(CRTM_SfcOptics_type)                         :: SfcOptics
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere(:)     ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)    ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)   ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)        ! M
    TYPE(CRTM_AncillaryInput_type) :: AncillaryInput
    TYPE(CRTM_GeometryInfo_type) :: GeometryInfo
    TYPE(CRTM_Predictor_type)    :: Predictor
    TYPE(CSvar_type)      :: CSvar  ! CloudScatter
    TYPE(AOvar_type)      :: AOvar  ! AtmOptics
    TYPE(RTV_type)        :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AOD'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: ln, H2O_idx
    REAL(fp), Allocatable :: Height(:), deltaZ(:)
    REAL(fp) :: transmittance1, transmittance2
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics
    TYPE(CRTM_PVar_type)  :: PVar   ! Predictor
    TYPE(ASVar_type) :: ASvar
    TYPE(CRTM_AAvar_type) :: AAvar  ! AtmAbsorption
!
    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS

    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check the number of channels
    IF ( SIZE(RTSolution,DIM=1) < n_Channels ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure array too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(RTSolution,DIM=2) /= n_Profiles ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for RTSolution argument.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF


    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles

    CALL CRTM_SfcOptics_Create( SfcOptics  , MAX_N_ANGLES, MAX_N_STOKES )

    Allocate( Height(0:Atmosphere(m)%n_Layers), deltaZ(Atmosphere(m)%n_Layers) )
      H2O_idx = 1  ! default
    CALL Geopotential_Height(Atmosphere(m)%Level_Pressure      , & ! Input
                            Atmosphere(m)%Temperature         , & ! Input
                            Atmosphere(m)%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
!    write(6,'(6f12.5)') Height
    deltaZ(:) = Height(0:Atmosphere(m)%n_Layers-1)-Height(1:Atmosphere(m)%n_Layers)

    CALL CRTM_RTSolution_Zero(RTSolution(:,m))
      ! ...Compute derived geometry
    CALL CRTM_GeometryInfo_SetValue( GeometryInfo, Geometry=Geometry(m) )
    CALL CRTM_GeometryInfo_Compute( GeometryInfo )
    CALL CRTM_Compute_SurfaceT( Surface(m), SfcOptics )
      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF
!
      CALL RTV_Create( RTV, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, Atmosphere(m)%n_Layers )
      
      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
        ! Check whether to skip this profile
        IF ( Options(m)%Skip_Profile ) CYCLE Profile_Loop
      END IF
      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF
      ! Check the RTSolution layer dimension
      IF ( ANY(RTSolution(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! Allocate AtmOptics based on Atmosphere dimension
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structure for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set the scattering switch
      AtmOptics%Include_Scattering = Options(m)%Include_Scattering
      ! ...Allocate the atmospheric optics internal structure
      CALL AOvar_Create( AOvar, Atmosphere(m)%n_Layers )

      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms = 16

      ! Allocate the scattering internal variables if necessary
      ! ...Cloud
      IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
        CALL CSvar_Create( CSvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Clouds    )
      END IF
      ! ...Aerosol
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Aerosols  )
      END IF
!
      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0
      Sensor_Loop: DO n = 1, n_Sensors

        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index

          CALL CRTM_Predictor_Create( &
                   Predictor, &
                   Atmosphere(m)%n_Layers,  &
                   SensorIndex    )
          IF ( .NOT. CRTM_Predictor_Associated(Predictor) ) THEN
            Error_Status=FAILURE
            WRITE( Message,'("Error allocating predictor structure for profile #",i0, &
                   &" and ",a," sensor.")' ) m, SC(SensorIndex)%Sensor_Id
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          END IF

          ! ...And now fill them
          CALL CRTM_Compute_Predictors( SensorIndex   , &  ! Input
                                        Atmosphere(m) , &  ! Input
                                        GeometryInfo  , &  ! Input
                                        AncillaryInput, &  ! Input
                                        Predictor , &  ! Output
                                        PVar        )  ! Internal variable output
        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
            ! Compute the gas absorption
          CALL CRTM_Compute_AtmAbsorption( SensorIndex   , &  ! Input
                                           ChannelIndex  , &  ! Input
                                           AncillaryInput, &  ! Input
                                           Predictor , &  ! Input
                                           AtmOptics , &  ! Output
                                           AAvar       )  ! Internal variable output

          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status = CRTM_Compute_AerosolScatter( Atmosphere(m), &  ! Input
                                                        SensorIndex  , &  ! Input
                                                        ChannelIndex , &  ! Input
                                                        AtmOptics    , &  ! In/Output
                                                        ASVar          )  ! Internal variable output
            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF
          END IF
          ! Compute the cloud particle absorption/scattering properties
          IF( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status = CRTM_Compute_CloudScatter( Atmosphere(m) , &  ! Input
                                                      SensorIndex , &  ! Input
                                                      ChannelIndex, &  ! Input
                                                      AtmOptics   , &  ! Output
                                                      CSvar         )  ! Internal variable output

            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing CloudScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            END IF
          END IF
!       
          IF ( Options(m)%Use_Emissivity ) THEN
            ! ...Cloudy/all-sky case
            SfcOptics%Compute = .FALSE.
            SfcOptics%Emissivity(1,1)       = Options(m)%Emissivity(ln)
            SfcOptics%Reflectivity(1,1,1,1) = ONE - Options(m)%Emissivity(ln)
            IF ( Options(m)%Use_Direct_Reflectivity ) THEN
              SfcOptics%Direct_Reflectivity(1,1) = Options(m)%Direct_Reflectivity(ln)
            ELSE
              SfcOptics%Direct_Reflectivity(1,1) = SfcOptics%Reflectivity(1,1,1,1)
            END IF
          END IF
          ! Save the nadir optical depth
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth
!
          Error_Status = Radar_Solution( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
        END DO Channel_Loop

      END DO Sensor_Loop
      ! Deallocate local sensor independent data structures
      CALL CRTM_Predictor_Destroy( Predictor )
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_SfcOptics_Destroy( SfcOptics )
      ! ...Internal variables
      CALL AOvar_Destroy( AOvar )
      CALL CSvar_Destroy( CSvar )
      CALL ASvar_Destroy( ASvar )
      CALL RTV_Destroy( RTV )
      DEALLOCATE( Height, deltaZ )
    END DO Profile_Loop

  END FUNCTION CRTM_ActiveSensor

END MODULE CRTM_ActiveSensor_Module
!
