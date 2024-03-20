!
! ActiveSensor_Model
!
! Module containing the CRTM radiative transfer solution for active sensors.
!
! CREATION HISTORY:  September 22, 2023
!       Contributed by:  Quanhua Liu
!                        Benjamin Johnson
!                        Isaac Moradi
!                        Yingtao Ma
!
MODULE ActiveSensor_Model

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds              , ONLY: fp
  USE Message_Handler         , ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters         , ONLY: ZERO, ONE, TWO, FOUR, FIVE, TEN, PI, MISSING_REFL
  USE Common_RTSolution       , ONLY: Assign_Common_Input, &
                                      Assign_Common_Output, &
                                      Assign_Common_Input_TL, &
                                      Assign_Common_Output_TL, &
                                      Assign_Common_Input_AD, &
                                      Assign_Common_Output_AD

  USE CRTM_SpcCoeff           , ONLY: SC
  USE CRTM_Atmosphere_Define  , ONLY: CRTM_Atmosphere_type
  USE CRTM_Surface_Define     , ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE CRTM_AtmOptics_Define   , ONLY: CRTM_AtmOptics_type
  USE CRTM_SfcOptics_Define   , ONLY: CRTM_SfcOptics_type
  USE CRTM_RTSolution_Define  , ONLY: CRTM_RTSolution_type
  USE CRTM_Utility
  USE RTV_Define
  USE Liu, ONLY: pVar_type => iVar_type, &
      Ocean_Permittivity => Liu_Ocean_Permittivity
  ! Disable all implicit typing
  IMPLICIT NONE

  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  ! RTSolution structure entities
  ! ...Datatypes
  ! Module procedures
  PUBLIC :: Radar_Solution
  REAL(fp), PARAMETER :: M6_MM6 = 1.0d18
  REAL(fp), PARAMETER :: REFLECTIVITY_THRESHOLD = TINY(REAL(fp))
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*),  PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_RTSolution.f90 60152 2015-08-13 19:19:13Z paul.vandelst@noaa.gov $'

CONTAINS

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!
! NAME:
!       Radar_Solution
!
! PURPOSE:
!       Function to solve radiative transfer equation for active sensors.
!
! CALLING SEQUENCE:
!       Error_Status = Radar_Solution( Atmosphere  , &  ! Input
!                                      Surface     , &  ! Input
!                                      AtmOptics   , &  ! Input
!                                      SfcOptics   , &  ! Input
!                                      GeometryInfo, &  ! Input
!                                      SensorIndex , &  ! Input
!                                      ChannelIndex, &  ! Input
!                                      deltaZ      , &  ! Input
!                                      wavenumber  , &  ! Input!
!                                      RTSolution  , &  ! Output
!                                      RTV           )  ! Internal variable output
!
! INPUT ARGUMENTS:
!       Atmosphere:     Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Surface:        Structure containing the surface state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Surface_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:      Structure containing the surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!                       UNITS:      N/A
!                       TYPE:       CRTM_SfcOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!                       See the SensorIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!      Delta_Z:         Atmospheric layer thickness in meters
!      Wavenumber:      Channel wavenumber
!
! OUTPUT ARGUMENTS:
!       RTSolution:     Structure containing the soluition to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       RTV:            Structure containing internal variables required for
!                       subsequent tangent-linear or adjoint model calls.
!                       The contents of this structure are NOT accessible
!                       outside of the CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
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
!       Note the INTENT on the output RTSolution argument is IN OUT rather than
!       just OUT. This is necessary because the argument is defined upon
!       input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!--------------------------------------------------------------------------------

  FUNCTION Radar_Solution( &
    Atmosphere  , &  ! Input
    Surface     , &  ! Input
    AtmOptics   , &  ! Input
    SfcOptics   , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    deltaZ      , &  ! Input
    wavenumber  , &  ! Input
    RTSolution  , &  ! Output
    RTV         ) &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments

    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics
    TYPE(CRTM_GeometryInfo_type), INTENT(IN OUT) :: GeometryInfo
    REAL(fp), INTENT(IN) :: deltaZ(:), wavenumber
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_RTSolution_type),   INTENT(IN OUT) :: RTSolution
    TYPE(RTV_type),               INTENT(IN OUT) :: RTV
    TYPE(pVar_type):: iVar

    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Radar_Solution'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: i, k, nZ, H2O_idx
    REAL(fp), PARAMETER :: Temperature0 = 273.16_fp, Salinity0 = 0.0_fp
    REAL(fp), PARAMETER :: Light_speed = 29.979246_fp
    REAL(fp) :: Frequency, Wavelength_meter,kw2
    REAL(fp) :: Level_Transmittance(0:size(deltaZ))
    REAL(fp), DIMENSION(size(deltaZ)) :: dZ, R, Reflectivity, Reflectivity_Attenuated
    COMPLEX(fp) :: Permittivity
    Error_Status = SUCCESS

    nZ = Atmosphere%n_Layers
    dZ = deltaZ * 1000.0_fp ! in meters
    dZ = dZ/GeometryInfo%Cosine_Sensor_Zenith
    ! Calculate the geometric heights of the pressure levels in km 
    Frequency = wavenumber * Light_speed  ! in GHz
    Wavelength_meter = 0.01_fp /wavenumber
!
    CALL Ocean_Permittivity( &
    Temperature0 , & ! Input
    Salinity0    , & ! Input
    Frequency   , & ! Input
    Permittivity, & ! Output
    iVar          ) ! Internal variable output   
!
    kw2 = ABS((Permittivity - ONE )/(Permittivity + TWO))**2
!   
    R = (M6_MM6 * Wavelength_meter**FOUR) / (PI**FIVE * Kw2)
    R = R / dZ  ! dZ_m to convert water_content to m/v or cloud water density
    ! Calculate level transmittance from top to layer k
    Level_Transmittance(0) = ONE
    DO k = 1, nZ
       Level_transmittance(k) = Level_transmittance(k-1)* &
       EXP(-TWO * AtmOptics%Optical_Depth(k))
    END DO
!    
    Reflectivity =  R * (AtmOptics%Back_Scattering) ! mm^6 m^-3
    Reflectivity_Attenuated = Reflectivity * Level_transmittance(0:nZ-1)

    ! Convert the unit to dBz
    WHERE (Reflectivity .GT.  REFLECTIVITY_THRESHOLD)
        RTSolution%Reflectivity = TEN * LOG10(Reflectivity) ! [dBZ]
    ELSE WHERE
        RTSolution%Reflectivity = MISSING_REFL
    END WHERE  
!
    WHERE (Reflectivity_Attenuated .GT.  REFLECTIVITY_THRESHOLD)
        RTSolution%Reflectivity_Attenuated = TEN * LOG10(Reflectivity_Attenuated)
    ELSE WHERE
        RTSolution%Reflectivity_Attenuated = MISSING_REFL
    END WHERE 
  END FUNCTION Radar_Solution
!
END MODULE ActiveSensor_Model
!
