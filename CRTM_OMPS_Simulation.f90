!
PROGRAM CRTM_OMPS_Simulation
! ------------------------------------------------------------------------------
!  Simulate OMPS radiance.
!
! This is a main program to drive the CRTM run at the channels of the OMPS NP and NM sensors aboard different mission satellites. 
! A composite data type is defined and implemented to control the CRTM IOs consistently thorughout the main program and
! other modules.  Users may use this composite data type to customize the IO paraemters of their individual runs.
!
! The surface reflectance can be either from surface reflectance databse or from CRTM derived reflectance by using OMPS measurements for
!    selected channels (typically 2 or 3 channels). The all channel surface reflectances are inter/extra-polated
!    from the surface reflectance at selected channels above.
! 
! Notes: opt%Derive_Surface_Refl=.true. and only RTV(nt)%mth_Azi==0 needed for deriving parameter by
!   call CRTM_SurfRef (in ADA_Module.f90) inside CRTM_Forward_Module needed for deriving surface reflectance.
!
! RTSolution(ln,m)%Surface_Planck_Radiance = alpha (Eq.(10a) doi: 10.1109/JSTARS.2022.3149767)
! RTSolution(ln,m)%Up_Radiance   = delta  (Eq.(10c)
! RTSolution(ln,m)%Down_Radiance = beta   (Eq.(10b)
!
!  References
!  Liu, Q., Yan, B., Garrett, K., Ma, Y., Liang, X., Huang, J., Wang, W., Cao, Y. (2022). 
!    Deriving Surface Reflectance from Visible/Near Infrared and Ultraviolet Satellite 
!    Observations Through the Community Radiative Transfer Model, in IEEE Journal of Selected 
!    Topics in Applied Earth Observations and Remote Sensing, doi: 10.1109/JSTARS.2022.3149767.
!
!  Liu, Q.; Cao, C.; Grassotti, C.; Liang, X.; Chen, Y. Experimental OMPS Radiance Assimilation 
!    through One-Dimensional Variational Analysis for Total Column Ozone in the Atmosphere. 
!    Remote Sens. 2021, 13, 3418. https://doi.org/10.3390/ rs13173418
!
!  HISTORY:
!       Written by:     Quanhua Liu, April-01-2022
!                       Quanhua.Liu@noaa.gov
!	04/05/2022	Jingfeng Huang	Integrated into the OMPS V-CRTM Interface for Multiple Granule and Global Processing
!
!	01/12/2024	Ming Chen	Implement the composite data type to control IOs, code refacoring and refinement
! ------------------------------------------------------------------------------
! 
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Message_Handler
  USE CRTM_Module
  USE Timing_Utility
  USE CRTM_SpcCoeff,           ONLY: SC
  USE Read_atms_inputs_module, ONLY: crtm_sat_inputs_type, Read_crtm_inputs, &
       Derive_surface_reflectance

  USE GOME2LER_Module, ONLY: GOME2LER_Atlas_Load, &
                                    GOME2LER_Atlas_Channels,  &
                                    GOME2LER_Mode, GOME2LER_Mode_pc, &
                                    GOME2LER_Atlas_Close
  USE OMPS_Namelist_Module, ONLY: OMPS_IO_NAMELIST, &
                                  ReadNameList, get_select_ch
  ! Disable all implicit typing
  IMPLICIT NONE
 
  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*),  PARAMETER :: PROGRAM_NAME = 'CRTM_OMPS_Simulation'
  CHARACTER(*),  PARAMETER :: PROGRAM_VERSION_ID = &
    '$Id: CRTM_OMPS_Simulation.f90  $'

  ! Variables
  ! ---------  
  CHARACTER(256 ) :: Namelist_File_Name 
  CHARACTER(256 ) :: crtm_output_file
  CHARACTER(256 ) :: Binary_File_Name
  CHARACTER(255)  :: CRTM_Coeff_Path, Sensor_Id
  CHARACTER(256)  :: err_msg, alloc_msg

  INTEGER, PARAMETER :: n_sensors = 1
  INTEGER :: n_aerosols, n_clouds, n_absorbers, n_layers
  INTEGER :: n_Channels, n_Stokes, n_Profiles

  TYPE(CRTM_Atmosphere_type) , ALLOCATABLE :: atm(:)
  TYPE(CRTM_Surface_type)    , ALLOCATABLE :: sfc(:)
  TYPE(CRTM_Geometry_type)   , ALLOCATABLE :: geo(:)  
  TYPE(CRTM_Options_type)    , ALLOCATABLE :: opt(:)
  TYPE(CRTM_RTSolution_type) , ALLOCATABLE :: RTSolution(:,:)

  TYPE(CRTM_ChannelInfo_type) :: chinfo(N_SENSORS)
  TYPE(crtm_sat_inputs_type)  :: crtm_sat_inputs
  TYPE(Timing_type) :: Timing

  REAL(fp), DIMENSION(:),   ALLOCATABLE :: CRTM_wavelength, CRTM_S_wavelength
  REAL(fp), DIMENSION(:,:), ALLOCATABLE :: mea_n_rad, CRTM_n_Rad
  REAL(fp), DIMENSION(:,:), ALLOCATABLE :: Surf_refl, CSurf_refl
  REAL(fp), DIMENSION(:,:), ALLOCATABLE :: m_Wavelength, LUT_wL

  REAL(fp), ALLOCATABLE :: wavelength(:), reflectivity(:)
  REAL(fp), ALLOCATABLE :: OMPS_mean_nR(:), CRTM_mean_nR(:)

!
!  Statistics
!  bias, std are absolute values; rel_bias, rel_std are relative values in % to mean

  REAL(fp), DIMENSION(:), ALLOCATABLE :: bias, std
  REAL(fp), DIMENSION(:), ALLOCATABLE :: rel_bias, rel_std
  REAL(fp), DIMENSION(:), ALLOCATABLE :: Obs_mean
  INTEGER,  DIMENSION(:), ALLOCATABLE :: QC_prof

  REAL(fp), ALLOCATABLE :: inter_lat(:), inter_lon(:)
  INTEGER,  ALLOCATABLE :: select_ch(:), save_FOV(:)
 
  INTEGER  :: n_Surf_ch, n_Surf_ch2 = 2
  REAL(fp) :: sur_ref_ch(2)
  REAL(fp) :: slope

 
  !There are 3 options fro surface reflectance under Lambertian Equivalent Reflectivity (LER) assumption.
  INTEGER, PARAMETER :: Refl_db = 1, Refl_edr = 2, Refl_ana = 3, Refl_crtm = 4
  INTEGER :: Refl_alg = 3

  INTEGER :: i, ic, n, k, m, m1, m2
  INTEGER :: FileID, Out_ID
  INTEGER :: IO_Status, err_stat, alloc_stat
  INTEGER :: Month
  INTEGER :: nsample 
 
  TYPE(OMPS_IO_NAMELIST) :: namelist
 
  ! *********  namelist Input File*************
  PRINT *,' read input file name '
  READ( *, '(a)' ) Namelist_File_Name
  namelist = ReadNameList(Namelist_File_Name)
  Sensor_Id = namelist%Sensor_Id
  CRTM_Coeff_Path = namelist%CRTM_Coeff_Path
  CRTM_Output_File = namelist%crtm_output_file
  CALL get_select_ch(namelist, n_Surf_ch, select_ch)
  namelist%Use_allFOV_Nadir = .TRUE.
  
  ! =========================================================
  !  STEP 1: CRTM initialization
  ! =========================================================
  n_Stokes = 1
  IF( n_Stokes > 1 ) THEN
     err_stat = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hence the (/../)
         chinfo  , &  ! Output
         CloudCoeff_File='Cloud_V3.bin', &
         AerosolCoeff_File='Aerosol_V3.bin', & 
         File_Path=trim(CRTM_Coeff_Path)) 
  ELSE
     err_stat = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hence the (/../)
         chinfo  , &  ! Output
         File_Path=trim(CRTM_Coeff_Path))    
  END IF

  n_Channels = sum(CRTM_ChannelInfo_n_Channels(chinfo(:)) )
  Allocate( OMPS_mean_nR(n_Channels), CRTM_mean_nR(n_Channels) )
  
  
  ! =========================================================
  !  STEP 2: Read inputs with file names in Namelist_File_Name 
  !          # 'CRTM_OMPS_input.txt'
  ! =========================================================
  crtm_sat_inputs%n_chan = n_Channels
  n_channels = crtm_sat_inputs%n_chan

  print *,' call Read_crtm_inputs ',namelist%nprof,namelist%nspec, namelist%max_FOV

  err_stat =  Read_crtm_inputs(namelist, Atm, Sfc, Geo, crtm_sat_inputs, &
    n_channels, Month,inter_lat,inter_lon,m_Wavelength)
 
  n_profiles = crtm_sat_inputs%n_profiles

  ALLOCATE(bias(n_channels), std(n_channels), rel_bias(n_channels), &
      rel_std(n_channels), Obs_mean(n_channels) )
  ALLOCATE( save_FOV(crtm_sat_inputs%n_profiles) )
  ALLOCATE( RTSolution(n_channels, crtm_sat_inputs%n_profiles), CRTM_n_Rad(n_channels,n_profiles))
  ALLOCATE( Opt(crtm_sat_inputs%n_profiles) )
  CALL CRTM_Options_Create( Opt, n_channels )
  DO k = 1, crtm_sat_inputs%n_profiles
      Opt(k)%Include_Scattering = .TRUE.
      Opt(k)%Use_n_Streams = .TRUE.
      Opt(k)%n_Streams = 2  ! two streams in downward and two strems in upwelling
      Opt(k)%n_Stokes = n_Stokes
      Opt(k)%RT_Algorithm_Id = RT_ADA
      Opt(k)%Use_Emissivity = .TRUE.
      Opt(k)%Use_Direct_Reflectivity = .TRUE.
      Opt(k)%Emissivity(:) = 0.3_fp !0.1_fp

      Opt(k)%Direct_Reflectivity(:) = (ONE-Opt(k)%Emissivity(:)) 
      IF( n_Stokes > 1 ) Opt(k)%RT_Algorithm_Id = RT_VMOM
      IF( Refl_alg == Refl_ana ) THEN
        Opt(k)%Derive_Surface_Refl = .TRUE.
      ELSE 
        Opt(k)%Derive_Surface_Refl = .false.
      END IF
  END DO

  ALLOCATE( CRTM_wavelength(n_channels), Surf_refl(n_channels,n_profiles), &
       mea_n_rad(n_Surf_ch, n_profiles), CSurf_refl(n_channels,n_profiles), &
       CRTM_S_wavelength(n_Surf_ch) )

  CRTM_wavelength(:) = 1.E07/SC(1)%wavenumber(:)
  CRTM_S_wavelength(:) = CRTM_wavelength(select_ch(:) )

  ! Note the AllWavenumber is slightly different from the
  ! m_Wavelength read from the standalone data file.
  ! m_Wavelength was used for the nadir-only coeff set
  ALLOCATE( LUT_wL( n_Channels, namelist%max_FOV) )
  IF (ALLOCATED(SC(1)%AllWavenumber)) THEN
     LUT_wL(:,:) = 1.E07/SC(1)%AllWavenumber(:,:)
  ELSE
     LUT_wL(:,:) = m_Wavelength 
  ENDIF

  ! Use the nadir FOV of the allFOV
  IF(namelist%Use_allFOV_Nadir) THEN
     opt(:)%nFOV = NINT(namelist%max_FOV/2.0) 
     save_FOV(:) = NINT(namelist%max_FOV/2.0)
  ELSE  
     opt(:)%nFOV = geo(:)%iFOV 
     save_FOV(:) = geo(:)%iFOV
  ENDIF

  m1 = 1 ; m2 = n_Profiles  ; n_Profiles = m2
  Allocate( QC_prof(n_profiles) )
  QC_prof(:) = 0
  
  ! =========================================================
  !  STEP 3: selectd algorithm for surface reflectance 
  ! =========================================================
  IF( Refl_alg == Refl_db ) THEN
     ! err_stat = GOME2LER_Atlas_Load('GOME-2_MetOp-AB_MSC_025x025_surface_LER_v3.2_Mar.nc', Month,wavelength,TRIM('./'))
     ! err_stat = GOME2LER_Atlas_Load('GOME-2_MetOp-ABC_PMD_025x025_surface_LER_v4.0.nc', Month,wavelength,TRIM('./'))
     err_stat = GOME2LER_Atlas_Load(trim(CRTM_Coeff_Path)//'/GOME-2_MetOp-ABC_MSC_025x025_surface_LER_v4.0.nc', Month,wavelength) !,TRIM('./'))
     Allocate( reflectivity(size(wavelength)) )
     QC_prof(:) = 0
      
    inter_lon(1:n_profiles) = inter_lon(1:n_profiles) -360.0

    DO k = 1, n_profiles   
      IF( inter_lat(k) < -100.0_fp .or. inter_lon(k) < -400.0 ) THEN
        QC_prof(k) = 10
        print *,' incorrect lat/lon ',k, inter_lat(k), inter_lon(k)
      ELSE
        err_stat = GOME2LER_Mode(inter_lat(k), inter_lon(k), month,  reflectivity)    
        !err_stat = GOME2LER_Mode_pc(inter_lat(k), inter_lon(k), month,  geo(k)%Sensor_Zenith_Angle, reflectivity)

        CALL myspline(size(wavelength(:)),wavelength(:),reflectivity(:),n_channels, &
          LUT_wL(:,geo(k)%iFOV),Surf_refl(:,k) )
        IF( maxval(Surf_refl(:,k)) > 1.0_fp .or. minval(Surf_refl(:,k)) < 0.0_fp) THEN
          print *,' Refl out of range ',k, maxval(Surf_refl(:,k)), minval(Surf_refl(:,k))
          !STOP
        END IF
      END If
    END DO

  ELSE IF( Refl_alg == Refl_edr ) THEN
     DO k = 1, n_profiles
       sur_ref_ch(1) =  crtm_sat_inputs%edrR331(k)
       sur_ref_ch(2) =  crtm_sat_inputs%edrR360(k)
       slope = (sur_ref_ch(2)-sur_ref_ch(1))/(360.0-331.0)
       Surf_refl(:,k) = sur_ref_ch(1) + slope*(CRTM_wavelength(:)-331.0)

       IF( minval(Surf_refl(:,k)) <= 0.001_fp .or. maxval(Surf_refl(:,k)) >= 0.999_fp ) QC_prof(k) = 1
       IF( Surf_refl(select_ch(1),k) <= 0.001_fp .or. Surf_refl(select_ch(1),k) >= 0.999_fp) QC_prof(k) = 2
       IF( Surf_refl(select_ch(2),k) <= 0.001_fp .or. Surf_refl(select_ch(2),k) >= 0.999_fp) QC_prof(k) = 3
       IF( Surf_refl(select_ch(1),k) <= 0.001_fp .and. Surf_refl(select_ch(2),k) <= 0.001_fp ) QC_prof(k) = 4
       IF( Surf_refl(select_ch(1),k) >= 0.999_fp .and. Surf_refl(select_ch(2),k) >= 0.999_fp ) QC_prof(k) = 4
     END DO
     print *,' edr complete ', n_profiles,n_channels
  
  ELSE IF( Refl_alg == Refl_ana ) THEN
     ! Determining surface reflectance assuming Lambertian reflectance
     print *,' Refl_ana ',n_profiles,n_channels
     chinfo(1)%Process_Channel(:) = .false.
     chinfo(1)%Process_Channel(select_ch(:)) = .true.
     !  interpolated measured normalized radiance to NR at FOV dependent CRTM wavelengths.
     DO k = 1, n_profiles
        CALL myspline(n_channels,m_Wavelength(:,opt(k)%nFOV),crtm_sat_inputs%Normalized_Radiance(:,k), &
             n_channels,LUT_wL(:,opt(k)%nFOV),CRTM_n_Rad(:,k) )
        !   print *,k,opt(k)%nFOV
        !   WRITE(6,'(10f9.2)') m_Wavelength(:,opt(k)%nFOV)
        !   WRITE(6,'(6E14.5)') crtm_sat_inputs%Normalized_Radiance(:,k)
        !   WRITE(6,'(6E14.5)') CRTM_n_Rad(:,k)
      
     END DO
     mea_n_rad(:,:) = CRTM_n_Rad(select_ch(:),:)
     ! Calculate needed parameters for calculating surface reflectance from selected channels
     err_stat = CRTM_Forward(atm(m1:m2)        , &
                             sfc(m1:m2)        , &
                             geo(m1:m2)        , &
                             chinfo            , &
                             RTSolution(1:n_channels,m1:m2) , &
                             Options = opt(m1:m2)  )

     ! In alighment with the allFOV coeff 
     CALL Derive_surface_reflectance(n_Surf_ch, save_FOV, n_channels, &
        n_Profiles, LUT_wL, RTSolution, mea_n_rad, select_ch, Surf_refl, QC_prof)
     ! In alighment with the nadir-only coeff, old version
     !CALL Derive_surface_reflectance(n_Surf_ch2, CRTM_S_wavelength(2:3), n_channels, &
     !  n_Profiles, CRTM_wavelength, RTSolution, mea_n_rad, select_ch(2:3), Surf_refl)
  END IF !  -----  completed surface reflectance calculations

  DO k = 1, n_profiles
      Opt(k)%Emissivity(:) = 1.0_fp - Surf_refl(:,k)  
      Opt(k)%Direct_Reflectivity(:) = (ONE-Opt(k)%Emissivity(:))
      Opt(k)%Derive_Surface_Refl = .false.   
  END DO

  ! =========================================================
  !  STEP 4: Simulate OMPS radiances
  ! =========================================================
  chinfo(1)%Process_Channel(:) = .true.
  print *,'Running CRTM_Forward .... '
  CALL CRTM_RTSolution_Zero( RTSolution )
  err_stat = CRTM_Forward( atm(m1:m2)          , &
                           sfc(m1:m2)          , &
                           geo(m1:m2)          , &
                           chinfo              , &
                           RTSolution(:,m1:m2) , &
                           Options = opt(m1:m2)  )
  print *,' forward ',err_stat


  ! =========================================================
  !  STEP 4.1: Calculate OMPS normalized CRTM_n_Rad
  !  Note: RTSolution(:,k)%SolarIrradiance depends on Sun and Earth distance 
  ! =========================================================
  Obs_mean(:) = ZERO
  nsample = 0
  DO k = 1, n_profiles
    !IF(QC_prof(k) == 0 .and. crtm_sat_inputs%edrCloudFraction(k) < 0.0001_fp) THEN
      CRTM_n_Rad(:,k) = RTSolution(:,k)%Radiance/RTSolution(:,k)%SolarIrradiance
      CALL myspline(n_channels,LUT_wL(:,opt(k)%nFOV),CRTM_n_Rad(:,k),n_channels, &
         m_Wavelength(:,opt(k)%nFOV),RTSolution(:,k)%Radiance )
       ! m_Wavelength(:,opt(k)%nFOV)-0.04,RTSolution(:,k)%Radiance )
      CRTM_n_Rad(:,k) = RTSolution(:,k)%Radiance
    
      Obs_mean(:) = Obs_mean(:) + crtm_sat_inputs%Normalized_Radiance(:,k)
      nsample = nsample + 1
    !END IF
  END DO
  Obs_mean(:) = Obs_mean(:)/float(nsample)

  ! =========================================================
  !  STEP 5: Write output files
  ! =========================================================
  CLOSE(10)
  OPEN( 10, FILE = crtm_output_file, &
            STATUS   = 'REPLACE', &                                      
            FORM     = 'UNFORMATTED', &                                  
            IOSTAT   = IO_Status )                                       
  IF ( IO_Status /= 0 ) THEN                                                 
     WRITE( *, * ) 'Error opening ', crtm_output_file, ' in Output' 
     STOP
  END IF 
  ! WRITE(10) n_profiles, n_channels
  WRITE(10) crtm_sat_inputs%Normalized_Radiance
  WRITE(10) CRTM_n_Rad
  WRITE(10) Surf_refl
  ! Added QC (based on derived surface ref whether in physics meaningful range [0.001, 0.999]
  ! QC_prof, 0 (normal), > 0 (abnormal
  WRITE(10) QC_prof(:)
  CLOSE(10)

  ! =========================================================
  !  STEP 6: Statistic analysis and output
  ! =========================================================
  nsample = 0
  bias(:) = ZERO
  std(:) = ZERO
  rel_bias(:) = ZERO
  rel_std(:) = ZERO    
  OMPS_mean_nR(:) = 0.0
  CRTM_mean_nR(:) = 0.0
   
  CLOSE(69)
  DO k = 1, n_profiles
     IF(QC_prof(k) == 0 .and. crtm_sat_inputs%edrCloudFraction(k) < 0.0001_fp) THEN
        bias(:) = bias(:) + crtm_sat_inputs%Normalized_Radiance(:,k)-CRTM_n_Rad(:,k)
        std(:) = std(:) + (crtm_sat_inputs%Normalized_Radiance(:,k)-CRTM_n_Rad(:,k))**2
        OMPS_mean_nR(:) = OMPS_mean_nR(:) + crtm_sat_inputs%Normalized_Radiance(:,k)
        CRTM_mean_nR(:) = CRTM_mean_nR(:) + CRTM_n_Rad(:,k)
     
        nsample = nsample + 1

        WRITE(69,'(I10)') k
        WRITE(69,'(10f9.3)') m_Wavelength(:,opt(k)%nFOV)
        WRITE(69,'(6E14.6)') crtm_sat_inputs%Normalized_Radiance(:,k)
        WRITE(69,'(6E14.6)') CRTM_n_Rad(:,k)

        rel_bias(:) = rel_bias(:) + (crtm_sat_inputs%Normalized_Radiance(:,k)-CRTM_n_Rad(:,k))/Obs_mean(:)
        rel_std(:) = rel_std(:) + ((crtm_sat_inputs%Normalized_Radiance(:,k)-CRTM_n_Rad(:,k))/Obs_mean(:))**2
     END IF
  END DO

  OMPS_mean_nR(:) = OMPS_mean_nR(:)/float(nsample)
  CRTM_mean_nR(:) = CRTM_mean_nR(:)/float(nsample)
     
  CLOSE(69)
  CLOSE(68)
  OPEN(68,file='OMPS_nR_mean.dat')
  WRITE(68,'(2I5)') n_profiles, nsample
  WRITE(68,'(10f9.3)') m_Wavelength(:,30)
  WRITE(68,'(6E14.6)') OMPS_mean_nR(:)
  WRITE(68,'(6E14.6)') CRTM_mean_nR(:)
  CLOSE(68)
     

  bias(:) = bias(:)/float(nsample)
  std(:) = sqrt(std(:)/float(nsample))
  rel_bias(:) = rel_bias(:)/float(nsample) * 100.0_fp    ! %
  rel_std(:) = sqrt(rel_std(:)/float(nsample))*100.0_fp  ! %

  CLOSE(66)
  DO k = 1, n_channels
     WRITE(6,'(I5,f9.3,3E13.4,2f10.3)') k,CRTM_wavelength(k),Obs_mean(k),bias(k),std(k),rel_bias(k),rel_std(k)
     WRITE(66,'(I5,f9.3,3E13.4,2f10.3)') k,CRTM_wavelength(k),Obs_mean(k),bias(k),std(k),rel_bias(k),rel_std(k)
  END DO    
  print *,' samples ',n_profiles, nsample, float(nsample)/float(n_profiles)*100.0
  WRITE(66,'(a,2I10,f10.4)') ' samples ',n_profiles, nsample, float(nsample)/float(n_profiles)*100.0     
  print *, maxval(geo(:)%iFOV), minval(geo(:)%iFOV)
  ! Cleanup
  err_stat = CRTM_Destroy(chinfo)
  DEALLOCATE( atm, sfc, geo )

END PROGRAM CRTM_OMPS_Simulation

