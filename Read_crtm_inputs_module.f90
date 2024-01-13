!
! Read_atms_inputs
!
! Module defining the needed CRTM inputs in the type: CRTM_user_inputs_type
! Parts of ICVS codes from Xin Jin are adopted in this module.
!
!  HISTORY:
!       Written by:     Xin Jin and Quanhua Liu, April-01-2022
!                       Quanhua.Liu@noaa.gov
!	07/07/2022 	Ding Liang	update nprof and nspec for SNPP and also update Convert_DOY_to_MonthDay
!
!	01/12/2024	Ming Chen	Implement the composite data type to control IOs, code refacoring and refinement
!					Implement the Derive_surface_reflectance interface
!
! 
 MODULE Read_atms_inputs_module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Intrinsic modules
  USE CRTM_Module
  USE Binary_File_Utility, ONLY: Open_Binary_File
  USE OMPS_Namelist_Module, ONLY: OMPS_IO_NAMELIST

  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE

  ! ...Structures
  PUBLIC :: crtm_sat_inputs_type
  ! ...Procedures
  PUBLIC :: Read_crtm_inputs
  PUBLIC :: Derive_surface_reflectance

 ! Atmosphere structure definition
  ! -------------------------------
  !:tdoc+:
  TYPE :: crtm_sat_inputs_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Dimension values
    INTEGER :: n_profiles               ! L
    INTEGER :: n_chan                   ! K
    ! EDR value
    Real(fp),ALLOCATABLE,DIMENSION(:) :: edrR331           ! L
    Real(fp),ALLOCATABLE,DIMENSION(:) :: edrR360           ! L
    Real(fp),ALLOCATABLE,DIMENSION(:) :: edrCloudFraction  ! L
    INTEGER,ALLOCATABLE,DIMENSION(:) ::  edrProfID         ! L    
!
    REAL(fp), ALLOCATABLE, DIMENSION(:,:) :: Normalized_Radiance
  END TYPE crtm_sat_inputs_type
  !:tdoc-:

  INTERFACE Derive_surface_reflectance
    ! In alighment with allFOV coeff set
    MODULE PROCEDURE Derive_surface_reflectance_a
    ! In alighment with  nadir-only coeff set
    MODULE PROCEDURE Derive_surface_reflectance_n 
  END INTERFACE Derive_surface_reflectance 

CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Atmosphere_Associated
!
! PURPOSE:
!       Elemental function to test the status of the allocatable components
!       of a CRTM Atmosphere object.
!
! CALLING SEQUENCE:
!       Status = CRTM_Atmosphere_Associated( Atm )
!
! OBJECTS:
!       Atm:       Atmosphere structure which is to have its member's
!                  status tested.
!                  UNITS:      N/A
!                  TYPE:       CRTM_Atmosphere_type
!                  DIMENSION:  Scalar or any rank
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:    The return value is a logical value indicating the
!                  status of the Atmosphere members.
!                    .TRUE.  - if the array components are allocated.
!                    .FALSE. - if the array components are not allocated.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Same as input
!
!:sdoc-:
!--------------------------------------------------------------------------------
  FUNCTION Read_crtm_inputs(namelist, Atm, Sfc, Geo, crtm_sat_inputs, &
    n_channels, Month, lat, lon, m_Wavelength) RESULT( Error_Status )

    ! Arguments
    USE CRTM_GeometryInfo, ONLY : Compute_AU_ratio2

    TYPE(OMPS_IO_NAMELIST) :: namelist
    TYPE( CRTM_Atmosphere_type ),  ALLOCATABLE,  DIMENSION( : ) :: Atm
    TYPE( CRTM_Surface_type ),     ALLOCATABLE,  DIMENSION( : ) :: Sfc
    TYPE( CRTM_Geometry_type ),    ALLOCATABLE,  DIMENSION( : ) :: Geo
    TYPE(crtm_sat_inputs_type), INTENT(INOUT) :: crtm_sat_inputs

    ! Function result
    INTEGER :: Error_Status

    CHARACTER(len=255) :: atm_file, sfc_file, angle_file, surf_ref_file, sdr_file
    CHARACTER(len=255) :: crtm_output_file, wL_file
    CHARACTER(len=255) :: err_msg, alloc_msg

    CHARACTER(*),  PARAMETER :: PROGRAM_NAME = 'Read_crtm_inputs'
    LOGICAL , PARAMETER :: QUIET = .TRUE.

    REAL (fp) :: Sensor_Zenith_Angle, Sensor_Azimuth_Angle, Earth_Sun_Distance_Ratio, & 
                 Source_Zenith_Angle, Source_Azimuth_Angle, Sensor_Scan_Angle
    REAL (fp) :: inter_fov, alat, alon
    REAL(fp), Allocatable :: lat(:), lon(:)
    INTEGER :: Year, Month, Day, DAYofYear
    INTEGER :: i, m, nkey
    INTEGER :: n_chan, n_channels
    INTEGER :: n_Scan_Angles,max_FOV
    INTEGER :: n_Profiles, nprof, nspec
    INTEGER :: FileID, Alloc_Stat, IO_Status

    !INTEGER, PARAMETER :: nprof = 175, nspec = 260
    !INTEGER, PARAMETER :: nprof = 7200, nspec = 260  ! 240*30 for J2 by Jingfeng
    !INTEGER, PARAMETER :: nprof = 1575, nspec = 260  ! J1 
    !INTEGER, PARAMETER :: nprof = 175, nspec = 196  ! 7/7/2022 by Ding Liang

    INTEGER(kind=2), ALLOCATABLE, DIMENSION(:) :: inter_year,inter_month,inter_day,inter_doy
    REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: inter_utc
    REAL(kind=4), ALLOCATABLE, DIMENSION(:,:) :: inter_robs, inter_sflx, N_Value
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: m_Wavelength
    REAL(kind=4), ALLOCATABLE, DIMENSION(:) :: inter_lat,inter_lon,inter_Reflectivity331
    REAL(kind=4), ALLOCATABLE, DIMENSION(:) :: inter_Reflectivity360,inter_CloudFraction,inter_satZenith
    REAL(kind=4), ALLOCATABLE, DIMENSION(:) :: inter_satAzimuth,inter_sunZenith,inter_sunAzimuth
    REAL(kind=4), ALLOCATABLE, DIMENSION(:) :: inter_tskin,inter_tpw,inter_ecmwf_flag,inter_water_flag
    
    nprof = namelist%nprof ; nspec = namelist%nspec
    max_FOV = namelist%max_FOV
    n_chan = crtm_sat_inputs%n_chan
    print *,'n_chan,nprof,nspec',n_chan,nprof,nspec

    ! allocate
    ALLOCATE( lat(nprof), lon(nprof) )
    ALLOCATE( m_Wavelength(n_channels,max_FOV) )
    ALLOCATE(inter_year(nprof),inter_month(nprof),inter_day(nprof),inter_doy(nprof) )
    ALLOCATE(inter_utc(nprof),inter_robs(nspec,nprof),inter_sflx(nspec,nprof),N_Value(nspec,nprof) )
    ALLOCATE(inter_lat(nprof),inter_lon(nprof),inter_Reflectivity331(nprof),inter_Reflectivity360(nprof) )
    ALLOCATE(inter_CloudFraction(nprof),inter_satZenith(nprof),inter_satAzimuth(nprof),inter_sunZenith(nprof) )
    ALLOCATE(inter_sunAzimuth(nprof),inter_tskin(nprof),inter_tpw(nprof),inter_ecmwf_flag(nprof),inter_water_flag(nprof) )

    ! io files
    atm_file = trim(namelist%atm_file)
    sfc_file = trim(namelist%sfc_file)
    angle_file = trim(namelist%angle_file)
    crtm_output_file = trim(namelist%crtm_output_file)
    surf_ref_file = trim(namelist%surf_ref_file)
    sdr_file = trim(namelist%sdr_file)                    
    wL_file = trim(namelist%wL_file)

    print *,trim(crtm_output_file)
    print *,trim(surf_ref_file)    
 
    ! read atmosphere
    Error_Status = CRTM_Atmosphere_InquireFile( TRIM(atm_file),   &  ! Input
                               n_Profiles=n_Profiles )
    ! ...Allocate input data structures
    crtm_sat_inputs%n_Profiles = n_Profiles
    ALLOCATE( Atm(n_profiles), Sfc(n_profiles), Geo(n_profiles), &
            STAT=alloc_stat, ERRMSG=alloc_msg )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message( PROGRAM_NAME, err_msg, alloc_stat ); STOP
    END IF  
    write(*,*) n_Profiles
    Error_Status = CRTM_Atmosphere_ReadFile( atm_file, atm, Quiet=QUIET )
    !print *,' atm%n_Clouds',atm(1)%n_Clouds, atm(2)%n_Clouds
    print *,' Absorber_ID ',atm(1)%Absorber_ID
    Error_Status = CRTM_Surface_ReadFile( sfc_file, sfc, Quiet=QUIET )

    !
    ! EDR surface reflectance at 331 nm and 360 nm, cloud fractiion
    ALLOCATE( crtm_sat_inputs%edrR331(n_profiles), crtm_sat_inputs%edrR360(n_profiles), &
        crtm_sat_inputs%edrCloudFraction(n_profiles),crtm_sat_inputs%edrprofID(n_profiles), &
        crtm_sat_inputs%Normalized_Radiance(crtm_sat_inputs%n_chan,n_profiles), &
        STAT=alloc_stat, ERRMSG=alloc_msg )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating edr reflectance arrays - '//TRIM(alloc_msg)
      CALL Display_Message( PROGRAM_NAME, err_msg, alloc_stat ); STOP
    END IF  


    FileID = Get_Lun()
    OPEN(FileID,file=TRIM(surf_ref_file))
    DO i= 1, n_Profiles
       READ(FileID,*) crtm_sat_inputs%edrR331(i), crtm_sat_inputs%edrR360(i), &
       crtm_sat_inputs%edrCloudFraction(i), crtm_sat_inputs%edrprofID(i)
    ENDDO
    CLOSE(FileID)
    !
    ! SDR data
    FileID = Get_Lun()
    open(FileID,file=trim(sdr_file),form='unformatted',convert='little_endian')
    read(FileID) inter_year      ,&
                     inter_doy       ,&
                     inter_utc       ,&
                     inter_robs      ,&
                     inter_sflx      ,&
                     inter_lat       ,&
                     inter_lon       ,&
                     inter_satZenith ,&
                     inter_satAzimuth,&
                     inter_sunZenith ,&
                     inter_sunAzimuth,&
                     inter_tskin     ,&
                     inter_tpw       ,&
                     inter_ecmwf_flag,&
                     inter_water_flag
                     !inter_fov 
                     !inter_for
    close(FileID)
    print *,' inter_year= ',inter_year(1)
    print *,' inter_doy= ',inter_doy(1)
    print *,' n_profiles=',n_profiles
  
    OPEN(FileID,file=TRIM(angle_file),form='unformatted')
    READ( FileID, IOSTAT = IO_Status ) nkey  ! should be  123456789
    READ( FileID, IOSTAT = IO_Status ) n_Scan_Angles 
    print *,' nS ',n_Scan_Angles
    IF(n_Profiles /= n_Scan_Angles ) THEN
      PRINT *,' error n_Profiles /= n_Scan_Angles '
      STOP
    END IF
     
    DO i = 1, n_Scan_Angles
      ! -------------------------------
      ! Read the scan_angle
      ! -------------------------------
      READ( FileID, IOSTAT = IO_Status )  alat, alon, Sensor_Zenith_Angle, &
                                          Sensor_Azimuth_Angle,  &
                                          Source_Zenith_Angle, &
                                          Source_Azimuth_Angle, &
                                          Sensor_Scan_Angle,&
                                          inter_fov 
      m = crtm_sat_inputs%edrprofID(i)+1
      YEAR = inter_year(m)
      DAYofYear = inter_doy(m)

      lat(i) = alat
      lon(i) = alon

      CALL Convert_DOY_to_MonthDay(YEAR, DAYofYear, MONTH, DAY)

      !print *,' leap ',i,m,YEAR, DAYofYear,MONTH, DAY
      CALL CRTM_Geometry_SetValue( Geo(i),   &
                               Longitude  = aLon,  &
                               Latitude   = aLat,  &         
                               Year       = Year, &         
                               Month      = Month,&         
                               Day        = Day,  &
                               Sensor_Scan_Angle   = Sensor_Scan_Angle ,   &
                               Sensor_Zenith_Angle = Sensor_Zenith_Angle , &
                               Sensor_Azimuth_Angle= Sensor_Azimuth_Angle, &
                               Source_Azimuth_Angle= Source_Azimuth_Angle, &
                               Source_Zenith_Angle = Source_Zenith_Angle,  &
                               iFOV   = int(inter_fov)+1  )
      IF (Geo(i)%Source_Azimuth_Angle .lt. 180) Geo(i)%Source_Azimuth_Angle = Geo(i)%Source_Azimuth_Angle + 180.0
      IF (Geo(i)%Source_Azimuth_Angle .gt. 180) Geo(i)%Source_Azimuth_Angle = Geo(i)%Source_Azimuth_Angle - 180.0
  
      !  7/4 has the longest and 1/3 has the shortest distance between Sun and Earth
      !  7/4 Earth_Sun_Distance_Ratio = 0.966589376173119; 1/3 = 1.03507737429997

      Earth_Sun_Distance_Ratio = Compute_AU_ratio2(YEAR, MONTH, DAY )
      crtm_sat_inputs%Normalized_Radiance(1:n_chan,i) = inter_robs(1:n_chan,m)/ &
        inter_sflx(1:n_chan,m)/Earth_Sun_Distance_Ratio

      !write(6,'(2I5,5E15.6)') i,m,crtm_sat_inputs%Normalized_Radiance(174,i), &
      !    inter_robs(174,m),inter_sflx(174,m),Earth_Sun_Distance_Ratio, &
      !    inter_robs(174,m)/inter_sflx(174,m)/Earth_Sun_Distance_Ratio
     
    END DO

    CLOSE(FileID)

    !FileID = Get_Lun()
    OPEN(FileID,file=TRIM(wL_file),form='unformatted')
    READ(FileID) m_Wavelength
    CLOSE(FileID)

    Deallocate(inter_year,inter_month,inter_day,inter_utc,inter_robs )
    Deallocate(inter_sflx,N_Value,inter_Reflectivity331 )
    Deallocate(inter_Reflectivity360,inter_CloudFraction,inter_satZenith)
    Deallocate(inter_satAzimuth,inter_sunZenith,inter_sunAzimuth,&
         inter_tskin,inter_tpw,inter_ecmwf_flag,inter_water_flag )

  END FUNCTION Read_crtm_inputs
 
  SUBROUTINE Convert_DOY_to_MonthDay(YEAR,DOY,MONTH,DAY)
  USE Date_Utility, ONLY: IsLeapYear
  INTEGER, INTENT(IN) :: YEAR, DOY
  INTEGER, INTENT(OUT) :: MONTH, DAY
  INTEGER :: i
  INTEGER, PARAMETER :: N_MONTHS = 12
  ! Days per Month in a non leap Year
  INTEGER, PARAMETER :: DAYS_PER_MONTH_IN_NONLEAP(N_MONTHS) = &
  (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  IF( IsLeapYear(YEAR) ) THEN
    IF(DOY <= 31 ) THEN
      MONTH = 1
      DAY = DOY
    ELSE IF( DOY <= 60 ) THEN
      MONTH = 2
      DAY = DOY - 31
    ELSE IF( DOY <= 91 ) THEN    
      MONTH = 3
      DAY = DOY - 60
    ELSE IF( DOY <= 121 ) THEN    
      MONTH = 4
      DAY = DOY - 91
    ELSE IF( DOY <= 152 ) THEN   ! replace 151 by 152 7/7/2022 Ding Liang 
      MONTH = 5
      DAY = DOY - 121      
    ELSE IF( DOY <= 182 ) THEN    
      MONTH = 6
      DAY = DOY - 152           !       
    ELSE IF( DOY <= 213 ) THEN    ! 
      MONTH = 7
      DAY = DOY - 182       
    ELSE IF( DOY <= 244 ) THEN    !
      MONTH = 8
      DAY = DOY - 213      
    ELSE IF( DOY <= 274 ) THEN    
      MONTH = 9
      DAY = DOY - 244
    ELSE IF( DOY <= 305 ) THEN    !
      MONTH = 10
      DAY = DOY - 274     
    ELSE IF( DOY <= 335 ) THEN    
      MONTH = 11
      DAY = DOY - 305             !
    ELSE
      MONTH = 12
      DAY = DOY - 335 
    END IF
  ELSE
    IF(DOY <= 31 ) THEN
      MONTH = 1
      DAY = DOY
    ELSE IF( DOY <= 59 ) THEN
      MONTH = 2
      DAY = DOY - 31
    ELSE IF( DOY <= 90 ) THEN    
      MONTH = 3
      DAY = DOY - 59
    ELSE IF( DOY <= 120 ) THEN    
      MONTH = 4
      DAY = DOY - 90
    ELSE IF( DOY <= 151 ) THEN    
      MONTH = 5
      DAY = DOY - 120
    ELSE IF( DOY <= 181 ) THEN    
      MONTH = 6
      DAY = DOY - 151   ! 181-151       7/7/2022 dliang
    ELSE IF( DOY <= 212 ) THEN    
      MONTH = 7
      DAY = DOY - 181       
    ELSE IF( DOY <= 243 ) THEN    
      MONTH = 8
      DAY = DOY - 212     
    ELSE IF( DOY <= 273 ) THEN    
      MONTH = 9
      DAY = DOY - 243
    ELSE IF( DOY <= 304 ) THEN    
      MONTH = 10
      DAY = DOY - 273    
    ELSE IF( DOY <= 334 ) THEN    
      MONTH = 11
      DAY = DOY - 304
    ELSE
      MONTH = 12
      DAY = DOY - 334
    END IF  
  
  END IF           
  END SUBROUTINE Convert_DOY_to_MonthDay
!

  SUBROUTINE Derive_surface_reflectance_a(n_Surf_ch, save_FOV, n_channels, &
       n_Profiles, LUT_wL, RTSolution, mea_n_rad, select_ch, Surf_refl,QC_prof)
  INTEGER :: n_Surf_ch, n_channels, n_Profiles
  INTEGER :: save_FOV(:)
  REAL(fp), DIMENSION(:,:) :: LUT_wL
  REAL(fp) :: CRTM_S_wavelength(n_Surf_ch), CRTM_wavelength(n_channels)
  REAL(fp), DIMENSION(:,:) ::  Surf_refl, mea_n_rad
  INTEGER, DIMENSION(:) :: select_ch
  REAL(fp), DIMENSION(n_Surf_ch) :: Rsphere, Rf0, Rf1, sur_ref_ch
  TYPE(CRTM_RTSolution_type) :: RTSolution(:,:)
  INTEGER :: i, L, k
  REAL(fp) :: slope
  INTEGER :: QC_prof(:)
!
    QC_Prof(:) = 0  
  DO L = 1, n_Profiles
     Rsphere(:) = RTSolution(1:n_Surf_ch,L)%Surface_Planck_Radiance
     Rf0(:) = (RTSolution(1:n_Surf_ch,L)%Radiance+RTSolution(1:n_Surf_ch,L)%Up_Radiance) &
     /RTSolution(1:n_Surf_ch,L)%SolarIrradiance
     Rf1(:) = RTSolution(1:n_Surf_ch,L)%Down_Radiance/RTSolution(1:n_Surf_ch,L)%SolarIrradiance  
!     
     Rf0(:) = mea_n_rad(:,L) - Rf0(:)     
! surface reflectance for selected channels
     sur_ref_ch(:) = Rf0(:)/(Rf1(:) + Rf0(:)*Rsphere(:))
!
     DO k = 1, n_Surf_ch
       IF( sur_ref_ch(k) < 0.001_fp ) THEN
         sur_ref_ch(k) = 0.001_fp
       ELSE IF( sur_ref_ch(k) > 0.999_fp ) THEN
         sur_ref_ch(k) = 0.999_fp
       END IF
     END DO
     IF( n_Surf_ch == 2 ) THEN
       slope = (sur_ref_ch(2)-sur_ref_ch(1))/(CRTM_S_wavelength(2)-CRTM_S_wavelength(1))
       Surf_refl(:,L) = sur_ref_ch(1) + slope*(CRTM_wavelength(:)-CRTM_S_wavelength(1))

       IF( minval(Surf_refl(:,L)) <= 0.001_fp .or. maxval(Surf_refl(:,L)) >= 0.999_fp ) QC_prof(L) = 1
       IF( Surf_refl(select_ch(1),L) <= 0.001_fp .or. Surf_refl(select_ch(1),L) >= 0.999_fp) QC_prof(L) = 2
       IF( Surf_refl(select_ch(2),L) <= 0.001_fp .or. Surf_refl(select_ch(2),L) >= 0.999_fp) QC_prof(L) = 3
       IF( Surf_refl(select_ch(1),L) <= 0.001_fp .and. Surf_refl(select_ch(2),L) <= 0.001_fp ) QC_prof(L) = 4
       IF( Surf_refl(select_ch(1),L) >= 0.999_fp .and. Surf_refl(select_ch(2),L) >= 0.999_fp ) QC_prof(L) = 4
       
       DO k = 1, n_channels
         IF( Surf_refl(k,L) < 0.001_fp ) THEN
           Surf_refl(k,L) = 0.001_fp
         ELSE IF( Surf_refl(k,L) > 0.999_fp ) THEN
           Surf_refl(k,L) = 0.999_fp
         END IF 
       END DO
! Three bands interpolation (TBD)
     ELSE  
      CRTM_S_wavelength(:) = LUT_wL(select_ch(:), save_FOV(L))
      CRTM_wavelength(:) = LUT_wL(:, save_FOV(L))
       DO k = 1, n_channels
         CALL Three_point_int(CRTM_S_wavelength, sur_ref_ch, CRTM_wavelength(k), Surf_refl(k,L))
         IF( Surf_refl(k,L) < 0.001_fp ) THEN
           Surf_refl(k,L) = 0.001_fp
         ELSE IF( Surf_refl(k,L) > 0.999_fp ) THEN
           Surf_refl(k,L) = 0.999_fp
         END IF       
       END DO     
!
       IF( minval(Surf_refl(:,L)) <= 0.001_fp .or. maxval(Surf_refl(:,L)) >= 0.999_fp ) QC_prof(L) = 1
       IF( Surf_refl(select_ch(1),L) <= 0.001_fp .or. Surf_refl(select_ch(1),L) >= 0.999_fp) QC_prof(L) = 2
       IF( Surf_refl(select_ch(2),L) <= 0.001_fp .or. Surf_refl(select_ch(2),L) >= 0.999_fp) QC_prof(L) = 3
       IF( Surf_refl(select_ch(3),L) <= 0.001_fp .and. Surf_refl(select_ch(3),L) <= 0.001_fp ) QC_prof(L) = 4
       IF( Surf_refl(select_ch(1),L) >= 0.999_fp .and. Surf_refl(select_ch(2),L) >= 0.999_fp ) QC_prof(L) = 5

     END IF

  END DO


  END SUBROUTINE Derive_surface_reflectance_a

  SUBROUTINE Derive_surface_reflectance_n(n_Surf_ch, CRTM_S_wavelength, n_channels, &
       n_Profiles, CRTM_wavelength, RTSolution, mea_n_rad, select_ch, Surf_refl)
  INTEGER :: n_Surf_ch, n_channels, n_Profiles
  REAL(fp), DIMENSION(:) :: CRTM_S_wavelength, CRTM_wavelength
  REAL(fp), DIMENSION(:,:) ::  Surf_refl, mea_n_rad
  INTEGER, DIMENSION(:) :: select_ch
  REAL(fp), DIMENSION(n_Surf_ch) :: Rsphere, Rf0, Rf1, sur_ref_ch
  TYPE(CRTM_RTSolution_type) :: RTSolution(:,:)
  INTEGER :: i, L, k
  REAL(fp) :: slope
  
  DO L = 1, n_Profiles
     Rsphere(:) = RTSolution(1:n_Surf_ch,L)%Surface_Planck_Radiance
     Rf0(:) = (RTSolution(1:n_Surf_ch,L)%Radiance+RTSolution(1:n_Surf_ch,L)%Up_Radiance) &
     /RTSolution(1:n_Surf_ch,L)%SolarIrradiance
     Rf1(:) = RTSolution(1:n_Surf_ch,L)%Down_Radiance/RTSolution(1:n_Surf_ch,L)%SolarIrradiance  
!     
!     print *,' Rf0 ',Rf0(:)
!     print *,' mea ', mea_n_rad(:,L)
     Rf0(:) = mea_n_rad(:,L) - Rf0(:)
!     print *,' Rsphere ',Rsphere(:)
!     print *,' Rf0 ',Rf0(:)
!     print *,' Rf1 ',Rf1(:)
          
! surface reflectance for selected channels
     sur_ref_ch(:) = Rf0(:)/(Rf1(:) + Rf0(:)*Rsphere(:))
!
     DO k = 1, n_Surf_ch
       IF( sur_ref_ch(k) < 0.01_fp ) THEN
         sur_ref_ch(k) = 0.01_fp
       ELSE IF( sur_ref_ch(k) > 0.99_fp ) THEN
         sur_ref_ch(k) = 0.99_fp
       END IF
     END DO
     IF( n_Surf_ch == 2 ) THEN
       slope = (sur_ref_ch(2)-sur_ref_ch(1))/(CRTM_S_wavelength(2)-CRTM_S_wavelength(1))
!         print *,slope,CRTM_S_wavelength(1:2),CRTM_wavelength(100:101)
       Surf_refl(:,L) = sur_ref_ch(1) + slope*(CRTM_wavelength(:)-CRTM_S_wavelength(1))
       DO k = 1, n_channels
         IF( Surf_refl(k,L) < 0.01_fp ) THEN
           Surf_refl(k,L) = 0.01_fp
         ELSE IF( Surf_refl(k,L) > 0.99_fp ) THEN
           Surf_refl(k,L) = 0.99_fp
         END IF       
!       write(6,'(I6,2f10.5)') L, sur_ref_ch(:)
!     write(6,'(10f9.4)') Surf_refl(:,L)
       END DO
! Three bands interpolation (TBD)
     ELSE  
       DO k = 1, n_channels
         CALL Three_point_int(CRTM_S_wavelength, sur_ref_ch, CRTM_wavelength(k), Surf_refl(k,L))
         IF( Surf_refl(k,L) < 0.01_fp ) THEN
           Surf_refl(k,L) = 0.01_fp
         ELSE IF( Surf_refl(k,L) > 0.99_fp ) THEN
           Surf_refl(k,L) = 0.99_fp
         END IF       
       END DO     
     END IF

  END DO


  END SUBROUTINE Derive_surface_reflectance_n
!
  SUBROUTINE Three_point_int(xx,yy,x,y)
  REAL(fp) :: xx(:),yy(:),x,y
  y = yy(1)*(x-xx(2))*(x-xx(3))/((xx(1)-xx(2))*(xx(1)-xx(3)))  &
    + yy(2)*(x-xx(1))*(x-xx(3))/((xx(2)-xx(1))*(xx(2)-xx(3)))  &
    + yy(3)*(x-xx(1))*(x-xx(2))/((xx(3)-xx(2))*(xx(3)-xx(1)))
  RETURN
  END SUBROUTINE Three_point_int
!

END MODULE Read_atms_inputs_module

