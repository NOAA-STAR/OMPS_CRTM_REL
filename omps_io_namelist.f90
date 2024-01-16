!----------------------------------------------------------------------------------
!:sdoc+:
!
!
! OMPS_NAMELIST_MODULE
!
! This module is used to define a composite data type with the
! parameters controling the CRTM IOs of the different OMPS sensors.
!
! The "static" IO controling parameters are set inside the moudle with
! the default values. The instance parameters are set "dynamically" by 
! reading in a namelist file. A static parameter value may be also 
! dynamically changed if users include it in the "dynamic" set. 
!
!
! CREATION HISTORY:
!       Written by:     Ming  Chen, 12-January-2024
!                       ming.chen@noaa.gov
!

MODULE OMPS_NAMELIST_MODULE 

  ! -----------------
  ! Environment setup
  ! -----------------

  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE

  PUBLIC :: OMPS_IO_NAMELIST
  PUBLIC :: ReadNameList
  PUBLIC :: get_select_ch

  INTEGER, PARAMETER :: n_max_select = 3
  TYPE :: OMPS_IO_NAMELIST
        PRIVATE 
        CHARACTER(len=255),PUBLIC :: namelist_file=""
        CHARACTER(len=255),PUBLIC :: path_dir="" 
        CHARACTER(len=255),PUBLIC :: Sensor_Id 
        CHARACTER(len=255),PUBLIC :: atm_file, sfc_file, angle_file
        CHARACTER(len=255),PUBLIC :: surf_ref_file
        CHARACTER(len=255),PUBLIC :: sdr_file
        CHARACTER(len=255),PUBLIC :: CRTM_Coeff_Path=""
        CHARACTER(len=255),PUBLIC :: crtm_output_file
        CHARACTER(len=255),PUBLIC :: wL_file
        INTEGER ,PUBLIC :: nprof, nspec
        INTEGER, PUBLIC :: max_FOV
        LOGICAL, PUBLIC :: Use_allFOV_Nadir = .FALSE.
        INTEGER, PUBLIC :: Refl_alg 

        INTEGER :: n_Surf_ch
        INTEGER :: select_ch(n_max_select)

  END TYPE OMPS_IO_NAMELIST

CONTAINS


  FUNCTION ReadNameList(namelist_file) Result(namelist)
    CHARACTER(len=255),OPTIONAL :: namelist_file

    TYPE(OMPS_IO_NAMELIST) :: namelist
    INTEGER :: FileID =  78
    ! FileID = Get_Lun()                                                         
    IF (PRESENT(namelist_file)) THEN
      namelist%namelist_file = namelist_file
    ENDIF
 
    CLOSE(FileID)
    OPEN(FileID,file = TRIM(namelist%namelist_file))
    READ(FileID,'( a )' ) namelist%Sensor_Id
    ! set the defaults
    CALL set_Defaults(namelist)
    ! read in the dynamic parameters
    READ(FileID,'( a )' ) namelist%path_dir
    READ(FileID,'( a )' ) namelist%atm_file
    READ(FileID,'( a )' ) namelist%sfc_file
    READ(FileID,'( a )' ) namelist%angle_file    
    READ(FileID,'( a )' ) namelist%crtm_output_file
    READ(FileID,'( a )' ) namelist%surf_ref_file
    READ(FileID,'( a )' ) namelist%sdr_file
    READ(FileID,'( a )' ) namelist%CRTM_Coeff_Path
    READ(FileID,'( a )' ) namelist%wL_file
    CLOSE(FileID)

  END FUNCTION ReadNameList
  

  SUBROUTINE set_Defaults(namelist)
      TYPE(OMPS_IO_NAMELIST) :: namelist
      INTEGER :: n_Surf_ch 
      CHARACTER(len=255) :: Sensor_Id 

      namelist%Refl_alg = 3
      namelist%n_Surf_ch = 3

      IF (namelist%n_Surf_ch > n_max_select) THEN
         WRITE(*,*) " Error: namelist%n_Surf_ch should be <= ", n_max_select
      ENDIF
      Sensor_Id = namelist%Sensor_Id
      n_Surf_ch = namelist%n_Surf_ch
      ! default values for different sensors
      namelist%nprof = 175
      namelist%nspec = 260
      namelist%select_ch(1:n_Surf_ch) = (/61,117,175/) 
      print*,namelist%select_ch(1:n_Surf_ch)
      namelist%max_FOV = 35 
      ! manually add or modify the default values
      IF(trim(Sensor_Id)=='u.omps-tcP18_npp') THEN
         namelist%nspec = 196
      END IF
   
  END SUBROUTINE set_Defaults
 
  SUBROUTINE get_select_ch(namelist, n_Surf_ch, select_ch ) 
    TYPE(OMPS_IO_NAMELIST) :: namelist
    INTEGER :: n_Surf_ch
    INTEGER,ALLOCATABLE :: select_ch(:)
    n_Surf_ch = namelist%n_Surf_ch
    ALLOCATE(select_ch(n_Surf_ch))
    select_ch(1:n_Surf_ch) = namelist%select_ch(1:n_Surf_ch)
  END SUBROUTINE get_select_ch
  
  #only used for the tuning purpose
  SUBROUTINE set_select_ch(namelist, select_ch ) 
    TYPE(OMPS_IO_NAMELIST) :: namelist
    INTEGER :: select_ch(:)
    INTEGER :: n_Surf_ch
    n_Surf_ch = size(select_ch)
    namelist%n_Surf_ch = n_Surf_ch 
    IF (n_Surf_ch > n_max_select) THEN
       WRITE(*,*) " Error: n_Surf_ch should be <= ", n_max_select
    ENDIF
    namelist%select_ch(1:n_Surf_ch) = select_ch(1:n_Surf_ch)
  END SUBROUTINE set_select_ch

END MODULE OMPS_NAMELIST_MODULE 

