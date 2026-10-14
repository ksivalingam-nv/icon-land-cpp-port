!> Contains structures and methods for the ICON-Land model configuration
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> Contains structures and methods for the ICON-Land model configuration
!>
!> Two configuration types are handled in this module: the basic model configuration ([[t_jsb_model_config]])
!> and [[t_jsb_config]], the part of the individual process configurations common to all ICON-Land processes.
!>
MODULE mo_jsb_config_class
#ifndef __NO_JSBACH__
  !>
  !> @TODO `filename_max` and explicit string length definitions should not be needed with
  !>       `CHARACTER(len=:), ALLOCATABLE :: ` syntax.
  !>
  USE mo_kind,               ONLY: wp
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_jsb_io_netcdf,      ONLY: t_input_file

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_config, t_jsb_config_p, t_jsb_model_config, new_model_config

  !>
  !> Type for basic model configuration
  !>
  !> @TODO Check if `grid_filename` and `grid_name` variables are needed at all.
  !>
  !> @TODO remove `relative_fractions_in_file` function. We agreed to always use relative fractions in the initial file, although absolute (grid box)
  !>      fractions are used elsewhere in the code.
  !>
  !>
  TYPE t_jsb_model_config

    ! TODO Check if `grid_filename` and `grid_name` variables are needed at all.
    CHARACTER(len=filename_max)    :: grid_filename       !< Model grid file name
    CHARACTER(len=filename_max)    :: fract_filename      !< Name of file containing tile fractions
    CHARACTER(len=30)              :: grid_name           !< Model grid name
    CHARACTER(len=40)              :: usecase             !< Model use case
    CHARACTER(len=20)              :: model_scheme_char   !< Model variant: 'jsbach'
    INTEGER                        :: model_scheme        !< Model variant identifier
    INTEGER                        :: npft                !< Number of PFTs
    LOGICAL                        :: use_lakes           !< Lakes are represented
    LOGICAL                        :: use_glacier         !< Glaciers are represented
    LOGICAL                        :: check_child_cf      !< Assert that cover fractions of children tiles sum up to the tile cf
    INTEGER                        :: enforce_water_budget!< How to deal with water balance errors?</br>
                                                          !< 'ignore' (0):  short, one-line warning,</br>
                                                          !< 'logging' (1): detailed debugging output,</br>
                                                          !< 'error' (2):   debugging output and model stop.
                                                          !< Note: HYDRO parameter [[t_hydro_config:enforce_water_budget]]
                                                          !< can override this setting for the HYDRO process.
    CHARACTER(len=SHORT_NAME_LEN)  :: tpe_scheme          !< Switch for Terraplanet setup: Handling of water
                                                          !< cycle: 'open' or 'closed'
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:)  !< List of the tiles for which output will
                                                                   !< be generated (default: 'box')
    LOGICAL                        :: relative_fractions_in_file   !< Tile fractions in [[t_jsb_model_config:fract_filename]]
                                                                   !< are relative to the parent tiles (true)
                                                                   !< or relative to the grid box (false).
    ! TODO remove `relative_fractions_in_file` function. We agreed to always use relative fractions in the initial file, although absolute (grid box)
    !      fractions are used elsewhere in the code.

    ! Initialization from ECMWF Integrated Forecast System (IFS) analysis data
    LOGICAL                        :: init_from_ifs            !< Initialize from IFS analysis
    CHARACTER(len=:), ALLOCATABLE  :: ifs_filename             !< IFS input data file name
    TYPE(t_input_file)             :: ifs_input_file           !< Structure for the IFS initial file

    ! Single site setup
    LOGICAL                        :: is_single_site        !< Single site model
    CHARACTER(len=:), ALLOCATABLE  :: single_site_name      !< Single site name
    CHARACTER(len=:), ALLOCATABLE  :: single_site_shortname !< Single site short name
    REAL(wp)                       :: single_site_lonlat(2) !< Lon/lat of site in degrees

  END TYPE t_jsb_model_config

  !> Abstract type for process configurations
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active                                          !< The process is active in the current setup.
    LOGICAL :: lrestart_cont                                   !< Continue although process variables are
                                                               !< missing in the restart file.
    CHARACTER(len=filename_max) :: ic_filename                 !< Name of the initial file
    CHARACTER(len=filename_max) :: bc_filename                 !< Name of the boundary conditions file
    CHARACTER(LEN=filename_max) :: namelist_filename           !< File containing the process namelist
                                                               !< **Note**: generally all ICON-Land namelists are
                                                               !< listed in one common namelist file
    TYPE(t_jsb_model_config), POINTER :: model_config => NULL() !< Initialize model configuration pointer
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:) !< List of tiles for which this process
                                                                !< will write output
  CONTAINS
    !> Each process needs to define an init function as extension of this abstract type.
    PROCEDURE (Init_config), DEFERRED, PASS(config) :: Init
  END TYPE t_jsb_config

  !> Pointer type fot a `t_jsb_config` type
  TYPE t_jsb_config_p
    CLASS(t_jsb_config), POINTER :: p
  END TYPE t_jsb_config_p

  ABSTRACT INTERFACE
    !> Definition of the interface of the init function (needed with deferred procedures).
    !> The users are enforced to use the identical calling structure for all processes.
    SUBROUTINE Init_config(config)
      IMPORT :: t_jsb_config
      CLASS(t_jsb_config), INTENT(inout) :: config
    END SUBROUTINE Init_config
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class' !< Name of the module

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Creates the model configuration
  !>
  !> The function reads the model configuration from namelist "jsb_model_nml" defined in the given namelist
  !> file and returns a pointer to the model configuration structure.
  !>
  !> @TODO Are `grid_filename`, `grid_name` and the `jsb_grid_nml` needed at all?
  !>
  !> @TODO Dissolve the dependency<br>
  !> `model_scheme` values are
  !> defined here as `INTEGER` replacing the ENUM `USE mo_jsb_model_class,    ONLY: MODEL_JSBACH`
  !> This is needed to avoid a dependency cycle, creating a compile error (tested with NAG):
  !>  stating that in `mo_jsb_memory_class` the `mo_jsb_config_class` is not available for `USE mo_jsb_config_class,  ONLY: t_jsb_config`
  FUNCTION new_model_config(namelist_filename) RESULT(model_config)

    USE mo_jsb_namelist_iface,  ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_impl_constants,  ONLY: WB_IGNORE, WB_LOGGING, WB_ERROR    ! handling of water balance errors
    USE mo_jsb_control,         ONLY: jsbach_runs_standalone
    USE mo_util,                ONLY: tolower

    CHARACTER(len=*),         INTENT(in) :: namelist_filename     !< File containing the model namelist
    TYPE(t_jsb_model_config), POINTER    :: model_config          !< Model configuration: structure
                                                                  !< containing the configurable parameters
    ! TODO Are `grid_filename`, `grid_name` and the `jsb_grid_nml` needed at all?
    CHARACTER(len=filename_max) :: grid_filename                  !< File with the model grid file
    CHARACTER(len=30)           :: grid_name                      !< Model grid name

    INTEGER                     :: MODEL_JSBACH                   !< Identifier for the model scheme

    NAMELIST /jsb_grid_nml/   &
      grid_filename,          &
      grid_name

    !> Variables defined in [[t_jsb_model_config]]
    CHARACTER(len=20)             :: model_scheme                 !< Model scheme: 'jsbach'
    LOGICAL                       :: use_lakes                    !< Represent lakes
    LOGICAL                       :: use_glacier                  !< Represent glaciers
    LOGICAL                       :: check_child_cf               !< Assert that cover fractions of children tiles sum up to the tile cf
    CHARACTER(len=10)             :: enforce_water_budget         !< Consequence of water balance errors
    CHARACTER(len=filename_max)   :: fract_filename               !< File with the tile fractions
    CHARACTER(len=40)             :: usecase                      !< ICON-Land use case
    INTEGER                       :: npft                         !< Number of PFTs
    CHARACTER(len=SHORT_NAME_LEN) :: tpe_scheme                   !< Terraplanet setup: Handling of water
                                                                  !< cycle: 'open' or 'closed'
    CHARACTER(len=SHORT_NAME_LEN) :: output_tiles(99)             !< List of the output tiles
    LOGICAL                       :: relative_fractions_in_file   !< Tile fractions in [[t_jsb_model_config:fract_filename]]
                                                                  !< are relative to the parent tiles.
    LOGICAL                       :: init_from_ifs                !< Initialize ICON-Land from IFS data
    CHARACTER(len=filename_max)   :: ifs_filename                 !< File containing IFS data
    LOGICAL                       :: is_single_site
    CHARACTER(len=filename_max)   :: single_site_name
    CHARACTER(len=SHORT_NAME_LEN) :: single_site_shortname
    REAL(wp)                      :: single_site_lonlat(2)

    ! Definition of the namelist
    NAMELIST /jsb_model_nml/  &
      model_scheme, &
      use_lakes, &
      use_glacier, &
      check_child_cf, &
      enforce_water_budget, &
      usecase, &
      npft, &
      tpe_scheme, &
      fract_filename, &
      output_tiles, &
      relative_fractions_in_file, &
      init_from_ifs, &
      ifs_filename, &
      is_single_site, &
      single_site_name, &
      single_site_shortname, &
      single_site_lonlat


    INTEGER ::       &
      & nml_handler, &     !< File handler for the namelist
      & nml_unit,    &     !< Unit number of the namelist file
      & istat,       &     !< Status reading the namelist
      & i                  !< Looping index

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_model_config' !< Name of the module

    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

    ALLOCATE(model_config)

    ! TODO Dissolve the dependency<br>
    ! `model_scheme` values are
    ! defined here as `INTEGER` replacing the ENUM `USE mo_jsb_model_class,    ONLY: MODEL_JSBACH`
    ! This is needed to avoid a dependency cycle, creating a compile error (tested with NAG):
    !  stating that in `mo_jsb_memory_class` the `mo_jsb_config_class` is not available for `USE mo_jsb_config_class,  ONLY: t_jsb_config`
    MODEL_JSBACH               = 1

    ! Set defaults
    model_scheme               = "jsbach"
    use_lakes                  = .TRUE.
    use_glacier                = .TRUE.
    check_child_cf             = .FALSE.
    enforce_water_budget       = 'ignore'
    grid_filename              = ''
    grid_name                  = ''
    usecase                    = ''
    npft                       = 0
    tpe_scheme                 = ''
    fract_filename             = 'bc_land_frac.nc'
    output_tiles(:)            = ''
    relative_fractions_in_file = .TRUE.
    init_from_ifs              = .FALSE.
    ifs_filename               = 'ifs2icon.nc'
    is_single_site             = .FALSE.
    single_site_name           = ''
    single_site_shortname      = ''
    single_site_lonlat         = [0._wp, 0._wp]

    ! Read namelist
    nml_handler = open_nml(TRIM(namelist_filename))

    ! jsb_grid_nml
    nml_unit = position_nml('jsb_grid_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_grid_nml)

    ! Write resulting values in config memory
    model_config%grid_filename = grid_filename
    model_config%grid_name     = grid_name

    ! jsb_model_nml
    nml_unit =  position_nml('jsb_model_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_model_nml)
    CALL close_nml(nml_handler)

    model_config%model_scheme_char = tolower(model_scheme)
    SELECT CASE (TRIM(model_config%model_scheme_char))
    CASE ("jsbach")
      model_config%model_scheme = MODEL_JSBACH
      CALL message(TRIM(routine), 'Use model_scheme JSBACH')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Selected model_scheme '//TRIM(model_scheme)//' not available.')
    END SELECT
    !>
    !> @TODO To be deleted: `relative_fractions_in_file` will be removed (clean up task)
    !>
    ! No absolute fraction files are available, therefore, assert accidental use of flag:
    IF(.NOT. relative_fractions_in_file) THEN
      CALL finish(TRIM(routine), 'relative_fractions_in_file was set to true, however, currently' &
      & // ' no absolute fraction files are available, please check.')
    ENDIF

    model_config%relative_fractions_in_file = relative_fractions_in_file
    ! set other namelist options
    model_config%use_lakes                  = use_lakes
    model_config%use_glacier                = use_glacier
    model_config%check_child_cf             = check_child_cf
    model_config%usecase                    = usecase
    model_config%npft                       = npft
    model_config%tpe_scheme                 = tpe_scheme
    model_config%fract_filename             = fract_filename

    !>
    !> @TODO To be deleted `relative_fractions_in_file` (clean up task).
    !>
    ! Currently, no absolute fraction files are available. Assert, the flag was not set accidentally.
    IF(.NOT. relative_fractions_in_file) THEN
      CALL finish(TRIM(routine), 'relative_fractions_in_file was set to false, however, currently' &
      & // ' no absolute fraction files are available, please check your settings.')
    ENDIF

    SELECT CASE (tolower(TRIM(enforce_water_budget)))
    CASE ("ignore")
      model_config%enforce_water_budget = WB_IGNORE
      CALL message(TRIM(routine), 'WARNING: Land surface water balance will not be checked during simulation '// &
        & '(unless set differently in hydro or hd config).')
    CASE ("logging")
      model_config%enforce_water_budget = WB_LOGGING
      CALL message(TRIM(routine), 'WARNING: Simulation will not stop due to any land water balance violation '// &
        & '(unless set differently in hydro or hd config) but information will be added to the log file')
    CASE ("error")
      model_config%enforce_water_budget = WB_ERROR
      CALL message(TRIM(routine), 'WARNING: Simulation will stop due to any land water balance violation '// &
        & '(unless set differently in hydro or hd config).')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'enforce_water_budget == '//tolower(TRIM(enforce_water_budget))//' not available.')
    END SELECT

    !>
    !> `model_config%init_from_ifs`, and `model_config%ifs_filename`.
    !>
    model_config%init_from_ifs              = init_from_ifs
    model_config%ifs_filename               = TRIM(ifs_filename)
    model_config%is_single_site             = is_single_site
    model_config%single_site_name           = TRIM(single_site_name)
    model_config%single_site_shortname      = TRIM(single_site_shortname)
    model_config%single_site_lonlat(1:2)    = single_site_lonlat(1:2)

    IF (is_single_site .AND. .NOT. jsbach_runs_standalone()) &
        & CALL finish(routine, 'Single site only allowed in standalone experiment')

    ! Get the number of valid entries in the vector
    !>
    !> @TODO use `SIZE(output_tiles)` instead of `99` to explore the
    !> `IF (output_tiles(i) == '') EXIT` cases.
    !>
    DO i=1,99
      IF (output_tiles(i) == '') EXIT
    END DO
    i = i - 1  ! Number of tile names in namelist variable output_tiles
    IF (i == 0) THEN             ! Default is 'box' if no tiles specified in namelist
      i = 1
      output_tiles(1) = 'box'
    END IF
    ALLOCATE(model_config%output_tiles(i))
    model_config%output_tiles(1:i) = output_tiles(1:i)

  END FUNCTION new_model_config

#endif
END MODULE mo_jsb_config_class
