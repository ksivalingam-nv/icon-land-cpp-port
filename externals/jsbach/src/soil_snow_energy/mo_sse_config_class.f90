!> Contains the configuration type and setup routines for the soil and snow energy process.
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
!>#### Configuration of the soil and snow energy process
!>
!> The module defines the data structure [[t_sse_config]] that contains the soil and snow energy
!> configuration, i.e. all soil and snow energy variables that are configurable via namelist.
!> In subroutine [[Init_sse_config]] we read the soil and snow energy namelist 'jsb_sse_nml'
!> and set the configuration data structure accordingly. We also perform some consistency checks.
!>
MODULE mo_sse_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string, int2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_config
  PUBLIC :: max_snow_layers

  !> The data structure contains the configurable parameters of the SSE process. It extends
  !> [[mo_jsb_config_class:t_jsb_config]], which provides the common process settings shared
  !> by all ICON-Land processes.
  !>
  !> The additional fields in this type control the snow discretization, the soil thermal
  !> parameterization, and the treatment of freezing, and soil-texture-based properties.
  TYPE, EXTENDS(t_jsb_config) :: t_sse_config
    INTEGER  :: nsnow           !< Maximum number of snow layers used in the vertical snow discretization
    LOGICAL  :: l_dynsnow       !< Switch for dynamic snow heat conductivity and heat capacity depending on snow density
    LOGICAL  :: l_heat_cap_dyn  !< Switch for dynamic soil heat capacity; .TRUE. uses dynamic calculation,
                                !< .FALSE. uses a static map or FAO-based estimate depending on l_heat_cap_map
    LOGICAL  :: l_heat_cond_dyn !< Switch for dynamic soil heat conductivity; .TRUE. uses dynamic calculation,
                                !< .FALSE. uses a static map or FAO-based estimate depending on l_heat_cond_map
    !
    LOGICAL  :: l_heat_cap_map  !< Use soil heat capacity from an input map instead of deriving it from other soil properties
    LOGICAL  :: l_heat_cond_map !< Use soil heat conductivity from an input map instead of
                                !< deriving it from heat capacity and FAO data
    LOGICAL  :: l_heat_cap_pond !< Account for the heat capacity of pond water in the top soil layer heat capacity
    LOGICAL  :: l_soil_texture  !< Derive mineral soil thermal parameters from soil texture information
    LOGICAL  :: l_freeze        !< Consider freezing and thawing of soil water
    LOGICAL  :: l_supercool     !< Allow for supercooled soil water
    LOGICAL  :: l_uniform_snow  !< Temporary workaround to use a uniform snow fraction for soil temperature computations
    REAL(wp) :: w_soil_critical !< Critical water/ice content in the upper soil layer used to correct
                                !< the surface temperature during freezing and melting [m water equivalent]

  CONTAINS
    PROCEDURE :: Init => Init_sse_config
  END TYPE t_sse_config

  INTEGER, PARAMETER :: max_snow_layers = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_config_class'

CONTAINS

  ! ============================================================================================== !
  !>
  !>#### Configuration of the soil and snow energy process
  !>
  !> Configure the soil and snow energy process from the settings in the namelist
  !> 'jsb_sse_nml'. Init_sse_config is a type-bound procedure of [[t_sse_config]] and is
  !> called during the model initialization phase from
  !> [[mo_jsb_model_class:Configure_model_processes]].
  !>
  !> The routine validates the requested configuration, creates the vertical grid for the
  !> soil-energy scheme when needed, and stores the resulting settings in the configuration
  !> object. For a more extensive description of the namelist parameters, see [[t_sse_config]].
  !>
  SUBROUTINE Init_sse_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_sse_constants,      ONLY: snow_depth_min
    USE mo_jsb_math_constants, ONLY: eps8

    CLASS(t_sse_config), INTENT(inout) :: config

    LOGICAL  :: active           !< Switch that activates the soil and snow energy process.
    LOGICAL  :: lrestart_cont    !< Continue even if SSE variables are missing in the restart data.
    LOGICAL  :: l_dynsnow        !< Switch for dynamic snow heat conductivity and heat capacity depending on snow density.
    LOGICAL  :: l_heat_cap_dyn   !< Switch for dynamic soil heat capacity.
    LOGICAL  :: l_heat_cond_dyn  !< Switch for dynamic soil heat conductivity.
    LOGICAL  :: l_heat_cap_map   !< Use soil heat capacity from an input map
                                 !< instead of deriving it from other soil properties.
    LOGICAL  :: l_heat_cond_map  !< Use soil heat conductivity from an input map
                                 !< instead of deriving it from heat capacity and FAO data.
    LOGICAL  :: l_heat_cap_pond  !< Account for the heat capacity of pond water in the top soil layer.
    LOGICAL  :: l_freeze         !< Consider freezing and thawing of soil water.
    LOGICAL  :: l_supercool      !< Allow for supercooled soil water.
    LOGICAL  :: l_soil_texture   !< Derive mineral soil thermal parameters from soil texture information.
    LOGICAL  :: l_uniform_snow   !< Temporary workaround to use a uniform snow fraction for soil temperature computations.
    REAL(wp) :: w_soil_critical  !< Critical water/ice content in the upper soil layer used to correct the surface temperature
                                 !< during freezing and melting [m water equivalent].
    INTEGER  :: nsnow            !< Maximum number of snow layers used in the vertical snow discretization.
    REAL(wp) :: dz_snow(max_snow_layers) !< Snow-layer thicknesses [m] used to define the vertical snow discretization.
    CHARACTER(len=filename_max) :: ic_filename !< Name of the input file containing the initial-condition data.
    CHARACTER(len=filename_max) :: bc_filename !< Name of the input file containing the boundary-condition data.

    NAMELIST /jsb_sse_nml/                      &
      & active,                                 &
      & lrestart_cont,                          &
      & nsnow, dz_snow, l_dynsnow,              &
      & l_heat_cap_dyn, l_heat_cond_dyn,        &
      & l_heat_cap_map, l_heat_cond_map,        &
      & l_heat_cap_pond,                        &
      & l_soil_texture, l_uniform_snow,         &
      & l_freeze, l_supercool, w_soil_critical, &
      & ic_filename,                            &
      & bc_filename

    INTEGER :: nml_handler                              !< Handle for the namelist input file.
    INTEGER :: nml_unit                                 !< IO unit number used to read the namelist.
    INTEGER :: istat                                    !< Status flag returned by namelist handling routines.
    INTEGER :: i                                        !< Loop index
    TYPE(t_input_file) :: input_file                    !< Input-file descriptor used for initial data.
    REAL(wp), POINTER  :: ptr_1D(:)                     !< Temporary pointer to a one-dimensional real field.
    REAL(wp), ALLOCATABLE :: depths(:)                  !< Depth levels used to construct the soil-layer profile.
    REAL(wp), ALLOCATABLE :: mids(:)                    !< Midpoint depths of the vertical layers.
    INTEGER :: nsoil                                    !< Number of soil layers used in the SSE process.

    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_e !< Vertical grid for the soil layers.
    TYPE(t_jsb_vgrid), POINTER :: vgrid_snow_e !< Vertical grid for the snow layers.

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting soil and snow energy configuration')

    ! Set defaults
    active           = .TRUE.
    lrestart_cont    = .FALSE.  ! TRUE: Continue although SSE variables are missing in restart file
    nsnow            = 5
    dz_snow(:)       = 0._wp
    dz_snow(1:nsnow) = 0.05_wp
    l_dynsnow        = .TRUE.
    l_heat_cap_dyn   = .TRUE.
    l_heat_cond_dyn  = .TRUE.
    l_heat_cap_map   = .FALSE.
    l_heat_cond_map  = .FALSE.
    l_heat_cap_pond  = .FALSE.
    l_soil_texture   = .FALSE.
    l_freeze         = .TRUE.
    l_supercool      = .TRUE.
    l_uniform_snow   = .TRUE.
    w_soil_critical  = 5.85036E-3_wp
    ic_filename      = 'ic_land_soil.nc'
    bc_filename      = 'bc_land_soil.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_sse_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_sse_nml)

    CALL close_nml(nml_handler)

    config%active        = active
    config%lrestart_cont = lrestart_cont
    config%ic_filename             = ic_filename
    config%bc_filename             = bc_filename

    config%nsnow         = nsnow
    IF (nsnow > max_snow_layers) THEN
      CALL finish(TRIM(routine), 'Too many snow layers, maximum is '//TRIM(int2string(max_snow_layers)))
    ELSE IF (nsnow < 3) THEN
      CALL finish(TRIM(routine), 'Number of snow layers must be larger than 2.')
    ELSE
      CALL message(TRIM(routine), 'Using '//TRIM(int2string(nsnow))//' snow layers.')
    END IF

    config%l_dynsnow       = l_dynsnow
    config%l_heat_cap_dyn  = l_heat_cap_dyn
    config%l_heat_cond_dyn = l_heat_cond_dyn
    config%l_heat_cap_map  = l_heat_cap_map
    config%l_heat_cond_map = l_heat_cond_map
    config%l_soil_texture  = l_soil_texture
    IF (l_soil_texture) THEN
      CALL message(TRIM(routine), 'Using soil texture to derive heat capacity and heat conductivity')
    END IF
    config%l_heat_cap_pond = l_heat_cap_pond
    IF (l_heat_cap_pond) THEN
      CALL message(TRIM(routine), 'The heat capacity of ponds is accounted for in top soil layer heat capacity.')
    END IF
    config%l_freeze        = l_freeze
    config%l_supercool     = l_supercool
    config%l_uniform_snow  = l_uniform_snow
    IF (l_uniform_snow) THEN
      CALL message(TRIM(routine), 'Assuming uniform snow cover for soil temperature computation')
    END IF

    ! @Todo: This parameter does not seem to be used anywhere in the code. Check and remove if not needed.
    config%w_soil_critical     = w_soil_critical
    CALL message(TRIM(routine), 'Critical water/ice content in upper soil layer for correction of '// &
      &                         'surface temperature for freezing/melting: '//TRIM(real2string(w_soil_critical)))

    IF (.NOT. active) RETURN

    !
    !> create vertical grid for soil energy calculations
    !
    ! init vgrid with values from 'ic_filename' (default MODEL_JSBACH)
    input_file = jsb_netcdf_open_input(ic_filename)

    ! Note: At the moment, the soil layers for the energy calculations are the same as for the hydrology
    ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
    nsoil = SIZE(ptr_1D)
    ALLOCATE(depths(nsoil+1))
    ALLOCATE(mids(nsoil))
    depths(1) = 0._wp
    depths(2:nsoil+1) = ptr_1D(1:nsoil)
    mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
    DEALLOCATE(ptr_1D)

    CALL input_file%Close()

    vgrid_soil_e  => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
      & levels    = mids                 (1:nsoil  ),                              &
      & lbounds   = depths               (1:nsoil  ),                              &
      & ubounds   = depths               (2:nsoil+1),                              &
      & units='m')
    CALL register_vgrid(vgrid_soil_e)

    DEALLOCATE(depths, mids)

    CALL message(TRIM(routine), 'Init vgrid soil_depth_energy with values from input file: '//TRIM(ic_filename))

    ! Info saved to logfile
    WRITE(message_text, *) 'Soil levels in soil energy (upper) [m]: ', vgrid_soil_e%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (mid)   [m]: ', vgrid_soil_e%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (lower) [m]: ', vgrid_soil_e%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil level depths in soil energy   [m]: ', vgrid_soil_e%dz
    CALL message(TRIM(routine), message_text)

    !
    !> create vertical grid for snow energy calculations
    !
    IF (dz_snow(1) < snow_depth_min) THEN
      CALL finish(TRIM(routine), 'Depth of first snow layer should be larger than snow_depth_min=' &
        &                         //real2string(snow_depth_min))
    END IF
    ALLOCATE(depths(nsnow+1))
    ALLOCATE(mids(nsnow))
    depths(1) = 0._wp
    DO i=1,nsnow
      depths(i+1) = depths(i) + dz_snow(i)
    END DO
    mids(1:nsnow) = (depths(1:nsnow) + depths(2:nsnow+1)) / 2._wp

    vgrid_snow_e  => new_vgrid('snow_depth_energy', ZAXIS_GENERIC, nsnow,     &
      & levels    = mids                 (1:nsnow  ),                         &
      & lbounds   = depths               (1:nsnow  ),                         &
      & ubounds   = depths               (2:nsnow+1),                         &
      & units='m')
    CALL register_vgrid(vgrid_snow_e)

    ! Info saved to logfile
    WRITE(message_text, *) 'Snow levels in soil energy (upper) [m]: ', vgrid_snow_e%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow levels in soil energy (mid)   [m]: ', vgrid_snow_e%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow levels in soil energy (lower) [m]: ', vgrid_snow_e%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow level depths in soil energy   [m]: ', vgrid_snow_e%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

  END SUBROUTINE Init_sse_config

#endif
END MODULE mo_sse_config_class
