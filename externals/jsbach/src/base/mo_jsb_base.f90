!> Orchestrates setting up the ICON-Land model structure
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
MODULE mo_jsb_base
#ifndef __NO_JSBACH__

  USE mo_jsb_class,   ONLY: jsbach
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message, finish
  USE mo_util,        ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_setup_models, jsbach_setup_tiles

CONTAINS

  ! ======================================================================================================= !
  !> Orchestrates setting up the ICON-Land model structure

  SUBROUTINE jsbach_setup_models(master_namelist_filename)
    USE mo_util,                ONLY: separator
    USE mo_jsb_control,         ONLY: jsb_models_nml, get_no_of_models, init_jsb_master_control, &
      &                               jsbach_runs_standalone, jsbach_is_restarted
    USE mo_jsb_parallel,        ONLY: init_parallel_iface
    USE mo_jsb_model_class,     ONLY: t_jsb_model, new_model, MODEL_JSBACH
    USE mo_jsb_config_class,    ONLY: new_model_config
    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process
    USE mo_jsb_time,            ONLY: init_time
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,            ONLY: Register_vgrid, Get_vgrid
    USE mo_jsb_io,              ONLY: ZAXIS_SURFACE
    USE mo_pheno_parameters,    ONLY: init_pheno_parameters
    !>
    !> 1. Setup the ICON-Land models
    !>
    !>     1a. Routine is called by both:
    !>      ICON-Land standalone: jsbach: src/drivers/mo_jsbach_model.f90: jsbach_model
    !>
    !>     1b. Host model (ICON):    icon: src/drivers/mo_atmo_model.f90: construct_atmo_model
    !>

    dsl4jsb_Use_processes SEB_, PHENO_, TURB_, SSE_, HYDRO_, RAD_

    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(HYDRO_)
    dsl4jsb_Use_config(RAD_)

    CHARACTER(len=*), INTENT(in) :: master_namelist_filename !< icon_master.namelist - which models to run, which calendar, ...

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(RAD_)

    TYPE(t_jsb_model), POINTER :: model => NULL()

    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w, snow_e       !< Vertical soil and snow grids
    TYPE(t_jsb_vgrid), POINTER :: surface                      !< Vertical grid for surface (only one layer)

    INTEGER :: no_of_models, model_idx, iproc                  !< IDs of models and process
    INTEGER :: nsoil_e, nsoil_w, nsnow_e                       !< Number of vertical soil and snow layers
    INTEGER :: model_scheme                                    !< Choice of model (JSBACH)

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_models' !< Name of this routine

    CALL message('', '')
    CALL message('', separator)
    CALL message('', 'ICON-Land model configuration:')
    CALL message('', '')

    CALL init_parallel_iface ! setup communicator for MPI

    CALL init_jsb_master_control(master_namelist_filename) ! read ICON-LAND part of master namelist and basic setup (timers, ...)

    !> 2. Register all model instances
    no_of_models = get_no_of_models()
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models

    !> 3. Create models according to master namelist and add to jsbach model registry
    DO model_idx=1,no_of_models
      NULLIFY(model)
      model => new_model(jsb_models_nml(model_idx)%model_id,                &
        &                  jsb_models_nml(model_idx)%model_name,              &
        &                  jsb_models_nml(model_idx)%model_shortname,         &
        &                  jsb_models_nml(model_idx)%model_description,       &
        &                  jsb_models_nml(model_idx)%model_namelist_filename)
      jsbach%models(jsb_models_nml(model_idx)%model_id)%m => model
    END DO

    !> 4. Configure models according to model namelist(s) and (grid namelist(s) if standalone)
    DO model_idx=1,no_of_models
      NULLIFY(model)
      model => jsbach%models(model_idx)%m
      model%config => new_model_config(model%namelist_filename)

      CALL message(TRIM(routine), 'New ICON-Land model: '//TRIM(model%shortname))
      CALL message('     ', '... ID: '//TRIM(int2string(model%id)))
      IF (model%config%is_single_site) &
        & CALL message('     ', '... Single site: '//TRIM(model%config%single_site_name))
      CALL message('     ', '... Usecase: '//TRIM(model%config%usecase))
      CALL message('     ', '... Using configuration from namelist: '//TRIM(model%namelist_filename))
      CALL message('     ', '... Using tile fractions from: '//TRIM(model%config%fract_filename))

      ! Require that only one model_scheme is used!
      IF (model_idx == 1) THEN
        model_scheme = model%config%model_scheme
      ELSE IF (model%config%model_scheme /= model_scheme) THEN
         CALL finish(routine, 'Only one model_scheme allowed!')
      END IF

      !> 6. Set parameter init
      !>
      IF (model_idx == 1) THEN
        ! assuming the shared/constants/ modules are available here already (impl, math, physical)
        ! model specific
        SELECT CASE (model_scheme)
        CASE (MODEL_JSBACH)
          CALL message(TRIM(routine), 'starting init of JSBACH process parameters')
          CALL init_pheno_parameters
        END SELECT
      END IF

      !> 7. Standalone simulation initialize all time infos here, coupled simulations get most of it from ICON
      IF (jsbach_runs_standalone()) THEN
        CALL init_time(TRIM(model%namelist_filename), TRIM(model%shortname), jsbach_is_restarted())
      ELSE
        CALL init_time()
      END IF

      model%use_multiple_queues = .TRUE.

      !>
      !> 8. Create and register surface grid, other grids are created by the processes that needs them.
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(surface)

      !>
      !> 9. Create model processes.
      !>    9.a Note, that the maximum number of processes are create, but processes
      !  beyond the used processes get a NULL pointer
      !>
      !>  @TODO We need a better error message in case the number of processes used is > max_no_of_processes
      !>        or a different way to the the maximum number of processes automatically.
      !>
      ALLOCATE(model%processes(max_no_of_processes))
      DO iproc=1,max_no_of_processes
        model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
      END DO

      ! Point from each associated process%config to the model_config
      DO iproc = 1,max_no_of_processes
        IF (ASSOCIATED(model%processes(iproc)%p)) THEN
          IF (ASSOCIATED(model%processes(iproc)%p%config)) THEN
            model%processes(iproc)%p%config%model_config => model%config
          ENDIF
        ENDIF
      ENDDO

      !> 10. Configure processes: run `Init_PROC_config` (read namelists and create vertical axis if necessary, method from mo_jsb_class)
      CALL model%Configure_processes()

      dsl4jsb_Get_config(RAD_)
      dsl4jsb_Get_config(SSE_)
      dsl4jsb_Get_config(HYDRO_)

      !>
      !>     10a. Make number of vertical layers available from the options type
      IF (model%Is_process_enabled(HYDRO_)) THEN
        soil_w => Get_vgrid('soil_depth_water')
        nsoil_w = soil_w%n_levels
      END IF
      IF (model%Is_process_enabled(SSE_)) THEN
        soil_e => Get_vgrid('soil_depth_energy')
        nsoil_e = soil_e%n_levels
        snow_e => Get_vgrid('snow_depth_energy')
        nsnow_e = snow_e%n_levels
      END IF
      CALL model%Set_options(nsoil_e=nsoil_e, nsoil_w=nsoil_w, nsnow_e=nsnow_e)

      CALL model%Set_options(npft=model%config%npft)

      !>
      !>     10b. Crosscheck process configs
      IF (model%Is_process_enabled(HYDRO_) .AND. model%Is_process_enabled(SSE_) &
        & .AND. (nsoil_e /= nsoil_w)) THEN
        CALL finish(routine, &
          & 'The number of soil layers used for water and energy currently need to be identical:'  &
          & //' nsoil_w = '//TRIM(int2string(nsoil_w))//', nsoil_e = '//TRIM(int2string(nsoil_e)))
      END IF

      IF (TRIM(model%config%tpe_scheme) /= '') THEN
        CALL message(TRIM(routine), 'Running terraplanet experiment with scheme '//TRIM(model%config%tpe_scheme))
        IF (model%config%use_lakes) THEN
          CALL message('    ', '... using lakes')
        ELSE
          CALL message('    ', '... not using lakes')
        END IF
        IF (model%config%use_glacier) THEN
          CALL message('    ', '... using glaciers')
        ELSE
          CALL message('    ', '... not using glaciers')
        END IF
      END IF

      !>
      !>     10d. Check process dependencies

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)

        CALL message(TRIM(routine), 'Using canopy albedo from ini file '//TRIM(dsl4jsb_Config(RAD_)%bc_filename)//&
          &                         ' ... therefore no influence of vegetation on albedo')

      END SELECT

      IF (model%config%init_from_ifs) THEN
        CALL message(TRIM(routine), 'Initialize JSBACH from IFS analysis file '//TRIM(model%config%ifs_filename))
        IF (dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
          CALL message(routine, 'Reading initial soil moisture from '//TRIM(dsl4jsb_Config(HYDRO_)%ic_filename))
        END IF
      ELSE
        dsl4jsb_Config(HYDRO_)%l_read_initial_moist = .FALSE.
      END IF

      IF (.NOT. dsl4jsb_Config(HYDRO_)%l_infil_subzero) THEN
        ! Print message if infiltration at sub-zero temperature is switched off in namelist
        CALL message(TRIM(routine), 'Not allowing infiltration at sub-zero soil temperature')
      END IF

    END DO

    CALL message('', '')
    CALL message('', separator)
    CALL message('', '')

  END SUBROUTINE jsbach_setup_models

  ! ======================================================================================================= !
  !>
  !> Set up tiles.
  !>
  !> Routine is called by both:
  !>
  !>  - ICON-Land standalone: jsbach: src/drivers/mo_jsbach_model.f90: jsbach_model
  !>
  !>  - Host model (ICON):    icon: src/drivers/mo_atmo_model.f90: construct_atmo_model
  !> after setting up the model(s).
  !>
  SUBROUTINE jsbach_setup_tiles(model_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_model_usecases, ONLY: init_usecase

    INTEGER, INTENT(in) :: model_id !< model id (id of the domain in ICON)

    TYPE(t_jsb_model), POINTER  :: model               !< Current instance of the model

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_tiles'

    model => Get_model(model_id)

    !>
    !> Initialize usecase
    !>
    IF (model%config%usecase /= '') THEN
      ! Shortcut model configuration through usecase
      CALL init_usecase(model)
      ! The tile tree is complete: collect its tiles once for all later loops over tiles
      CALL model%Build_tile_lists()
    ELSE
      ! Explicit model configuration through namelist
      !>
      !> @TODO  At least add a finish in the if structure to initialize the model configuration thorugh the usecase or the namelist.
      !>
    END IF

  END SUBROUTINE jsbach_setup_tiles

#endif
END MODULE mo_jsb_base
