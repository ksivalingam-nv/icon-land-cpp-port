!> Contains methods for ICON-Land master control.
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
!> Methods for ICON-Land master control
!>
!> The module contains the subroutines to read the master namelist and set different times and debug options
!> for the ICON-Land model instances.
MODULE mo_jsb_control
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_io_units,  ONLY: filename_max
  USE mo_exception, ONLY: finish
  USE mo_util,      ONLY: int2string
  USE mo_timer,     ONLY: new_timer    ! ICON module

  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN ! Number of characters for short names
  USE mo_jsb_time_iface,     ONLY: l_timer_host   ! ICON namelist flag indicating if timers are enabled in the ICON host
                                                  ! This flag is independent of the ICON-Land timer_level.
  USE mo_jsb_io_iface,       ONLY: ldebugio
  USE mo_jsb_thermo_iface,   ONLY: init_jsb_thermo_iface, iaes, inwp

#ifdef _OPENACC
  USE openacc, ONLY: acc_async_sync
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_jsb_master_control
  PUBLIC :: jsbach_runs_standalone, jsbach_is_restarted, force_from_observations
  PUBLIC :: debug_on, timer_on, l_timer_host, get_debug_memory_level
  PUBLIC :: model_base_dir
  ! t_jsb_model_nml needs to be public for the NEC compiler
  PUBLIC :: t_jsb_model_nml, jsb_models_nml, get_no_of_models
  PUBLIC :: timer_jsbach, timer_aggregate, timer_integrate, timer_process_task, timer_forcing
  PUBLIC :: max_no_of_models, acc_stream, acc_async_sync

  ! Needed for inlining timer_on and debug_on on NEC.
  PUBLIC :: debug_level, timer_level
  PROTECTED :: debug_level, timer_level

  LOGICAL,           SAVE :: is_standalone         !< True if ICON-Land runs standalone
  LOGICAL,           SAVE :: restart_jsbach        !< True if ICON-Land is restarted from a previous run
  CHARACTER(len=99), SAVE :: model_base_dir        !< Base directory where the model runs
  INTEGER,           SAVE :: debug_level           !< Debug level for ICON-Land (0 (no output) to 3)
                                                   !< The higher the number the more detail
  INTEGER,           SAVE :: debug_io              !< Debug level for IO - deprecated - only used for ECHAM, TODO: should be deleted
  INTEGER,           SAVE :: debug_memory_level    !< Debug level for memory use (0 (no output) to 3)
                                                   !< Number defines the report interval (1: monthly, 2: daily, 3: each timestep)
  INTEGER,           SAVE :: timer_level           !< Timer level (0 to 3). The higher the more detail (1: basic, 2: detail, 3: all)
                                                   !< 0 disables (additional) icon land timers (or all if icon-land standalone),
                                                   !< but land timers can still be enabled by the icon host.
  LOGICAL,           SAVE :: l_force_from_obs      !< Always true for standalone runs, TODO: should be deleted
  INTEGER,           SAVE :: acc_stream = 1        !< Used for CUDA graphs

#ifndef _OPENACC
  INTEGER, PARAMETER      :: acc_async_sync = -42  !< Use for CUDA graphs if not using OpenACC
#endif

  !> Structure for JSBACH model namelist
  !>
  !> **Note**: while in ICON `timer_jsbach`, `timer_aggregate`,
  !>`timer_integrate`, `timer_process_task`, and
  !>`timer_forcing` are hardcoded, they are not hardcoded in ICON-Land because
  !> we do not know in advance which processes will be active.
  !>
  TYPE t_jsb_model_nml
    INTEGER            :: model_id = 0                       !< Model ID
    CHARACTER(len=30)  :: model_name                         !< Model name
    CHARACTER(len=SHORT_NAME_LEN) :: model_shortname         !< Model short name
    CHARACTER(len=132) :: model_description                  !< Model description
    CHARACTER(len=filename_max) :: model_namelist_filename   !< ICON-Land namelist filename
  END TYPE t_jsb_model_nml

  INTEGER, PARAMETER    :: max_no_of_models = 10                  !< Maximum number of allowed models
  INTEGER :: no_of_models = 0                                     !< Number of actual model
  TYPE(t_jsb_model_nml), SAVE :: jsb_models_nml(max_no_of_models) !< Array of namelist structures for each model


  ! **Note**: while in ICON `timer_jsbach`, `timer_aggregate`,
  ! `timer_integrate`, `timer_process_task`, and
  ! `timer_forcing` are hardcoded, they are not hardcoded in ICON-Land because
  ! we do not know in advance which processes will be active.
  INTEGER, SAVE, DIMENSION(max_no_of_models) :: &
    & timer_jsbach,       &  !< Timer for the whole land model
    & timer_aggregate,    &  !< Timer for aggregation
    & timer_integrate,    &  !< Timer for update routines
    & timer_process_task, &  !< Timer for processing tasks on the tile tree
    & timer_forcing          !< Time for forcing

  CHARACTER(len=filename_max) :: master_namelist_filename = 'jsbach_master.namelist' !< Name of the master namelist file

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_control'  !< Module name for error messages

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initializes JSBACH master control.
  !>
  !> Read the master namelist and set timers for each model instance.
  !>
  SUBROUTINE init_jsb_master_control(namelist_filename)

    USE mo_jsb_io,      ONLY: init_jsb_io

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: namelist_filename  !< Name of the master namelist file to read

    INTEGER :: i !< Loop index for number of models

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_jsb_master_control' !< Routine name for error messages

    IF (PRESENT(namelist_filename)) master_namelist_filename = TRIM(namelist_filename)

    !>
    !> 1. Read the control namelist and set the `no_of_models` variable.
    !>
    CALL read_jsb_control_namelist(TRIM(master_namelist_filename)) ! This sets *no_of_models* and reads model namelists

    CALL init_jsb_io(TRIM(master_namelist_filename))

    !>
    !> 2. Set up the ICON-Land timers
    !>
    IF (timer_on() .OR. l_timer_host) THEN
      IF (no_of_models > 1) THEN
        DO i=1,no_of_models
          timer_jsbach(i)  = new_timer('jsbach_'//int2string(i))
          timer_forcing(i) = new_timer('jsb:forcing_'//int2string(i))
        END DO
      ELSE
        timer_jsbach(1)  = new_timer('jsbach')
        timer_forcing(1) = new_timer('jsb:forcing')
      END IF
    END IF
    IF (timer_on()) THEN
    !> @TODO Right now, checking if timers are enabled for ICON-Land specific timers is not done according to `timer_level`. This could be improved!
      IF (no_of_models > 1) THEN
        DO i=1,no_of_models
          timer_aggregate(i)    = new_timer('jsb:aggregate_'//int2string(i))
          timer_integrate(i)    = new_timer('jsb:integrate_'//int2string(i))
          timer_process_task(i) = new_timer('jsb:process_task_'//int2string(i))
        END DO
      ELSE
        timer_aggregate(1)    = new_timer('jsb:aggregate')
        timer_integrate(1)    = new_timer('jsb:integrate')
        timer_process_task(1) = new_timer('jsb:process_task')
      END IF
    END IF

  END SUBROUTINE init_jsb_master_control

  ! ======================================================================================================= !
  !>
  !> Reads JSBACH control parameters.
  !>
  SUBROUTINE read_jsb_control_namelist(namelist_filename)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CHARACTER(len=*), INTENT(in) :: namelist_filename       !< Name of the master namelist file to read

    INTEGER            :: model_id                          !< Model ID
    CHARACTER(len=30)  :: model_name                        !< Model name
    CHARACTER(len=SHORT_NAME_LEN) :: model_shortname        !< Model short name
    CHARACTER(len=132) :: model_description                 !< Model description
    CHARACTER(len=filename_max) :: model_namelist_filename  !< ICON-Land namelist filename

    ! Variables for namelist handling
    INTEGER :: nml_handler, &  !< Handle to deal with namelist file
      & nml_unit,           &  !< IO unit of the namelist file
      & istat,              &  !< Status flag
      & i                      !< Loop index for number of models
    LOGICAL :: rwnd            !< True: rewind namelist unit for reading again

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_jsb_control_namelist'  !< Routine name for error messages

    ! These two namelists are created in the master namelist file (i.e. separately from the process namelists)
    NAMELIST /jsb_control_nml/ &
      is_standalone,           &
      restart_jsbach,          &
      model_base_dir,          &
      timer_level,             &
      debug_level,             &
      debug_io,                &
      debug_memory_level,      &
      l_force_from_obs

    NAMELIST /jsb_model_nml/   &
      model_id,                &
      model_name,              &
      model_shortname,         &
      model_description,       &
      model_namelist_filename

    ! Namelist defaults
    is_standalone      = .TRUE.   ! Set in the runscript
    restart_jsbach     = .FALSE.  ! Only relevant for standalone model
    model_base_dir     = '.'
    timer_level        = 0        ! I.e. no (additional) ICON-Land timers
    debug_level        = 0        ! I.e. no debugging output
    debug_io           = 0
    debug_memory_level = 0        ! I.e. no debugging output for memory use
    l_force_from_obs   = .TRUE.

    nml_handler = open_nml(TRIM(namelist_filename))

    !>
    !> 1. Reads namelist `jsb_control_nml` and namelists `jsb_model_nml` from master
    !> namelist file for each model instance.
    !>
    nml_unit = position_nml('jsb_control_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_control_nml)

    ! Get number of models first by counting the number of jsb_model_nml sections in the namelist
    rwnd = .TRUE. ! start reading at the beginning of the namelist file (as opposed to continuing at the current position)
    no_of_models = 0
    DO
      nml_unit = position_nml('jsb_model_nml', nml_handler, lrewind=rwnd, STATUS=istat)
      IF (istat /= POSITIONED) EXIT
      rwnd=.FALSE. ! do not rewind again
      no_of_models = no_of_models + 1
      READ(nml_unit, jsb_model_nml)
    END DO

    IF (no_of_models == 0 ) THEN
      CALL finish(routine, 'Namelist jsb_model_nml missing in '//TRIM(namelist_filename)// &
        & '. One jsb_model_nml section is needed for each model incidence.')
    END IF
    IF (no_of_models > max_no_of_models) THEN
      CALL finish(routine, 'no_of_models > max_no_of_models')
    ENDIF

    !>
    !> @TODO Delete `IF (is_standalone .AND. .NOT. l_force_from_obs) THEN`
    !> (seems to be obsolete and such a test outside of a loop over models
    !> would not make sense)
    !>
    IF (is_standalone .AND. .NOT. l_force_from_obs) THEN
      CALL finish(routine, 'l_force_from_obs needs to be true - also with forcing data generated from '// &
        & 'model output. Forcing from atmosphere interface data currently not implemented.')
    END IF

    !>
    !> @TODO `ldebugio = debug_io > 0` this is a bit arbitrary since `debug_io`
    !> will be the one of the last model.
    !> It could also be taken from the first
    !> and should then be compared in the loop below and safeguarded that its
    !> the same for all models.
    !>
    ldebugio = debug_io > 0

    !>
    !> 2. Initialize the JSBACH thermodynamic interface in icon:mo_util_jsbach for different saturation
    !> vapor pressure formulations depending on the atmospheric physics used in the host model (ICON)
    !> or as specified in the namelist via `iforcing = 2/3` in `run_nml` for 2=AES, 3=NWP (standalone; default: NWP).
    !>
    CALL init_jsb_thermo_iface(is_standalone)

    !> 3. Now read namelists again and set parameters.
    !> **Note**: `no_of_models` needs to be known here (therefore extra loop above) in order to test for a valid model_id below
    !>
    rwnd = .TRUE.
    DO i=1,no_of_models
      nml_unit = position_nml('jsb_model_nml', nml_handler, lrewind=rwnd, STATUS=istat)
      IF (istat /= POSITIONED) THEN
        CALL finish(TRIM(routine), 'jsb_model_nml not found in namelist')
      END IF
      rwnd=.FALSE.

      model_id = i
      model_name = ''
      model_shortname = ''
      model_description = ''
      model_namelist_filename = ''

      READ(nml_unit, jsb_model_nml)

      IF (jsb_models_nml(model_id)%model_id > 0) THEN
        CALL finish(TRIM(routine), 'Model id already taken.')
      END IF

      IF (model_id > no_of_models) THEN
        CALL finish(TRIM(routine), 'Model id larger than number of models.')
      END IF

      jsb_models_nml(i)%model_id                = model_id
      jsb_models_nml(i)%model_name              = model_name
      jsb_models_nml(i)%model_shortname         = model_shortname
      jsb_models_nml(i)%model_description       = model_description
      jsb_models_nml(i)%model_namelist_filename = model_namelist_filename
    END DO

    CALL close_nml(nml_handler)

  END SUBROUTINE read_jsb_control_namelist

  !> Return the value of the `is_standalone` variable
  LOGICAL FUNCTION jsbach_runs_standalone()

    jsbach_runs_standalone = is_standalone

  END FUNCTION jsbach_runs_standalone

  !> Returns true if the model is forced from observations
  LOGICAL FUNCTION force_from_observations()
    force_from_observations = l_force_from_obs

  END FUNCTION force_from_observations

  !> Returns true if the model was restarted
  LOGICAL FUNCTION jsbach_is_restarted()

    !> @TODO Get `restart_jsbach` from parent model if not standalone
    jsbach_is_restarted = restart_jsbach

  END FUNCTION jsbach_is_restarted

  !> Returns the number of models
  INTEGER FUNCTION get_no_of_models()

    get_no_of_models = no_of_models

  END FUNCTION get_no_of_models

  ! ======================================================================================================= !
  !>
  !> Returns true if the debug level has been set to a level exceeding the level passed as argument
  LOGICAL FUNCTION debug_on(level)

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: level !< Debugging level

    CHARACTER(len=*), PARAMETER :: routine = modname//':debug_on' !< Name of the function

    IF (PRESENT(level)) THEN
      SELECT CASE (TRIM(level))
      CASE ('basic')
        debug_on = debug_level >= 1
      CASE ('detail')
        debug_on = debug_level >= 2
      CASE ('tasks')
        debug_on = debug_level >= 3
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown debug level')
      END SELECT
    ELSE
      debug_on = debug_level > 0
    END IF

  END FUNCTION debug_on

  !> Returns the debug level for memory
  INTEGER FUNCTION get_debug_memory_level()

    get_debug_memory_level = debug_memory_level

  END FUNCTION get_debug_memory_level


  ! ======================================================================================================= !
  !>
  !> Returns true if the timer level exceeds the level passed as argument
  LOGICAL FUNCTION timer_on(level)

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: level !< Debugging level

    CHARACTER(len=*), PARAMETER :: routine = modname//':timer_on' !< Name of the function

    IF (PRESENT(level)) THEN
      SELECT CASE (TRIM(level))
      CASE ('basic')
        timer_on = timer_level >= 1
      CASE ('detail')
        timer_on = timer_level >= 2
      CASE ('all')
        timer_on = timer_level >= 3
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown timer level')
      END SELECT
    ELSE
      timer_on = timer_level > 0
    END IF

  END FUNCTION timer_on

#endif
END MODULE mo_jsb_control
