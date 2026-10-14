!> Contains time control for JSBACH
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
MODULE mo_jsb_time
#ifndef __NO_JSBACH__

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_jsb_control,    ONLY: jsbach_runs_standalone
  USE mo_jsb_time_iface, ONLY: & ! Uses the mtime library
    & t_datetime,              & ! Basic mtime type
                                 ! Object containing year, month, day, hour, minutes, seconds, milliseconds
                                 ! and routines to work with it
    & deallocateDatetime,      & ! Deallocate datetime (mtime)
    & get_time_next,           & ! Returns an object of t_datetime containing the next timestep - jsb util routine (not mtime)
    & get_time_previous,       & ! Returns an object of t_datetime containing the previous timestep - jsb util routine
    & get_date_components,     & ! jsb util routine: returns present components of this datetime object
    & get_time_start,          &
    & get_time_dt,             &
    & is_time_experiment_start,  & ! jsb util routine: test if this is first timestep (i.e. start) of experiment
    & get_time_experiment_start, &
    & is_time_restart,           & ! jsb util routine: test if this is first timestep after restart
    & get_year_length, get_month_length, get_day_length, get_year_day, &
    & is_time_ltrig_rad_m1,      &      ! jsb util routine: check if this is the timestep before radiation calculation
    & get_time_interpolation_weights, & ! jsb util routine: provides month index and weights for interpolation monthly values
    & get_asselin_coef, &
    & cache_time_interpolation_weights

  IMPLICIT NONE
  PRIVATE

!!$  PUBLIC :: init_jsb_time
  PUBLIC :: init_time
  PUBLIC :: &
    & t_datetime, &         ! Basic mtime type provided for further use in ICON-Land
    & deallocateDatetime, & ! mtime routine
    & get_time_next         ! jsb util routine provided for further use in ICON-Land
  PUBLIC :: get_date_components, get_time_start, get_time_experiment_start, &
    & get_time_dt,              & ! jsb util routine
    & is_time_experiment_start, & ! jsb util routine
    & is_time_restart,          & ! jsb util routine
    & is_time_ltrig_rad_m1,     & ! jsb util routine
            get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: is_newday, is_newmonth, is_newyear, get_year_at_model_start, get_year_at_experiment_start, &
            get_year, get_month, get_day, get_secs_of_day, get_previous_month_length, get_previous_year_length
  PUBLIC :: cache_time_interpolation_weights
  PUBLIC :: get_time_interpolation_weights ! jsb util routine
  PUBLIC :: get_asselin_coef
  PUBLIC :: timestep_in_days, timesteps_per_day

  TYPE(t_datetime), POINTER, SAVE :: start_date
  TYPE(t_datetime), POINTER, SAVE :: final_date

  LOGICAL, SAVE :: time_initialized = .FALSE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time'

CONTAINS

  SUBROUTINE init_time(namelist_filename, model_shortname, restart_jsbach)

    USE mo_jsb_time_iface, ONLY: read_time_namelist, & !, configure_time
                                 configure_time_and_events, &
                                 get_time_stop

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: namelist_filename !< Jsbach namelist file name
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname   !< Jsbach model
    LOGICAL, OPTIONAL, INTENT(IN) :: restart_jsbach             !< Restart state

    REAL(wp) :: eps

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_time'

    ! NOTE this routine does not seem to be necessary for JSBACH within ICON anymore and is just
    ! returning start_date and final_date. This should be solved in a more simple way.

    !TODO? currently the time is initalised from the namelist of the first processed model
    !      but then used globally for all models!
    !      If different models should have different time information this would need to change
    !      -> Attention: it is not possible to specify this global information before knowing the names
    !         of the models: knowing the models is a precondition for restart with echam, because echam restart files
    !        contain the short name of a model. The time-information of a restart file in turn needs
    !        to be processed before times and events can be initialised.
    IF (.NOT. time_initialized) THEN
      CALL message(TRIM(routine), 'Initializing JSBACH time')

      IF (PRESENT(restart_jsbach)) CONTINUE    ! Avoid compiler warnings about dummy arguments not being used
      IF (PRESENT(namelist_filename)) CONTINUE ! Avoid compiler warnings about dummy arguments not being used
      ! IF (PRESENT(namelist_filename)) THEN                   ! Present only for standalone model
      !   CALL read_time_namelist(TRIM(namelist_filename))     ! Already done from mo_util_jsbach.f90 in icon
      ! END IF

      IF (PRESENT(model_shortname)) THEN
        CALL configure_time_and_events(model_shortname)
      ELSE
        CALL configure_time_and_events()
      ENDIF

      start_date => get_time_start()
      final_date => get_time_stop()

      ! Just to avoid compiler warnings about variable set but never referenced
      IF (.FALSE.) THEN
        PRINT*, get_year_day(start_date)
        PRINT*, get_year_day(final_date)
      END IF

      ! Check if Asselin time filter is used and print message
      eps = get_asselin_coef()
      IF (eps > 0._wp) THEN
        WRITE(message_text,*) 'Using Asselin time filter for surface temperature: eps =', eps
      ELSE
        WRITE(message_text,*) 'Not using Asselin time filter for surface temperature'
      END IF
      CALL message(TRIM(routine), TRIM(message_text))

      time_initialized = .TRUE.
    ENDIF

  END SUBROUTINE init_time

  ! ======================================================================================================= !
  !>
  !> Returns true if the current and previous time step are from different days
  LOGICAL FUNCTION is_newday(current, dt)

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_day, previous_day

    previous => get_time_previous(current, dt)

    current_day = INT(get_year_day(current))
    previous_day = INT(get_year_day(previous))

    is_newday = current_day /= previous_day

    CALL deallocateDatetime(previous)

  END FUNCTION is_newday

  ! ======================================================================================================= !
  !>
  !> Returns true if the current and previous time step are from different months
  LOGICAL FUNCTION is_newmonth(current, dt)

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_month, previous_month

    previous => get_time_previous(current, dt)

    current_month  = INT(get_month(current))
    previous_month = INT(get_month(previous))

    is_newmonth = current_month /= previous_month

    CALL deallocateDatetime(previous)

  END FUNCTION is_newmonth

  ! ======================================================================================================= !
  !>
  !> Returns true if the current and previous time step are from different years
  LOGICAL FUNCTION is_newyear(current, dt)

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_year, previous_year

    previous => get_time_previous(current, dt)

    current_year  = get_year(current)
    previous_year = get_year(previous)

    is_newyear = current_year /= previous_year

    CALL deallocateDatetime(previous)

  END FUNCTION is_newyear

  ! ======================================================================================================= !
  !>
  !> Returns the year of the provided datetime variable
  INTEGER FUNCTION get_year(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: year

    CALL get_date_components(this_datetime, year=year)

    get_year = year

  END FUNCTION get_year

  ! ======================================================================================================= !
  !>
  !> Returns the month of the provided datetime variable
  INTEGER FUNCTION get_month(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: month

    CALL get_date_components(this_datetime, month=month)

    get_month = month

  END FUNCTION get_month

  ! ======================================================================================================= !
  !>
  !> Returns the day of the provided datetime variable
  INTEGER FUNCTION get_day(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: day

    CALL get_date_components(this_datetime, day=day)

    get_day = day

  END FUNCTION get_day

  ! ======================================================================================================= !
  !>
  !> Returns the number of seconds elapsed between start of a day and the provided datetime variable
  REAL(wp) FUNCTION get_secs_of_day(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: hour, minute, second, ms

    CALL get_date_components(this_datetime, hour=hour, minute=minute, second=second, ms=ms)

    get_secs_of_day = 3600._wp * REAL(hour,wp) + 60._wp * REAL(minute,wp) + REAL(second,wp) + REAL(ms)/1000._wp

  END FUNCTION get_secs_of_day

  ! ======================================================================================================= !
  !>
  !> Returns the number of seconds of the previous month
  !
  INTEGER FUNCTION get_previous_month_length(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: month, prev_month, year

    CALL get_date_components(this_datetime, month=month, year=year)
    IF (month == 1) THEN
      get_previous_month_length = NINT(31._wp * REAL(get_day_length(),wp))
    ELSE
      prev_month = month-1
      get_previous_month_length = NINT(REAL(get_month_length(year,prev_month),wp) * REAL(get_day_length(),wp))
    END IF

  END FUNCTION get_previous_month_length

  ! ======================================================================================================= !
  !>
  !> Returns the number of seconds of the previous year
  !
  INTEGER FUNCTION get_previous_year_length(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: year, prev_year

    CALL get_date_components(this_datetime, year=year)
    prev_year = year - 1
    get_previous_year_length = NINT(REAL(get_year_length(prev_year),wp) * REAL(get_day_length(),wp))

  END FUNCTION get_previous_year_length

  ! ======================================================================================================= !
  !>
  !> Returns the starting year of the current simulation
  !! e.g. 1900 for a run (re)started 1900 even though the full experiment might have started earlier
  INTEGER FUNCTION get_year_at_model_start()

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_year_at_model_start'

    IF (.NOT. time_initialized) &
      CALL finish(TRIM(routine), 'Time not initialized yet.')

    get_year_at_model_start = get_year(start_date)

  END FUNCTION get_year_at_model_start

  ! ======================================================================================================= !
  !>
  !> Returns the starting year of an experiment including time passed during earlier simulation periods
  !! e.g. 1850 for an experiment which started 1850 even if the current run (re)started 1900
  INTEGER FUNCTION get_year_at_experiment_start()

    TYPE(t_datetime), POINTER :: temp_date

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_year_at_experiment_start'

    IF (.NOT. time_initialized) &
      CALL finish(TRIM(routine), 'Time not initialized yet.')

    temp_date => get_time_experiment_start()
    get_year_at_experiment_start = get_year(temp_date)
    CALL deallocateDatetime(temp_date)

  END FUNCTION get_year_at_experiment_start

  ! ======================================================================================================= !
  !>
  !> Returns the time step in days (as fraction of day)
  !> e.g. two and a half days would yield 2.5
  REAL(wp) FUNCTION timestep_in_days(dt, model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), OPTIONAL, INTENT(in) :: dt        !< Timestep in seconds
    INTEGER,  OPTIONAL, INTENT(in) :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: dt_loc
    CHARACTER(len=*), PARAMETER :: routine = modname//':timestep_in_days'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (PRESENT(dt)) THEN
      dt_loc = dt
    ELSE IF (PRESENT(model_id)) THEN
      dt_loc = get_time_dt(model_id, jsbach_runs_standalone())
    ELSE
      CALL finish(routine, 'Must provide either dt oder model_id')
    END IF

    timestep_in_days = dt_loc / 86400._wp ! dt by number of seconds in a day

  END FUNCTION timestep_in_days

  ! ======================================================================================================= !
  !>
  !> Returns the number of time steps per day
  INTEGER FUNCTION timesteps_per_day(dt, model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), OPTIONAL, INTENT(in) :: dt        !< Timestep in seconds
    INTEGER,  OPTIONAL, INTENT(in) :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: dt_loc
    CHARACTER(len=*), PARAMETER :: routine = modname//':timesteps_per_day'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (PRESENT(dt)) THEN
      dt_loc = dt
    ELSE IF (PRESENT(model_id)) THEN
      dt_loc = get_time_dt(model_id, jsbach_runs_standalone())
    ELSE
      CALL finish(routine, 'Must provide either dt oder model_id')
    END IF

    timesteps_per_day = NINT(86400._wp/dt_loc) ! Assert that day length is an integer multiple of the time step
    IF (ABS(timesteps_per_day * dt_loc - 86400._wp) > EPSILON(1._wp)) &
      & CALL finish(routine, 'Time step must fit into day length')

  END FUNCTION timesteps_per_day

#endif
END MODULE mo_jsb_time
