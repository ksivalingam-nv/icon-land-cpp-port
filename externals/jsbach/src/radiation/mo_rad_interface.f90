!> Contains the interfaces to the radiation processes
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/radiation/mo_rad_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_rad_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,         ONLY: debug_on, acc_stream
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish

  USE mo_jsb_control,         ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, RAD_, SEB_, HYDRO_

  ! Use process configurations
  dsl4jsb_Use_config(RAD_)
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  ! Use of processes in this module
  dsl4jsb_Use_processes PHENO_

  ! Use of process memories
  dsl4jsb_Use_memory(PHENO_)
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_rad_tasks

  !> Type definition for radiation_surface task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_radiation
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_radiation
  END TYPE tsk_surface_radiation

  INTERFACE tsk_surface_radiation
    PROCEDURE Create_task_surface_radiation
  END INTERFACE tsk_surface_radiation

  !> Type definition for albedo task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_albedo
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_albedo
    PROCEDURE, NOPASS :: Aggregate => aggregate_albedo
  END TYPE tsk_albedo

  INTERFACE tsk_albedo
    PROCEDURE Create_task_albedo
  END INTERFACE tsk_albedo

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_radiation task
  !>
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_surface_radiation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_radiation::return_ptr)
    CALL return_ptr%Construct(name='surface_radiation', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for albedo task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "albedo"
  !!
  FUNCTION Create_task_albedo(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_albedo::return_ptr)
    CALL return_ptr%Construct(name='albedo', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_albedo

  ! ======================================================================================================= !
  !> Register tasks for radiation process
  !>
  SUBROUTINE Register_rad_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_surface_radiation       (model_id))
    CALL this%Register_task(tsk_albedo                  (model_id))

  END SUBROUTINE Register_rad_tasks

  ! ================================================================================================================================
  !>
  !> Implementation to update the surface radiation
  !! Task "surface_radiation" calculates the net radiations at the surface.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_radiation(tile, options)

    USE mo_rad_process,             ONLY: calc_radiation_surface_net
    USE mo_jsb_physical_constants,  ONLY: zemiss_def, stbo

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    ! Declare pointers for variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & swvis_srf_net,        &
      & swnir_srf_net,        &
      & sw_srf_net,           &
      & lw_srf_net,           &
      & rad_srf_net,          &
      & lw_srf_down,          &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
      & alb_vis_lnd,          &
      & alb_nir_lnd,          &
      & fract_lice,           &
      & rad_net_lwtr,         &
      & rad_net_lice,         &
      & albedo_lwtr,          &
      & albedo_lice,          &
      & t,                    &
      & t_old,                &
      & t_lwtr,               &
      & t_lice

    REAL(wp) ::                  &
      & lw_net_lwtr(options%nc), &
      & lw_net_lice(options%nc), &
      & sw_net_lwtr(options%nc), &
      & sw_net_lice(options%nc)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: iblk , ics, ice, nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_radiation'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(RAD_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)    ! in

    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, lw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)    ! out

    ! ---------------------------
    ! Go

    !> 1.0 lakes
    !!
    IF (tile%is_lake) THEN

      ! Water
      dsl4jsb_Get_var2D_onChunk(SEB_, t_lwtr)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lwtr)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lwtr) ! out

      !$ACC DATA CREATE(sw_net_lwtr, lw_net_lwtr, sw_net_lice, lw_net_lice) ASYNC(acc_stream)

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),           & ! in
        & swnir_srf_down(:),           & ! in
        & albedo_lwtr(:),              & ! in
        & albedo_lwtr(:),              & ! in
        & lw_srf_down(:),              & ! in
        & t_lwtr(:),                   & ! in
        & rad_net_lwtr(:),             & ! out
        & sw_net=sw_net_lwtr(:),       & ! out
        & lw_net=lw_net_lwtr(:)        & ! out
        & )

      ! Ice on lake
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(SEB_, t_lice)       ! in
        dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice)   ! in
        dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lice)  ! in
        dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lice) ! out

        CALL calc_radiation_surface_net( &
          & swvis_srf_down(:),           & ! in
          & swnir_srf_down(:),           & ! in
          & albedo_lice(:),              & ! in
          & albedo_lice(:),              & ! in
          & lw_srf_down(:),              & ! in
          & t_lice(:),                   & ! in
          & rad_net_lice(:),             & ! out
          & sw_net=sw_net_lice(:),       & ! out
          & lw_net=lw_net_lice(:)        & ! out
          & )

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          rad_srf_net(ic) = (1._wp - fract_lice(ic)) * rad_net_lwtr(ic) + fract_lice(ic) * rad_net_lice(ic)
          sw_srf_net (ic) = (1._wp - fract_lice(ic)) * sw_net_lwtr(ic)  + fract_lice(ic) * sw_net_lice(ic)
          lw_srf_net (ic) = (1._wp - fract_lice(ic)) * lw_net_lwtr(ic)  + fract_lice(ic) * lw_net_lice(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          rad_srf_net(ic) = rad_net_lwtr(ic)
          sw_srf_net (ic) = sw_net_lwtr (ic)
          lw_srf_net (ic) = lw_net_lwtr (ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC END DATA

    !> 2.0 all but lakes
    !!
    !!
    ELSE
      dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_lnd)    ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_lnd)    ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      swvis_srf_net)  ! out
      dsl4jsb_Get_var2D_onChunk(RAD_,      swnir_srf_net)  ! out

      CALL calc_radiation_surface_net(  &
        & swvis_srf_down(:),            & ! in
        & swnir_srf_down(:),            & ! in
        & alb_vis_lnd(:),               & ! in
        & alb_nir_lnd(:),               & ! in
        & lw_srf_down(:),               & ! in
        & t(:),                         & ! in
        & rad_srf_net(:),               & ! out
        & swvis_net = swvis_srf_net(:), & ! out
        & swnir_net = swnir_srf_net(:), & ! out
        & sw_net    = sw_srf_net(:),    & ! out
        & lw_net    = lw_srf_net(:)     & ! out
        & )

    END IF  ! NOT LAKES

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_radiation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_radiation(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(RAD_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_radiation'

    INTEGER :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(RAD_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(RAD_, rad_srf_net,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, sw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, swvis_srf_net,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, swnir_srf_net,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, lw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lwtr,   weighted_by_fract)

    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lice,   weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_radiation

  ! ================================================================================================================================
  !>
  !> Implementation to update the albedo
  !! Task "albedo" calculates several albedos.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_albedo(tile, options)

    ! Use declarations
    USE mo_rad_process,         ONLY: &
      & get_surface_albedo_simple, calc_pond_albedo, calc_alb_lwtr, calc_alb_lice, &
      & calc_glacier_albedo, calc_snow_albedo, &
      & Has_minimal_radiation, calc_snow_age
    USE mo_jsb_time,            ONLY: is_time_experiment_start, is_time_ltrig_rad_m1
    USE mo_jsb_math_constants,  ONLY: eps8, eps1, seconds_per_day
    USE mo_jsb_control,         ONLY: jsbach_runs_standalone
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER      :: model
    INTEGER                         :: iblk, ic, is
    INTEGER                         :: ics, ice, nc
    INTEGER                         :: model_scheme
    LOGICAL, ALLOCATABLE            :: l_day(:)
    REAL(wp)                        :: dtime                     !< timestep
    REAL(wp), ALLOCATABLE           :: temperature_surface(:)    !< use for jsbach & quincy, as helper variable
    CHARACTER(len=*), PARAMETER     :: routine = modname//':update_albedo'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(PHENO_)
    ! Declare configuration variables
    REAL(wp) ::            &
      & albedo_age_weight       ! 0 < albedo_age_weight < 1: snow albedo is calculated by linearly weighting
                                !                            the snow albedo resulting from both schemes

    LOGICAL ::                     &
      & use_alb_veg_simple

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  &
      & cos_zenith_angle,      &
      & t_unfilt,              &
      & snow,                  &
      & alb_vis,               &
      & alb_nir,               &
      & alb_vis_lnd,           &
      & alb_nir_lnd,           &
      & alb_vis_soil,          &
      & alb_nir_soil,          &
      & alb_vis_soil_and_snow, &
      & alb_nir_soil_and_snow, &
      & alb_vis_snow,          &
      & alb_nir_snow,          &
      & alb_vis_pond,          &
      & alb_nir_pond,          &
      & alb_background,        &
      & snow_age,              &
      & lai,                   &
      & fract_forest,          &
      & fract_snow_can,        &
      & fract_snow_soil,       &
      & fract_pond,            &
      & wtr_pond,              &
      & ice_pond,              &
      & swvis_srf_down,        &
      & swnir_srf_down,        &
      & fract_par_diffuse,     &
      & sai
    ! Lakes
    dsl4jsb_Real2D_onChunk ::  &
      & albedo_lwtr,           &
      & fract_lice,            &
      & t_lice,                &
      & weq_snow_lice,         &
      & albedo_lice

    dsl4jsb_Real2D_onChunk :: &
      & t

    LOGICAL :: &
      & l_ponds_param

    ! Set local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set model
    model => Get_model(tile%owner_model_id)
    model_scheme = model%config%model_scheme

    ! Update albedo only at experiment start, if model runs standalone or one time step before the radiation time step
    IF (.NOT. is_time_experiment_start(options%current_datetime) .AND. &
      & .NOT. is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id, jsbach_runs_standalone()) .AND. &
      & .NOT. jsbach_runs_standalone()) THEN
      IF (debug_on() .AND. iblk == 1) THEN
        WRITE (message_text,*) 'Update albedo only at experiment start, ', &
          &                    'if model runs standalone or one time step before the radiation time step ', &
          &                    'RETURN'
        CALL message(TRIM(routine), message_text)
      END IF
      RETURN
    END IF

    ! Enable for debugging:
    ! IF (iblk == 1 .AND. .NOT. ASSOCIATED(tile%next_sibling_tile)) THEN
    !   IF (is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id, jsbach_runs_standalone()) .OR. &
    !     & is_time_experiment_start(options%current_datetime)) THEN
    !     CALL message('update_albedo ', 'compute surface albedo')
    !   END IF
    ! END IF

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    IF (tile%contains_vegetation) THEN
      SELECT CASE (model_scheme)
      CASE (MODEL_JSBACH)
        dsl4jsb_Get_memory(PHENO_)
      END SELECT
    END IF

    ! Set pointers to variables in config
    use_alb_veg_simple        = dsl4jsb_Config(RAD_) %use_alb_veg_simple
    albedo_age_weight         = dsl4jsb_Config(RAD_) %albedo_age_weight

    l_ponds_param             = dsl4jsb_Config(HYDRO_)%l_ponds

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   cos_zenith_angle)      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   swvis_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   swnir_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   fract_par_diffuse)     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)                  ! in

    ALLOCATE(temperature_surface(nc))
    !$ACC ENTER DATA CREATE(temperature_surface) ASYNC(acc_stream)
    SELECT CASE (model_scheme)
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)              ! in
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1,nc
        temperature_surface(ic) = t_unfilt(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END SELECT

    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_background)        ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_snow)          ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_snow)          ! inout
    IF (.NOT. tile%is_lake) THEN
      dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_lnd)         ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_lnd)         ! inout
      IF (.NOT. tile%is_glacier) THEN
        dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_pond)      ! inout
        dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_pond)      ! inout
      END IF
    END IF

    ! ---------------------------
    ! Go

    !> 0.9  Set day-or-night-switch
    !!
    ALLOCATE(l_day(nc))
    !$ACC ENTER DATA CREATE(l_day) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      l_day(ic) = Has_minimal_radiation(swvis_srf_down(ic), swnir_srf_down(ic))
    END DO
    !$ACC END PARALLEL LOOP

    !> 1.0 snow and pont albedo
    !!
    !! Calculate snow albedos. Used later only for bare and vegetation tile types
    !! Note, for vegetation tiles this snow albedos are used only for snow on soil
    IF (tile%is_bare .OR. tile%is_vegetation) THEN

      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)  ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_pond)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)         ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,   snow_age)         ! in

      ! Update (non-dimensional) snow age
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_snow_age(          &
          & dtime,                   & ! in
          & temperature_surface(ic), & ! in
          & snow(ic),                & ! in
          & fract_snow_soil(ic),     & ! in
          & snow_age(ic)             & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_snow_albedo(    &
          & l_day(ic),            & ! in
          & fract_snow_soil(ic),  & ! in
          & snow_age(ic),         & ! in
          & cos_zenith_angle(ic), & ! in
          & albedo_age_weight,    & ! in
          & temperature_surface(ic), & ! in
          & alb_vis_snow(ic),     & ! inout
          & alb_nir_snow(ic)      & ! inout
          & )

        IF (l_ponds_param) THEN
          CALL calc_pond_albedo( &
            & l_day(ic),            & ! in
            & cos_zenith_angle(ic), & ! in
            & temperature_surface(ic), & ! in
            & wtr_pond(ic),         & ! in
            & ice_pond(ic),         & ! in
            & alb_vis_pond(ic),     & ! inout
            & alb_nir_pond(ic)      & ! inout
            & )
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END IF  ! IF tile%is_bare .OR. tile%is_vegetation

    !> 2.0 lakes & lakes with ice
    !!
    !!
    IF (tile%is_lake) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,   t_lice)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_lwtr)  ! inout

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_alb_lwtr( &
          & l_day(ic),       & ! in
          & albedo_lwtr(ic)  & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

      ! Ice on lake
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_lice)  ! in
        dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)     ! in
        dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_lice)    ! inout

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          CALL calc_alb_lice(    &
            & l_day(ic),         & ! in
            & t_lice(ic),        & ! in
            & weq_snow_lice(ic), & ! in
            & albedo_lice(ic)    & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          IF(l_day(ic)) THEN
            alb_vis(ic) = (1._wp - fract_lice(ic)) * albedo_lwtr(ic) + fract_lice(ic) * albedo_lice(ic)
            alb_nir(ic) = alb_vis(ic)  ! TODO
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          IF(l_day(ic)) THEN
            alb_vis(ic) = albedo_lwtr(ic)
            alb_nir(ic) = alb_vis(ic)  ! TODO
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    !> 3.0 glacier
    !!
    ELSE IF (tile%is_glacier) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_glacier_albedo(  &
          & l_day(ic),             & ! in
          & temperature_surface(ic), & ! in
          & alb_vis_lnd(ic),       & ! inout
          & alb_nir_lnd(ic)        & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

    !> 4.0 soil
    !! currently jsbach4 does not consider "bare tiles"
    !!
    ELSE IF (tile%is_bare) THEN
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_soil)          ! inout

      ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from bc_land_phy.nc file
      ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
      ! the canopy above.
      IF (debug_on() .AND. iblk == 1) THEN
        CALL message(TRIM(routine), 'Soil albedo scheme not used for bare soil tiles but taken from bc_land_phy.nc file')
      END IF

      ! Merge the albedo of soil with the snow on it for the overall bare soil tile albedo
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        alb_vis_lnd(ic) = fract_snow_soil(ic)                                       * alb_vis_snow(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * fract_pond(ic)            * alb_vis_pond(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * (1.0_wp - fract_pond(ic)) * alb_vis_soil(ic)

        alb_nir_lnd(ic) = fract_snow_soil(ic)                                       * alb_nir_snow(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * fract_pond(ic)            * alb_nir_pond(ic) &
                     & + (1.0_wp - fract_snow_soil(ic)) * (1.0_wp - fract_pond(ic)) * alb_nir_soil(ic)
      END DO
      !$ACC END PARALLEL LOOP

    !> 5.0 vegetation / canopy
    !!
    ELSE IF (tile%is_vegetation) THEN

      ! TODO: correct indentation >>>
    SELECT CASE (model_scheme)
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)                   ! in
      IF (use_alb_veg_simple) THEN
        dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_forest)          ! in
      ENDIF
    END SELECT
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)       ! in

      !> 5.1  simple calculation of alb_vis_lnd / alb_nir_lnd
      !!
      !! TRUE for jsbach_lite, FALSE for jsbach_pfts and quincy use cases
      IF (use_alb_veg_simple) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          alb_vis_lnd(ic) = get_surface_albedo_simple( &
            & t_unfilt(ic),                            &
            & fract_forest(ic),                        &
            & fract_snow_soil(ic),                     &
            & fract_snow_can(ic),                      &
            & fract_pond(ic),                          &
            & lai(ic),                                 &
            & alb_background(ic)                       & ! Note, alb_background <-- bc_land_phy.nc
                        )
          alb_nir_lnd(ic) = alb_vis_lnd(ic)
        END DO
        !$ACC END PARALLEL LOOP
      !> 5.2  detailed calculation of alb_vis_lnd & alb_nir_lnd
      !!
      !! merging albedo of snow, soil and canopy
      ELSE

        !> 5.2.1 Get bare soil albedo for vegetation tiles
        ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from the bc_land_phy.nc file
        ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
        ! the canopy above.

        IF (debug_on() .AND. iblk == 1) THEN
          CALL message(TRIM(routine), 'Soil albedo scheme not used for veg tiles but taken from bc_land_phy.nc file')
        END IF

        SELECT CASE (model_scheme)
        ! TODO: adapt indent for jsbach case
        CASE (MODEL_JSBACH)
        CALL finish(TRIM(routine), 'lctlib not available on tile')
        END SELECT
      END IF      ! NOT use_alb_veg_simple
    END IF      ! IF (tile%is_vegetation)

    !> 6.0 pass the "land albedo" to the "tile albedo"
    !! this is needed because for lake tiles, the "tile albedo" is calculated directly
    !! but for vegetation, bare and glacier, the "land albedo" is calculated
    !!

    IF (tile%is_bare .OR. tile%is_glacier .OR. tile%is_vegetation) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        ! for lake tile, alb_vis/nur are already set above to albedo_lwtr
        alb_vis(ic) = alb_vis_lnd(ic)
        alb_nir(ic) = alb_nir_lnd(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC EXIT DATA DELETE(l_day, temperature_surface) ASYNC(acc_stream)
    DEALLOCATE(l_day, temperature_surface)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_albedo

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "albedo"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_albedo(tile, options)

    USE mo_jsb_time,          ONLY: is_time_experiment_start, is_time_ltrig_rad_m1
    USE mo_jsb_control,       ONLY: jsbach_runs_standalone
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_JSBACH
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(RAD_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    TYPE(t_jsb_model),       POINTER :: model
    INTEGER                          :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_albedo'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Get_memory(RAD_)

    CALL weighted_by_fract%BeginAggregate()

    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
    ! Update albedo only at experiment start, if model runs standalone or one time step before the radiation time step
    IF (.NOT. is_time_experiment_start(options%current_datetime) .AND. &
      & .NOT. is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id, jsbach_runs_standalone()) .AND. &
      & .NOT. jsbach_runs_standalone()) RETURN

      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_pond,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_pond,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, snow_age,            weighted_by_fract)
      IF (.NOT. tile%is_vegetation) THEN
        dsl4jsb_Aggregate_onChunk(RAD_, albedo_lwtr,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(RAD_, albedo_lice,         weighted_by_fract)
      ENDIF

    END SELECT

    CALL weighted_by_fract%EndAggregate()


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_albedo

#endif
END MODULE mo_rad_interface
