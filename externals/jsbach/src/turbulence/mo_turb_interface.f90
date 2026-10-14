!> Contains the interfaces to the turbulence process
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
!NEC$ options "-finline-file=externals/jsbach/src/soil_snow_energy/mo_sse_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"

MODULE mo_turb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,         ONLY: debug_on, jsbach_runs_standalone, acc_stream
  USE mo_kind,                ONLY: wp
  USE mo_jsb_math_constants,  ONLY: eps8
  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model, MODEL_JSBACH
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes TURB_, A2L_, SEB_, HYDRO_, SSE_

  ! Use process configurations
  dsl4jsb_Use_config(TURB_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SEB_)

  ! Use process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  ! Use of processes in this module
  dsl4jsb_Use_processes PHENO_

  ! Use process memories
  dsl4jsb_Use_memory(PHENO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_turb_tasks !,t_turb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_interface'

  !> Type definition for surface_roughness task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_roughness
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_roughness     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_roughness  !< Aggregates computed task variables
  END TYPE tsk_surface_roughness

  !> Constructor interface for surface_roughness task
  INTERFACE tsk_surface_roughness
    PROCEDURE Create_task_surface_roughness                        !< Constructor function for task
  END INTERFACE tsk_surface_roughness

  !> Type definition for humidity_scaling task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_humidity_scaling
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_humidity_scaling     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_humidity_scaling  !< Aggregates computed task variables
  END TYPE tsk_humidity_scaling

  !> Constructor interface for humidity scaling task
  INTERFACE tsk_humidity_scaling
    PROCEDURE Create_task_humidity_scaling                        !< Constructor function for task
  END INTERFACE tsk_humidity_scaling

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_roughness task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_surface_roughness"
  !!
  FUNCTION Create_task_surface_roughness(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_roughness::return_ptr)
    CALL return_ptr%Construct(name='surface_roughness', process_id=TURB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_roughness

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for humidity_scaling task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_humidity_scaling"
  !!
  FUNCTION Create_task_humidity_scaling(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_humidity_scaling::return_ptr)
    CALL return_ptr%Construct(name='humidity_scaling', process_id=TURB_, owner_model_id=model_id)

  END FUNCTION Create_task_humidity_scaling

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for turbulence process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of TURB process class
  !!
  SUBROUTINE Register_turb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_roughness(model_id))
    CALL this%Register_task(tsk_humidity_scaling(model_id))

  END SUBROUTINE Register_turb_Tasks

  ! ======================================================================================================= !
  !>
  !> Implementation of "update" for task "surface_roughness"
  !>
  !> the variables rough_m_star & rough_h_star are aggregated, and
  !> in the aggregate_roughness() rough_m and rough_h are calculated from these var
  !>
  SUBROUTINE update_roughness(tile, options)

    USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_config(SEB_)

    ! Declare process memories
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: rough_m
    dsl4jsb_Real2D_onChunk :: rough_h
    dsl4jsb_Real2D_onChunk :: rough_m_star
    dsl4jsb_Real2D_onChunk :: rough_h_star
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_lice

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp), DIMENSION(options%nc) :: &
      & rough_bare_star, rough_snow_star
    REAL(wp) :: &
      & config_rough_bare, config_rough_snow, config_rough_water, config_blending_height, &
      & config_roughness_momentum_to_heat
    LOGICAL                         :: l_veg_bare_glac_tile      !< use for jsbach & quincy, as helper variable
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_roughness'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(TURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(TURB_)
    dsl4jsb_Get_config(SEB_)

    !$ACC DATA CREATE(rough_bare_star, rough_snow_star) ASYNC(acc_stream)

    ! Get process memories
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m)         ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h)         ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m_star)    ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h_star)    ! out

    config_blending_height = dsl4jsb_Config(TURB_)%blending_height
    config_rough_snow      = dsl4jsb_Config(TURB_)%roughness_snow
    config_rough_water     = dsl4jsb_Config(TURB_)%roughness_water


    ! set logical (vegetated, bare and glacier tiles) model specific to TRUE/FALSE
    l_veg_bare_glac_tile = .FALSE.
    ! icon-land
    IF (tile%Has_process_memory(SSE_)) THEN  ! Should catch all vegetated, bare and glacier tiles
      l_veg_bare_glac_tile = .TRUE.
    END IF

    ! all vegetated, bare and glacier tiles
    IF (l_veg_bare_glac_tile) THEN

      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil) ! in

      config_rough_bare                 = dsl4jsb_Config(TURB_)%roughness_bare
      config_roughness_momentum_to_heat = dsl4jsb_Config(TURB_)%roughness_momentum_to_heat

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        rough_bare_star(ic) = 1._wp  / LOG(config_blending_height / config_rough_bare) ** 2
        rough_snow_star(ic) = 1._wp  / LOG(config_blending_height / config_rough_snow) ** 2
      END DO
      !$ACC END PARALLEL LOOP

      ! Vegetated tile without PFT, e.g., jsbach_lite
      IF (tile%is_vegetation) THEN

        ! rough_m and rough_m_star are taken from values from init
        !
        ! We don't have PFTs so use rough_m from ini file
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          rough_h_star(ic) = (1._wp - fract_snow_soil(ic)) * 1._wp / LOG(config_blending_height / MIN(1._wp, rough_m(ic)))**2 &
            &               +        fract_snow_soil(ic)  * rough_snow_star(ic)
        END DO
        !$ACC END PARALLEL LOOP

      ELSE  ! bare or glacier tile
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          rough_m_star(ic) = rough_bare_star(ic)
          ! Modify roughness length heat for bare part by snow cover
          rough_h_star(ic) = (1._wp - fract_snow_soil(ic)) * rough_bare_star(ic) + fract_snow_soil(ic) * rough_snow_star(ic)
        END DO
        !$ACC END PARALLEL LOOP

      END IF  ! vegetated or bare tile

      IF (.NOT. tile%is_vegetation) THEN ! bare/glacier
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          rough_m(ic) = config_blending_height * EXP(-1._wp / SQRT(rough_m_star(ic)))
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        rough_h(ic) = config_blending_height * EXP(-1._wp / SQRT(rough_h_star(ic))) / config_roughness_momentum_to_heat
      END DO
      !$ACC END PARALLEL LOOP

    ELSE IF (tile%contains_lake) THEN

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_memory(SEB_)
        dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice) ! IN

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          rough_m_star(ic) =            fract_lice(ic)  / (LOG(config_blending_height / config_rough_snow)  ** 2)  &
            &               + (1._wp - fract_lice(ic)) / (LOG(config_blending_height / config_rough_water) ** 2)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          rough_m_star(ic) = 1._wp / (LOG(config_blending_height / config_rough_water) ** 2)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        rough_h_star(ic) = rough_m_star(ic)
        rough_m(ic) = config_blending_height * EXP(-1._wp / SQRT(rough_m_star(ic)))
        rough_h(ic) = rough_m(ic)
      END DO
      !$ACC END PARALLEL LOOP

    ELSE
      CALL finish(TRIM(routine), 'Don"t know what to do.')
    END IF

    ! Limit roughness length for heat by 1 m (see below Eq.5.6 in ECHAM5 manual)
    ! $noACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    ! DO ic = 1, nc
    !   rough_h(ic) = MIN(1._wp, blending_height * EXP(-1._wp / SQRT(rough_h(ic))))
    ! END DO
    ! $noACC END PARALLEL LOOP

    !$ACC END DATA

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_roughness

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_roughness"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_roughness(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_memory(TURB_)

    dsl4jsb_Real2D_onChunk :: rough_m
    dsl4jsb_Real2D_onChunk :: rough_h
    dsl4jsb_Real2D_onChunk :: rough_m_star
    dsl4jsb_Real2D_onChunk :: rough_h_star

    REAL(wp) :: config_blending_height, config_roughness_momentum_to_heat

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_roughness'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc, ic

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Get_memory(TURB_)

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(TURB_, rough_m_star,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, rough_h_star,   weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(TURB_)

    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m)
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h)
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m_star)
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h_star)

    config_blending_height            = dsl4jsb_Config(TURB_)%blending_height
    config_roughness_momentum_to_heat = dsl4jsb_Config(TURB_)%roughness_momentum_to_heat

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      rough_m(ic) = config_blending_height * EXP(-1._wp / SQRT(rough_m_star(ic)))
      rough_h(ic) = config_blending_height * EXP(-1._wp / SQRT(rough_h_star(ic))) / config_roughness_momentum_to_heat
    END DO
    !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_roughness

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "humidity_scaling"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_humidity_scaling(tile, options)

    !USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_sse_process,            ONLY: relative_humidity_soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(HYDRO_)

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::   &
      & fract_fpc,              &
      & fract_snow,             &
      & fract_wet,              &
      & wtr_rootzone_avail_max, &
      & wtr_rootzone_avail,     &
      & fact_q_air,             &
      & fact_qsat_srf,          &
      & fact_qsat_trans_srf,    &
      & qsat_star,              &
      & canopy_cond_limited,    &
      & pch,                    &
      & wind_air,               &
      & q_air

    dsl4jsb_Real3D_onChunk :: &
      & wtr_soil_sl,          &
      & soil_depth_sl,        &
      & vol_porosity_sl,      &
      & vol_wres_sl

    ! Local variables
    REAL(wp) ::                    &
      fact_qsat_veg(options%nc),   &
      fact_qsat_trans(options%nc), &
      fact_qsat_soil(options%nc),  &
      fact_qair_veg(options%nc),   &
      fact_qair_soil(options%nc),  &
      zfract_fpc(options%nc),      &
      w_soil_wilt_fract_config
    REAL(wp) ::                    &
      rel_hum,                     &
      weq_wsat_sl1,                &
      weq_wres_sl1
    REAL(wp) ::                    &
      ch_tmp(options%nc)                      !< local variable for pch / kh from jsbach, quincy, A2L_, TURB_
    LOGICAL :: &
      mask(options%nc)
    REAL(wp), ALLOCATABLE      :: canopy_conductance(:)       !< use for jsbach & quincy, as helper variable
    LOGICAL                    :: jsb_standalone              !< model runs standalone?

    INTEGER :: iblk, ics, ice, nc, ic
    INTEGER :: model_scheme

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_humidity_scaling'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(TURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    model_scheme = model%config%model_scheme
    jsb_standalone = jsbach_runs_standalone()

    dsl4jsb_Get_config(HYDRO_)

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(fact_qsat_veg, fact_qsat_trans, fact_qsat_soil) &
    !$ACC   CREATE(fact_qair_veg, fact_qair_soil, zfract_fpc, mask, ch_tmp)

    ! Get pointers to process configs and memory
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    SELECT CASE (model_scheme)
    CASE (MODEL_JSBACH)
      IF (tile%contains_vegetation) THEN
        dsl4jsb_Get_memory(PHENO_)
      END IF
    END SELECT

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   pch)               ! in
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ch_tmp(ic) = pch(ic)
    END DO
    !$ACC END PARALLEL LOOP
    dsl4jsb_Get_var2D_onChunk(A2L_,   wind_air)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)          ! in
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_q_air)          ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_qsat_srf)       ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_qsat_trans_srf) ! OUT
    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,    qsat_star)              ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  fract_wet)              ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  wtr_soil_sl)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  wtr_rootzone_avail_max) ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  wtr_rootzone_avail)     ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  soil_depth_sl)          ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  vol_porosity_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  vol_wres_sl)            ! in
    END IF

    ! ----------------------------------------------------------------------------------------------------- !

    IF (tile%contains_vegetation) THEN
      ALLOCATE(canopy_conductance(nc))
      !$ACC ENTER DATA CREATE(canopy_conductance) ASYNC(acc_stream)
      SELECT CASE (model_scheme)
      CASE (MODEL_JSBACH)
        dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc)           ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_limited)         ! in
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1,nc
          canopy_conductance(ic) = canopy_cond_limited(ic)
        END DO
        !$ACC END PARALLEL LOOP
        w_soil_wilt_fract_config = dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          mask(ic) = wtr_rootzone_avail(ic) > w_soil_wilt_fract_config * wtr_rootzone_avail_max(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END SELECT

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF (mask(ic)) THEN
          fact_qsat_veg(ic) = fract_snow(ic) + (1._wp - fract_snow(ic)) *  &
            & (fract_wet(ic) + ( 1._wp - fract_wet(ic)) /                  &
            & (1._wp + ch_tmp(ic) * MAX(1._wp, wind_air(ic)) / MAX(1.E-20_wp, canopy_conductance(ic))))
          fact_qsat_trans(ic) = (1._wp - fract_snow(ic)) *  &
            & (1._wp - fract_wet(ic)) /                     &
            & (1._wp + ch_tmp(ic) * MAX(1._wp, wind_air(ic)) / MAX(1.E-20_wp, canopy_conductance(ic)))
        ELSE
          fact_qsat_veg(ic)   = fract_snow(ic) + (1._wp - fract_snow(ic)) * fract_wet(ic)
          fact_qsat_trans(ic) = 0._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        fact_qair_veg(ic) = fact_qsat_veg(ic)
        zfract_fpc(ic)    = fract_fpc(ic)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC EXIT DATA DELETE(canopy_conductance) ASYNC(acc_stream)
      DEALLOCATE(canopy_conductance)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        fact_qsat_veg(ic)   = 0._wp
        fact_qsat_trans(ic) = 0._wp
        fact_qair_veg(ic)   = 0._wp
        zfract_fpc(ic)      = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (tile%contains_soil) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) PRIVATE(rel_hum)
      DO ic = 1, nc
        weq_wres_sl1 = vol_wres_sl(ic,1) * soil_depth_sl(ic,1)
        weq_wsat_sl1 = vol_porosity_sl(ic,1) * soil_depth_sl(ic,1)
        rel_hum      = relative_humidity_soil(wtr_soil_sl(ic,1), weq_wres_sl1, weq_wsat_sl1)
        IF (qsat_star(ic) * rel_hum > q_air(ic) .AND. rel_hum > 1.e-10_wp) THEN
          fact_qsat_soil(ic) = fract_snow(ic) + (1._wp - fract_snow(ic)) * &
            & (fract_wet(ic) + (1._wp - fract_wet(ic)) * rel_hum)
          fact_qair_soil(ic) = 1._wp
        ELSE
          fact_qsat_soil(ic) = fract_snow(ic) + (1._wp - fract_snow(ic)) * fract_wet(ic)
          fact_qair_soil(ic) = fact_qsat_soil(ic)
        END IF
        IF (q_air(ic) > qsat_star(ic)) THEN
          fact_qsat_soil(ic) = 1._wp
          fact_qair_soil(ic) = 1._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP

    ELSE IF (tile%contains_glacier .OR. tile%is_lake) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        fact_qsat_soil(ic) = 1._wp
        fact_qair_soil(ic) = 1._wp
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        fact_qsat_soil(ic) = 0._wp
        fact_qair_soil(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      fact_qsat_srf(ic) = zfract_fpc(ic) * fact_qsat_veg(ic) + (1._wp - zfract_fpc(ic)) * fact_qsat_soil(ic)
      fact_q_air(ic)    = zfract_fpc(ic) * fact_qair_veg(ic) + (1._wp - zfract_fpc(ic)) * fact_qair_soil(ic)
      fact_qsat_trans_srf(ic) = zfract_fpc(ic) * fact_qsat_trans(ic)
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_humidity_scaling

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "humidity_scaling"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_humidity_scaling(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(TURB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_humidity_scaling'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(TURB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(TURB_, fact_q_air,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, fact_qsat_srf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, fact_qsat_trans_srf, weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_humidity_scaling

  ! ================================================================================================================================
  !>
  !> Another task for "<PROCESS_NAME_LOWER_CASE>"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !  ....

#endif
END MODULE mo_turb_interface
