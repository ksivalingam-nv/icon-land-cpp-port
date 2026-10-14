!> Contains the task interfaces for the soil and snow energy process.
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
!>#### Task interfaces for the SSE process
!>
!> The module defines the process tasks that compute snow and soil layer temperatures
!> and update soil and snow properties. The tasks are registered within the ICON-LAND
!> process framework [[mo_jsb_model_init:jsbach_init]].
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/soil_snow_energy/mo_sse_process.pp-jsb.f90"

MODULE mo_sse_interface
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_control,     ONLY: debug_on, acc_stream
  USE mo_jsb_time,        ONLY: is_time_experiment_start

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes, configuration and memory in this module
  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_sse_tasks

  !> Type definition for the soil-and-snow temperature task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_and_snow_temp
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_and_snow_temperature
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_and_snow_temperature
  END TYPE tsk_soil_and_snow_temp

  !> Constructor interface for the soil-and-snow temperature task
  INTERFACE tsk_soil_and_snow_temp
    PROCEDURE Create_task_soil_and_snow_temperature
  END INTERFACE tsk_soil_and_snow_temp

  !> Type definition for the soil-and-snow property-update task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_and_snow_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_and_snow_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_and_snow_properties
  END TYPE tsk_soil_and_snow_properties

  !> Constructor interface for the soil-and-snow property-update task
  INTERFACE tsk_soil_and_snow_properties
    PROCEDURE Create_task_soil_and_snow_properties
  END INTERFACE tsk_soil_and_snow_properties

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_interface'

CONTAINS

  ! ===================================================================================================== !
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Construct the task that advances the soil-and-snow temperature state.
  !>
  FUNCTION Create_task_soil_and_snow_temperature(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id    !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "soil and snow temperature"

    ALLOCATE(tsk_soil_and_snow_temp::return_ptr)
    CALL return_ptr%Construct(name='soil_and_snow_temperature', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_and_snow_temperature

  ! -------------------------------------------------------------------------------------------------------
  !> Construct the task that updates the soil-and-snow properties.
  !>
  FUNCTION Create_task_soil_and_snow_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id    !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "soil and snow properties"

    ALLOCATE(tsk_soil_and_snow_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_and_snow_properties', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_and_snow_properties

  ! ===================================================================================================== !
  !>
  !>#### Register the SSE tasks within the ICON-LAND process queue
  !>
  !> Register the temperature and soil property update tasks.
  !>
  SUBROUTINE Register_sse_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process       !< Process
    INTEGER,              INTENT(in)    :: model_id      !< Model ID

    CALL process%Register_task(tsk_soil_and_snow_temp(model_id))
    CALL process%Register_task(tsk_soil_and_snow_properties(model_id))

  END SUBROUTINE Register_sse_tasks

  ! ===================================================================================================== !
  !>
  !>#### Calculate new snow and soil temperature states
  !>
  !> Using the current surface temperature, a zero flux condition at the bottom of the soil column and
  !> last timesteps Richtmeyer-Morton coefficients, we compute the new snow and soil layer temperatures,
  !> phase changes within the soil column as well as the ground heat flux and capacity.
  !>
  SUBROUTINE update_soil_and_snow_temperature(tile, options)

    USE mo_sse_process,            ONLY: calc_soil_temperature, calc_snow_temperature, calc_snow_abcoeff
    USE mo_jsb_physical_constants, ONLY: rhoh2o, dens_snow
    USE mo_sse_constants,          ONLY: snow_depth_min
    USE mo_jsb_time,               ONLY: is_newyear

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    TYPE(t_jsb_model), POINTER :: model                    !< Current instance of the model

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(SSE_)     !< Configurable soil and snow energy parameters
    dsl4jsb_Def_config(HYDRO_)   !< Configurable hydrology parameters
    dsl4jsb_Def_memory(SSE_)     !< Memory of the soil and snow energy process
    dsl4jsb_Def_memory(SEB_)     !< Memory of the surface energy balance process
    dsl4jsb_Def_memory(HYDRO_)   !< Memory of the hydrology process

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_srf                 !< Surface temperature [K]
    dsl4jsb_Real2D_onChunk :: snowmelt              !< Snowmelt rate [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: t_soil_sl             !< Soil layer temperature [K]
    dsl4jsb_Real3D_onChunk :: t_soil_acoef          !< Soil A coefficient of the Richtmeyer-Morton scheme
    dsl4jsb_Real3D_onChunk :: t_soil_bcoef          !< Soil B coefficient of the Richtmeyer-Morton scheme
    dsl4jsb_Real3D_onChunk :: vol_heat_cap_sl       !< Volumetric heat capacity of soil layers [J m-3 K-1]
    dsl4jsb_Real3D_onChunk :: heat_cond_sl          !< Thermal conductivity of soil layers [J m-1 s-1 K-1]
    dsl4jsb_Real3D_onChunk :: t_snow                !< Snow layer temperature [K]
    dsl4jsb_Real3D_onChunk :: t_snow_acoef          !< Snow A coefficient of the Richtmeyer-Morton scheme
    dsl4jsb_Real3D_onChunk :: t_snow_bcoef          !< Snow B coefficient of the Richtmeyer-Morton scheme
    dsl4jsb_Real2D_onChunk :: vol_heat_cap_snow     !< Volumetric heat capacity of snow [J m-3 K-1]
    dsl4jsb_Real2D_onChunk :: heat_cond_snow        !< Thermal conductivity of snow [J m-1 s-1 K-1]
    dsl4jsb_Real2D_onChunk :: hcap_grnd             !< Heat capacity of the uppermost surface layer [J m-2 K-1]
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old         !< Previous heat capacity of the uppermost surface layer [J m-2 K-1]
    dsl4jsb_Real2D_onChunk :: grnd_hflx             !< Ground heat flux [J m-2 s-1]
    dsl4jsb_Real2D_onChunk :: grnd_hflx_old         !< Previous ground heat flux [J m-2 s-1]
    dsl4jsb_Real3D_onChunk :: snow_depth_sl         !< Depth of snow layers [m]
    dsl4jsb_Real3D_onChunk :: soil_depth_sl         !< Depth of soil layers [m]
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl           !< Liquid water content of soil layers [m]
    dsl4jsb_Real3D_onChunk :: ice_soil_sl           !< Ice content of soil layers [m]
    dsl4jsb_Real3D_onChunk :: wtr_freeze_sl         !< Water frozen in soil layers [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: ice_melt_sl           !< Ice melted in soil layers [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl !< Potentially supercooled water in soil layers [m]
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl       !< Volumetric porosity of soil layers []
    dsl4jsb_Real3D_onChunk :: matric_pot_sl         !< Soil matric potential [m]
    dsl4jsb_Real3D_onChunk :: bclapp_sl             !< Clapp-Hornberger parameter of soil layers []
    dsl4jsb_Real2D_onChunk :: thaw_depth            !< Current thaw depth [m]
    dsl4jsb_Real2D_onChunk :: thaw_depth_max        !< Maximum thaw depth reached during this year [m]
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1    !< Maximum thaw depth of the previous year [m]
    dsl4jsb_Real2D_onChunk :: heat_cond             !< Effective thermal conductivity [J m-1 s-1 K-1]
    dsl4jsb_Real2D_onChunk :: vol_heat_cap          !< Effective volumetric heat capacity [J m-3 K-1]
    dsl4jsb_Real2D_onChunk :: weq_snow_soil         !< Snow water equivalent on the ground [m]
    dsl4jsb_Real2D_onChunk :: snow_soil_dens        !< Snow density [kg m-3]
    dsl4jsb_Real2D_onChunk :: fract_snow_soil       !< Fraction of the ground covered by snow []
    dsl4jsb_Real2D_onChunk :: le_pc_remain          !< Latent energy remaining for phase change [J m-2]

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & snow_depth,                    & !< Snow depth after conversion from water equivalent [m]
      & snow_fract,                    & !< Fraction of the ground covered by snow []
      & t_soil_top,                    & !< Temperature at the top of the soil column [K]
      & hcap_grnd_snow,                & !< Ground heat capacity contributed by snow [J m-2 K-1]
      & hcap_grnd_soil,                & !< Ground heat capacity contributed by soil [J m-2 K-1]
      & grnd_hflx_snow,                & !< Ground heat flux contributed by snow [J m-2 s-1]
      & grnd_hflx_soil                   !< Ground heat flux contributed by soil [J m-2 s-1]
    LOGICAL, DIMENSION(options%nc)  :: &
      & ldglac,                        &    !< Glacier-tile mask
      & is_snow_covered                     !< T: soil layer is snow covered
    REAL(wp), ALLOCATABLE :: ws_max_sl(:,:) !< Maximum soil water storage by layer [m]

    INTEGER, DIMENSION(options%nc) ::  &
      & itop, itop_old                   !< Index of the top snow layer at the current and previous time step

    LOGICAL  :: lstart, l_uniform_snow       !< Experiment start flag and uniform-snow configuration flag
    INTEGER  :: iblk, ics, ice, nc, ic       !< Block bounds and horizontal loop indices
    REAL(wp) :: dtime                        !< Time step length [s]

    TYPE(t_jsb_vgrid), POINTER :: soil_e, snow_e !< Soil and snow energy vertical grids
    INTEGER                    :: nsoil, is, nsnow, isnow !< Vertical grid sizes and layer indices

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_temperature'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(SSE_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    snow_e => Get_vgrid('snow_depth_energy')
    nsnow = snow_e%n_levels

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, le_pc_remain)    ! inout

    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)         ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_acoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_bcoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_snow)            ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_snow_acoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_snow_bcoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)   ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)      ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap_snow) ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond_snow)    ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, hcap_grnd)         ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, hcap_grnd_old)     ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, grnd_hflx)         ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, grnd_hflx_old)     ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, snow_depth_sl)     ! out
    IF (.NOT. tile%is_glacier) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, matric_pot_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_freeze_sl)    ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_melt_sl)      ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_scool_sl) ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth)          ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_ym1)  ! out

    dsl4jsb_Get_var2D_onChunk(SSE_,   heat_cond)
    dsl4jsb_Get_var2D_onChunk(SSE_,   vol_heat_cap)

    dsl4jsb_Get_var2D_onChunk(SEB_,   t_srf)           ! in

    lstart = is_time_experiment_start(options%current_datetime)
    l_uniform_snow = dsl4jsb_Config(SSE_)%l_uniform_snow

    !$ACC ENTER DATA ASYNC(acc_stream) &
    !$ACC   CREATE(snow_depth, snow_fract, t_soil_top, is_snow_covered) &
    !$ACC   CREATE(hcap_grnd_snow, hcap_grnd_soil, grnd_hflx_snow, grnd_hflx_soil, ldglac, itop, itop_old)

#ifdef _CRAYFTN
    ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
    !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
    DO ic = 1, nc
      ! Save old ground heat flux and heat capacity before updating
      ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
      ! given back to the atmosphere so that the surface energy balance can be computed there
      grnd_hflx_old(ic) = grnd_hflx(ic)
      hcap_grnd_old(ic) = hcap_grnd(ic)
    END DO
    !$ACC END PARALLEL LOOP

    IF (.NOT. tile%is_glacier) THEN

      ! Calculate snow depth in m (no longer in m water equivalent)
      IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
#ifdef _CRAYFTN
        ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
        !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
        DO ic = 1,nc
          snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / snow_soil_dens(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
#ifdef _CRAYFTN
        ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
        !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
        DO ic = 1,nc
          snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / dens_snow
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Calculate top snow layer for previous time step
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        itop_old(ic) = nsnow + 1 - COUNT(snow_depth_sl(ic,:) > 0._wp)
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate depths of snow layers
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
#endif
      DO isnow=1,nsnow
        DO ic = 1,nc
          snow_depth_sl(ic,isnow) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        snow_fract(ic) = fract_snow_soil(ic)
        IF (snow_depth(ic) > snow_depth_min .AND. snow_fract(ic) > EPSILON(1._wp)) THEN
          IF (.NOT. l_uniform_snow) THEN
            snow_depth(ic) = snow_depth(ic) / snow_fract(ic)
          END IF
        ELSE
          snow_depth(ic) = 0._wp
          snow_fract(ic) = 0._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL PRESENT(sse__mem, seb__mem, hydro__mem) ASYNC(acc_stream)
#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
#endif
      !$ACC LOOP SEQ
      DO isnow=nsnow,1,-1
        !$ACC LOOP GANG VECTOR
        DO ic = 1,nc
          IF (snow_depth(ic) > 0._wp) THEN
            snow_depth_sl(ic,isnow) = MIN(snow_depth(ic), snow_e%dz(isnow))
            snow_depth(ic) = snow_depth(ic) - snow_depth_sl(ic,isnow)
          END IF
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        snow_depth_sl(ic,nsnow) = snow_depth_sl(ic,nsnow) + snow_depth(ic)
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate top snow layer for current time step

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        itop(ic) = nsnow + 1 - COUNT(snow_depth_sl(ic,:) > 0._wp)
        is_snow_covered(ic) = (weq_snow_soil(ic) * rhoh2o / dens_snow > snow_depth_min)
      END DO
      !$ACC END PARALLEL LOOP

      CALL calc_snow_temperature( &
        ! in
        & nc,                     &
        & nsnow,                  &
        & lstart,                 &
        & itop_old           (:), &
        & itop               (:), &
        & t_srf              (:), &
        & t_soil_acoef     (:,1), &
        & t_soil_bcoef     (:,1), &
        ! inout
        & t_snow           (:,:), &
        & t_snow_acoef     (:,:), &
        & t_snow_bcoef     (:,:), &
        ! out
        & t_soil_top         (:)  &
        & )

      IF (.NOT. l_uniform_snow) THEN
#ifdef _CRAYFTN
        ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
        !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
        DO ic = 1, nc
          IF (snow_depth_sl(ic, nsnow) > 0._wp) THEN
            ! At least one snow layer above soil
            t_soil_top(ic) = snow_fract(ic) * t_soil_top(ic) + (1._wp - snow_fract(ic)) * t_srf(ic)
          ELSE
            ! No snow
            t_soil_top(ic) = t_srf(ic)
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    ELSE
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
#endif
      DO isnow = 1,nsnow
        DO ic = 1, nc
            snow_depth_sl(ic,isnow) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1, nc
        snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / dens_snow
        snow_fract(ic) = 0._wp
        t_soil_top(ic) = t_srf(ic)
        grnd_hflx_snow(ic) = 0._wp
        hcap_grnd_snow(ic) = 0._wp
        is_snow_covered(ic) = .FALSE.
      END DO
      !$ACC END PARALLEL LOOP

    END IF

    ! @todo: move water/ice calculations to HYDRO_ process
    IF (tile%is_glacier) THEN
      CALL calc_soil_temperature( &
        ! in
        & nc                                       , &
        & nsoil                                    , &
        & soil_e%dz(:)                             , &
        & soil_e%levels(:)                         , &
        & soil_e%ubounds(:)                        , &
        & dtime                                    , &
        & lstart                                   , &
        & dsl4jsb_Config(SSE_)%l_freeze            , &
        & dsl4jsb_Config(SSE_)%l_supercool         , &
        & dsl4jsb_Config(HYDRO_)%enforce_water_budget, & ! Water balance check setting
        & t_soil_top                          (:  ), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        ! inout
        & t_soil_sl                           (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd                           (:  ), &
        & grnd_hflx                           (:  ), &
        & thaw_depth                          (:  ), &
        & is_snow_covered                     (:  )  &
        & )
    ELSE
      ALLOCATE(ws_max_sl(nc,nsoil))
      !$ACC ENTER DATA CREATE(ws_max_sl) ASYNC(acc_stream)
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
#endif
      DO is = 1,nsoil
        DO ic = 1, nc
          ws_max_sl(ic,is) = soil_depth_sl(ic,is) * vol_porosity_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      CALL calc_soil_temperature( &
        ! in
        & nc                                       , &
        & nsoil                                    , &
        & soil_e%dz(:)                             , &
        & soil_e%levels(:)                         , &
        & soil_e%ubounds(:)                        , &
        & dtime                                    , &
        & lstart                                   , &
        & dsl4jsb_Config(SSE_)%l_freeze            , &
        & dsl4jsb_Config(SSE_)%l_supercool         , &
        & dsl4jsb_Config(HYDRO_)%enforce_water_budget, & ! Water balance check setting
        & t_soil_top                          (:  ), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        ! inout
        & t_soil_sl                           (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  ), &
        & is_snow_covered                     (:  ), &
        ! optional, in
        & soil_depth_sl                       (:,:), &
        & ws_max_sl                           (:,:), &
        & matric_pot_sl                       (:,:), &
        & bclapp_sl                           (:,:), &
        ! optional, inout
        & wtr_soil_sl                         (:,:), &
        & ice_soil_sl                         (:,:), &
        & le_pc_remain                        (:  ), &
        & wtr_soil_pot_scool_sl               (:,:), &
        ! optional, out
        & wtr_freeze_sl                       (:,:), &
        & ice_melt_sl                         (:,:)  &
        & )
      !$ACC EXIT DATA DELETE(ws_max_sl) ASYNC(acc_stream)
      DEALLOCATE(ws_max_sl)

      CALL calc_snow_abcoeff( &
        ! in
        & nc,                     &
        & nsnow,                  &
        & soil_e%dz(1:2),         &
        & dtime,                  &
        & snow_depth_sl    (:,:), &
        & itop               (:), &
        & vol_heat_cap_snow  (:), &
        & heat_cond_snow     (:), &
        & vol_heat_cap_sl  (:,1), &
        & heat_cond_sl     (:,1), &
        & t_soil_sl        (:,:), &
        & t_soil_acoef     (:,:), &
        & t_soil_bcoef     (:,:), &
        ! inout
        & t_snow           (:,:), &
        & t_snow_acoef     (:,:), &
        & t_snow_bcoef     (:,:), &
        ! out
        & hcap_grnd_snow     (:), &
        & grnd_hflx_snow     (:)  &
        & )

    END IF

    IF ((.NOT. lstart) .AND. is_newyear(options%current_datetime, options%dtime)) THEN
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        thaw_depth_max_ym1(ic) = thaw_depth_max(ic)
        thaw_depth_max(ic)     = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        thaw_depth_max(ic) = MAX(thaw_depth_max(ic),thaw_depth(ic))
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (.NOT. tile%is_glacier) THEN
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(acc_stream)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
#endif
      DO ic = 1,nc
        IF (.NOT. l_uniform_snow .AND. snow_depth_sl(ic,nsnow) > 0._wp) THEN
          grnd_hflx(ic) = snow_fract(ic) * grnd_hflx_snow(ic) + (1._wp - snow_fract(ic)) * grnd_hflx_soil(ic)
          hcap_grnd(ic) = snow_fract(ic) * hcap_grnd_snow(ic) + (1._wp - snow_fract(ic)) * hcap_grnd_soil(ic)
        ELSE IF (l_uniform_snow .AND. snow_depth(ic) > snow_depth_min) THEN
          grnd_hflx(ic) = grnd_hflx_snow(ic)
          hcap_grnd(ic) = hcap_grnd_snow(ic)
        ELSE
          grnd_hflx(ic) = grnd_hflx_soil(ic)
          hcap_grnd(ic) = hcap_grnd_soil(ic)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC EXIT DATA ASYNC(acc_stream) &
    !$ACC   DELETE(snow_depth, snow_fract, t_soil_top, hcap_grnd_snow, is_snow_covered) &
    !$ACC   DELETE(hcap_grnd_soil, grnd_hflx_snow, grnd_hflx_soil, ldglac, itop, itop_old)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_and_snow_temperature

  SUBROUTINE aggregate_soil_and_snow_temperature(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    TYPE(t_jsb_model), POINTER :: model                    !< Current instance of the model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract !< Aggregator weighted by tile fraction

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_temperature'

    INTEGER :: iblk, ics, ice                             !< Block bounds

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_acoef,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_bcoef,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_sl,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   grnd_hflx,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   grnd_hflx_old,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   hcap_grnd,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   hcap_grnd_old,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   snow_depth_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_snow_acoef,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_snow_bcoef,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_snow,                weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_freeze_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_melt_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pot_scool_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, le_pc_remain,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_ym1,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   vol_heat_cap,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   heat_cond,             weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_temperature

  ! ===================================================================================================== !
  !>
  !>#### Update soil and snow properties for one time step
  !>
  SUBROUTINE update_soil_and_snow_properties(tile, options)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min ! , cs

    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, vol_hcap_org_top, vol_hcap_org_below, &
                                       & hcond_ice, hcond_snow, hcond_snow_min, &
                                       & hcond_org_top, hcond_dry_org_top, hcond_org_below, &
                                       & hcond_dry_org_below !, vol_hcap_snow_min
    USE mo_hydro_constants,        ONLY: vol_porosity_org_top, vol_porosity_org_below
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity, add_heat_capacity_pond

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    ! Local variables
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & snow_soil_dens,       & !< Snow density [kg m-3]
      & vol_heat_cap,         & !< Effective volumetric heat capacity [J m-3 K-1]
      & heat_cond,            & !< Effective thermal conductivity [J m-1 s-1 K-1]
      & vol_heat_cap_snow,    & !< Volumetric heat capacity of snow [J m-3 K-1]
      & heat_cond_snow,       & !< Thermal conductivity of snow [J m-1 s-1 K-1]
      & wtr_pond,             & !< Liquid water stored in ponds [m]
      & ice_pond                !< Ice stored in ponds [m]
    dsl4jsb_Real3D_onChunk :: &
      & soil_depth_sl,        & !< Depth of soil layers [m]
      & fract_org_sl,         & !< Organic matter fraction of soil layers []
      & vol_heat_cap_sl,      & !< Volumetric heat capacity of soil layers [J m-3 K-1]
      & heat_cond_sl,         & !< Thermal conductivity of soil layers [J m-1 s-1 K-1]
      & wtr_soil_sl,          & !< Liquid water content of soil layers [m]
      & ice_soil_sl,          & !< Ice content of soil layers [m]
      & vol_porosity_sl         !< Volumetric porosity of soil layers []

    ! Locally allocated vectors
    !
    dsl4jsb_Real3D_onChunk :: &
      & fract_water,          & !< Fraction of soil water in each layer []
      & fract_ice               !< Fraction of soil ice in each layer []

    TYPE(t_jsb_model), POINTER :: model    !< Current instance of the model
    TYPE(t_jsb_vgrid), POINTER :: soil_e   !< Soil energy vertical grid

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & vol_heat_cap_tmp,                    & !< Temporary volumetric heat capacity by soil layer [J m-3 K-1]
      & vol_heat_cap_org_tmp,                & !< Temporary organic volumetric heat capacity by soil layer [J m-3 K-1]
      & vol_porosity_org_tmp,                & !< Temporary organic soil porosity by soil layer []
      & heat_cond_tmp,                       & !< Temporary thermal conductivity by soil layer [J m-1 s-1 K-1]
      & hcond_org_tmp,                       & !< Temporary organic thermal conductivity [J m-1 s-1 K-1]
      & hcond_dry_org_tmp                      !< Temporary dry organic thermal conductivity [J m-1 s-1 K-1]

    INTEGER  :: iblk, ics, ice, nc, nsoil, ic, is !< Block bounds and horizontal and vertical loop indices
    REAL(wp), POINTER :: soil_e_dz(:)             !< Thicknesses of soil energy layers [m]

    REAL(wp) :: hlp1                              !< Temporary scalar helper variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    soil_e_dz => soil_e%dz

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap_snow) ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond_snow)    ! out


    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          vol_heat_cap_sl(ic,is) = vol_hcap_ice
          heat_cond_sl   (ic,is) = hcond_ice
        END DO
      END DO
      !$ACC END LOOP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        vol_heat_cap_snow(ic) = vol_hcap_snow
        heat_cond_snow   (ic) = hcond_snow
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL
      IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
      RETURN
    END IF

    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap)
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)

    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nc, nsoil))
      !$ACC ENTER DATA CREATE(fract_org_sl) ASYNC(acc_stream)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          fract_org_sl(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nc, nsoil), fract_ice(nc, nsoil))
      !$ACC ENTER DATA CREATE(fract_water, fract_ice) ASYNC(acc_stream)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (soil_depth_sl(ic,is) > 0._wp) THEN
            fract_water(ic,is) = wtr_soil_sl(ic,is) / soil_depth_sl(ic,is)
            fract_ice  (ic,is) = ice_soil_sl(ic,is) / soil_depth_sl(ic,is)
          ELSE
            fract_water(ic,is) = 0._wp
            fract_ice  (ic,is) = 0._wp
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ALLOCATE(vol_heat_cap_tmp(nc,nsoil), vol_heat_cap_org_tmp(nc,nsoil))
      !$ACC ENTER DATA CREATE(vol_heat_cap_tmp, vol_heat_cap_org_tmp) ASYNC(acc_stream)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        vol_heat_cap_tmp    (ic,1) = vol_heat_cap(ic)
        vol_heat_cap_org_tmp(ic,1) = vol_hcap_org_top
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 2, nsoil
        DO ic = 1, nc
          vol_heat_cap_tmp    (ic,is) = vol_heat_cap(ic)
          vol_heat_cap_org_tmp(ic,is) = vol_hcap_org_below
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          vol_heat_cap_sl(ic,is) = calc_vol_heat_capacity(          &
            &                         vol_heat_cap_tmp(ic,is)     , &
            &                         vol_heat_cap_org_tmp(ic,is) , &
            &                         vol_porosity_sl(ic,is)      , &
            &                         fract_org_sl(ic,is)         , &
            &                         fract_water(ic,is)          , &
            &                         fract_ice(ic,is)            , &
            &                         soil_depth_sl(ic,is)        , &
            &                         soil_e_dz(is)                 &
            & )
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC EXIT DATA DELETE(fract_water, fract_ice, vol_heat_cap_tmp, vol_heat_cap_org_tmp) ASYNC(acc_stream)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_tmp, vol_heat_cap_org_tmp)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          vol_heat_cap_sl(ic,is) = vol_heat_cap(ic)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ENDIF

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_pond) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
      !$ACC   PRIVATE(hlp1)
      DO ic = 1, nc
        hlp1 = vol_heat_cap_sl(ic,1)
        vol_heat_cap_sl(ic,1) = add_heat_capacity_pond( &
          & soil_e_dz(1),                               &
          & wtr_pond(ic),                               &
          & ice_pond(ic),                               &
          & hlp1                                        &
          & )
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN

      ALLOCATE(heat_cond_tmp(nc,nsoil), vol_porosity_org_tmp(nc,nsoil), &
             & hcond_org_tmp(nc,nsoil), hcond_dry_org_tmp(nc,nsoil))
      !$ACC ENTER DATA CREATE(heat_cond_tmp, vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp) ASYNC(acc_stream)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        heat_cond_tmp(ic,1)         = heat_cond(ic)
        vol_porosity_org_tmp(ic,1)  = vol_porosity_org_top
        hcond_org_tmp(ic,1)         = hcond_org_top
        hcond_dry_org_tmp(ic,1)     = hcond_dry_org_top
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 2, nsoil
        DO ic = 1, nc
          heat_cond_tmp(ic,is)        = heat_cond(ic)
          vol_porosity_org_tmp(ic,is) = vol_porosity_org_below
          hcond_org_tmp(ic,is)        = hcond_org_below
          hcond_dry_org_tmp(ic,is)    = hcond_dry_org_below
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          heat_cond_sl (ic,is) = calc_thermal_conductivity( &
            &                                heat_cond_tmp(ic,is)          , &
            &                                wtr_soil_sl(ic,is)            , &
            &                                ice_soil_sl(ic,is)            , &
            &                                vol_porosity_sl(ic,is)        , &
            &                                fract_org_sl(ic,is)           , &
            &                                vol_porosity_org_tmp(ic,is)   , &
            &                                hcond_org_tmp(ic,is)          , &
            &                                hcond_dry_org_tmp(ic,is)      , &
            &                                soil_depth_sl(ic,is)          , &
            &                                soil_e_dz(is)                   &
            & )
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC EXIT DATA DELETE(heat_cond_tmp, vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp) ASYNC(acc_stream)
      DEALLOCATE(heat_cond_tmp, vol_porosity_org_tmp,hcond_org_tmp, hcond_dry_org_tmp)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          heat_cond_sl(ic,is) = heat_cond(ic)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$ACC EXIT DATA DELETE(fract_org_sl) ASYNC(acc_stream)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF (snow_soil_dens(ic) > dens_snow_min) THEN         ! "old" or mixture of old and fresh snow
          vol_heat_cap_snow(ic) = ci * snow_soil_dens(ic)    ! Use specific heat of ice (Verseghy 1991)
          ! heat_cond_snow   (:) = 2.576e-6_wp * snow_soil_dens(:) ** 2._wp + 0.074_wp                      ! Verseghy (1991)
          ! heat_cond_snow   (:) = 2.9e-6_wp*snow_soil_dens(:)**2._wp                                       ! Goodrich (1982)
          ! Calonne et al. (2011):
          heat_cond_snow   (ic) = 2.5e-6_wp*snow_soil_dens(ic)**2._wp - 0.123e-3_wp*snow_soil_dens(ic) + 0.024_wp
        ELSE                                         ! fresh snow
          ! @todo: these minimum values from mo_sse_constants.f90 don't match with the above equations and minimum snow density
          ! vol_heat_cap_snow(:) = vol_hcap_snow_min
          heat_cond_snow   (ic) = hcond_snow_min
          vol_heat_cap_snow(ic) = ci * dens_snow_min
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        vol_heat_cap_snow(ic) = vol_hcap_snow  ! Use specific heat of snow (cs * dens_snow)
        heat_cond_snow   (ic) = hcond_snow
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_and_snow_properties
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_and_snow_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregator weighted by tile fraction

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_properties'

    INTEGER :: iblk, ics, ice                              !< Block bounds

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(SSE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SSE_, vol_heat_cap_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, heat_cond_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, vol_heat_cap_snow, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, heat_cond_snow,    weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_properties

#endif
END MODULE mo_sse_interface
