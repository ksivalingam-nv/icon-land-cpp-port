!> Interface to run ICON-Land for one time step
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
!>#### Interface and routines to run a ICON-Land time step
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish !, message_text
  ! USE mo_util,         ONLY: int2string

  USE mo_jsb_control,         ONLY: debug_on, timer_on, l_timer_host, timer_jsbach, jsbach_runs_standalone
  USE mo_timer,               ONLY: timer_start, timer_stop
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
!  USE mo_jsb_controller,      ONLY: Run_process_task
  USE mo_atmland_interface,   ONLY: update_atm2land, update_land2atm
  USE mo_jsb_parallel,        ONLY: Get_omp_thread, Is_omp_inside_serial
  USE mo_jsb_time,            ONLY: t_datetime, is_newyear, is_newday, &
      &                             is_time_experiment_start, is_time_restart
  USE mo_jsb_physical_constants, ONLY: cvdifts
  USE mo_jsb_time_iface,      ONLY: cache_time_interpolation_weights

  dsl4jsb_Use_processes SEB_, HYDRO_, PHENO_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_interface, jsbach_start_timestep, jsbach_finish_timestep
  PUBLIC :: jsbach_get_var
  PUBLIC :: invalidate_cuda_graphs

  INTERFACE jsbach_interface
    MODULE PROCEDURE interface_full
    ! MODULE PROCEDURE interface_inquire
  END INTERFACE jsbach_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_interface'
  LOGICAL :: invalidate_cuda_graphs ! If true, all the existing graphs will be deleted
                                    ! For example, after reading new data or at the start of a new day

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Called by the host model (i.e. ICON) in simulations with atmosphere
  !>
  !> @NOTE In standalone simulations 'mo_jsbach_model:jsbach_model' is called
  !>       and variables are provided via 'mo_jsbach_model:jsbach_model:get_standalone_driver'
  !>
  SUBROUTINE interface_full( &
    & model_id,         &
    & iblk, ics, ice,   &
    & current_datetime, & !< in -- object containing date and time of the timestep we are working on
    & dtime,            &
    & t_air,            &
    & q_air,            &
    & rain,             &
    & snow,             &
    & wind_air,         &
    & wind_10m,         &
    & lw_srf_down,      &
    & swvis_srf_down,   &
    & swnir_srf_down,   &
    & swpar_srf_down,   &
    & fract_par_diffuse, &
    & press_srf,        &
    & drag_srf,         &
    & t_acoef,          &
    & t_bcoef,          &
    & q_acoef,          &
    & q_bcoef,          &
    & pch,              &
    & cos_zenith_angle, &
    & CO2_air,          &
    & t_srf,            &
    & t_eff_srf,        &
    & qsat_srf,         &
    & s_srf,            &
    & fact_q_air,       &
    & fact_qsat_srf,    &
    & evapotrans,       &
    & latent_hflx,      &
    & sensible_hflx,    &
    & grnd_hflx,        &
    & grnd_hcap,        &
    & rough_h_srf,      &
    & rough_m_srf,      &
    & q_snocpymlt,      &
    & alb_vis_dir,      &
    & alb_nir_dir,      &
    & alb_vis_dif,      &
    & alb_nir_dif,      &
    & CO2_flux,         &
    ! OPTIONAL arguments
    & DEBUG_VAR,        &
    ! For lakes
    & drag_wtr, drag_ice,                                               & !< in
    & t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr,               & !< in
    & t_acoef_ice, t_bcoef_ice, q_acoef_ice, q_bcoef_ice,               & !< in
    & t_lwtr, qsat_lwtr, s_lwtr,                                        & !< out
    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                    & !< out
    & albedo_lwtr,                                                      & !< out
    & t_lice, qsat_lice, s_lice,                                        &
    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                    & !< out
    & albedo_lice,                                                      & !< out
    & ice_fract_lake,                                                   & !< out
    ! For ICON-Land standalone
    & evapopot                                                          & !< out
    & )

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: &
      & model_id,          &
      & iblk,              &
      & ics,               &
      & ice
    TYPE(t_datetime), POINTER, INTENT(in) :: &
      & current_datetime
    REAL(wp), INTENT(in) :: &
      & dtime
    REAL(wp), INTENT(in) ::    &
      & t_air             (:), &
      & q_air             (:), &
      & rain              (:), &
      & snow              (:), &
      & wind_air          (:), &
      & wind_10m          (:), &
      & lw_srf_down       (:), &
      & swvis_srf_down    (:), &
      & swnir_srf_down    (:), &
      & swpar_srf_down    (:), &
      & fract_par_diffuse  (:), &
      & press_srf         (:), &
      & drag_srf          (:), &
      & t_acoef           (:), &
      & t_bcoef           (:), &
      & q_acoef           (:), &
      & q_bcoef           (:), &
      & pch               (:), &
      & cos_zenith_angle  (:), &
      & CO2_air           (:)
    REAL(wp), INTENT(in), OPTIONAL :: &
      & DEBUG_VAR         (:), &
      ! For lakes:
      & drag_wtr          (:), &
      & drag_ice          (:), &
      & t_acoef_wtr       (:), &
      & t_bcoef_wtr       (:), &
      & q_acoef_wtr       (:), &
      & q_bcoef_wtr       (:), &
      & t_acoef_ice       (:), &
      & t_bcoef_ice       (:), &
      & q_acoef_ice       (:), &
      & q_bcoef_ice       (:)
    REAL(wp), INTENT(out) ::   &
      & t_srf             (:), &
      & t_eff_srf         (:), &
      & qsat_srf          (:), &
      & s_srf             (:), &
      & fact_q_air        (:), &
      & fact_qsat_srf     (:), &
      & evapotrans        (:), &
      & latent_hflx       (:), &
      & sensible_hflx     (:), &
      & grnd_hflx         (:), &
      & grnd_hcap         (:), &
      & rough_h_srf       (:), &
      & rough_m_srf       (:), &
      & q_snocpymlt       (:), &
      & alb_vis_dir       (:), &
      & alb_nir_dir       (:), &
      & alb_vis_dif       (:), &
      & alb_nir_dif       (:), &
      & CO2_flux          (:)
    ! For lakes:
    REAL(wp), INTENT(out), OPTIONAL ::   &
      & t_lwtr            (:), &
      & qsat_lwtr         (:), &
      & s_lwtr            (:), &
      & evapo_wtr         (:), &
      & latent_hflx_wtr   (:), &
      & sensible_hflx_wtr (:), &
      & albedo_lwtr       (:), &
      & t_lice            (:), &
      & qsat_lice         (:), &
      & s_lice            (:), &
      & evapo_ice         (:), &
      & latent_hflx_ice   (:), &
      & sensible_hflx_ice (:), &
      & albedo_lice       (:), &
      & ice_fract_lake    (:)
    ! For standalone jsbach:
    REAL(wp), INTENT(out), OPTIONAL :: &
      & evapopot(:)

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_model),         POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                             :: nc
    INTEGER                             :: no_omp_thread
    ! INTEGER                           :: p_pe
    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_full'
    ! ----------------------------------------------------------------------------------------------------- !

    IF (timer_on() .OR. l_timer_host) CALL timer_start(timer_jsbach(model_id))

    no_omp_thread = Get_omp_thread()

    nc = ice - ics + 1
    model => Get_model(model_id)

    IF (debug_on('basic') .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.0 Update options
    !>
    CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, current_datetime=current_datetime, &
      & dtime=dtime, alpha=MERGE(1._wp, cvdifts, jsbach_runs_standalone()) )

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.0 Put atmospheric forcing into the a2l memory
    !>

    CALL model%Get_top_tile(tile)
      CALL update_atm2land(                           &
        & tile,                                       &
        & model%options(no_omp_thread),               &
        & t_air              = t_air(:),              &
        & q_air              = q_air(:),              &
        & rain               = rain(:),               &
        & snow               = snow(:),               &
        & wind_air           = wind_air(:),           &
        & wind_10m           = wind_10m(:),           &
        & lw_srf_down        = lw_srf_down(:),        &
        & swvis_srf_down     = swvis_srf_down(:),     &
        & swnir_srf_down     = swnir_srf_down(:),     &
        & swpar_srf_down     = swpar_srf_down(:),     &
        & fract_par_diffuse   = fract_par_diffuse(:), &
        & press_srf          = press_srf(:),          &
        & drag_srf           = drag_srf(:),           &
        & t_acoef            = t_acoef(:),            &
        & t_bcoef            = t_bcoef(:),            &
        & q_acoef            = q_acoef(:),            &
        & q_bcoef            = q_bcoef(:),            &
        & pch                = pch(:),                &
        & cos_zenith_angle   = cos_zenith_angle(:),   &
        & CO2_air            = CO2_air(:),            &
        & DEBUG_VAR          = DEBUG_VAR,             &
        ! For lakes
        & drag_wtr           = drag_wtr,              &
        & drag_ice           = drag_ice,              &
        & t_acoef_wtr        = t_acoef_wtr,           &
        & t_bcoef_wtr        = t_bcoef_wtr,           &
        & q_acoef_wtr        = q_acoef_wtr,           &
        & q_bcoef_wtr        = q_bcoef_wtr,           &
        & t_acoef_ice        = t_acoef_ice,           &
        & t_bcoef_ice        = t_bcoef_ice,           &
        & q_acoef_ice        = q_acoef_ice,           &
        & q_bcoef_ice        = q_bcoef_ice            &
        & )

    NULLIFY(tile)


    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.0 Run tasks in queue
    !>
    !>  Run tasks of processes as defined in 'mo_jsb_model_init:jsbach_init'
    !>
    CALL model%Run_tasks(debug_on('tasks') .AND. iblk == 1)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.0 Give the required variables back to the atmosphere
    !>
    CALL model%Get_top_tile(tile)
      CALL update_land2atm(                        &
        & tile,                                    &
        & model%options(no_omp_thread),            &
        & t_srf              = t_srf,              &
        & t_eff_srf          = t_eff_srf,          &
        & qsat_srf           = qsat_srf,           &
        & s_srf              = s_srf,              &
        & fact_q_air         = fact_q_air,         &
        & fact_qsat_srf      = fact_qsat_srf,      &
        & evapopot           = evapopot,           &
        & evapotrans         = evapotrans,         &
        & latent_hflx        = latent_hflx,        &
        & sensible_hflx      = sensible_hflx,      &
        & grnd_hflx          = grnd_hflx,          &
        & grnd_hcap          = grnd_hcap,          &
        & rough_h_srf        = rough_h_srf,        &
        & rough_m_srf        = rough_m_srf,        &
        & q_snocpymlt        = q_snocpymlt,        &
        & alb_vis_dir        = alb_vis_dir,        &
        & alb_nir_dir        = alb_nir_dir,        &
        & alb_vis_dif        = alb_vis_dif,        &
        & alb_nir_dif        = alb_nir_dif,        &
        & CO2_flux           = CO2_flux,           &
        & t_lwtr             = t_lwtr,             &
        & qsat_lwtr          = qsat_lwtr,          &
        & s_lwtr             = s_lwtr,             &
        & evapo_wtr          = evapo_wtr,          &
        & latent_hflx_wtr    = latent_hflx_wtr,    &
        & sensible_hflx_wtr  = sensible_hflx_wtr,  &
        & t_lice             = t_lice,             &
        & qsat_lice          = qsat_lice,          &
        & s_lice             = s_lice,             &
        & evapo_ice          = evapo_ice,          &
        & latent_hflx_ice    = latent_hflx_ice,    &
        & sensible_hflx_ice  = sensible_hflx_ice,  &
        & ice_fract_lake     = ice_fract_lake,     &
        & albedo_lwtr        = albedo_lwtr,        &
        & albedo_lice        = albedo_lice         &
        & )
    NULLIFY(tile)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.0 messages
    !>
    IF (debug_on('basic') .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')

    ! ----------------------------------------------------------------------------------------------------- !
    !> 6.0 update timer
    !>
    IF (timer_on() .OR. l_timer_host) CALL timer_stop(timer_jsbach(model_id))

  END SUBROUTINE interface_full

  ! ======================================================================================================= !
  !> @TODO not used at the moment
  !>
  ! SUBROUTINE interface_inquire( model_id, iblk, ics, ice,         &
  !   & t_srf, t_eff_srf, qsat_srf, zh_srf, zm_srf, alb_vis, alb_nir)
  !
  !   CHARACTER(len=*), PARAMETER :: routine = modname//':interface_inquire'
  !
  !   INTEGER, INTENT(in) :: &
  !     & model_id,          &
  !     & iblk,              &
  !     & ics,               &
  !     & ice
  !
  !   REAL(wp), INTENT(inout) :: &
  !     & t_srf    (:), &
  !     & t_eff_srf(:), &
  !     & qsat_srf (:), &
  !     & zh_srf   (:), &
  !     & zm_srf   (:), &
  !     & alb_vis  (:), &
  !     & alb_nir  (:)

  ! END SUBROUTINE interface_inquire

  ! ======================================================================================================= !
  !>
  !> Start ICON-Land timestep
  !>
  !>   is called by both:
  !>      ICON-Land standalone: 'mo_jsbach_model:run_one_timestep'
  !>      host model (ICON):    './atm_phy_echam/mo_interface_iconam_echam.f90'
  !>
  SUBROUTINE jsbach_start_timestep(model_id, current_datetime, dtime)

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    REAL(wp),                  INTENT(in) :: dtime
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),  POINTER :: model
    INTEGER                     :: no_omp_thread
    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_start_timestep'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    no_omp_thread = Get_omp_thread()   ! Should always be 1
    IF (no_omp_thread /= 1) THEN
      CALL finish(routine, 'should not be called inside OMP parallel region!')
    END IF
    model => Get_model(model_id)

    CALL model%Set_options(current_datetime=current_datetime, dtime=dtime)


    ! If this is .true., old graphs will be invalidated and the new ones will be captured
    invalidate_cuda_graphs = .FALSE.
    IF (     is_time_experiment_start(current_datetime) .OR. is_time_restart(current_datetime)  &
        .OR. is_newday(current_datetime, dtime) .OR. is_newyear(current_datetime, dtime) ) THEN
      ! Invalidate graphs when the workflow is changing
      invalidate_cuda_graphs = .TRUE.
    END IF

    ! Save time interpolation weights for further usage
    ! This simplifies access of the weights in case of CUDA graphs
    CALL cache_time_interpolation_weights()

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE jsbach_start_timestep

  ! ======================================================================================================= !
  !>
  !> Finish ICON-land timestep
  !>
  !>   is called by both:
  !>      ICON-Land standalone: 'mo_jsbach_model:run_one_timestep'
  !>      host model (ICON):    './atm_phy_echam/mo_interface_iconam_echam.f90'
  !>
  SUBROUTINE jsbach_finish_timestep(model_id, current_datetime, dtime, init)
    USE mo_util,        ONLY: report_memory_usage
    USE mo_jsb_control, ONLY: get_debug_memory_level
    USE mo_jsb_time,    ONLY: is_newmonth

    USE mo_seb_interface,     ONLY: seb_check_temperature_range
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    REAL(wp),                  INTENT(in) :: dtime
    LOGICAL, OPTIONAL,         INTENT(in) :: init
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    LOGICAL                             :: l_trace_memory
    INTEGER                             :: no_omp_thread
    LOGICAL                             :: linit
    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_finish_timestep'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')
    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF
    no_omp_thread = Get_omp_thread()

    model => Get_model(model_id)

    linit = .FALSE.
    IF (PRESENT(init)) linit = init
    ! ----------------------------------------------------------------------------------------------------- !


    ! Update options
    CALL model%Set_options(current_datetime=current_datetime, dtime=dtime)

    CALL model%Get_top_tile(tile)

    ! Make sure temperatures are within an exceptable range
    IF (model%processes(SEB_)%p%config%active) THEN
      CALL seb_check_temperature_range(model_id)
    END IF


! #ifdef _OPENACC
!    CALL model%gpu_to_cpu()
! #endif

    ! Global diagnostics
    CALL global_diagnostics(tile)

    l_trace_memory = .FALSE.
    SELECT CASE (get_debug_memory_level())
    CASE (1)
      l_trace_memory = is_newmonth(current_datetime, dtime)
    CASE (2)
      l_trace_memory = is_newday(current_datetime, dtime)
    CASE (3)
      l_trace_memory = .TRUE.
    END SELECT

    IF (l_trace_memory) CALL report_memory_usage()

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')
  END SUBROUTINE jsbach_finish_timestep


  ! ======================================================================================================== !
  !>
  !>#### Routine encapsulating the call of global land diagnostic routines on the different processes
  !>
  !> Calculates area scaling factors and calls global diagnostic routines for active processes which
  !> provide such a routine
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_diagnostics(tile)

    USE mo_jsb_control,         ONLY: acc_stream
    USE mo_sync,                ONLY: global_sum_array
    USE mo_jsb_grid,            ONLY: Get_grid
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid
    USE mo_jsb_process_factory, ONLY: max_no_of_processes

    USE mo_seb_interface,      ONLY: global_seb_diagnostics
    USE mo_hydro_interface,    ONLY: global_hydrology_diagnostics
    USE mo_pheno_interface,    ONLY: global_phenology_diagnostics
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile  !< Current tile
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_diagnostics'

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    INTEGER                :: iproc, jc, blk
    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain(:,:)    ! 1._wp: cell in domain, 0._wp: halo cell
                                                ! - to mask all halo cells for global sums (otherwise these cells are counted twice)
    REAL(wp), ALLOCATABLE  :: area_scaling(:,:)
    REAL(wp), ALLOCATABLE  :: land_area_scaling(:,:)
    REAL(wp)               :: global_land_area
    ! ----------------------------------------------------------------------------------------------------- !
    notsea => tile%fract(:,:)   ! fraction of the box tile: notsea
    model  => Get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ALLOCATE (in_domain(grid%nproma,grid%nblks), area_scaling(grid%nproma,grid%nblks), land_area_scaling(grid%nproma,grid%nblks))
    !$ACC ENTER DATA CREATE(in_domain, area_scaling, land_area_scaling) ASYNC(acc_stream)

    ! Calculate global land area
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
    DO blk = 1, grid%nblks
      DO jc = 1, grid%nproma
        IF (is_in_domain(jc,blk)) THEN
          in_domain(jc,blk) = 1._wp
        ELSE
          in_domain(jc,blk) = 0._wp
        END IF

        area_scaling(jc,blk) = area(jc,blk) * in_domain(jc,blk)
        land_area_scaling(jc,blk) = area_scaling(jc,blk) * notsea(jc,blk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    global_land_area = global_sum_array(land_area_scaling, lacc=.TRUE.)

    DO iproc = 1, max_no_of_processes
      IF (ASSOCIATED(model%processes(iproc)%p)) THEN
        IF (ASSOCIATED(model%processes(iproc)%p%config)) THEN
          IF (model%processes(iproc)%p%config%active) THEN
            SELECT CASE (iproc)
            CASE (SEB_)
              CALL global_seb_diagnostics(tile, global_land_area, land_area_scaling)
            CASE (HYDRO_)
              CALL global_hydrology_diagnostics(tile, global_land_area, in_domain, land_area_scaling)
            CASE (PHENO_)
              CALL global_phenology_diagnostics(tile, global_land_area, land_area_scaling)
            END SELECT
          END IF
        END IF
      END IF
    END DO

    !$ACC WAIT(acc_stream)
    !$ACC EXIT DATA DELETE(in_domain, area_scaling, land_area_scaling) ASYNC(acc_stream)
    DEALLOCATE (area_scaling, in_domain, land_area_scaling)

  END SUBROUTINE global_diagnostics


  ! ====================================================================================================== !
  !
  !> Routine to be called from host model in order to inquire a variable (section) on a specified tile
  !
  SUBROUTINE jsbach_get_var( &
    & var_name,              &
    & model_id,              &
    & tile,                  &
    & ics, ice, iblk,        &
    & arr1d, arr2d, arr3d,   &
    & ptr1d, ptr2d, ptr3d,   &
    & lacc, opt_acc_async_queue)

    USE mo_jsb_var_class,     ONLY: t_jsb_var, REAL2D, REAL3D
    USE mo_jsb_process_class, ONLY: Get_process_id
    USE mo_fortran_tools,     ONLY: copy, set_acc_host_or_device, set_acc_async_queue

    CHARACTER(len=*),  INTENT(in)          :: var_name
                                              !< Name of requested variable including the leading process name
    INTEGER,           INTENT(in)          :: model_id
                                              !< Model id
    ! Optional arguments, exactly one of the ptr?d or arr?d arguments is required
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: tile
                                              !< Name of tile from which variable is taken; box tile if not given
    INTEGER, OPTIONAL, INTENT(in)          :: ics, ice, iblk
                                              !< If given, chunk (ics:ice,iblk) is returned, otherwise whole domain
    REAL(wp), OPTIONAL, INTENT(inout)      :: arr1d(:), arr2d(:,:), arr3d(:,:,:)
                                              !< Requested variable (section) is copied into these arrays
    REAL(wp), OPTIONAL, POINTER            :: ptr1d(:), ptr2d(:,:), ptr3d(:,:,:)
                                              !< Pointers to the requested variable (section); must not be associated
    LOGICAL, OPTIONAL, INTENT(in)          :: lacc
                                              !< Switch whether the variables should be copied to GPU
    INTEGER, OPTIONAL, INTENT(in)          :: opt_acc_async_queue
                                              !< OpenACC async queue to use, default is 1
    INTEGER :: acc_async_queue

    TYPE(t_jsb_model),           POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile_ptr
    CLASS(t_jsb_var),            POINTER :: var_ptr
    CHARACTER(len=:), ALLOCATABLE :: process_name, variable_name
    INTEGER :: process_id, pos, ic, nc, icount, itile
    LOGICAL :: l_on_chunk, l_found
    INTEGER :: no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_get_var'

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)
    CALL set_acc_async_queue(acc_async_queue, opt_acc_async_queue)

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    icount = 0
    IF (PRESENT(ptr1d)) THEN
      IF (ASSOCIATED(ptr1d)) CALL finish(routine, 'ptr1d must not be ASSOCIATED for var '//var_name)
      icount = icount + 1
    END IF
    IF (PRESENT(ptr2d)) THEN
      IF (ASSOCIATED(ptr2d)) CALL finish(routine, 'ptr2d must not be ASSOCIATED for var '//var_name)
      icount = icount + 1
    END IF
    IF (PRESENT(ptr3d)) THEN
      IF (ASSOCIATED(ptr3d)) CALL finish(routine, 'ptr3d must not be ASSOCIATED for var '//var_name)
      icount = icount + 1
    END IF
    IF (PRESENT(arr1d)) icount = icount + 1
    IF (PRESENT(arr2d)) icount = icount + 1
    IF (PRESENT(arr3d)) icount = icount + 1
    IF (icount == 0) CALL finish(routine, 'no target specified')
    IF (icount > 1) CALL finish(routine, 'more than one target specified')

    model => Get_model(model_id)

    CALL model%Get_top_tile(tile_ptr)
    ! If `tile` is not present the variable on the box tile is used, otherwise search for tile with given name
    l_found = .TRUE.
    IF (PRESENT(tile)) THEN
      l_found = .FALSE.
      DO itile=1,SIZE(model%tiles)
        IF (TRIM(model%tiles(itile)%p%name) == tile) THEN
          tile_ptr => model%tiles(itile)%p
          l_found = .TRUE.
          EXIT
        END IF
      END DO
    END IF
    IF (.NOT. l_found) CALL finish(routine, 'tile '//tile//' not found!')

    pos = INDEX(var_name, '_')
    IF (pos < 1) CALL finish(routine, 'var_name must have format <process name>_<variable name>')
    process_name = var_name(1:pos-1)
    variable_name = var_name(pos+1:LEN(TRIM(var_name)))
    process_id = Get_process_id(process_name)
    var_ptr => tile_ptr%mem(process_id)%p%Get_var(variable_name)
    IF (.NOT. ASSOCIATED(var_ptr)) CALL finish(routine, 'variable '//var_name//'_'//TRIM(tile_ptr%name)//' not found')

    IF (PRESENT(ics) .OR. PRESENT(ice) .OR. PRESENT(iblk)) THEN
      IF (.NOT. (PRESENT(ics) .AND. PRESENT(ice) .AND. PRESENT(iblk))) THEN
        CALL finish(routine, 'Must provide either none or all of ics,ice,iblk arguments!')
      ELSE
        l_on_chunk = .TRUE.
        nc = ice - ics + 1
        IF (var_ptr%type == REAL2d) THEN
          IF (.NOT. PRESENT(ptr1d) .AND. .NOT. PRESENT(arr1d)) CALL finish(routine, 'Must provide argument ptr1d or arr1d')
        ELSE
          IF (.NOT. PRESENT(ptr2d) .AND. .NOT. PRESENT(arr2d)) CALL finish(routine, 'Must provide argument ptr2d or arr2d')
        END IF
      END IF
    ELSE
      IF (.NOT. Is_omp_inside_serial()) THEN
        CALL finish(TRIM(routine), 'Should not be called inside parallel OMP region')
      END IF
      l_on_chunk = .FALSE.
      IF (var_ptr%type == REAL2d) THEN
        IF (.NOT. PRESENT(ptr2d) .AND. .NOT. PRESENT(arr2d)) CALL finish(routine, 'Must provide argument ptr2d or arr2d')
      ELSE
        IF (.NOT. PRESENT(ptr3d) .AND. .NOT. PRESENT(arr3d)) CALL finish(routine, 'Must provide argument ptr3d or arr3d')
      END IF
    END IF

    IF (l_on_chunk) THEN
      ! For chunks, it is assumed that this routine is already called within a `!$OMP PARALLEL` region for the call to `copy`
      IF (var_ptr%type == REAL2d) THEN
        IF (PRESENT(ptr1d)) THEN
          ptr1d => var_ptr%ptr2d(ics:ice,iblk)
        ELSE IF (PRESENT(arr1d)) THEN
          IF (SIZE(arr1d,1) /= nc) THEN
            CALL finish(routine, 'arr1d dimension mismatch')
          ELSE
            CALL copy(var_ptr%ptr2d(ics:ice,iblk), arr1d(:), lacc=lzacc, opt_acc_async_queue=acc_async_queue)
          END IF
        ELSE
          CALL finish(routine, 'argument ptr1d or arr1d needed for variable '//var_name//'_'//TRIM(tile_ptr%name))
        END IF
      ELSE IF (var_ptr%type == REAL3d) THEN
        IF (PRESENT(ptr2d)) THEN
          ptr2d => var_ptr%ptr3d(ics:ice,:,iblk)
        ELSE IF (PRESENT(arr2d)) THEN
          IF (SIZE(arr2d,1) /= nc .OR. SIZE(arr2d,2) /= SIZE(var_ptr%ptr3d,2)) THEN
            CALL finish(routine, 'arr2d dimension mismatch')
          ELSE
            CALL copy(var_ptr%ptr3d(ics:ice,:,iblk), arr2d(:,:), lacc=lzacc, opt_acc_async_queue=acc_async_queue)
          END IF
        ELSE
          CALL finish(routine, 'argument ptr2d or arr2d needed for variable '//var_name//'_'//TRIM(tile_ptr%name))
        END IF
      ELSE
        CALL finish(routine, 'not implemented - '//var_name//'_'//TRIM(tile_ptr%name)//' is not 2D or 3D')
      END IF
    ELSE
      IF (var_ptr%type == REAL2d) THEN
        IF (PRESENT(ptr2d)) THEN
          ptr2d => var_ptr%ptr2d(:,:)
        ELSE IF (PRESENT(arr2d)) THEN
          IF (SIZE(arr2d,1) /= SIZE(var_ptr%ptr2d,1) .OR. SIZE(arr2d,2) /= SIZE(var_ptr%ptr2d,2)) THEN
            CALL finish(routine, 'arr2d dimension mismatch')
          ELSE
!$OMP PARALLEL
            CALL copy(var_ptr%ptr2d(:,:), arr2d(:,:), lacc=lzacc, opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL
          END IF
        ELSE
          CALL finish(routine, 'argument ptr2d or arr2d needed for variable '//var_name//'_'//TRIM(tile_ptr%name))
        END IF
      ELSE IF (var_ptr%type == REAL3d) THEN
        IF (PRESENT(ptr3d)) THEN
          ptr3d => var_ptr%ptr3d(:,:,:)
        ELSE IF (PRESENT(arr3d)) THEN
          IF (SIZE(arr3d,1) /= SIZE(var_ptr%ptr3d,1) .OR. SIZE(arr3d,2) /= SIZE(var_ptr%ptr3d,2) .OR. &
            & SIZE(arr3d,3) /= SIZE(var_ptr%ptr3d,3)) THEN
            CALL finish(routine, 'arr3d dimension mismatch')
          ELSE
!$OMP PARALLEL
            CALL copy(var_ptr%ptr3d(:,:,:), arr3d(:,:,:), lacc=lzacc, opt_acc_async_queue=acc_async_queue)
!$OMP END PARALLEL
          END IF
        ELSE
          CALL finish(routine, 'argument ptr3d or arr3d needed for variable '//var_name//'_'//TRIM(tile_ptr%name))
        END IF
      ELSE
        CALL finish(routine, 'not implemented - '//var_name//'_'//TRIM(tile_ptr%name)//' is not 2D or 3D')
      END IF
    END IF

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE jsbach_get_var

#endif
END MODULE mo_jsb_interface
