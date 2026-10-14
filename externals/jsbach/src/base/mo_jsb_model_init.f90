!> Methods for JSBACH initializations
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
!>#### Initialize JSBACH instance with grid, task queue, tile structure, initial/boundary conditions etc.
!>
MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

  dsl4jsb_Use_processes SEB_, TURB_, SSE_, HYDRO_, RAD_, PHENO_, A2L_

  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(HYDRO_)

  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_control,      ONLY: get_no_of_models, debug_on
  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_grid,         ONLY: Get_grid
  USE mo_jsb_grid_class,   ONLY: t_jsb_grid

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish
  USE mo_jsb_parallel,     ONLY: Is_omp_inside_serial
  USE mo_jsb_time,         ONLY: is_time_restart, get_time_start

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart

  !> Interface to set up grid of a model instance
  !
  INTERFACE jsbach_setup_grid
    MODULE PROCEDURE setup_grid_standalone
    MODULE PROCEDURE setup_grid_patch
  END INTERFACE jsbach_setup_grid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'


CONTAINS


  ! ====================================================================================================== !
  !
  !> Set up grid for standalone model if using ECHAM infrastructure
  !
  SUBROUTINE setup_grid_standalone(nproca, nprocb, npedim, model_id, model_name)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)                    :: nproca !< Number of processors in 1st dimension
                                                     !< of 2-d decomposition
    INTEGER, INTENT(in)                    :: nprocb
        !< Number of processors in 2nd dimension of 2-d decomposition
    INTEGER, INTENT(in)                    :: npedim
        !< Size of blocks on each domain
    INTEGER, INTENT(in)                    :: model_id
        !< ID of model instance
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name
        !< Name of model instance
    ! -------------------------------------------------------------------------------------------------- !

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_grid_standalone'
    ! ---------------------------------------------------------------------------------------------------!

    ! avoid compiler warnings about dummy arguments not being used
    IF (.FALSE.) WRITE(*,*) nproca, nprocb, npedim, model_id, model_name

    CALL finish(TRIM(routine),'Not implemented for ICON')

  END SUBROUTINE setup_grid_standalone


  ! ====================================================================================================== !
  !
  !> Set up grid provided by host model (ICON coupled/standalone or ECHAM coupled)
  !
  SUBROUTINE setup_grid_patch( model_id, patch, type)

    USE mo_jsb_domain_iface,   ONLY: t_patch
    USE mo_jsb_grid_class,     ONLY: new_grid
    USE mo_jsb_grid,           ONLY: Register_grid
    USE mo_jsb_math_constants, ONLY: deg2rad

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)                  :: model_id !< ID of model instance
    TYPE(t_patch), TARGET, INTENT(inout) :: patch    !< Grid information provided by host model
    CHARACTER(LEN=*), INTENT(in)         :: type     !< Type of host model ('icon' or 'echam')
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_grid_patch'
    ! -------------------------------------------------------------------------------------------------- !

    IF (model_id > get_no_of_models()) &
      & CALL finish(routine, 'Host request to more jsbach models than are defined')

    model => Get_model(model_id)

    grid => new_grid(patch, type)
    IF (grid%ntotal_g == 1) THEN
      IF (.NOT. model%config%is_single_site) &
          & CALL finish(routine, 'Number of global cells == 1 but config not for single site')
      ! Overwrite lon/lat from single cell grid file with sitelevel model-specific namelist location
      patch%cells%center(:,:)%lon = deg2rad * model%config%single_site_lonlat(1)
      patch%cells%center(:,:)%lat = deg2rad * model%config%single_site_lonlat(2)
      grid%lon(:,:) = model%config%single_site_lonlat(1)
      grid%lat(:,:) = model%config%single_site_lonlat(2)
      grid%coslat(:,:) = COS(deg2rad * grid%lat(:,:))
      grid%sinlat(:,:) = SIN(deg2rad * grid%lat(:,:))
      grid%coslon(:,:) = COS(deg2rad * grid%lon(:,:))
      grid%sinlon(:,:) = SIN(deg2rad * grid%lon(:,:))
      grid%area(:,:) = 1._wp
    END IF

    CALL Register_grid(grid)
    model%grid_id = grid%id

    ! CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE setup_grid_patch


  ! ====================================================================================================== !
  !>
  !> Initializes an ICON-Land model instance
  !>
  !>   model_scheme: MODEL_JSBACH
  !>
  SUBROUTINE jsbach_init(model_id)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_model_class,     ONLY: print_model, MODEL_JSBACH
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, &
      &                               t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_lct_class,       ONLY: Contains_lct
    USE mo_jsb_io_netcdf,       ONLY: jsb_netcdf_open_input
    USE mo_seb_init,            ONLY: seb_init
    USE mo_turb_init,           ONLY: turb_init
    USE mo_sse_init,            ONLY: sse_init, init_soil_properties
    USE mo_hydro_init,          ONLY: hydro_init
    USE mo_pheno_init,          ONLY: pheno_init
    USE mo_rad_init,            ONLY: rad_init

    USE mo_jsb_io,              ONLY: cdiDefMissval, missval
    USE mo_util,                ONLY: one_of

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: model_id !< ID of the model instance to be initialized
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid),  POINTER          :: grid
    CLASS(t_jsb_tile_abstract), POINTER :: tile, child
    INTEGER                             :: ilct, jlct, i, no_children, lct_start, blks, blke, blk, jcs, jce
    REAL(wp), ALLOCATABLE               :: fract_tmp(:,:), tile_lct_fract_tmp(:,:)
    INTEGER                             :: itile
    LOGICAL                             :: l_fixed_fractions
    CLASS(t_jsb_aggregator), POINTER    :: weighted_by_fract

    dsl4jsb_Def_config(PHENO_)

    CHARACTER(len=*), PARAMETER         :: routine = modname//':jsbach_init'
    ! -------------------------------------------------------------------------------------------------- !

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL cdiDefMissval(missval)

    model => Get_model(model_id)

    dsl4jsb_Get_config(PHENO_)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 1. Populate task queue
    !>
    !>     @todo Should be in separate subroutine
    !>
    !>
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      CALL model%Queue_task(RAD_,    "surface_radiation")
      CALL model%Queue_task(PHENO_,  "phenology")
        ! JN moved phenology up - for allom it needs to happen before C_NPP_pot_allocation and
        ! with area changes it should happen before assimi
        ! --- before it was between humidity_scaling and foliage_projected_cover
        ! + prev comment: can be placed between the task "surface_radiation" and "asselin_filter"
      CALL model%Queue_task(PHENO_,  "foliage_projected_cover")
      ! The first call to task `snow_and_wet_fraction` currently needs to be before `surface_hydrology`,
      ! `soil_hydrology` and `update_albedo` because of pond scheme
      ! TODO: add cross-check
      CALL model%Queue_task(HYDRO_,  "snow_and_wet_fraction")
      CALL model%Queue_task(HYDRO_,  "soil_properties")
      CALL model%Queue_task(SEB_,    "surface_energy")
      CALL model%Queue_task(HYDRO_,  "evaporation")
      CALL model%Queue_task(SEB_,    "surface_fluxes")
      CALL model%Queue_task(HYDRO_,  "surface_hydrology")
      CALL model%Queue_task(HYDRO_,  "snow_and_wet_fraction")
      CALL model%Queue_task(SSE_,    "soil_and_snow_properties")
      CALL model%Queue_task(SSE_,    "soil_and_snow_temperature")
      CALL model%Queue_task(HYDRO_,  "soil_hydrology")
      CALL model%Queue_task(HYDRO_,  "canopy_cond_unstressed")
      CALL model%Queue_task(HYDRO_, "water_stress")
      CALL model%Queue_task(HYDRO_,  "canopy_cond_stressed")
      ! CALL model%Queue_task(HYDRO_,  "snow_and_wet_fraction")
      CALL model%Queue_task(TURB_,   "humidity_scaling")
      CALL model%Queue_task(TURB_,   "surface_roughness")
      CALL model%Queue_task(RAD_,    "albedo")
          ! Albedo is updated at the end of the time step, one time step before the radiation is updated in
          ! the atmosphere. Otherwise, in between, the albedo changes would be consindered in the jsbach
          ! radiation but where not considered in the atmospheric radiation calculations. This would be
          ! inconsistent. Furthermore, atmospheric cloud changes affect radiation also only at radiation
          ! time step and jsbach should be also consistent with this.
          ! In the albedo update time step the albedo should be updated after the SEB process, so that also
          ! the SEB remains consistent within this time step.
      CALL model%Queue_task(HYDRO_,  "water_balance")
      CALL model%Queue_task(SEB_,    "asselin_filter")      ! Should always be the last task
      CALL model%Queue_task(SEB_,    "energy_balance")

      ! Fuel is calculated last in a timestep to be right available in the next timestep

    END SELECT
    ! Task Queue created
    !-------------------------------------------------------------------------------------------- !

    ! Debug message after creating task queue
    IF (debug_on()) THEN
      CALL print_model(model)
    END IF

    ! At this point, the model configuration and tile structure is set up.

    ! We continue with initializing the memory structure and different model processes.

    ! ------------------------------------------------------------------------------------------- !
    !>
    !> 2. Iterate over tiles and call process initialization if necessary
    !>
    l_fixed_fractions = .NOT. model%Do_fractions_change()
    IF (l_fixed_fractions) THEN
      CALL message(routine, 'Tile fractions are fixed over time')
    ELSE
      CALL message(routine, 'Tile fractions vary over time')
    END IF

    IF (.NOT. ALLOCATED(model%tiles)) CALL finish(TRIM(routine), 'Tile tree not set up')

    ! Loop over the tiles in pre-order, i.e. parent tiles are initialized before their children
    DO itile=1,SIZE(model%tiles)
      tile => model%tiles(itile)%p

      ! ---------------------------------------------------------------------------------------------- !
      !>
      !>     2.1 Initialize each tile
      !>
      ! Initialize tile process memories and read tile fractions
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name),                &
        & grid_id = model%grid_id,                                              &
        & in_var_groups=one_of(TRIM(tile%name), model%config%output_tiles) > 0)

      ! ---------------------------------------------------------------------------------------------- !
      !>
      !>     2.2 Initialize fractions of each tile
      !>
      CALL tile%Init_fractions(                                           &
        & TRIM(model%shortname), '', TRIM(tile%name),                     &
        & l_fixed_fractions=l_fixed_fractions,                            &
        & l_rel_fractions=model%config%relative_fractions_in_file,        &
        & is_sitelevel=model%config%is_single_site )

      ! ---------------------------------------------------------------------------------------------- !
      !>
      !>     2.3 Initialize processes on each tile
      !>
      IF (model%config%init_from_ifs) &
        & model%config%ifs_input_file = jsb_netcdf_open_input(model%config%ifs_filename, model%grid_id)

      ! TBD: Loop over processes instead of explicit

      IF (tile%Is_process_active(RAD_)) THEN
          CALL rad_init(tile)
      END IF

      IF (tile%Is_process_active(SEB_)) THEN
        CALL seb_init(tile)
      END IF

      IF (tile%Is_process_active(TURB_)) THEN
        CALL turb_init(tile)
      END IF

      ! Initialize PHENO_ before HYDRO_
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(PHENO_)) THEN
          CALL pheno_init(tile)
        END IF
      END SELECT

      ! Initialize SSE_ before HYDRO_
      IF (tile%Is_process_active(SSE_)) THEN
        CALL sse_init(tile)
      END IF
      IF (tile%Is_process_active(HYDRO_)) THEN
        CALL hydro_init(tile)
      END IF
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
      ! Needs both HYDRO_ and SSE_ to continue soil initialization from IFS
        IF (tile%Is_process_active(SSE_)) THEN
          CALL init_soil_properties(tile)
        END IF
      END SELECT

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
      END SELECT

      ! ---------------------------------------------------------------------------------------------- !
    END DO

    IF (model%config%ifs_input_file%is_open) CALL model%config%ifs_input_file%Close()

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 3. Iterate over tiles again and
    !>
    !>     3.0 if demanded: correct tile fractions to assert consistency
    !>
    !>     3.1 set initial tile fractions for aggregator(s)
    !>
    !>     3.2 set LCT fractions of the specific tile
    !>
    ! It is important to only do this after all tile fractions have been set above
    !
    CALL model%Get_top_tile(tile)
    ALLOCATE(fract_tmp(SIZE(tile%lcts(1)%fract,1), SIZE(tile%lcts(1)%fract,2)))
    ALLOCATE(tile_lct_fract_tmp(SIZE(tile%lcts(1)%fract,1), SIZE(tile%lcts(1)%fract,2)))

    ! If demanded, check (and correct) the cover fractions of the direct children of each tile,
    ! starting from the top tile (i.e. box)
    IF (model%config%check_child_cf) THEN
      DO itile=1,SIZE(model%tiles)
        tile => model%tiles(itile)%p
        IF (tile%Has_children()) CALL tile%Check_child_cover_fractions()
      END DO
    END IF

    ! Loop over the tiles in post-order, i.e. the children of a tile are handled before the tile itself
    DO itile=1,SIZE(model%tiles_postorder)
      tile => model%tiles_postorder(itile)%p

      ! Set initial tile fractions for 'weighted_by_fract' aggregator
      IF (tile%Has_children()) THEN
        weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
        SELECT TYPE (weighted_by_fract)
        TYPE IS (t_jsb_aggregator_weighted_by_fract)
          CALL weighted_by_fract%Set_fractions(tile)
        END SELECT
      END IF

      ! Set (relative) land cover type fractions of the current tile

      ! The rel. fraction of the primary LCT is 1: The tile is fully covered by this LCT.
      ! Todo: Are these arrays (i.e. fract_glac_glac, fract_veg_veg, ...) of any interest, or
      ! should we remove the definition here and also remove the corresponding Add_var calls?
      tile%lcts(1)%fract(:,:) = 1._wp

      ! We start with second LCT, as the tile's primary LCT fraction has just been set.
      lct_start=2
      ! There is however an exception: The primary LCT of the box tile is the LCT of the first child.
      IF (.NOT. ASSOCIATED(tile%parent_tile)) lct_start=1
      DO ilct=lct_start,SIZE(tile%lcts)

        ! Leaf tiles should have only the primary LCT and should not enter this loop.
        IF (.NOT. tile%Has_children()) CALL finish(TRIM(routine), 'Leave tiles should only have one LCT')

        ! Tile has children, which have already been handled.
        ! Loop over the children and add their contributions to the specific LCT fraction.
        tile_lct_fract_tmp(:,:) = 0._wp
        no_children = tile%Get_no_of_children()
        child => tile%Get_first_child_tile()
        DO i=1,no_children
          IF (Contains_lct(child%lcts, tile%lcts(ilct)%id)) THEN
            ! The child contains the specific land cover type
            DO jlct=1,SIZE(child%lcts)
              IF (child%lcts(jlct)%id == tile%lcts(ilct)%id) THEN
                ! Get the (absolute) child tile fraction
                CALL child%Get_fraction(fract=fract_tmp(:,:))
                ! In case the child has only one LTC, LCT fraction and tile fraction are
                ! identical and we are done. If it has several LCTs, the fraction of the
                ! current LCT still needs to be calculated:
                IF (jlct > 1) THEN
                  fract_tmp(:,:) = fract_tmp(:,:) * child%lcts(jlct)%fract(:,:)
                END IF
                ! Sum up the children's contributions
                tile_lct_fract_tmp(:,:) = tile_lct_fract_tmp(:,:) + fract_tmp(:,:)
                EXIT
              END IF
            END DO
            ! Get the (absolute) fraction of the current tile
            CALL tile%Get_fraction(fract=fract_tmp(:,:))
            ! Calculate the (relative) LCT fraction of the tile
            WHERE (tile_lct_fract_tmp(:,:) == 0._wp)
              tile%lcts(ilct)%fract(:,:) = 0._wp
            ELSE WHERE (fract_tmp(:,:) > 0._wp)
              tile%lcts(ilct)%fract(:,:) = tile_lct_fract_tmp(:,:) / fract_tmp(:,:)
            ELSE WHERE
              tile%lcts(ilct)%fract(:,:) = 1._wp
            END WHERE
          END IF
          child => child%Get_next_sibling_tile()
        END DO
        !$ACC UPDATE DEVICE(tile%lcts(ilct)%fract) ASYNC(1)
      END DO

    END DO

    !$ACC WAIT(1)

    DEALLOCATE(fract_tmp)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 4. Iterate over tiles again and initialize hierarchical structure for memories and variables.
    !>
    ! Loop over the non-leaf tiles in post-order
    DO itile=1,SIZE(model%tiles_postorder)
      tile => model%tiles_postorder(itile)%p
      IF (.NOT. tile%Has_children()) CYCLE

      CALL tile%Init_vars()
      CALL tile%Cache_GPU_pointers()
    END DO

    ! Iterate over tiles, crosscheck tile structure, and collect conserved quantities
    DO itile=1,SIZE(model%tiles)
      tile => model%tiles(itile)%p

      IF (tile%Is_process_calculated(PHENO_)) THEN
        dsl4jsb_Get_config(PHENO_)
        IF (        tile%is_vegetation                                   &
            & .AND. TRIM(dsl4jsb_Config(PHENO_)%scheme) /= 'climatology' &
            & .AND. tile%lcts(1)%lib_id == 0)                            &
          & CALL finish(TRIM(routine), &
          &             'Must use phenology climatology on general, non-PFT VEG tile '//TRIM(tile%name))
      END IF

    END DO

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 5. Iterate over tiles again and print tile information (if debug is ON)
    !>
    IF (debug_on()) THEN
      DO itile=1,SIZE(model%tiles)
        CALL model%tiles(itile)%p%Print()
      END DO
    END IF

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 6. Build the schedule of the tasks on the tile tree, which is used in each time step.
    !>    This needs to be done after the tile initialization above, since it changes the process actions.
    !>
    CALL model%Build_task_schedule()

! #ifdef _OPENACC
!     IF (.NOT. is_time_restart(get_time_start())) CALL model%cpu_to_gpu()
! #endif

  END SUBROUTINE jsbach_init


  ! ====================================================================================================== !
  !
  !> Do some more initialization after restart files have been read
  !
  SUBROUTINE jsbach_init_after_restart(model_id)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,     ONLY: MODEL_JSBACH
    USE mo_hydro_init,          ONLY: hydro_sanitize_state
    USE mo_util,                ONLY: one_of
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: model_id !< ID of model instance to be initialized
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                             :: itile, blks, blke, blk, jcs, jce

    dsl4jsb_Def_config(HYDRO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_init_after_restart'
    ! -------------------------------------------------------------------------------------------------- !

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    model => Get_model(model_id)

    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_config(HYDRO_)
    END SELECT

    ! Iterate over tiles and call initialization if necessary
    DO itile=1,SIZE(model%tiles)
      tile => model%tiles(itile)%p

      ! Indicate to the weighted average aggregator that tile fractions might have changed after restart
      tile%l_fract_children_changed(:,:) = model%Do_fractions_change()

      ! TBD: Loop over processes instead of explicit
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (dsl4jsb_Config(HYDRO_)%sanitize_restart .AND. &
          & tile%Is_process_active(HYDRO_)) THEN
          CALL hydro_sanitize_state(tile)
        END IF
      END SELECT

    END DO

#ifdef _OPENACC
    CALL model%cpu_to_gpu()
#endif

  END SUBROUTINE jsbach_init_after_restart
  ! ====================================================================================================== !

#endif
END MODULE mo_jsb_model_init
