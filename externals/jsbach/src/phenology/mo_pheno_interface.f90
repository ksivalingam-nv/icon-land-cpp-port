!> Contains the interfaces to the phenology processes
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
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_pheno_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,        ONLY: debug_on, acc_stream
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes PHENO_

  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pheno_tasks, global_phenology_diagnostics

  !> Type definition for update_phenology task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_phenology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_phenology    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_phenology !< Aggregates computed task variables
  END TYPE   tsk_phenology

  !> Constructor interface for update_phenology task
  !>
  INTERFACE   tsk_phenology
    PROCEDURE Create_task_phenology         !< Constructor function for task
  END INTERFACE   tsk_phenology

  !> Type definition for foliage_projected_cover task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_foliage_projected_cover
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_fpc
    PROCEDURE, NOPASS :: Aggregate => aggregate_fpc
  END TYPE tsk_foliage_projected_cover

  !> Constructor interface for foliage_projected_cover task
  !>
  INTERFACE tsk_foliage_projected_cover
    PROCEDURE Create_task_foliage_projected_cover
  END INTERFACE tsk_foliage_projected_cover

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_interface'

CONTAINS

  ! ======================================================================================================= !
  !> Constructor for phenology task
  !>
  FUNCTION Create_task_phenology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_phenology::return_ptr)
    CALL return_ptr%Construct(name='phenology', process_id=PHENO_, owner_model_id=model_id)

  END FUNCTION Create_task_phenology

  ! ======================================================================================================= !
  !> Constructor for foliage_projected_cover task
  !>
  FUNCTION Create_task_foliage_projected_cover(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_foliage_projected_cover::return_ptr)
    CALL return_ptr%Construct(name='foliage_projected_cover', process_id=PHENO_, owner_model_id=model_id)

  END FUNCTION Create_task_foliage_projected_cover

  ! ======================================================================================================= !
  !> Register tasks for pheno process
  !>
  SUBROUTINE Register_pheno_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in) :: model_id

    CALL this%Register_task(tsk_phenology(model_id))
    CALL this%Register_task(tsk_foliage_projected_cover(model_id))
  END SUBROUTINE Register_pheno_tasks

  ! ======================================================================================================= !
  !>Implementation of "update" of the "phenology" task
  !>
  !>  calls a model specific routine
  !>
  SUBROUTINE update_phenology(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    CALL update_phenology_jsbach(tile, options)
  END SUBROUTINE update_phenology

  ! ======================================================================================================= !
  !>
  !> Implementation to update the phenology state
  !! Task "phenology" calculates the phenology state of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE update_phenology_jsbach(tile, options)

    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER             :: model
    dsl4jsb_Def_config(PHENO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_phenology_jsbach'

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(PHENO_)

    IF (TRIM(dsl4jsb_Config(PHENO_)%scheme) /= 'climatology') THEN
      CALL message(TRIM(routine), 'Vegetation tile has no PFTs, using climatology for LAI')
      CALL update_phenology_climatology(tile, options)
    ELSE
      SELECT CASE (TRIM(dsl4jsb_Config(PHENO_)%scheme))
      CASE ('climatology')
        CALL update_phenology_climatology(tile, options)
      CASE DEFAULT
        CALL finish('Invalid scheme for phenology')
      END SELECT
    END IF
  END SUBROUTINE update_phenology_jsbach

  ! ================================================================================================================================
  !>
  !> Implementation of 'climatology' phenology scheme
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  !
  SUBROUTINE update_phenology_climatology(tile, options)

    USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: lai
    REAL(wp), POINTER      :: lai_mon(:,:)

    CHARACTER(len=*), PARAMETER            :: routine = modname//':update_phenology_climatology'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                              ! out
    lai_mon(1:,0:) => dsl4jsb_memory(PHENO_)%lai_mon (ics:ice, iblk, :) ! in

    !$ACC DATA CREATE(wgt1, wgt2, nmw1, nmw2) ASYNC(acc_stream)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2, acc_stream)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRESENT(wgt1, nmw1, wgt2, nmw2) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      lai(ic) = REAL(wgt1,wp) * lai_mon(ic,nmw1) + REAL(wgt2,wp) * lai_mon(ic,nmw2)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_phenology_climatology

  ! ======================================================================================================= !
  !> Implementation of "aggregate" of the "phenology" task
  !>
  SUBROUTINE aggregate_phenology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(PHENO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_phenology'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(PHENO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(PHENO_, lai,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(PHENO_, lai_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc_max, weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    !X Implementation: Start your process scheme here...

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_phenology

  ! ======================================================================================================= !
  !>update foliage projected cover
  !>
  SUBROUTINE update_fpc(tile, options)

    USE mo_jsb_time,          ONLY: get_time_interpolation_weights

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables

    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: fract_fpc
    REAL(wp), POINTER      :: fract_fpc_mon(:,:)


    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_fpc'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc)      ! out

    ! Update foliage projected cover
    !
    fract_fpc_mon(1:,0:) => dsl4jsb_memory(PHENO_)%fract_fpc_mon (ics:ice, iblk, :) ! in
    !$ACC DATA CREATE(wgt1, wgt2, nmw1, nmw2) ASYNC(acc_stream)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2, acc_stream)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRESENT(wgt1, nmw1, wgt2, nmw2) GANG VECTOR ASYNC(acc_stream)

    DO ic = 1, nc
      fract_fpc(ic) = REAL(wgt1,wp) * fract_fpc_mon(ic,nmw1) + REAL(wgt2,wp) * fract_fpc_mon(ic,nmw2)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fpc

  SUBROUTINE aggregate_fpc(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(PHENO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_fpc'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(PHENO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc, weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fpc

  !-----------------------------------------------------------------------------------------------------
  !> Calculations of diagnostic global land mean phenology output
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_phenology_diagnostics(tile, global_land_area, land_area_scaling)

    USE mo_util,                    ONLY: scaled_global_sum_2D
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile                   !< Current tile
    REAL(wp),                   INTENT(in) :: global_land_area       !< Global land area
    REAL(wp),                   INTENT(in) :: land_area_scaling(:,:) !< Land area scaling factor (accounting for halo cells)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(PHENO_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_phenology_diagnostics'

    dsl4jsb_Real2D_onDomain :: lai_ta
    dsl4jsb_Real2D_onDomain :: fract_fpc
    REAL(wp), POINTER       :: lai_ta_gmean(:)
    REAL(wp), POINTER       :: fract_fpc_gmean(:)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_var2D_onDomain(PHENO_,  lai_ta)             ! in
    dsl4jsb_Get_var2D_onDomain(PHENO_,  fract_fpc)          ! in

    lai_ta_gmean     => PHENO__mem%lai_ta_gmean%ptr(:)    ! out
    fract_fpc_gmean  => PHENO__mem%fract_fpc_gmean%ptr(:) ! out

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Calculate global land means phenology variables, if requested for output
    IF (PHENO__mem%lai_ta_gmean%is_in_output) THEN
      CALL scaled_global_sum_2D(lai_ta_gmean, lai_ta, land_area_scaling, 1.0_wp / global_land_area, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(PHENO__mem%lai_ta_gmean%ptr(:))
    END IF

    IF (PHENO__mem%fract_fpc_gmean%is_in_output) THEN
      CALL scaled_global_sum_2D(fract_fpc_gmean, fract_fpc, land_area_scaling, 1.0_wp / global_land_area, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(PHENO__mem%fract_fpc_gmean%ptr(:))
    END IF

  END SUBROUTINE global_phenology_diagnostics

#endif
END MODULE mo_pheno_interface
