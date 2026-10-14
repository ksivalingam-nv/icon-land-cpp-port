!> Contains concrete tile class that contains the surface structure and memory state.
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
!>### Definition of concrete class and methods for tiles
!> This module implements 't_jsb_tile' as an extension of [[mo_jsb_tile_class:t_jsb_tile_abstract]] and
!>   implements, among others, the methods
!>
!> - Init => [[Init_tile]]
!> - Init_fractions => [[Init_fractions_on_tile]]
!> - Init_vars => [[Init_vars_on_tile]]
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_tile
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_control,         ONLY: debug_on, get_no_of_models

  USE mo_jsb_tile_tree_class, ONLY: t_jsb_tile_tree_node
  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  USE mo_jsb_process_factory, ONLY: max_no_of_processes
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task

#ifdef _OPENACC
  USE openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile, t_jsb_tile_p

  TYPE, EXTENDS(t_jsb_tile_abstract) :: t_jsb_tile           !< provides the components defined in tile_class abstract object
    TYPE(t_jsb_model),         POINTER :: model => NULL()    !< connects tile to the model object
    LOGICAL :: in_var_groups                                 !< defines whether tile is part of a output group
  CONTAINS
    ! These deferred routines must match the interfaces provided in mo_jsb_tile_class
    PROCEDURE :: Init                             => Init_tile
                 !< Initialize a tile
    PROCEDURE :: Init_fractions                   => Init_fractions_on_tile
                 !< Initialize tile fractions
    PROCEDURE :: Init_vars                        => Init_vars_on_tile
                 !< Store variable positions in memory of child tiles
    PROCEDURE :: Is_process_calculated            => Is_process_calculated_on_tile
    PROCEDURE :: Is_process_active                => Is_process_active_on_tile
    PROCEDURE :: Has_process_memory               => Has_process_memory_on_tile
    PROCEDURE :: Check_child_cover_fractions      => Check_child_cover_fractions_on_tile
    PROCEDURE :: Print                            => Print_tile
  END TYPE t_jsb_tile

  !> Used to create a vector of tiles
  ! (Define a type only containing a pointer to enable referring to objects of t_jsb_tile in vectors)
  TYPE t_jsb_tile_p
    TYPE(t_jsb_tile), POINTER :: p
  END TYPE t_jsb_tile_p

  !> Derived-type constructor used to create the tile instance
  INTERFACE t_jsb_tile
    PROCEDURE Construct_tile
  END INTERFACE t_jsb_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Constructs the tile object, as part of the tile tree, and returns the pointer to the tile
  !>
  FUNCTION Construct_tile(name, description, processes, lct, &
    & process_actions, parent, model_id, fract_filename, fract_varname) RESULT(return_ptr)

    USE mo_jsb_class,           ONLY: Get_model
    USE mo_util,                ONLY: one_of
    USE mo_jsb_lct_class,       ONLY: t_jsb_lct, LAND_TYPE, LAKE_TYPE
    USE mo_jsb_process_class,   ONLY: SEB_, A2L_, L2A_
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*),  INTENT(in)           :: name
    CHARACTER(len=*),  INTENT(in)           :: description
    INTEGER,           INTENT(in), OPTIONAL :: processes(:)
    TYPE(t_jsb_lct),   INTENT(in), OPTIONAL :: lct                !< primary lct of this tile
    INTEGER,           INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_jsb_tile), POINTER,    OPTIONAL :: parent
    INTEGER,           INTENT(in), OPTIONAL :: model_id           !< Only needed if parent is not present
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_filename     !< Name of file with tile fractions
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_varname      !< Name of fraction variable in fract_filename
    CLASS(t_jsb_tile), POINTER              :: return_ptr
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),            POINTER :: model
    INTEGER                               :: i, j, iproc
    CLASS(t_jsb_tile_tree_node),  POINTER :: parent_tmp   !< Temporary pointer to the parent tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':Construct_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting construction of tile '//TRIM(name))

    ! Allocate and initialize base structure (=tile-tree node) for tile
    ALLOCATE(t_jsb_tile::return_ptr)

    ! Add the parent tile to the actual tile (for all but the box tile)
    ! init the parent with information (properties) from the tile-tree structure only
    IF (PRESENT(parent)) THEN
      parent_tmp => parent
      CALL return_ptr%Init_node(name, description, parent_tmp)
    ELSE
      CALL return_ptr%Init_node(name, description)
    END IF

    ! Initialize further components of t_jsb_tile_abstract extended type
    return_ptr%self_tile => return_ptr

    ! Allocate parent tile
    ALLOCATE(return_ptr%process_action(max_no_of_processes))
    ! If defined in usecase, check if each process has one action (e.g., ON_TILE_ etc.)
    IF (PRESENT(processes) .AND. PRESENT(process_actions)) THEN
      IF (SIZE(processes) /= SIZE(process_actions)) &
        & CALL finish(TRIM(routine), 'Size of processes and process_actions must be equal!')
    END IF
    ! Give owner_model_id to child tiles from parent tile, hence, only box tile needs to get a owner_model_id
    IF (PRESENT(parent)) THEN
      SELECT TYPE (parent)
      CLASS IS (t_jsb_tile)
        return_ptr%parent_tile => parent
        return_ptr%owner_model_id = parent%owner_model_id
      END SELECT
    ELSE
      ! box tile already has an owner_model_id
      IF (PRESENT(model_id)) THEN
        return_ptr%owner_model_id = model_id
      ELSE
        CALL finish(TRIM(routine), 'Need model_id')
      END IF
    END IF

    ! Get model from model_id
    model => Get_model(return_ptr%owner_model_id)
    return_ptr%model => model
    ! ..
    IF (ASSOCIATED(return_ptr%parent_tile)) THEN
      IF (.NOT. ASSOCIATED(return_ptr%parent_tile%first_child_tile)) THEN
          return_ptr%parent_tile%first_child_tile => return_ptr
      END IF
    END IF
    IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
      ! to avoid further SELECT TYPE, just define the CLASS/TYPE here
      SELECT TYPE (prev_sibling => return_ptr%prev_sibling)
      CLASS IS (t_jsb_tile)
        prev_sibling%next_sibling_tile => return_ptr
      END SELECT
    END IF
    ! Define lct - land cover type
    ! Each tile can have multiple lct, the 1st one (1st element of the below vector) is the important one
    ! Other lct are "inherited" from leave tiles
    ! If an lct is not provided for a tile, this tile gets the same lct as the parent
    IF (PRESENT(lct)) THEN
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        IF (ASSOCIATED(return_ptr%parent_tile%parent_tile) .AND. (lct%id == LAND_TYPE .OR. lct%id == LAKE_TYPE)) &
          & CALL finish(TRIM(routine), 'Tiles with lct LAND_TYPE or LAKE_TYPE must be direct children of root tile.')
      END IF
      ! Add lct
      CALL return_ptr%Add_lct(lct)
      ! Add the lct of this tile to all parent tiles -> needed for aggregation of variables
      CALL return_ptr%Add_lct_to_parents(lct)
    ELSE
      ! If no lct is passed upon tile construction the primary lct of the parent tile is used
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        CALL return_ptr%Add_lct(return_ptr%parent_tile%lcts(1))
      END IF
    END IF
    ! Filename of the fract file
    IF (PRESENT(fract_filename)) THEN
      ! Use this filename
      return_ptr%fract_filename = TRIM(fract_filename)
    ELSE
      ! Use filename from namelist
      return_ptr%fract_filename = TRIM(model%config%fract_filename)
    END IF
    ! Which variable name to read the fraction from
    ! 'fract_varname' is defined in the usecase
    ! If not defined in the usecase, it is set to "empty" here, but changed later
    IF (PRESENT(fract_varname)) THEN
      return_ptr%fract_varname = TRIM(fract_varname)
    ELSE
      return_ptr%fract_varname = ''
    END IF

    ! Loop over all processes and set actions for this tile
    ! (init is going from box to bottom tile)
    ! process action of each process, as defined in the usecase
    DO i=1,max_no_of_processes
      ! Default is SKIP, i.e., update_*() routine may not do anything
      return_ptr%process_action(i) = SKIP_
      ! If process and config is associated
      IF (ASSOCIATED(model%processes(i)%p)) THEN
        IF (ASSOCIATED(model%processes(i)%p%config)) THEN
          ! Go on only if process is set ACTIVE, the default is FALSE; needs namelist setting to activate
          IF (.NOT. model%processes(i)%p%config%active) CYCLE
        END IF
      END IF
      iproc = -1
      j = 0
      IF (PRESENT(processes)) THEN
        IF (SIZE(processes) > 0) THEN
          ! j is the process that was explicitly defined in the usecase for this tile
          j = one_of(i, processes) ! i counts through the number of all processes; j is the index of the i-process within
                                   ! the explicitly given input argument `processes` for this tile (e.g. defined in init_usecase)
        END IF
      END IF
      ! The i-process was explicitly given in the input argument `processes` at position j
      IF (j > 0) THEN
        ! Get the ID of the current process i from input argument `processes`
        iproc = processes(j)
        SELECT CASE (model%config%model_scheme)
        CASE (MODEL_JSBACH)
          IF (iproc == SEB_ .AND. ALLOCATED(return_ptr%lcts)) THEN
            ! TODO hard coded checks of process dependencies - should not be handled here
            IF (return_ptr%lcts(1)%id /= LAND_TYPE .AND. return_ptr%lcts(1)%id /= LAKE_TYPE) THEN
              CALL finish(TRIM(routine), 'SEB_ process only allowed on tiles with lct LAND_TYPE or LAKE_TYPE.')
            END IF
          END IF
        END SELECT
        IF (PRESENT(process_actions)) THEN
          ! Use the given `process_actions(j)` for this tile and process `iproc` ...
          CALL return_ptr%Set_process_action(iproc, process_actions(j), action_to_parents=AGGREGATE_)
        ELSE          ! use defaults if process_actions are not given
          ! ... otherwise set process action for this tile and process to ON_LEAFS (except to ON_TILE_ for A2L_ and L2A_ processes)
          IF (iproc == A2L_) THEN
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)
          ELSE IF (iproc == L2A_) THEN
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)
          ELSE
            CALL return_ptr%Set_process_action(iproc, ON_LEAFS_, action_to_parents=AGGREGATE_)
          END IF
        END IF
      ! If the i-process was NOT explicitly defined for this tile and this tile is not the top tile,
      ! then use the process action from the parent ... except the process action of the parent is ON_TILE_ in
      ! which case we set the action on this tile to INHERIT_ , resp. SKIP for process L2A_
      ELSE IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        iproc = i
        SELECT CASE (return_ptr%parent_tile%process_action(iproc))
        CASE (SKIP_, INHERIT_, ON_SUBTREE_, ON_LEAFS_)
          CALL return_ptr%Set_process_action(iproc, return_ptr%parent_tile%process_action(iproc))
        CASE (ON_TILE_)
          IF (iproc == L2A_) THEN
            ! For the L2A_ process we want to do nothing on all descendant tiles of the box tile
            CALL return_ptr%Set_process_action(iproc, SKIP_)
          ELSE
            CALL return_ptr%Set_process_action(iproc, INHERIT_)
          END IF
        END SELECT
      END IF

    END DO

    IF (debug_on()) CALL message(TRIM(routine), 'Finished construction of tile '//TRIM(name))

  END FUNCTION Construct_tile

  ! ======================================================================================================= !
  !>
  !> Init tile - in particular initialise memory for processes on the tile
  !>
  ! Note: It is important that Init_tile is only called for the tiles after the hierarchical tile tree has
  !       been completely constructed. Init_tile is called for each tile from jsbach_init (CALL tile%Init,
  !       in mo_jsb_model_init.f90). The corresponding Init_interface is defined in mo_jsb_tile_class.
  SUBROUTINE Init_tile(this, varlist_name, prefix, suffix, grid_id, in_var_groups)

    USE mo_jsb_process_factory, ONLY: Create_process_memory
    USE mo_jsb_process_class,   ONLY: Get_process_name
    USE mo_jsb_memory_class,    ONLY: t_jsb_memory

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this              !< Current tile
    CHARACTER(len=*),  INTENT(in)    :: varlist_name      !< Model shortname
    CHARACTER(len=*),  INTENT(in)    :: prefix            !< Prefix for variable names (generally empty)
    CHARACTER(len=*),  INTENT(in)    :: suffix            !< Suffix for variable names
    INTEGER,           INTENT(in)    :: grid_id           !< Grid id
    LOGICAL,           INTENT(in)    :: in_var_groups     !< Should output be written for this tile?
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER                               :: iproc
    TYPE(t_jsb_model), POINTER            :: model
    CLASS(t_jsb_memory), POINTER          :: tmp_ptr
    CLASS(t_jsb_tile_abstract), POINTER   :: parent_tile => NULL()
    CHARACTER(LEN=:), ALLOCATABLE         :: process_name

    CHARACTER(len=:), ALLOCATABLE :: varlist_name_loc     ! Local varlist name
    CHARACTER(len=3)              :: model_id_suffix      ! Suffix for domain if more than one model
    INTEGER                       :: no_of_children, child_no

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile '//TRIM(this%name))

    WRITE(model_id_suffix, '(A1,I2.2)') 'D', this%owner_model_id

    model => this%model
    this%grid_id = grid_id                ! Defined for convenience

    ! Allocate all memory components for new tile and initialize their data part with NULL()
    ALLOCATE(this%mem(max_no_of_processes))
    !$ACC ENTER DATA CREATE(this%mem)

    ! Loop over all processes and initialize memory for this tile
    parent_tile => this%parent_tile       ! Defined for convenience
    DO iproc=1,max_no_of_processes

      this%mem(iproc)%p => NULL()

      IF (this%process_action(iproc) == SKIP_) CYCLE      ! Nothing needs to be done

      ! If current tile is a leaf, we replace ON_LEAF by ON_TILE and set the parent tile to AGGREGATE
      IF (this%process_action(iproc) == ON_LEAFS_) THEN
        IF (ASSOCIATED(parent_tile)) parent_tile%process_action(iproc) = AGGREGATE_
        IF (.NOT. this%Has_children()) this%process_action(iproc) = ON_TILE_
      END IF

      ! Not all processes need memory (compare process_factory):
      ! Only create / link memory for processes that have a memory
      IF (ASSOCIATED(model%processes(iproc)%p)) THEN
        IF (.NOT. model%processes(iproc)%p%has_memory) CYCLE
      END IF

      process_name = Get_process_name(iproc)

      IF (this%process_action(iproc) == INHERIT_) THEN
        ! We do not need to generate memory, just use the memory of the parent.
        this%mem(iproc)%p => parent_tile%mem(iproc)%p
        ! No output is needed for the tile, as the parent output is identical (???).
        this%mem(iproc)%p%in_var_groups = in_var_groups
        __acc_attach(this%mem(iproc)%p)
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Linked memory for process '//TRIM(process_name)//&
          &                             ' on tile '//TRIM(this%name)//' to parent')
      ELSE
        tmp_ptr => Create_process_memory(iproc) ! Create memory for the specific process
        __acc_attach(tmp_ptr)
        this%mem(iproc)%p => tmp_ptr            !TODO: Still needed? Same as this%mem(iproc)%p%owner_proc_id  (see below)
        __acc_attach(this%mem(iproc)%p)         !      (Needs to be tested on gpu)
        this%mem(iproc)%p%id = iproc

        ! Preparations needed for the Init process memory function called below
        ! Generate variable list name. Each process has its own variable list; used in add_var
        IF (get_no_of_models() > 1) THEN
          varlist_name_loc = TRIM(varlist_name)//'_'//TRIM(process_name)//'_'//model_id_suffix
        ELSE
          varlist_name_loc = TRIM(varlist_name)//'_'//TRIM(process_name)
        END IF
        this%mem(iproc)%p%varlist_name = TRIM(varlist_name_loc)

        ! Also the memory of the parent tile is available on the current tile (if it has a parent).
        IF (ASSOCIATED(parent_tile)) this%mem(iproc)%p%parent => parent_tile%mem(iproc)%p

        this%mem(iproc)%p%grid_id = grid_id
        this%mem(iproc)%p%owner_model_id = this%owner_model_id
        this%mem(iproc)%p%owner_proc_id = iproc
        this%mem(iproc)%p%owner_proc_name = process_name
        ALLOCATE(this%mem(iproc)%p%owner_tile_path(SIZE(this%path)))
        this%mem(iproc)%p%owner_tile_path(:) = this%path(:)
        this%mem(iproc)%p%owner_tile_name = TRIM(this%name)
        this%mem(iproc)%p%in_var_groups = in_var_groups    ! Must be before CALL to %Init !
        IF (ASSOCIATED(model%processes(iproc)%p%config)) THEN
          IF (this%Has_process_memory(iproc)) THEN
            this%mem(iproc)%p%lrestart_cont = model%processes(iproc)%p%config%lrestart_cont
          END IF
        END IF
        !$ACC UPDATE DEVICE(this%mem(iproc)%p) ASYNC(1)
        !$ACC ENTER DATA CREATE(this%mem(iproc)%p%vars)
        IF (prefix /= '') THEN
          !TODO: check if prefix is still needed ...
          CALL this%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), &
            &     TRIM(suffix), this%lcts(:)%id, this%lcts(1)%lib_id, this%owner_model_id)
        ELSE
          CALL this%mem(iproc)%p%Init(TRIM(process_name), &
            &     TRIM(suffix), this%lcts(:)%id, this%lcts(1)%lib_id, this%owner_model_id)
        END IF
        !$ACC UPDATE DEVICE(this%mem(iproc)%p%vars) ASYNC(1)
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(this%name))

        ! Make the current tile's memory available from the parent tile
        IF (ASSOCIATED(parent_tile)) THEN                               ! Tile has a parent (i.e. this is not the box tile)
          no_of_children = parent_tile%Get_no_of_children()             ! Total number of siblings
          IF (.NOT. ASSOCIATED(parent_tile%mem(iproc)%p%children) .AND. no_of_children > 0) THEN
            parent_tile%mem(iproc)%p%no_of_children = no_of_children    !TODO: this line is probably not needed
            ALLOCATE(parent_tile%mem(iproc)%p%children(no_of_children)) ! Allocate memory for all siblings
            !$ACC ENTER DATA CREATE(parent_tile%mem(iproc)%p%children)
          END IF
          ! Index of the child (i.e. the current tile; comes from tile-tree structure)
          child_no = this%path(this%level)
          parent_tile%mem(iproc)%p%children(child_no)%p => this%mem(iproc)%p
          !$ACC UPDATE DEVICE(parent_tile%mem(iproc)%p%children(child_no)%p) ASYNC(1)
          __acc_attach(parent_tile%mem(iproc)%p%children(child_no)%p)
        END IF

      END IF

    END DO
    !$ACC UPDATE DEVICE(this%mem) ASYNC(1)
    !$ACC WAIT(1)

  END SUBROUTINE Init_tile

  ! ======================================================================================================= !
  !>
  !> Add and initialise fraction variables for this tile - either from a file or according to specification
  !>
  SUBROUTINE Init_fractions_on_tile(this, varlist_name, prefix, suffix, l_fixed_fractions, l_rel_fractions, is_sitelevel)

    USE mo_jsb_tile_class,      ONLY: t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_varlist,         ONLY: jsb_add_var, VARNAME_LEN
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_INSTANT, TSTEP_CONSTANT
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this               !< Current tile
    CHARACTER(len=*),  INTENT(in)    :: varlist_name       !< Model shortname
    CHARACTER(len=*),  INTENT(in)    :: prefix             !< Prefix for variable names (generally empty)
    CHARACTER(len=*),  INTENT(in)    :: suffix             !< Suffix for variable names
    LOGICAL,           INTENT(in)    :: l_fixed_fractions
      !< If tile fractions can change (active lcc process). I.a. determines if fracts go to restart file
    LOGICAL,           INTENT(in)    :: l_rel_fractions
      !< Tile fractions in input file are usually relative with respect to parent (historical reasons)
      !< while the fractions in the code are absolute (with respect to the box tile)
    LOGICAL,           INTENT(in)    :: is_sitelevel       !< Is this a site-level model configuration?
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator_weighted_by_fract
    CLASS(t_jsb_aggregator),                  POINTER :: aggregator
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER                    :: table
    TYPE(t_input_file)         :: input_file
    REAL(wp),          POINTER :: return_pointer(:,:), tmp_fract(:,:)
    CHARACTER(len=50)          :: name_loc
    CHARACTER(len=:), ALLOCATABLE :: varlist_name_loc
    CHARACTER(len=3)           :: model_id_suffix                    ! Suffix for domain if more than one model

    INTEGER :: isteptype, i, ilct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fractions_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing fractions on tile '//TRIM(this%name))

    hgrid   => Get_grid(this%grid_id)
    surface => Get_vgrid('surface')

    table = tables(1)

    ! -------------------------------------------------------------------------------------------------- !
    ! Add fraction(s) of this tile to varlist

    ! isteptype us used in the add_var below
    ! - if constant cdi knows that its not time-dependent
    IF (l_fixed_fractions) THEN
      isteptype = TSTEP_CONSTANT
    ELSE
      isteptype = TSTEP_INSTANT
    END IF

    WRITE(model_id_suffix, '(A1,I2.2)') 'D', this%owner_model_id
    IF (get_no_of_models() > 1) THEN
      varlist_name_loc = TRIM(varlist_name)//'_'//model_id_suffix
    ELSE
      varlist_name_loc = TRIM(varlist_name)
    END IF

    ! Note: uses a much simpler add_var routine than process variables do
    ! Further, the fractions are direct components of the tile class and do not 'live' on a process memory
    ! Therefore, as opposed to the process variables a temporal pointer (tmp_fract) is used here,
    ! where variables on the process memory already use the corresponding pointer on the memory
    ! In addition, there is also no process prefix on this variable in the output
    CALL jsb_add_var(varlist_name_loc, 'fract', tmp_fract,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('tile_fraction', '',                                                                             &
      &      'fraction of grid box covered by tile '//TRIM(this%name)),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix, lrestart=.NOT. l_fixed_fractions, lrestart_cont=.TRUE.,                               &
      & initval_r=0.0_wp, isteptype=isteptype, groups=[character(len=VARNAME_LEN) :: 'jsb_tile_fractions'])

    ! Now we need to set the tile component fract pointer to the just allocated memory
    CALL this%Set_fract_ptr(tmp_fract)
    NULLIFY(tmp_fract)

    CALL jsb_add_var(varlist_name_loc, 'fract_old', tmp_fract,                                 &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_old', '',                                                                          &
      &      'old fraction of grid box covered by tile '//TRIM(this%name)),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix, loutput=.NOT. l_fixed_fractions, lrestart=.NOT. l_fixed_fractions, lrestart_cont=.TRUE., &
      & initval_r=0.0_wp, isteptype=isteptype, groups=[character(len=VARNAME_LEN) :: 'jsb_tile_fractions'])
    CALL this%Set_fract_old_ptr(tmp_fract)
    NULLIFY(tmp_fract)


    ! -------------------------------------------------------------------------------------------------- !
    ! Initialize fraction(s) of this tile from input file

    ALLOCATE(this%l_fract_children_changed(hgrid%nproma, hgrid%nblks))
    this%l_fract_children_changed(:,:) = .TRUE. ! default: fraction of tile is changing

    ALLOCATE(tmp_fract(hgrid%nproma, hgrid%nblks))
    IF (TRIM(this%fract_varname) == 'equal') THEN ! fract_varname == 'equal' ==> each child tile gets the same fraction
      IF (ASSOCIATED(this%parent_tile)) THEN ! The tile is not the box tile.
        ! The tile fraction is set to the fraction of the parent tile divided by
        ! the number of children of the parent tile (i.e. number of sibling tiles + 1)
        CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        tmp_fract(:,:) = tmp_fract(:,:) / REAL(this%parent_tile%Get_no_of_children(), wp)
      ELSE ! tile is the box tile, fraction of the tile is set to 1
        tmp_fract(:,:) = 1._wp
      END IF
      ! Set the fraction of the tile
      CALL this%Set_fraction(fract=tmp_fract(:,:))
    ELSE IF (this%fract_varname == 'allFirst') THEN
      IF(ASSOCIATED(this%parent_tile)) THEN
        CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        IF (.NOT. this%Is_first_child()) THEN
          tmp_fract(:,:) = 0._wp
        END IF
      ELSE
        tmp_fract(:,:) = 1._wp
      END IF
      CALL this%Set_fraction(fract=tmp_fract(:,:))
    ELSE
      input_file = jsb_netcdf_open_input(TRIM(this%fract_filename), this%grid_id)
      ! Determine the name of the tile fraction variable in the input file
      IF (TRIM(this%fract_varname) == '') THEN
        name_loc = 'fract'
        IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
        IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)
      ELSE
        name_loc = TRIM(this%fract_varname)
      END IF
      ! Read fraction of the tile
      IF (input_file%Has_var(TRIM(name_loc))) THEN
        return_pointer => input_file%Read_2d(variable_name=TRIM(name_loc))
      ELSE
        CALL finish(TRIM(routine), 'Fraction '//TRIM(name_loc)//' for tile '//TRIM(this%name)// &
          &                        ' not found in '//TRIM(this%fract_filename))
      END IF
      ! Set fraction of the tile
      IF (ASSOCIATED(this%parent_tile) .AND. l_rel_fractions) THEN ! l_rel_fractions == .true. is default, fraction read from
                                                                   ! input file is multiplied by fraction of parent tile, i.e.
                                                                   ! fraction in the model is relative to box tile, but
                                                                   ! fraction in the input file is relative to parent tile
        CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        CALL this%Set_fraction(fract=MAX(0._wp, return_pointer(:,:)*tmp_fract(:,:)))
      ELSE
        IF (is_sitelevel .AND. .NOT. ASSOCIATED(this%parent_tile)) THEN
          ! For site-level, set fraction of box tile to 1
          tmp_fract(:,:) = 1._wp
          CALL this%Set_fraction(fract=tmp_fract(:,:))
        ELSE
          CALL this%Set_fraction(fract=return_pointer(:,:))
        END IF
        CALL this%Get_fraction(fract=hgrid%lsf(:,:))
        hgrid%lsm(:,:) = hgrid%lsf(:,:) > 0._wp
        !$ACC UPDATE DEVICE(hgrid%lsf, hgrid%lsm) ASYNC(1)
      END IF
      ! Close the input file
      CALL input_file%Close()

    END IF

    ! If tile has children, create and register aggregators
    ! Presenty only the "weighted_by_fract" aggregator is used and registered
    IF (this%Has_children()) THEN
      aggregator_weighted_by_fract => t_jsb_aggregator_weighted_by_fract(this%grid_id, this%Get_no_of_children())
      CALL this%Register_aggregator(aggregator_weighted_by_fract)
    END IF

    ! If tile has parent, put fractions of current tile into corresponding aggregator of parent
    IF (ASSOCIATED(this%parent)) THEN ! tile is not the box tile
      i = this%Get_pos_of_child() ! index i is used to put fractions of the current tile into aggregator of parent
                                  ! at the right place with respect ot the siblings of the current tile
      ! Parent must be of the correct class (should always be the case since only the abstract class tile is "above")
      SELECT TYPE (parent => this%parent)
      CLASS IS (t_jsb_tile)
        aggregator => parent%Get_aggregator('weighted_by_fract')
        ! If parent tile's "weighted_by_fract" aggregator found and of the right type (NOT class?!) set the fraction of
        !   the correct sibling of this aggregator to the tile's current fraction
        IF (ASSOCIATED(aggregator)) THEN
          SELECT TYPE (aggregator)
          TYPE IS (t_jsb_aggregator_weighted_by_fract)
            ! Get fractions for aggregation (i: index of the child)
            ! These fraction might be updated later on, e.g. in case of land use change.
            CALL this%Get_fraction(aggregator%fractions(:,:,i))
          END SELECT
        END IF
      END SELECT
    END IF

    ! Loop over lcts of this tile and add variables for fractions for each lct,
    ! i.e. allocate memory for the current tile and all its children
    ! ( lcts(:) are a part of TYPE(t_jsb_tile_abstract) defined in mo_jsb_tile_class.f90 (s. comments there)
    ! and lcts(:) are TYPE(t_jsb_lct) defined in mo_jsb_lct_class.f90 )
    !$ACC ENTER DATA CREATE(this%lcts)
    DO ilct=1,SIZE(this%lcts)

      !! !$acc enter data create(this%lcts(ilct)%fract)
      CALL jsb_add_var(varlist_name_loc, 'fract_'//TRIM(this%lcts(ilct)%name), this%lcts(ilct)%fract,             &
        & hgrid, surface,                                                                                         &
        & t_cf('lct_fraction', '',                                                                                &
        &      'fraction of '//TRIM(this%name)//' covered by LCT '//TRIM(this%lcts(ilct)%name)),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix, lrestart=.NOT. l_fixed_fractions, lrestart_cont=.TRUE.,                                 &
        & initval_r=0.0_wp, isteptype=isteptype, groups=[character(len=VARNAME_LEN) :: 'jsb_tile_fractions'])
      __acc_attach(this%lcts(ilct)%fract)

    END DO

    !$ACC WAIT(1)
    DEALLOCATE(tmp_fract)

  END SUBROUTINE Init_fractions_on_tile

  ! ======================================================================================================= !
  !>
  !> Initialise the child_idx vector for each variable on each process running on child tiles
  !> - this vector then collects the process memory position of the according var on each child of the tile
  !>
  !@TODO: this routine should be renamed!
  !
  SUBROUTINE Init_vars_on_tile(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< the current tile
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER  :: mem  !< temporary: pointer to current memory
    CHARACTER(len=VARNAME_LEN)    :: name !< temporary: name of the current variable
    INTEGER                       :: iproc, i, no_of_children, child_no

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_vars_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile structure on vars '//TRIM(this%name))

    no_of_children = this%Get_no_of_children()
    ! Assertion: this function should not be called for a leaf
    IF (no_of_children == 0) CALL finish(routine, 'This should not happen')

    ! Loop over all processes and set-up the bookkeeping
    DO iproc=1,max_no_of_processes

      ! If the current process is not on leaves of current tile exit DO loop
      SELECT CASE (this%process_action(iproc))
      CASE (SKIP_, ON_TILE_, INHERIT_)
        CYCLE
      END SELECT

      ! If the current process does not have memory we can skip it
      IF (.NOT. ASSOCIATED(this%mem(iproc)%p)) CYCLE

      ! Do the bookkeeping for all variables on this process (find positions of the vars on the children)
      mem => this%mem(iproc)%p
      __acc_attach(mem)
      DO i=1,mem%no_of_vars
        name = mem%vars(i)%p%name
        ! child_idx gets allocated on first use
        IF (.NOT. ALLOCATED(mem%vars(i)%p%child_idx)) THEN
          ALLOCATE(mem%vars(i)%p%child_idx(no_of_children))
          !$ACC ENTER DATA CREATE(mem%vars(i)%p%child_idx)
        END IF
        DO child_no=1,no_of_children
          ! Check if this child has this memory
          IF (ASSOCIATED(mem%children(child_no)%p)) THEN
            ! Get the position from the linked list
            mem%vars(i)%p%child_idx(child_no) = mem%children(child_no)%p%Get_var_position(TRIM(name))
          ELSE
            ! The variable does not exist on this child
            mem%vars(i)%p%child_idx(child_no) = 0
          END IF
        END DO
        !$ACC UPDATE DEVICE(mem%vars(i)%p%child_idx) ASYNC(1)
      END DO

    END DO
    !$ACC WAIT(1)

  END SUBROUTINE Init_vars_on_tile

  ! ====================================================================================================== !
  !>
  !> Check if a given process is to be integrated on the tile:
  !>  i.e. if the process action on this tile is one of ON_TILE_ or ON_SUBTREE_
  !>
  LOGICAL FUNCTION Is_process_calculated_on_tile(this, iproc)

    USE mo_util,                ONLY: one_of ! one_of(arg1,arg2(:)), returns i, if arg2(i) is equal to arg1, and -1 otherwise
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this  !< tile on which the process activity is to be checked
    INTEGER,           INTENT(in) :: iproc !< process for which the activity is to be checked
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Is_process_calculated_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    ! Note: It's not necessary to test for ON_LEAFS_ because, after the tile tree has been completely
    !       set up, ON_LEAFS_ on a tile has been replaced by AGGREGATE_ and ON_TILE_ on the sub-tree.
    Is_process_calculated_on_tile = (one_of(this%process_action(iproc), (/ON_TILE_, ON_SUBTREE_/)) > 0)

  END FUNCTION Is_process_calculated_on_tile

  ! ====================================================================================================== !
  !>
  !> Check if a given process is to be active (integrated or aggregated) on the tile:
  !>  i.e. if the process action on this tile is not SKIP_ or INHERIT_
  !>
  LOGICAL FUNCTION Is_process_active_on_tile(this, iproc)

    USE mo_util,                ONLY: one_of ! one_of(arg1,arg2(:)), returns i, if arg2(i) is equal to arg1, and -1 otherwise
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this  !< tile on which the process activity is to be checked
    INTEGER,           INTENT(in) :: iproc !< process for which the activity is to be checked
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Is_process_active_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    Is_process_active_on_tile = (one_of(this%process_action(iproc), (/SKIP_, INHERIT_/)) < 0)

  END FUNCTION Is_process_active_on_tile

  ! ====================================================================================================== !
  !
  !> Check if the tile has memory for the process
  !
  LOGICAL FUNCTION Has_process_memory_on_tile(this, iproc)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this  !< tile for which the process memory is to be checked
    INTEGER,           INTENT(in) :: iproc !< process for which it is checked if memory is available
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Has_process_memory_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    Has_process_memory_on_tile = ASSOCIATED(this%mem(iproc)%p)

  END FUNCTION Has_process_memory_on_tile

  ! ======================================================================================================= !
  !>
  !> Checks and if required corrects the cover fractions of the children of a tile
  !>
  SUBROUTINE Check_child_cover_fractions_on_tile(this)
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid
    USE mo_jsb_grid,            ONLY: Get_grid
    USE mo_util,                ONLY: real2string
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< Tile for which the child cover fractions are checked
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),  POINTER :: current_child
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    INTEGER :: nproma, nblks, i_tile
    REAL(wp),           POINTER :: tmp_fract(:,:)
    REAL(wp),       ALLOCATABLE :: sum_fract(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':Check_child_cover_fractions_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    ! Only relevant for tiles with children
    IF (.NOT. this%Has_children()) RETURN

    IF (debug_on()) CALL message(TRIM(routine), TRIM(this%name))

    hgrid   => Get_grid(this%grid_id)
    nproma = hgrid%nproma
    nblks  = hgrid%nblks

    ALLOCATE(tmp_fract(hgrid%nproma, hgrid%nblks))
    ALLOCATE(sum_fract(hgrid%nproma, hgrid%nblks))
    sum_fract(:,:) = 0.0_wp

    current_child => this%Get_first_child_tile()

    i_tile = 0
    DO WHILE (ASSOCIATED(current_child))
      i_tile = i_tile + 1
      CALL current_child%Get_fraction(fract=tmp_fract(:,:))

      tmp_fract(:,:) = MIN(1.0_wp, tmp_fract(:,:))
      tmp_fract(:,:) = MAX(0.0_wp, tmp_fract(:,:))

      CALL current_child%Set_fraction(fract=tmp_fract(:,:))

      sum_fract(:,:) = sum_fract(:,:) + tmp_fract(:,:)
      current_child => current_child%Get_next_sibling_tile()
    END DO

    ! Assert that sum is "equal" to tile fraction (whereby the threshold is arbitraryly chosen here)
    CALL this%Get_fraction(fract=tmp_fract(:,:))
    IF (ANY(ABS((tmp_fract(:,:) - sum_fract(:,:))) > 1.e-15_wp)) THEN

      CALL finish(TRIM(routine), 'Fractions of '//TRIM(this%name)//' not equal to sum of child cf.'// &
        &  ' MAX diff: ' //real2string(ABS(MAXVAL(tmp_fract(:,:) - sum_fract(:,:)))))
    END IF

    DEALLOCATE(tmp_fract, sum_fract)

  END SUBROUTINE Check_child_cover_fractions_on_tile


  ! ======================================================================================================= !
  !>

  ! ======================================================================================================= !
  !>
  !> Routine that prints information about the tile to the log-file
  !>
  SUBROUTINE Print_tile(this)

    USE mo_jsb_process_class,  ONLY: Get_process_name, Get_action_name
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this !< tile to be printed
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: ilct, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Print_tile'
    ! -------------------------------------------------------------------------------------------------- !
    ! Print for basic type
    CALL this%Print_node() ! print position of the tile in the whole hierarchical tree of tiles as a series of integers

    IF (this%fract_filename /= '') THEN
      CALL message('     fract_filename', TRIM(this%fract_filename))
    END IF

    ! Additional prints for tile type
    IF (ALLOCATED(this%lcts)) THEN
      DO ilct=1,SIZE(this%lcts)
        IF (this%lcts(ilct)%lib_id > 0) THEN
          WRITE(message_text,*) TRIM(this%lcts(ilct)%Get_name()),' ... ',TRIM(this%lcts(ilct)%Get_longname()), &
            &                   ' (',this%lcts(ilct)%lib_id,')'
        ELSE
          WRITE(message_text,*) TRIM(this%lcts(ilct)%Get_name())//' ... '//TRIM(this%lcts(ilct)%Get_longname())
        END IF
        CALL message('     LCT', message_text) ! print names of lcts of the tile
      END DO
    END IF

    DO iproc=1,max_no_of_processes
      IF (this%process_action(iproc) > 0) THEN
        ! Print names of processes running on the tile
        CALL message('     Process '//Get_process_name(iproc), 'Action '//Get_action_name(this%process_action(iproc)))
      END IF
    END DO

  END SUBROUTINE Print_tile

#endif
END MODULE mo_jsb_tile
