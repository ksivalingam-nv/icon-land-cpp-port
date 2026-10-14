!> Contains main JSBACH model class.
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_var_class.pp-jsb.f90"

MODULE mo_jsb_model_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_kind,              ONLY: wp
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_parallel,      ONLY: Get_omp_no_of_threads, Get_omp_thread, Is_omp_inside_serial

  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_model_config
  USE mo_jsb_process_class, ONLY: t_jsb_process, t_jsb_process_p, Get_process_name
  USE mo_jsb_task_class,    ONLY: t_Task_queue, t_jsb_task_options, t_jsb_task_schedule, t_jsb_task_step, &
    &                             TASK_INTEGRATE, TASK_AGGREGATE, Get_task_action_name
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract, t_jsb_tile_abstract_p
  USE mo_jsb_impl_constants,ONLY: SHORT_NAME_LEN
  !USE mo_util,              ONLY: report_memory_usage, toupper

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_model_m, t_jsb_model
  PUBLIC :: new_model, print_model, delete_model
  PUBLIC :: MODEL_JSBACH

  !>
  !! Type for model structure
  !!
  !! This derived type defines one instance of the JSBACH model
  !!
  TYPE :: t_jsb_model
    CHARACTER(len=30)                  :: name              !< Model name
    INTEGER                            :: id                !< Model id
    CHARACTER(len=SHORT_NAME_LEN)      :: shortname         !< Model short name
    CHARACTER(len=132)                 :: description       !< Model description
    CHARACTER(len=filename_max)        :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model_config),  POINTER :: config            !< Configuration
    TYPE(t_jsb_task_options),  ALLOCATABLE :: options(:)    !< Options needed by tasks
    INTEGER                            :: grid_id
    TYPE(t_jsb_process_p),     POINTER :: processes(:)
    TYPE(t_Task_queue),        POINTER :: queue        => NULL()
      !< Task queue: the tasks of the processes in the order in which they are executed in a time step, as
      !< queued in `jsbach_init`. It only lists the tasks, without any information on tiles.
    CLASS(t_jsb_tile_abstract), POINTER :: top => NULL()          !< Top (root) tile of the tile tree
    TYPE(t_jsb_tile_abstract_p), ALLOCATABLE :: tiles(:)
      !< All tiles of the tile tree in pre-order (parents before their children)
    TYPE(t_jsb_tile_abstract_p), ALLOCATABLE :: tiles_postorder(:)
      !< All tiles of the tile tree in post-order (children before their parents)
    TYPE(t_jsb_task_schedule), ALLOCATABLE :: schedule(:)
      !< Task schedule, derived once from the task queue and the tile tree: for each queued task of an active
      !< process, the steps (tile and integrate/aggregate action) to execute the task on the tile tree
    LOGICAL                            :: use_multiple_queues = .FALSE.
      !< Use multiple OpenACC async queues for the integration of a task on different tiles
  CONTAINS
    PROCEDURE :: Configure_processes => Configure_model_processes
    PROCEDURE :: Do_fractions_change => Any_process_changes_fractions
    PROCEDURE :: Get_shortname       => Get_model_shortname
    PROCEDURE :: Get_top_tile
    PROCEDURE :: Get_process
    PROCEDURE :: Is_process_enabled
    PROCEDURE :: Queue_task
    PROCEDURE :: Build_tile_lists
    PROCEDURE :: Build_task_schedule
    PROCEDURE :: Print_task_schedule
    PROCEDURE :: Run_tasks
    PROCEDURE :: Set_options
#ifdef _OPENACC
    PROCEDURE :: gpu_to_cpu => Sync_vars_from_gpu_to_cpu
    PROCEDURE :: cpu_to_gpu => Sync_vars_from_cpu_to_gpu
#endif
  END TYPE t_jsb_model

  !>
  !! Container CLASS for JSBACH model
  !!
  !! Contains pointer to an instance of type t_jsb_model in order to allow arrays
  !! of pointers to models.
  !!
  TYPE t_jsb_model_m
    TYPE(t_jsb_model), POINTER :: m
  END TYPE t_jsb_model_m


  !>model schemes
  !!  model_scheme namlist option in mo_jsb_config_class:t_jsb_model_config
  !!  NOTE in mo_jsb_config_class:new_model_config this ENUM is not used, but replaced by INTEGER, to avoid a dependency cycle
  ENUM, BIND(C)
    ENUMERATOR :: MODEL_JSBACH=1
  END ENUM


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_class'

CONTAINS

  !>
  !! Constructor method for model instance.
  !!
  !! Subroutine to construct new instance of basic model structure.
  !!
  FUNCTION new_model(id, name, shortname, description, namelist_filename) RESULT(return_ptr)

!!$    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process

    INTEGER,                    INTENT(in) :: id
    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model), POINTER             :: return_ptr

    INTEGER :: nthreads

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:new_model'

    IF (debug_on()) CALL message(TRIM(routine), 'starting construction of new ICON-Land model instance.')

    ALLOCATE(return_ptr)

    return_ptr%name              = TRIM(name)
    return_ptr%id                = id
    return_ptr%shortname         = TRIM(shortname)
    return_ptr%description       = TRIM(description)
    return_ptr%namelist_filename = TRIM(namelist_filename)
    return_ptr%config => NULL()

    nthreads = Get_omp_no_of_threads()
    ALLOCATE(return_ptr%options(nthreads))

!!$    ALLOCATE(return_ptr%processes(max_no_of_processes))
!!$    DO iproc=1,max_no_of_processes
!!$      return_ptr%processes(iproc)%p => Create_process(iproc, return_ptr%id, TRIM(namelist_filename))
!!$    END DO

    ALLOCATE(return_ptr%queue)

    IF (debug_on()) CALL message(TRIM(routine), 'construction of ICON-Land model completed.')

  END FUNCTION new_model

  FUNCTION Get_process(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_process), POINTER :: return_ptr

    return_ptr => this%processes(iproc)%p

  END FUNCTION Get_process

  LOGICAL FUNCTION Is_process_enabled(this, iproc)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER, INTENT(in) :: iproc

    Is_process_enabled = .FALSE.
    IF (ASSOCIATED(this%processes(iproc)%p)) THEN
      IF (ASSOCIATED(this%processes(iproc)%p%config)) THEN
        Is_process_enabled = this%processes(iproc)%p%config%active
      ELSE
        ! If the process is associated but the configuration is not associated
        ! we assume here that the process is active on the model -- currently this is the case for l2a and a2l
        Is_process_enabled = .TRUE.
      ENDIF
    ENDIF

  END FUNCTION Is_process_enabled

  SUBROUTINE Queue_task(this, iproc, name)

    CLASS(t_jsb_model),        INTENT(inout) :: this
    INTEGER,                   INTENT(in)    :: iproc
    CHARACTER(LEN=*),          INTENT(in)    :: name

    CALL this%queue%Append(this%processes(iproc)%p%Get_task(TRIM(name)), Get_process_name(iproc))

  END SUBROUTINE Queue_task

  SUBROUTINE Configure_model_processes(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    INTEGER :: iproc

    DO iproc=1,SIZE(this%processes)
      IF (ASSOCIATED(this%processes(iproc)%p)) &
        & CALL this%processes(iproc)%p%Configure()
    END DO

  END SUBROUTINE Configure_model_processes

  ! ======================================================================================================= !
  !>
  !> Checks if there is an active process that changes fractions on this tile
  !>
  FUNCTION Any_process_changes_fractions(this) RESULT(l_fractions_change)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_model), INTENT(inout) :: this                 !< tile to check
    LOGICAL                           :: l_fractions_change
      !< true if there is any active land cover change (lcc) process
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: iproc
    ! -------------------------------------------------------------------------------------------------- !
    l_fractions_change = .FALSE.
    DO iproc=1,SIZE(this%processes)
      IF (ASSOCIATED(this%processes(iproc)%p)) THEN
        IF (ASSOCIATED(this%processes(iproc)%p%config)) THEN
          IF (this%processes(iproc)%p%config%active .AND. this%processes(iproc)%p%l_changes_fractions) THEN
            l_fractions_change = .TRUE.
            EXIT
          END IF
        END IF
      END IF
    END DO

  END FUNCTION Any_process_changes_fractions


  ! ====================================================================================================== !
  !>
  !> Collect all tiles of the tile tree below `this%top` in pre-order (`this%tiles`) and post-order
  !> (`this%tiles_postorder`)
  !>
  !> Must be called once after the tile tree has been completely set up (see
  !> [[mo_jsb_base:jsbach_setup_tiles]]). Loops over all tiles then simply iterate over these lists.
  !>
  SUBROUTINE Build_tile_lists(this)

    USE mo_jsb_tile_tree_class, ONLY: t_jsb_tile_tree_node_p

    CLASS(t_jsb_model), INTENT(inout) :: this

    TYPE(t_jsb_tile_tree_node_p), ALLOCATABLE :: preorder(:), postorder(:)
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Build_tile_lists'

    IF (.NOT. ASSOCIATED(this%top)) CALL finish(routine, 'Top tile not set')

    CALL this%top%Get_subtree_nodes(preorder, postorder)

    IF (ALLOCATED(this%tiles))           DEALLOCATE(this%tiles)
    IF (ALLOCATED(this%tiles_postorder)) DEALLOCATE(this%tiles_postorder)
    ALLOCATE(this%tiles(SIZE(preorder)), this%tiles_postorder(SIZE(postorder)))
    DO i=1,SIZE(preorder)
      this%tiles(i)%p           => As_tile(preorder(i)%p)
      this%tiles_postorder(i)%p => As_tile(postorder(i)%p)
    END DO

  CONTAINS

    FUNCTION As_tile(node) RESULT(tile)

      USE mo_jsb_tile_tree_class, ONLY: t_jsb_tile_tree_node

      CLASS(t_jsb_tile_tree_node), POINTER, INTENT(in) :: node
      CLASS(t_jsb_tile_abstract),  POINTER             :: tile

      SELECT TYPE (node)
      CLASS IS (t_jsb_tile_abstract)
        tile => node
      CLASS DEFAULT
        CALL finish(routine, 'Node '//TRIM(node%name)//' of tile tree is not a tile')
      END SELECT

    END FUNCTION As_tile

  END SUBROUTINE Build_tile_lists

  ! ====================================================================================================== !
  !>
  !> Build the schedule of all tasks in the task queue
  !>
  !> For each task whose process is active, the steps (tile, action) in which the task is executed on the
  !> tile tree are determined once, so that [[Run_tasks]] only needs to iterate over them. The tile tree
  !> is walked depth-first starting from the top tile; for each tile, depending on the action of the
  !> task's process on the tile (`tile%process_action`):
  !>
  !> - `SKIP_`, `INHERIT_`: the tile and its sub-tree are skipped
  !> - `ON_TILE_`, `ON_SUBTREE_`, or tile without children: the task is integrated on the tile
  !> - otherwise (`AGGREGATE_`, `ON_LEAFS_`): the children are processed first, then the task is
  !>   aggregated on the tile
  !>
  !> Must be called after the task queue has been set up and the tiles have been initialized, since the
  !> tile initialization modifies the process actions (see [[mo_jsb_model_init:jsbach_init]]).
  !>
  SUBROUTINE Build_task_schedule(this)

    USE mo_jsb_process_class, ONLY: AGGREGATE_, ON_TILE_, ON_SUBTREE_

    CLASS(t_jsb_model), INTENT(inout) :: this

    TYPE(t_Task_queue),    POINTER     :: queue_element
    TYPE(t_jsb_task_step), ALLOCATABLE :: steps(:)
    INTEGER :: ntasks, nsteps, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Build_task_schedule'

    IF (.NOT. ALLOCATED(this%tiles)) CALL finish(routine, 'Tile lists not set up')

    ! Each tile appears at most once in the steps of a task
    ALLOCATE(steps(SIZE(this%tiles)))

    ! First pass counts the scheduled tasks, second pass fills the schedule
    IF (ALLOCATED(this%schedule)) DEALLOCATE(this%schedule)
    ntasks = Scan_queue(.FALSE.)
    ALLOCATE(this%schedule(ntasks))
    ntasks = Scan_queue(.TRUE.)

    DEALLOCATE(steps)

    IF (debug_on()) CALL this%Print_task_schedule()

  CONTAINS

    !> Walk through the task queue, determine the tasks to be scheduled and their steps on the tile
    !> tree, and (if `l_fill`) store them in `this%schedule`; returns the number of scheduled tasks
    !>
    !> For each task in the queue whose process is active, the steps of the task on the tile tree are
    !> determined (in `steps`, see `Add_steps`). Tasks of inactive processes and tasks without any step
    !> (process not active on any tile) are not scheduled. The scheduled tasks keep the order of the queue.
    !>
    !> The function is called twice: first with `l_fill = .FALSE.`, only to get the number of scheduled
    !> tasks needed to allocate `this%schedule`, then with `l_fill = .TRUE.` to actually fill the schedule.
    INTEGER FUNCTION Scan_queue(l_fill)

      LOGICAL, INTENT(in) :: l_fill !< .FALSE.: only count the tasks, since `this%schedule` is not
                                    !< allocated yet; .TRUE.: also store the tasks and their steps in it

      Scan_queue = 0
      queue_element => this%queue
      DO WHILE (ASSOCIATED(queue_element))
        IF (ASSOCIATED(queue_element%task)) THEN
          iproc = queue_element%task%process_id
          IF (this%processes(iproc)%p%config%active) THEN   ! Only if process this task belongs to is active
            nsteps = 0
            CALL Add_steps(this%top)
            IF (nsteps > 0) THEN
              Scan_queue = Scan_queue + 1
              IF (l_fill) THEN
                this%schedule(Scan_queue)%task => queue_element%task
                ALLOCATE(this%schedule(Scan_queue)%steps(nsteps))
                this%schedule(Scan_queue)%steps(:) = steps(1:nsteps)
              END IF
            END IF
          END IF
        END IF
        queue_element => queue_element%next
      END DO

    END FUNCTION Scan_queue

    !> Append the steps of the current task (process `iproc`) on `tile` and its sub-tree to `steps`
    !>
    !> Recursive depth-first walk: a tile on which the process is not active is skipped together with
    !> its sub-tree; a tile on which the process is calculated (`ON_TILE_`, `ON_SUBTREE_`, or a leaf)
    !> gets an integrate step and its sub-tree is not entered; for any other tile the steps of all
    !> children are added first (in the order of the siblings), followed by an aggregate step for the
    !> tile itself. The steps are thus ordered leaves first, so that each aggregation finds its
    !> children already processed.
    RECURSIVE SUBROUTINE Add_steps(tile)

      CLASS(t_jsb_tile_abstract), POINTER, INTENT(in) :: tile !< Root tile of the sub-tree to be added

      CLASS(t_jsb_tile_abstract), POINTER :: child  !< Child tile currently processed
      INTEGER                             :: action !< Action of the process on `tile` (`SKIP_`, ...)

      action = tile%process_action(iproc)

      ! SKIP_ and INHERIT_: nothing to do on this tile and its sub-tree
      IF (action < AGGREGATE_) RETURN

      IF (tile%Has_children() .AND. action /= ON_TILE_ .AND. action /= ON_SUBTREE_) THEN
        child => tile%first_child_tile
        DO WHILE (ASSOCIATED(child))
          CALL Add_steps(child)
          child => child%next_sibling_tile
        END DO
        CALL Add_step(tile, TASK_AGGREGATE)
      ELSE
        CALL Add_step(tile, TASK_INTEGRATE)
      END IF

    END SUBROUTINE Add_steps

    !> Append one step (execute `action` of the current task on `tile`) to `steps`
    !>
    !> Aborts if `steps` is full, which cannot happen as long as each tile gets at most one step per task.
    SUBROUTINE Add_step(tile, action)

      CLASS(t_jsb_tile_abstract), POINTER, INTENT(in) :: tile   !< Tile on which the task is to be executed
      INTEGER,                             INTENT(in) :: action !< `TASK_INTEGRATE` or `TASK_AGGREGATE`

      nsteps = nsteps + 1
      IF (nsteps > SIZE(steps)) &
        & CALL finish(routine, 'Too many steps for task '//TRIM(queue_element%task%name))
      steps(nsteps)%tile   => tile
      steps(nsteps)%action =  action

    END SUBROUTINE Add_step

  END SUBROUTINE Build_task_schedule

  ! ====================================================================================================== !
  !>
  !> Print the task schedule
  !>
  SUBROUTINE Print_task_schedule(this)

    CLASS(t_jsb_model), INTENT(in) :: this

    INTEGER :: itask, istep

    CHARACTER(len=*), PARAMETER :: routine = modname//':Print_task_schedule'

    IF (.NOT. ALLOCATED(this%schedule)) RETURN

    CALL message(routine, 'Task schedule for model '//TRIM(this%name)//':')
    DO itask=1,SIZE(this%schedule)
      ASSOCIATE (task => this%schedule(itask)%task, steps => this%schedule(itask)%steps)
        CALL message('', '  Task '//TRIM(Get_process_name(task%process_id))//':'//TRIM(task%name))
        DO istep=1,SIZE(steps)
          CALL message('', '    '//Get_task_action_name(steps(istep)%action)//' on '// &
            &              TRIM(steps(istep)%tile%name))
        END DO
      END ASSOCIATE
    END DO

  END SUBROUTINE Print_task_schedule

  ! ====================================================================================================== !
  !>
  !> Run all tasks on the current chunk by executing the steps of the task schedule
  !>
  !> With OpenACC and `use_multiple_queues`, consecutive `TASK_INTEGRATE` steps of a task (on different
  !> tiles) are launched on separate async queues, which are joined into the base queue before the next
  !> `TASK_AGGREGATE` step and at the end of the task.
  !>
  SUBROUTINE Run_tasks(this, debug)

    USE mo_jsb_control, ONLY: timer_on, timer_process_task
    USE mo_timer,       ONLY: timer_start, timer_stop
#ifdef _OPENACC
    USE mo_jsb_control, ONLY: acc_stream
#endif

    CLASS(t_jsb_model), INTENT(inout) :: this
    LOGICAL, OPTIONAL,  INTENT(in)    :: debug

    LOGICAL :: l_debug, l_timer
    INTEGER :: no_omp_thread, itask, istep
#ifdef _OPENACC
    INTEGER :: base_stream, n_forked
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':Run_tasks'

    no_omp_thread = Get_omp_thread()

    l_debug = .FALSE.
    IF (PRESENT(debug)) l_debug = debug
    l_timer = timer_on('all') .OR. (timer_on() .AND. l_debug)

    IF (.NOT. ALLOCATED(this%schedule)) CALL finish(routine, 'Task schedule not set up')

    DO itask=1,SIZE(this%schedule)
      ASSOCIATE (task => this%schedule(itask)%task, steps => this%schedule(itask)%steps)

        IF (l_timer) CALL timer_start(timer_process_task(this%id))
#ifdef _OPENACC
        base_stream = acc_stream
        n_forked = 0
#endif

        DO istep=1,SIZE(steps)
          IF (l_debug) CALL message(routine, 'Task '//TRIM(task%name)//': '// &
            &   Get_task_action_name(steps(istep)%action)//' on tile '//TRIM(steps(istep)%tile%name))
#ifdef _OPENACC
          IF (this%use_multiple_queues) THEN
            IF (steps(istep)%action == TASK_INTEGRATE) THEN
              ! Launch integration on its own queue, after all previous work on the base queue
              n_forked = n_forked + 1
              acc_stream = base_stream + n_forked
              !$ACC WAIT(base_stream) ASYNC(acc_stream)
            ELSE
              CALL Join_queues()
            END IF
          END IF
#endif
          CALL task%Do_it(steps(istep)%action, steps(istep)%tile, this%options(no_omp_thread))
        END DO

#ifdef _OPENACC
        IF (this%use_multiple_queues) CALL Join_queues()
#endif
        IF (l_timer) CALL timer_stop(timer_process_task(this%id))

      END ASSOCIATE
    END DO

#ifdef _OPENACC
  CONTAINS

    !> Join all forked queues into the base queue and continue on the base queue
    SUBROUTINE Join_queues()

      INTEGER :: stream, dummy

      IF (n_forked == 0) RETURN
      DO stream = base_stream+1, base_stream+n_forked
        !$ACC WAIT(stream) ASYNC(base_stream)
      END DO
      ! Dummy kernel for the faster graph dependency resolution
      !$ACC PARALLEL ASYNC(base_stream)
      dummy = 1
      !$ACC END PARALLEL
      acc_stream = base_stream
      n_forked = 0

    END SUBROUTINE Join_queues
#endif

  END SUBROUTINE Run_tasks

  FUNCTION Get_model_name(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%name

  END FUNCTION Get_model_name

  FUNCTION Get_model_shortname(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%shortname

  END FUNCTION Get_model_shortname

  SUBROUTINE Get_top_tile(this, return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    CLASS(t_jsb_tile_abstract), POINTER, INTENT(inout) :: return_ptr

    return_ptr => this%top

  END SUBROUTINE Get_top_tile

  SUBROUTINE Set_options(this, current_datetime, dtime, alpha, iblk, ics, ice, nc, &
    &                    nsoil_e, nsoil_w, nsnow_e, npft)

    USE mo_fortran_tools,   ONLY: assign_if_present
    USE mo_jsb_time,        ONLY: t_datetime

    CLASS(t_jsb_model),                  INTENT(inout) :: this
    TYPE(t_datetime), OPTIONAL, POINTER, INTENT(in)    :: current_datetime
    REAL(wp),         OPTIONAL,          INTENT(in)    :: dtime
    REAL(wp),         OPTIONAL,          INTENT(in)    :: alpha
    INTEGER,          OPTIONAL,          INTENT(in)    :: iblk
    INTEGER,          OPTIONAL,          INTENT(in)    :: ics
    INTEGER,          OPTIONAL,          INTENT(in)    :: ice
    INTEGER,          OPTIONAL,          INTENT(in)    :: nc
    INTEGER,          OPTIONAL,          INTENT(in)    :: nsoil_e
    INTEGER,          OPTIONAL,          INTENT(in)    :: nsoil_w
    INTEGER,          OPTIONAL,          INTENT(in)    :: nsnow_e
    INTEGER,          OPTIONAL,          INTENT(in)    :: npft

    INTEGER :: no_omp_thread
    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_options'

    no_omp_thread = Get_omp_thread()

    CALL assign_if_present(this%options(no_omp_thread)%dtime,   dtime)
    CALL assign_if_present(this%options(no_omp_thread)%alpha,   alpha)
    CALL assign_if_present(this%options(no_omp_thread)%iblk,    iblk)
    CALL assign_if_present(this%options(no_omp_thread)%ics,     ics)
    CALL assign_if_present(this%options(no_omp_thread)%ice,     ice)
    CALL assign_if_present(this%options(no_omp_thread)%nc,      nc)
    IF (PRESENT(current_datetime)) this%options(no_omp_thread)%current_datetime => current_datetime
    ! Set number of soil and snow layers to options on all threads
    ! Note: The assign_if_present interface does not support assigning a value to a vector.
    IF (PRESENT(npft)) THEN
      IF (Is_omp_inside_serial()) THEN
        this%options(:)%npft = npft
      ELSE
        CALL finish(routine, 'npft may only be set outside OpenMP parallel regions')
      END IF
    END IF

    IF (PRESENT(nsoil_e)) THEN
      IF (Is_omp_inside_serial()) THEN
        this%options(:)%nsoil_e = nsoil_e
      ELSE
        CALL finish(routine, 'nsoil_e may only be set outside OpenMP parallel regions')
      END IF
    END IF

    IF (PRESENT(nsoil_w)) THEN
      IF (Is_omp_inside_serial()) THEN
        this%options(:)%nsoil_w = nsoil_w
      ELSE
        CALL finish(routine, 'nsoil_w may only be set outside OpenMP parallel regions')
      END IF
    END IF

    IF (PRESENT(nsnow_e)) THEN
      IF (Is_omp_inside_serial()) THEN
        this%options(:)%nsnow_e = nsnow_e
      ELSE
        CALL finish(routine, 'nsnow_e may only be set outside OpenMP parallel regions')
      END IF
    END IF
  END SUBROUTINE Set_options

  ! ====================================================================================================== !
  !>
  !! Destructor method for model instance.
  !!
  !! Subroutine to delete model instance from memory.
  !!
  SUBROUTINE delete_model(self)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:delete_model'

    CALL message(TRIM(routine), 'starting destruction of JSBACH model instance')

    DEALLOCATE(self)

    CALL message(TRIM(routine), 'destruction of JSBACH model instance completed.')

  END SUBROUTINE delete_model

  !>
  !! Print method for model instance.
  !!
  !! Subroutine to print description of model instance.
  !!
  SUBROUTINE print_model(self)

    USE mo_jsb_control,    ONLY: model_base_dir
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    TYPE(t_jsb_model), INTENT(in) :: self

    TYPE(t_jsb_grid), POINTER :: grid
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:print_model'

    CALL message(TRIM(routine), 'print JSBACH model description')

    WRITE(message_text,'(A,I2,A,A,A,A,A)') 'Model id: ', self%id, ', model name: ', TRIM(self%name), &
                                     ' (', TRIM(self%shortname), ')'
    CALL message('', message_text)
    CALL message('', '... '//TRIM(self%description))

    CALL message('', 'Namelist file: '//TRIM(model_base_dir)//'/'//TRIM(self%namelist_filename))

    grid => Get_grid(self%grid_id)
    IF (ASSOCIATED(grid)) CALL grid%Print

  END SUBROUTINE print_model

#ifdef _OPENACC
  SUBROUTINE Sync_vars_from_gpu_to_cpu(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_

    CLASS(t_jsb_model), INTENT(inout) :: this

    CLASS(t_jsb_memory), POINTER  :: mem
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                       :: iproc, i, itile

    CHARACTER(len=*), PARAMETER :: routine = modname//':Sync_vars_from_gpu_to_cpu'

    IF (debug_on()) CALL message(TRIM(routine), 'Syncing variables from GPU to CPU '//TRIM(this%name))

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    DO itile=1,SIZE(this%tiles)
      tile => this%tiles(itile)%p

      ! Loop over all processes and initialize memory for this tile
      DO iproc=1,SIZE(tile%mem)

        ! if the current process is not on leafs of current tile exit DO loop
        SELECT CASE (tile%process_action(iproc))
        CASE (SKIP_, INHERIT_)
          CYCLE
        END SELECT

        IF (.NOT. tile%Has_process_memory(iproc)) CYCLE

        mem => tile%mem(iproc)%p
        DO i=1,mem%no_of_vars
          IF (ASSOCIATED(mem%vars(i)%p%ptr2d)) THEN
            !$ACC UPDATE HOST(mem%vars(i)%p%ptr2d) ASYNC(1) IF(mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart)
            ! IF (mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart) CALL message(routine//': gpu2cpu', TRIM(mem%vars(i)%p%full_name))
            CONTINUE
          ELSE IF (ASSOCIATED(mem%vars(i)%p%ptr3d)) THEN
            !$ACC UPDATE HOST(mem%vars(i)%p%ptr3d) ASYNC(1) IF(mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart)
            ! IF (mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart) CALL message(routine//': gpu2cpu', TRIM(mem%vars(i)%p%full_name))
            CONTINUE
          END IF
        END DO

      END DO

    END DO

    !$ACC WAIT(1)

  END SUBROUTINE Sync_vars_from_gpu_to_cpu

  SUBROUTINE Sync_vars_from_cpu_to_gpu(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_

    CLASS(t_jsb_model), INTENT(inout) :: this

    CLASS(t_jsb_memory), POINTER  :: mem
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                       :: iproc, i, itile

    CHARACTER(len=*), PARAMETER :: routine = modname//':Sync_vars_from_cpu_to_gpu'

    IF (debug_on()) CALL message(TRIM(routine), 'Syncing variables from CPU to GPU '//TRIM(this%name))

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    DO itile=1,SIZE(this%tiles)
      tile => this%tiles(itile)%p

      ! Loop over all processes and initialize memory for this tile
      DO iproc=1,SIZE(tile%mem)

        ! if the current process is not on leafs of current tile exit DO loop
        SELECT CASE (tile%process_action(iproc))
        CASE (SKIP_, INHERIT_)
          CYCLE
        END SELECT

        IF (.NOT. tile%Has_process_memory(iproc)) CYCLE

        mem => tile%mem(iproc)%p
        DO i=1,mem%no_of_vars
          ! IF (mem%vars(i)%p%is_in_restart) CALL message(routine//': cpu2gpu', TRIM(mem%vars(i)%p%full_name))
          IF (ASSOCIATED(mem%vars(i)%p%ptr2d)) THEN
            !$ACC UPDATE DEVICE(mem%vars(i)%p%ptr2d) ASYNC(1) IF(mem%vars(i)%p%is_in_restart)
            CONTINUE
          ELSE IF (ASSOCIATED(mem%vars(i)%p%ptr3d)) THEN
            !$ACC UPDATE DEVICE(mem%vars(i)%p%ptr3d) ASYNC(1) IF(mem%vars(i)%p%is_in_restart)
            CONTINUE
          END IF
        END DO

      END DO

    END DO

    !$ACC WAIT(1)

  END SUBROUTINE Sync_vars_from_cpu_to_gpu
#endif

#endif
END MODULE mo_jsb_model_class
