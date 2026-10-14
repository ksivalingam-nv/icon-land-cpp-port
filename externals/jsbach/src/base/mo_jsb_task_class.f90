!> Contains process task class
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
!>### Definition of abstract class and methods for ICON-Land process tasks
!>
!> The module defines the abstract process task class [[t_jsb_process_task]], which is the base class
!> for all process tasks in ICON-Land. Each task belongs to the task list of a process [[t_jsb_process]]
!> and can have timers for the `Integrate` and `Aggregate` methods.
!>
!> It provides the interface for task-specific implementations of the `Integrate` and `Aggregate`
!> methods, as well as the method `Do_it` to execute a task on a tile and adding the task to a task
!> queue (`Append`).
!>
!> The module also defines the task schedule [[t_jsb_task_schedule]]: for each task of the queue the
!> sequence of steps (tile and action `TASK_INTEGRATE` or `TASK_AGGREGATE`) in which the task is executed
!> on the tile tree. The schedule is computed once during initialization (see
!> [[mo_jsb_model_class:Build_task_schedule]]).
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_task_class
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish
  USE mo_timer,       ONLY: new_timer, timer_start, timer_stop    ! Functions from timer infrastructure

  USE mo_jsb_control, ONLY: timer_aggregate, timer_integrate, timer_on, debug_on ! Objects used for timing

  USE mo_jsb_tile_class,  ONLY: t_jsb_tile_abstract               ! Tile abstract class
  USE mo_jsb_time,        ONLY: t_datetime                        ! DateTime class from mtime

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process_task, t_jsb_process_task_p, t_Task_queue, t_jsb_task_options
  PUBLIC :: t_jsb_task_step, t_jsb_task_schedule, TASK_INTEGRATE, TASK_AGGREGATE, Get_task_action_name

  !> Actions of a task on a tile
  ENUM, BIND(C)
    ENUMERATOR :: TASK_INTEGRATE = 1, TASK_AGGREGATE
  END ENUM

  TYPE t_jsb_task_options
    ! Set of parameters that are needed by tasks during their execution
    REAL(wp)                  :: dtime   = -1._wp !< Time step length [s]
    REAL(wp)                  :: alpha   = -1._wp !< Implicitness factor
    TYPE(t_datetime), POINTER :: current_datetime !< Date and time of current time step
    INTEGER                   :: iblk    = -1     !< Number of current block (chunk)
    INTEGER                   :: ics     = -1     !< Index of chunk start
    INTEGER                   :: ice     = -1     !< Index of chunk end
    INTEGER                   :: nc      = -1     !< Length of chunk
    INTEGER                   :: nsoil_e = -1     !< Number of soil layers for energy, set by ICON
    INTEGER                   :: nsoil_w = -1     !< Number of soil layers for water, set by ICON
    INTEGER                   :: nsnow_e = -1     !< Number of snow layers for energy, set by ICON
    INTEGER                   :: npft    =  0     !< Number of PFTs
  END TYPE t_jsb_task_options

  TYPE, ABSTRACT :: t_jsb_process_task
    ! Abstract base class for a process task containing the interface for specific task implementations
    CHARACTER(len=50)              :: name                !< Name of the task
    INTEGER                        :: process_id     = 0  !< ID of the process the task belongs to
    INTEGER                        :: owner_model_id = 0  !< ID of the model the task belongs to
    INTEGER                        :: timer_integrate = 0 !< ID of integrate timer for the task
    INTEGER                        :: timer_aggregate = 0 !< ID of aggregate timer for the task
  CONTAINS
    ! Precedures that must be implemented by specific task types
    PROCEDURE(Interface_task), DEFERRED, NOPASS :: Integrate                   !< Update method for a specific tasks
    PROCEDURE(Interface_task), DEFERRED, NOPASS :: Aggregate                   !< Aggregate method for a specific tasks
    PROCEDURE                                   :: Do_it                       !< Executer of the task on a tile
    PROCEDURE                                   :: Construct => Construct_task !< Initialize a process task
  END TYPE t_jsb_process_task

  ABSTRACT INTERFACE
    SUBROUTINE Interface_task(tile, options)
      ! Each call needs to always pass the tile and options objects
      ! options: e.g. number of layers, array sizes, time step length, current time, etc.
      IMPORT :: t_jsb_tile_abstract, t_jsb_task_options !  Equivalent to a USE statement, but used in interfaces
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Tile on which the task is processed
      TYPE(t_jsb_task_options),   INTENT(in) :: options !< Container for infos typically needed for the update task
    END SUBROUTINE Interface_task
  END INTERFACE

  TYPE t_jsb_process_task_p
    ! Pointer type for process task queue
    CLASS(t_jsb_process_task), POINTER :: p
  END TYPE t_jsb_process_task_p

  TYPE t_Task_queue
    ! Type for a sequential line of task entries: only the tasks in the order of their execution, without
    ! any information on which tiles they are executed (see t_jsb_task_schedule for that)
    CLASS(t_jsb_process_task), POINTER :: task    => NULL() !< Element of the task queue
    TYPE(t_Task_queue),        POINTER :: next    => NULL() !< Pointer to next element in the queue
  CONTAINS
    PROCEDURE :: Append => Append_task                      !< Method to append a task to the queue
  END TYPE t_Task_queue

  !> One step of a task schedule: execute the `action` of the task on `tile`
  TYPE t_jsb_task_step
    CLASS(t_jsb_tile_abstract), POINTER :: tile   => NULL() !< Tile on which the task is executed
    INTEGER                             :: action = 0       !< TASK_INTEGRATE or TASK_AGGREGATE
  END TYPE t_jsb_task_step

  !> Schedule of a task: the steps in which the task is executed on the tile tree, in order
  !>
  !> While the task queue ([[t_Task_queue]]) only defines the order of the tasks, the schedule also
  !> contains, for each task, on which tiles it is to be integrated or aggregated and in which order.
  TYPE t_jsb_task_schedule
    CLASS(t_jsb_process_task), POINTER  :: task => NULL()  !< Task
    TYPE(t_jsb_task_step),  ALLOCATABLE :: steps(:)        !< Steps of the task in order of execution
  END TYPE t_jsb_task_schedule

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_task_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !>#### Initialize a process task
  !>
  !> This routine is used to initialize the components of a process task, such as its name, the
  !> ID of the process and model it belongs to, and timers for the `Integrate` and `Aggregate` methods.
  SUBROUTINE Construct_task(this, name, process_id, owner_model_id)

    CLASS(t_jsb_process_task), INTENT(inout) :: this           !< Process instance
    CHARACTER(len=*),          INTENT(in)    :: name           !< Name of the task
    INTEGER,                   INTENT(in)    :: process_id     !< ID of the process the task belongs to
    INTEGER,                   INTENT(in)    :: owner_model_id !< ID of the model the task belongs to

    ! Set attributes for the task
    this%name = name
    this%process_id = process_id
    this%owner_model_id = owner_model_id

    IF (timer_on('all')) THEN
      ! Create timers for this task
      this%timer_integrate = new_timer('jsb:update_'//TRIM(name)) ! appears in the log file of the simulation
      this%timer_aggregate = new_timer('jsb:aggregate_'//TRIM(name))
    END IF

  END SUBROUTINE Construct_task

  ! ======================================================================================================= !
  !>
  !>#### Executor for a process task
  !>
  !> This subroutine calls the appropriate function depending on the action (`TASK_INTEGRATE` or
  !> `TASK_AGGREGATE`) for a specific task. It also starts and stops timers if timing is enabled.
  !>
  SUBROUTINE Do_it(this, action, tile, options)

    CLASS(t_jsb_process_task),  INTENT(in)    :: this    !< Process instance
    INTEGER,                    INTENT(in)    :: action  !< Action to be performed
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile on which the task is processed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Container for infos typically needed for the task

    CHARACTER(len=*), PARAMETER :: routine = modname//':Do_it'

    ! Note: The task schedule only contains steps on tiles where the process of the task is active, so
    !       this is not checked here.

    ! Provide the compiler with information about the particular instance of the task
    SELECT TYPE (this)
    CLASS IS (t_jsb_process_task)

      SELECT CASE (action)
      CASE (TASK_INTEGRATE)
        IF (timer_on('detail')) THEN
          CALL timer_start(timer_integrate(this%owner_model_id))                ! Start the general integrate timer
          IF (this%timer_integrate > 0) CALL timer_start(this%timer_integrate)  ! Start the timer for this task
          ! Note: this%timer_integrate > 0 will only be true if timer_level is 3 ('all')
        END IF
        CALL this%Integrate(tile, options) ! Call the update method of the specific task
        IF (timer_on('detail')) THEN
          IF (this%timer_integrate > 0) CALL timer_stop(this%timer_integrate)
          CALL timer_stop(timer_integrate(this%owner_model_id))
        END IF
      CASE (TASK_AGGREGATE)
        IF (timer_on('detail')) THEN
          CALL timer_start(timer_aggregate(this%owner_model_id)) ! Start the general aggregate timer
          IF (this%timer_aggregate > 0) CALL timer_start(this%timer_aggregate) ! Start the timer for this task
        END IF
        CALL this%Aggregate(tile, options) ! Call the aggregate method of the specific task
        IF (timer_on('detail')) THEN
          IF (this%timer_aggregate > 0) CALL timer_stop(this%timer_aggregate)
          CALL timer_stop(timer_aggregate(this%owner_model_id))
        END IF
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown action for task'//TRIM(this%name))
      END SELECT

    END SELECT

  END SUBROUTINE Do_it

  !> Returns name of task action
  FUNCTION Get_task_action_name(action) RESULT(return_value)

    INTEGER, INTENT(in) :: action           !< TASK_INTEGRATE or TASK_AGGREGATE
    CHARACTER(len=9)    :: return_value

    SELECT CASE (action)
    CASE (TASK_INTEGRATE)
      return_value = 'INTEGRATE'
    CASE (TASK_AGGREGATE)
      return_value = 'AGGREGATE'
    CASE DEFAULT
      return_value = 'UNKNOWN'
    END SELECT

  END FUNCTION Get_task_action_name

  ! ======================================================================================================= !
  !>
  !>#### Append task to task queue
  !>
  !> This subroutine appends a specific process task to a task queue.
  !>
  SUBROUTINE Append_task(this, task, process_name)

    CLASS(t_Task_queue), TARGET, INTENT(inout)  :: this          !< The task queue object itself
    CLASS(t_jsb_process_task), TARGET           :: task          !< Process task to append
    CHARACTER(len=*),            INTENT(in)     :: process_name  !< Name of the process the task belongs to

    CLASS(t_Task_queue), POINTER :: current, next                !< Pointers to navigate and create new queue elements

    CHARACTER(len=*), PARAMETER :: routine = modname//':Append_task'

    IF (.NOT. ASSOCIATED(this%task)) THEN
      ! If queue is empty just add task as first element
      this%task => task
    ELSE
      ! else find last task in queue (go through the list until there is no next element pointer)
      current => this
      DO WHILE (ASSOCIATED(current%next))
        current => current%next
      END DO

      ! Allocate and append a new element to the end of the queue
      ALLOCATE(next)
      next%task => task      ! Next element is set to be the task to append
      current%next => next   ! Set 'next' element (the newly added entry) as next of current last element
    END IF

    IF (debug_on()) CALL message(TRIM(routine), &
      & 'Appended task "'//TRIM(process_name)//':'//TRIM(task%name)//'" to model queue')

  END SUBROUTINE Append_task

#endif
END MODULE mo_jsb_task_class
