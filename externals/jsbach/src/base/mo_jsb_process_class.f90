!> Contains process class
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
!>#### Provides functionality for process objects
!>
!> Scientific processes are usually capsuled within an own sub-folder within the ICON-LAND code structure,
!> e.g. hydrology, phenology, srf_energy_bal, and more. A process is realized by a couple of obligatory
!> and optional process modules (depending on the type of process and process requirements).
!> A usecase specifies which process is running on which tiles and thereby determines the process actions,
!> which further can be modified (processes switched on and off) via the namelist.
!> When implementing a new process in ICON-LAND, mandatory process modules need to be added, and the process
!> needs to be made known in several modules of the ICON-LAND infrastructure. This functionality is provided
!> by the jsb_process_class module.
!>
MODULE mo_jsb_process_class
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish
  USE mo_util,                ONLY: tolower
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_process_task_p
  USE mo_jsb_config_class,    ONLY: t_jsb_config
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process, t_jsb_process_p, Get_process_name, Get_process_id
  PUBLIC :: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, PHENO_
  PUBLIC :: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_, Get_action_name

  TYPE t_jsb_process
    INTEGER                       :: id = 0                         !< Process id
    CHARACTER(len=SHORT_NAME_LEN) :: name = ''                      !< Process name
    CLASS(t_jsb_config), POINTER  :: config => NULL()               !< Pointer to process configuration
    TYPE(t_jsb_process_task_p), ALLOCATABLE :: task_list(:)         !< List of tasks to be done for a given process
    INTEGER                       :: owner_model_id                 !< Model instance which uses a process
    LOGICAL                       :: l_changes_fractions = .FALSE.  !< Does this process change tile fractions?
    LOGICAL                       :: has_memory = .TRUE.            !< Does this process have an own memory module
  CONTAINS
    PROCEDURE :: Configure => Configure_process                     !< Call process initialization
    PROCEDURE :: Register_task                                      !< Add a task to the process task list
    PROCEDURE :: Get_task                                           !< Get the pointer to a process task
  END TYPE t_jsb_process

  ! The constructor for t_jsb_process is Create_process in mo_jsb_process_factory.

  ! -------------------------------------------------------------------------------------------------- !
  ! Create a vector of pointers for processes
  TYPE t_jsb_process_p
    CLASS(t_jsb_process), POINTER :: p
  END TYPE t_jsb_process_p

  ! -------------------------------------------------------------------------------------------------- !
  ! Assign identifiers to processes
  ENUM, BIND(C)
    ENUMERATOR :: A2L_=1, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, PHENO_
  END ENUM

  ! -------------------------------------------------------------------------------------------------- !
  ! Assign identifiers for process actions, e.g where to run or what to do for a process
  ! It is important that AGGREGATE comes before the other actions
  ! INHERIT and AGGREGATE are usually set automatically.
  !   SKIP:       default, do not do anything for this process but allow for own memory
  !   INHERIT:    give access to the memory of the parent tile (automatically
  !               for childs of a ON_TILE process)
  !   AGGREGATE:  aggregate the child tiles variables to the parent tiles
  !   ON_LEAFS:   runs process on all sub-tiles of the current tile that don't
  !               have any children
  !   ON_TILE:    runs process on the current tile
  !   ON-SUBTREE: runs process on current tile but can modify memory on children tiles
  ENUM, BIND(C)
    ENUMERATOR :: SKIP_=0, INHERIT_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_class'

CONTAINS

  ! ====================================================================================================== !
  !>
  !>#### Assign a process short name to a given process
  !>
  FUNCTION Get_process_name(iproc) RESULT(return_value)

    INTEGER, INTENT(in) :: iproc                   !< Process identifier
    CHARACTER(len=:), ALLOCATABLE :: return_value  !< Process name

    SELECT CASE(iproc)
    CASE (A2L_)
      return_value = 'a2l'
    CASE (L2A_)
      return_value = 'l2a'
    CASE (SEB_)
      return_value = 'seb'
    CASE (TURB_)
      return_value = 'turb'
    CASE (SSE_)
      return_value = 'sse'
    CASE (HYDRO_)
      return_value = 'hydro'
    CASE (RAD_)
      return_value = 'rad'
    CASE (PHENO_)
      return_value = 'pheno'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION Get_process_name

  ! ====================================================================================================== !
  !>
  !>#### Return process identifier from process name
  !>
  FUNCTION Get_process_id(name) RESULT(return_value)

    CHARACTER(len=*), INTENT(in) :: name           !< process name
    INTEGER                      :: return_value   !< process identifier

    SELECT CASE(tolower(TRIM(name)))
    CASE ('a2l')
      return_value = A2L_
    CASE ('l2a')
      return_value = L2A_
    CASE ('seb')
      return_value = SEB_
    CASE ('turb')
      return_value = TURB_
    CASE ('sse')
      return_value = SSE_
    CASE ('hydro')
      return_value = HYDRO_
    CASE ('rad')
      return_value = RAD_
    CASE ('pheno')
      return_value = PHENO_
    CASE DEFAULT
      return_value = 0
    END SELECT

  END FUNCTION Get_process_id

  ! ====================================================================================================== !
  !>
  !>#### Add a process task to the task list for a given process
  !>
  SUBROUTINE Register_task(this, task)

    CLASS(t_jsb_process),     INTENT(inout)       :: this !< Process task list
    CLASS(t_jsb_process_task), TARGET, INTENT(in) :: task !< Process task

    TYPE(t_jsb_process_task_p), ALLOCATABLE :: task_temp(:) !< Temporary task list
    INTEGER :: n !< Number of already existing tasks in the task list

    CHARACTER(len=*), PARAMETER :: routine = modname//':Register_task'

    ! Check whether the task belongs to the actual process
    IF (this%id /= task%process_id) &
      CALL finish(TRIM(routine), 'Task is for another process')

    ! Add task as first element of the task list
    IF (.NOT. ALLOCATED(this%task_list)) THEN
      ALLOCATE(this%task_list(1))
      this%task_list(1) = t_jsb_process_task_p(task)
    ! Enlarge task list by one and add the task to it
    ELSE
      n = SIZE(this%task_list)
      ALLOCATE(task_temp(n+1))
      task_temp(1:n) = this%task_list
      task_temp(n+1) = t_jsb_process_task_p(task)
      CALL MOVE_ALLOC(task_temp, this%task_list)
    END IF

    IF (debug_on()) CALL message(TRIM(routine), 'Registered task "'//TRIM(task%name)// &
      &                                         '" for process "'//TRIM(this%name)//'".')

  END SUBROUTINE Register_task

  ! ====================================================================================================== !
  !>
  !>#### Return a pointer to a task in the task list
  !>
  FUNCTION Get_task(this, name) RESULT(return_ptr)

    CLASS(t_jsb_process), INTENT(inout)          :: this          !< Process object
    CHARACTER(LEN=*),              INTENT(in)    :: name          !< Process task name
    CLASS(t_jsb_process_task), POINTER           :: return_ptr    !< Pointer to the task

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_task'

    return_ptr => NULL()
    DO i=1,SIZE(this%task_list)
      !IF (.NOT. ASSOCIATED())
      IF (TRIM(this%task_list(i)%p%name) == TRIM(name) ) THEN
        return_ptr => this%task_list(i)%p
        EXIT
      END IF
    END DO

    IF (.NOT.(ASSOCIATED(return_ptr))) &
      & CALL finish(TRIM(routine), 'Task "'//TRIM(name)//'" not found for process "'// &
      &                            TRIM(this%name)//'".')

  END FUNCTION Get_task

  ! ====================================================================================================== !
  !>
  !>#### Call the process init function
  !>
  SUBROUTINE Configure_process(this)

    CLASS(t_jsb_process), INTENT(inout) :: this      !< Process object

    CHARACTER(len=*), PARAMETER :: routine = modname//':Configure_process'

    ! Test whether a config exists for this process
    IF (.NOT. ASSOCIATED(this%config)) RETURN

    CALL this%config%Init()

  END SUBROUTINE Configure_process

  ! ====================================================================================================== !
  !>
  !>#### Return name for a process action identifier
  !>
  FUNCTION Get_action_name(iact) RESULT(return_value)

    INTEGER, INTENT(in) :: iact                    !< Process action identifier
    CHARACTER(len=:), ALLOCATABLE :: return_value  !< Process action name

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_action_name'

    SELECT CASE(iact)
    CASE (SKIP_)
      return_value = 'SKIP'
    CASE (AGGREGATE_)
      return_value = 'AGGREGATE'
    CASE (ON_LEAFS_)
      return_value = 'ON_LEAFS'
    CASE(ON_TILE_)
      return_value = 'ON_TILE'
    CASE (ON_SUBTREE_)
      return_value = 'ON_SUBTREE'
    CASE (INHERIT_)
      return_value = 'INHERIT'
    END SELECT

  END FUNCTION Get_action_name
  !
  ! ====================================================================================================== !

#endif
END MODULE mo_jsb_process_class
