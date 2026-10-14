!> Contains basic definitions for variable class
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
MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: finish, message
  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p, REAL1D, REAL2D, REAL3D

  ENUM, BIND(C)
    ENUMERATOR :: REAL1D=1, REAL2D, REAL3D
  END ENUM

  !================================================================================================================================
  ! t_jsb_var
  !================================================================================================================================
  TYPE, ABSTRACT :: t_jsb_var
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr1d(:) => NULL()
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr2d(:,:) => NULL()
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr3d(:,:,:) => NULL()
    CHARACTER(len=VARNAME_LEN)    :: name                      ! Name of variable
    CHARACTER(len=VARNAME_LEN)    :: standard_name             ! Standard name of the variable
    CHARACTER(len=VARNAME_LEN)    :: full_name                 ! Full name of variable in varlist, including prefix and suffix
    CHARACTER(len=:), ALLOCATABLE :: unit
    INTEGER                       :: type                      ! REAL1D / REAL2D / REAL3D
    CHARACTER(len=:), ALLOCATABLE :: element_name              !
    CHARACTER(len=:), ALLOCATABLE :: element_shortname         ! Short name of element
    INTEGER                       :: element_id                !
    INTEGER                       :: owner_model_id            !
    INTEGER                       :: owner_proc_id = -1        !
    INTEGER, ALLOCATABLE          :: owner_tile_path(:)        !
    INTEGER, ALLOCATABLE          :: child_idx(:)
      !< Position-numbers of the same var (as this var) on the children of the current tile (i.e. within mem%vars(:) of that child).
    REAL(dp)                      :: missval
    LOGICAL                       :: l_aggregate_all = .FALSE. ! When aggregating, consider fractions of all child tiles even
                                                               ! of those where variable doesn't exist
    LOGICAL                       :: is_in_output  = .FALSE.   ! Whether variable is requested to be written to output
    LOGICAL                       :: is_in_restart = .FALSE.   ! Whether variable is requested to be written to restart file
  END TYPE t_jsb_var

!#pragma GCC diagnostic push
!#pragma GCC diagnostic ignored "-Wsurprising"

  !================================================================================================================================
  ! t_jsb_var_real1d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real1d
    REAL(wp), POINTER &
      :: ptr(:) => NULL()
  END TYPE t_jsb_var_real1d

  !================================================================================================================================
  ! t_jsb_var_real2d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
      :: ptr(:,:) => NULL()    ! associated with ptr2d(:,:)
  END TYPE t_jsb_var_real2d

  !================================================================================================================================
  ! t_jsb_var_real3d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()     ! associated with ptr3d(:,:,:)
  END TYPE t_jsb_var_real3d

!#pragma GCC diagnostic pop

  !================================================================================================================================
  ! t_jsb_var_p
  !================================================================================================================================
  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p => NULL()
  END TYPE t_jsb_var_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'

#endif
END MODULE mo_jsb_var_class
