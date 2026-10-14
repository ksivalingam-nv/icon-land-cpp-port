! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
!/// @file
!/// Fortran counterpart of icon_f2c_basic_types.hpp.
!/// Defines the working precision kind, C-compatible math types (TangVec),
!/// and the f2c type-kind integer constants used for array descriptor tagging.

MODULE mo_ragnarok_f2c_basic_types
  USE ISO_C_BINDING, ONLY: c_int, c_float, c_double
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wp
  PUBLIC :: t_f2c_kind
  PUBLIC :: f2c_kind_values
  PUBLIC :: f2c_index_shift

  ! Working precision - mirrors the Real alias in common/real.hpp.
#ifdef __SINGLE_PRECISION
  INTEGER, PARAMETER :: wp = c_float
#else
  INTEGER, PARAMETER :: wp = c_double
#endif

  ! The ICON math types are not redefined here;
  ! we rely on their C representations for interoperability.
  ! We never use the Fortran type definitions in ragnarok.
  !
  ! cartesian coordinates
  !TYPE, BIND(C) :: t_cartesian_coordinates
  !  REAL(wp) :: x(3)
  !END TYPE t_cartesian_coordinates
  !
  ! geographical coordinates
  !TYPE, BIND(C) :: t_geographical_coordinates
  !  REAL(wp) :: lon
  !  REAL(wp) :: lat
  !END TYPE t_geographical_coordinates
  !
  ! the two coordinates on the tangent plane
  !TYPE, BIND(C) :: t_tangent_vectors
  !  REAL(wp) :: v1
  !  REAL(wp) :: v2
  !END TYPE t_tangent_vectors
  !

  ! Integer type-kind identifiers - mirrors the f2c_typekind<T> specialisations
  ! in icon_f2c_basic_types.hpp.
  TYPE, BIND(c) :: t_f2c_kind
    INTEGER(c_int) :: undef = 0
    INTEGER(c_int) :: c_bool = 1
    INTEGER(c_int) :: c_int = 2
    INTEGER(c_int) :: c_float = 3
    INTEGER(c_int) :: c_double = 4
    INTEGER(c_int) :: cart_c_float = 5
    INTEGER(c_int) :: geoc_c_float = 6
    INTEGER(c_int) :: tanv_c_float = 7
    INTEGER(c_int) :: cart_c_double = 8
    INTEGER(c_int) :: geoc_c_double = 9
    INTEGER(c_int) :: tanv_c_double = 10
  END TYPE t_f2c_kind
  TYPE(t_f2c_kind), PARAMETER :: f2c_kind_values = t_f2c_kind()

CONTAINS

  ! Fortran 1-based index -> ragnarok 0-based index
  PURE ELEMENTAL FUNCTION f2c_index_shift(idx) RESULT(c_idx)
    INTEGER, INTENT(IN) :: idx
    INTEGER(c_int) :: c_idx
    c_idx = INT(idx, c_int) - 1_c_int
  END FUNCTION f2c_index_shift

END MODULE mo_ragnarok_f2c_basic_types
