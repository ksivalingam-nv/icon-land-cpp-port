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

PROGRAM test_add_var_fail3

  USE mo_impl_constants, ONLY: REAL_T
  USE mo_var, ONLY: t_var
  USE mo_var_list_register, ONLY: vlr_add
  USE mo_var_list, ONLY: add_var, t_var_list_ptr
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_test_common, ONLY: test_pass, test_fail
  USE memman, ONLY: mm_initialize

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: testname = "test_add_var_fail3"
  CHARACTER(*), PARAMETER :: dname = "clash"
  TYPE(t_var), POINTER :: vl_elem => NULL()
  TYPE(t_var_metadata) :: info
  TYPE(t_var_list_ptr) :: list_a, list_b

  IF (mm_initialize() /= 0) CALL test_fail(testname, "mm_initialize() failed")

  CALL vlr_add(list_a, "foo_a", model_type="test")
  CALL vlr_add(list_b, "foo_b", model_type="test")

  ! Both lists are MemMan-backed, so the same variable name (same domain, same
  ! timelevel) maps to the same MemMan descriptor key. The second add_var must be
  ! rejected by the descriptor-collision detector: MemMan would otherwise register a
  ! second entry under the same key, and every key-based lookup silently resolves to
  ! the first one, aliasing the two variables.
  CALL add_var(REAL_T, list_a, dname, info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)

  ! FAIL 3: expected to abort here
  CALL add_var(REAL_T, list_b, dname, info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)

  CALL test_pass()

END PROGRAM test_add_var_fail3
