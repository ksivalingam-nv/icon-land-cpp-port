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

PROGRAM test_add_var

  USE mo_kind, ONLY: dp, sp
  USE mo_impl_constants, ONLY: REAL_T, SINGLE_T, BOOL_T, INT_T
  USE mo_var, ONLY: t_var
  USE mo_var_list_register, ONLY: vlr_add
  USE mo_var_list, ONLY: add_var, t_var_list_ptr
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_test_common, ONLY: test_pass, test_fail
  USE memman, ONLY: mm_initialize, mm_finalize

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: testname = "test_add_var"
  CHARACTER(*), PARAMETER :: listname = "foo"
  CHARACTER(*), PARAMETER :: dname = "bar"
  TYPE(t_var), POINTER :: vl_elem => NULL()
  TYPE(t_var_metadata) :: info
  TYPE(t_var_list_ptr) :: var_list_ptr, plain_list_ptr

  ! MemMan must be initialized before add_var uses it (device allocation on GPU)
  IF (mm_initialize() /= 0) CALL test_fail(testname, "mm_initialize() failed")

  ! Register var_list_ptr
  CALL vlr_add(var_list_ptr, TRIM(listname), model_type="test")

  ! PASS: data_type=REAL_T, and passing rvals. info%ndims defaults to 0, so this is a
  ! rank-0 variable: not representable in MemMan and must take the plain-ALLOCATE
  ! fallback.
  CALL add_var(REAL_T, var_list_ptr, dname//"1", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, info%used_dimensions(1:info%ndims), vl_elem, &
      & initval_r=info%initval%rval, resetval_r=info%resetval%rval, &
      & missval_r=info%missval%rval)
  IF (vl_elem%info%lmemman) &
    & CALL test_fail(testname, "rank-0 REAL_T must not be MemMan-backed")

  ! PASS: data_type=SINGLE_T, and passing svals (rank 0, plain fallback as above)
  CALL add_var(SINGLE_T, var_list_ptr, dname//"2", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, info%used_dimensions(1:info%ndims), vl_elem, &
      & initval_s=info%initval%sval, resetval_s=info%resetval%sval, &
      & missval_s=info%missval%sval)
  IF (vl_elem%info%lmemman) &
    & CALL test_fail(testname, "rank-0 SINGLE_T must not be MemMan-backed")

  ! PASS: data_type=REAL_T with real dimensions -> allocated through MemMan.
  CALL add_var(REAL_T, var_list_ptr, dname//"3", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)
  IF (.NOT. vl_elem%info%lmemman) &
    & CALL test_fail(testname, "rank-3 REAL_T must be MemMan-backed")
  IF (.NOT. ASSOCIATED(vl_elem%r_ptr)) &
    & CALL test_fail(testname, "rank-3 REAL_T: r_ptr not associated")
  IF (SIZE(vl_elem%r_ptr) /= 24) &
    & CALL test_fail(testname, "rank-3 REAL_T: MemMan storage has the wrong size")

  ! PASS: data_type=SINGLE_T with real dimensions -> MemMan.
  CALL add_var(SINGLE_T, var_list_ptr, dname//"4", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)
  IF (.NOT. vl_elem%info%lmemman) &
    & CALL test_fail(testname, "rank-3 SINGLE_T must be MemMan-backed")

  ! PASS: data_type=INT_T with real dimensions -> MemMan.
  CALL add_var(INT_T, var_list_ptr, dname//"5", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem, &
      & initval_i=info%initval%ival, resetval_i=info%resetval%ival, &
      & missval_i=info%missval%ival)
  IF (.NOT. vl_elem%info%lmemman) &
    & CALL test_fail(testname, "rank-3 INT_T must be MemMan-backed")

  ! PASS: data_type=BOOL_T. LOGICAL arrays are not MemMan-managed (MemMan's Bool is
  ! c_bool, ICON's is LOGICAL) and must use a plain Fortran allocation.
  CALL add_var(BOOL_T, var_list_ptr, dname//"6", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem, &
      & initval_l=info%initval%lval, resetval_l=info%resetval%lval, &
      & missval_l=info%missval%lval)
  IF (vl_elem%info%lmemman) &
    & CALL test_fail(testname, "BOOL_T must not be MemMan-backed")

  ! PASS: data_type=REAL_T with a degenerate (zero-length) extent. ext_data creates
  ! placeholder fields sized by dimension counts (nclass_lu/nhori/nmonths_ext) that are
  ! 0 until the external parameter data is read. Such a field must not reach MemMan
  ! (get_var_data would leave the field pointer disassociated, and a negative extent
  ! folds into a huge unsigned size in the binding); a plain ALLOCATE yields a valid
  ! zero-size array.
  CALL add_var(REAL_T, var_list_ptr, dname//"7", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 0 /), vl_elem)
  IF (vl_elem%info%lmemman) &
    & CALL test_fail(testname, "degenerate extents must not be MemMan-backed")
  IF (.NOT. ASSOCIATED(vl_elem%r_ptr)) &
    & CALL test_fail(testname, "degenerate extents: r_ptr not associated")
  IF (SIZE(vl_elem%r_ptr) /= 0) &
    & CALL test_fail(testname, "degenerate extents: expected a zero-size array")

  ! PASS: as above but with lopenacc=.TRUE. -> in GPU builds this also exercises the
  ! device sub-path with a zero size (ext_data placeholder fields such as horizon and
  ! emi_* are created with lopenacc=.TRUE.).
  CALL add_var(REAL_T, var_list_ptr, dname//"8", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 0 /), vl_elem, lopenacc=.TRUE.)
  IF (vl_elem%info%lmemman) &
    & CALL test_fail(testname, "degenerate extents (lopenacc) must not be MemMan-backed")

  ! PASS: two time levels of the same base name. The MemMan descriptor encodes the
  ! parsed timelevel separately from the stripped name, so the two must register as
  ! distinct variables rather than as a descriptor collision.
  CALL add_var(REAL_T, var_list_ptr, dname//"9.TL1", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)
  IF (.NOT. vl_elem%info%lmemman) &
    & CALL test_fail(testname, "timelevel 1 variable must be MemMan-backed")
  CALL add_var(REAL_T, var_list_ptr, dname//"9.TL2", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)
  IF (.NOT. vl_elem%info%lmemman) &
    & CALL test_fail(testname, "timelevel 2 variable must be MemMan-backed")

  ! PASS: a var list registered with lmemman=.FALSE. opts out of MemMan for all its
  ! variables (used by the derived post-processing lists, whose variables reuse their
  ! source variables' names). The same variable name as in the MemMan-backed list above
  ! must be accepted here without a descriptor collision precisely because the copy
  ! never registers with MemMan.
  CALL vlr_add(plain_list_ptr, TRIM(listname)//"_plain", model_type="test", &
      & lmemman=.FALSE.)
  CALL add_var(REAL_T, plain_list_ptr, dname//"3", info%hgrid, info%vgrid, info%cf, &
      & info%grib2, (/ 4, 3, 2 /), vl_elem)
  IF (vl_elem%info%lmemman) &
    & CALL test_fail(testname, "variable in an opted-out list must not be MemMan-backed")

  IF (mm_finalize() /= 0) CALL test_fail(testname, "mm_finalize() failed")

  CALL test_pass()

END PROGRAM test_add_var
