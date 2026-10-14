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
!/// Fortran-C interface types and helpers for ICON-Ragnarok calls.
!/// Defines C-compatible data structures and routines for exchanging grid, patch, and communication info between Fortran and C++.

MODULE mo_ragnarok_f2c
  USE mo_ragnarok_f2c_basic_types
  USE ISO_C_BINDING, ONLY: c_int, c_float, c_double, c_ptr, c_null_ptr,&
      & c_loc, c_funptr, c_null_funptr, c_bool, c_char
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_f2c_ftable
  PUBLIC :: t_dom_info
  PUBLIC :: t_f2c_patch_descr
  PUBLIC :: t_patch_info
  PUBLIC :: init_ragnarok_f2c
  PUBLIC :: t_comm_patch
  PUBLIC :: t_comm_pattern_cdescr
  PUBLIC :: ragnarok_sync_patch_array_r3_sync_c
  PUBLIC :: t_process_info
  PUBLIC :: t_f2c_grid_cells
  PUBLIC :: t_f2c_grid_edges
  PUBLIC :: t_f2c_grid_vertices
  PUBLIC :: gen_repr5, make_repr5
  PUBLIC :: f2c_kind_values
  PUBLIC :: t_f2c_repr5
  PUBLIC :: wp, t_f2c_kind, f2c_index_shift

  TYPE, BIND(c) :: t_f2c_ftable
    TYPE(c_funptr) :: get_domain_info = c_null_funptr
    TYPE(c_funptr) :: get_patch_info = c_null_funptr
    TYPE(c_funptr) :: get_mo_model_domain_p_patch_descr = c_null_funptr
    TYPE(c_funptr) :: get_comm_patch = c_null_funptr
    TYPE(c_funptr) :: get_grid_cells = c_null_funptr
    TYPE(c_funptr) :: get_grid_edges = c_null_funptr
    TYPE(c_funptr) :: get_grid_vertices = c_null_funptr
    TYPE(c_funptr) :: f2c_message = c_null_funptr
    TYPE(c_funptr) :: f2c_finish = c_null_funptr
    TYPE(c_funptr) :: f2c_exchange_data_r3d = c_null_funptr
    TYPE(c_funptr) :: get_process_info = c_null_funptr
    TYPE(c_funptr) :: f2c_new_timer = c_null_funptr
    TYPE(c_funptr) :: f2c_timer_start = c_null_funptr
    TYPE(c_funptr) :: f2c_timer_stop = c_null_funptr
    TYPE(c_funptr) :: f2c_timer_value = c_null_funptr
    TYPE(c_funptr) :: f2c_get_timers_level = c_null_funptr
  END TYPE t_f2c_ftable

  ! domain info:
  TYPE, BIND(c) :: t_dom_info
    INTEGER(c_int) :: nproma = 0
    INTEGER(c_int) :: id_min = -1
    INTEGER(c_int) :: id_max = -2
  END TYPE t_dom_info

  TYPE, BIND(c) :: t_f2c_patch_descr
    TYPE(c_ptr) :: cptr
  END TYPE t_f2c_patch_descr

  ! portable part of top level t_patch components
  TYPE, BIND(c) :: t_patch_info
    INTEGER(c_int) :: id = -1
    INTEGER(c_int) :: nlev = 0
    INTEGER(c_int) :: nblks_c = 0
    INTEGER(c_int) :: npromz_c = 0
    INTEGER(c_int) :: nblks_e = 0
    INTEGER(c_int) :: npromz_e = 0
    INTEGER(c_int) :: nblks_v = 0
    INTEGER(c_int) :: npromz_v = 0
  END TYPE t_patch_info

  TYPE, BIND(c) :: t_comm_pattern_cdescr
    TYPE(c_ptr) :: cptr
  END TYPE t_comm_pattern_cdescr

  TYPE, BIND(c) :: t_comm_patch
    TYPE(t_comm_pattern_cdescr) :: comm_pat_c
    TYPE(t_comm_pattern_cdescr) :: comm_pat_e
    TYPE(t_comm_pattern_cdescr) :: comm_pat_v
    TYPE(t_comm_pattern_cdescr) :: comm_pat_c1
  END TYPE t_comm_patch

  TYPE, BIND(c) :: t_f2c_repr5
    TYPE(c_ptr) :: cptr = c_null_ptr
    INTEGER(c_int) :: f2c_kind = f2c_kind_values%undef
    INTEGER(c_int) :: rank = 0
    INTEGER(c_int) :: shape_vec5(5) = 0
    INTEGER(c_int) :: lbound_vec5(5) = 1
  END TYPE t_f2c_repr5


  TYPE, BIND(c) :: t_f2c_grid_cells
    TYPE(t_f2c_repr5) :: num_edges
    TYPE(t_f2c_repr5) :: area
    TYPE(t_f2c_repr5) :: center
    TYPE(t_f2c_repr5) :: edge_idx_f
    TYPE(t_f2c_repr5) :: edge_blk_f
    TYPE(t_f2c_repr5) :: start_block_f
    TYPE(t_f2c_repr5) :: end_block_f
    TYPE(t_f2c_repr5) :: start_index_f
    TYPE(t_f2c_repr5) :: end_index_f
  END TYPE t_f2c_grid_cells

  TYPE, BIND(c) :: t_f2c_grid_edges
    TYPE(t_f2c_repr5) :: primal_normal_cell
    TYPE(t_f2c_repr5) :: primal_normal_vert
    TYPE(t_f2c_repr5) :: dual_normal_cell
    TYPE(t_f2c_repr5) :: dual_normal_vert
    TYPE(t_f2c_repr5) :: tangent_orientation
    TYPE(t_f2c_repr5) :: inv_primal_edge_length
    TYPE(t_f2c_repr5) :: inv_dual_edge_length
    TYPE(t_f2c_repr5) :: inv_vert_vert_length
    TYPE(t_f2c_repr5) :: edge_cell_length
    TYPE(t_f2c_repr5) :: cell_idx_f
    TYPE(t_f2c_repr5) :: cell_blk_f
    TYPE(t_f2c_repr5) :: vertex_idx_f
    TYPE(t_f2c_repr5) :: vertex_blk_f
    TYPE(t_f2c_repr5) :: start_block_f
    TYPE(t_f2c_repr5) :: end_block_f
    TYPE(t_f2c_repr5) :: start_index_f
    TYPE(t_f2c_repr5) :: end_index_f
  END TYPE t_f2c_grid_edges

  TYPE, BIND(c) :: t_f2c_grid_vertices
    TYPE(t_f2c_repr5) :: vertex
    TYPE(t_f2c_repr5) :: cell_idx_f
    TYPE(t_f2c_repr5) :: cell_blk_f
    TYPE(t_f2c_repr5) :: start_block_f
    TYPE(t_f2c_repr5) :: end_block_f
    TYPE(t_f2c_repr5) :: start_index_f
    TYPE(t_f2c_repr5) :: end_index_f
  END TYPE t_f2c_grid_vertices

  TYPE, BIND(c) :: t_process_info
    LOGICAL(c_bool) :: is_mpi_parallel
  END TYPE t_process_info

  INTERFACE
    SUBROUTINE init_ragnarok_f2c(funtab) BIND(c)
      IMPORT t_f2c_ftable
      TYPE(t_f2c_ftable), INTENT(IN) :: funtab
    END SUBROUTINE init_ragnarok_f2c


    SUBROUTINE ragnarok_sync_patch_array_r3_sync_c(pid, arr, arr_shape) BIND(c)
      IMPORT :: c_double, c_int
      INTEGER(c_int), VALUE :: pid
      REAL(c_double) :: arr(*)
      INTEGER(c_int), INTENT(IN) :: arr_shape(3)
    END SUBROUTINE ragnarok_sync_patch_array_r3_sync_c
  END INTERFACE

  INTERFACE gen_repr5
    MODULE PROCEDURE gen_repr5_from_int_1d
    MODULE PROCEDURE gen_repr5_from_int_2d
    MODULE PROCEDURE gen_repr5_from_int_3d
    MODULE PROCEDURE gen_repr5_from_sp_3d
    MODULE PROCEDURE gen_repr5_from_dp_2d
    MODULE PROCEDURE gen_repr5_from_dp_3d
    MODULE PROCEDURE gen_repr5_from_dp_4d
    MODULE PROCEDURE gen_repr5_from_dp_3d_pointer
  END INTERFACE gen_repr5

CONTAINS

  ! helper function to fill a shape vector of size 5,
  ! e.g., a given shape of [2,3] gets normalized to [2,3,1,1,1]
  ! works also for offsets in the same way
  PURE FUNCTION gen_normalized_vec5(vec) RESULT(res)
    INTEGER, INTENT(IN) :: vec(:)
    INTEGER(c_int) :: res(5)
    INTEGER :: n
    n = SIZE(vec)
    res(1:n) = INT(vec, c_int)
    res(n + 1:) = 1
  END FUNCTION gen_normalized_vec5


  ! helper function to generate a t_repr5 value
  PURE FUNCTION make_repr5(ptr, f2c_kind, sh, lb) RESULT(r)
    USE ISO_C_BINDING, ONLY: c_ptr, c_int
    TYPE(c_ptr), INTENT(IN) :: ptr
    INTEGER(c_int), INTENT(IN) :: f2c_kind
    INTEGER, INTENT(IN) :: sh(:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    r%cptr = ptr
    r%f2c_kind = f2c_kind
    r%shape_vec5 = gen_normalized_vec5(sh)
    r%lbound_vec5 = gen_normalized_vec5(lb)
    r%rank = SIZE(sh)
  END FUNCTION make_repr5


  PURE FUNCTION gen_repr5_from_int_1d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_int
    INTEGER(c_int), TARGET, INTENT(IN) :: arr(:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0], [1])
    END IF
  END FUNCTION gen_repr5_from_int_1d


  PURE FUNCTION gen_repr5_from_int_2d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_int
    INTEGER(c_int), TARGET, INTENT(IN) :: arr(:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0], [1,1])
    END IF
  END FUNCTION gen_repr5_from_int_2d


  PURE FUNCTION gen_repr5_from_int_3d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_int
    INTEGER(c_int), TARGET, INTENT(IN) :: arr(:,:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0,0], [1,1,1])
    END IF
  END FUNCTION gen_repr5_from_int_3d


  PURE FUNCTION gen_repr5_from_sp_3d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_float
    REAL(c_float), TARGET, INTENT(IN) :: arr(:,:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0,0], [1,1,1])
    END IF
  END FUNCTION gen_repr5_from_sp_3d


  PURE FUNCTION gen_repr5_from_dp_3d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_double
    REAL(c_double), TARGET, INTENT(IN) :: arr(:,:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0,0], [1,1,1])
    END IF
  END FUNCTION gen_repr5_from_dp_3d


  ! Variant for POINTER components whose association status is only conditionally
  ! established (e.g. add_var'd under some physics config, otherwise left disassociated).
  ! Precondition (caller's responsibility): arr's association status must already be
  ! DEFINED (associated or disassociated) when this is called - never undefined. This
  ! avoids ever calling SIZE/SHAPE/LBOUND/C_LOC on a disassociated pointer, which is
  ! non-conforming Fortran regardless of how "well-defined" the disassociation is.
  PURE FUNCTION gen_repr5_from_dp_3d_pointer(arr) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_double
    REAL(c_double), POINTER, INTENT(IN) :: arr(:,:,:)
    TYPE(t_f2c_repr5) :: r
    IF (ASSOCIATED(arr)) THEN
      IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
        r = make_repr5(C_LOC(arr), tk, SHAPE(arr), LBOUND(arr))
      ELSE
        r = make_repr5(c_null_ptr, tk, [0,0,0], [1,1,1])
      END IF
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0,0], [1,1,1])
    END IF
  END FUNCTION gen_repr5_from_dp_3d_pointer


  PURE FUNCTION gen_repr5_from_dp_2d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_double
    REAL(c_double), TARGET, INTENT(IN) :: arr(:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0], [1,1])
    END IF
  END FUNCTION gen_repr5_from_dp_2d


  PURE FUNCTION gen_repr5_from_dp_4d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%c_double
    REAL(c_double), TARGET, INTENT(IN) :: arr(:, :, :, :)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0,0,0], [1,1,1,1])
    END IF
  END FUNCTION gen_repr5_from_dp_4d

END MODULE mo_ragnarok_f2c
