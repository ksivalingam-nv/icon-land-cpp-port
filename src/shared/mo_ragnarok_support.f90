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

! provides access to icon functionlity for ragnarok
#ifndef __NO_RAGNAROK__

MODULE mo_ragnarok_support
  USE ISO_C_BINDING, ONLY: c_int, c_ptr, c_null_ptr, c_loc, c_funptr, c_FUNLOC, &
       & c_char, c_f_pointer, c_bool, c_double
  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message, message_text, finish, warning
  USE mo_parallel_config, ONLY: nproma
  USE mo_model_domain, ONLY: t_patch, p_patch
  USE mo_io_units, ONLY: filename_max
  USE mo_ragnarok_f2c, ONLY: t_f2c_ftable, t_dom_info, t_f2c_patch_descr, &
       & t_patch_info, init_ragnarok_f2c, &
       & t_comm_patch, t_comm_pattern_cdescr, &
       & t_process_info, t_f2c_grid_cells, t_f2c_grid_edges, t_f2c_grid_vertices, &
       & gen_repr5, make_repr5, f2c_kind_values, t_f2c_repr5
  USE mo_ragnarok_timer, only: init_ragnarok_timer
  USE mo_communication, ONLY: t_comm_pattern,  exchange_data
  USE mo_communication_types, ONLY: t_comm_pattern_descr
  USE mo_timer, ONLY: new_timer, timer_start, timer_stop, ltimer, timers_level
  USE mo_real_timer, ONLY: timer_val
  USE mo_math_types, ONLY: t_tangent_vectors, t_geographical_coordinates
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_ragnarok_support

  CHARACTER(len=*), PARAMETER :: module_name = 'mo_ragnarok_support'
  LOGICAL, PRIVATE :: is_initialized = .FALSE.

  ! extending imported interface:
  INTERFACE gen_repr5
    MODULE PROCEDURE gen_repr5_from_tanv_dp_3d
    MODULE PROCEDURE gen_repr5_from_geoc_dp_2d
  END INTERFACE gen_repr5

  ! X, LBOUND(X) as a single comma-expanding actual-argument pair for gen_repr5(arr, lb).
  ! Needed because gen_repr5's array dummy is a plain (non-pointer) assumed-shape argument,
  ! which can't recover the actual's real bounds on its own - the caller must capture
  ! LBOUND(x) itself, right next to the array expression, to avoid the two getting out of sync.
#define X_COMMA_LBOUND_X(x) x, LBOUND(x)

CONTAINS

  ! Call once p_patch exists and its communication patterns are set up.
  SUBROUTINE init_ragnarok_support()
    CHARACTER(len=*), PARAMETER :: context = module_name //'::init_ragnarok_support'
    IF (is_initialized) RETURN
    CALL init_ragnarok_f2c( t_f2c_ftable( &
         & get_domain_info = c_funloc(get_domain_info), &
         & get_patch_info = c_funloc(get_patch_info), &
         & get_mo_model_domain_p_patch_descr = c_funloc(get_mo_model_domain_p_patch_descr), &
         & get_comm_patch = c_funloc(get_comm_patch), &
         & get_grid_cells = c_funloc(get_f2c_grid_cells), &
         & get_grid_edges = c_funloc(get_f2c_grid_edges), &
         & get_grid_vertices = c_funloc(get_f2c_grid_vertices), &
         & f2c_message = c_funloc(f2c_message), &
         & f2c_finish = c_funloc(f2c_finish), &
         & f2c_exchange_data_r3d = c_funloc(f2c_exchange_data_r3d), &
         & get_process_info = c_FUNLOC(get_process_info), &
         & f2c_new_timer = c_FUNLOC(f2c_new_timer), &
         & f2c_timer_start = c_FUNLOC(f2c_timer_start), &
         & f2c_timer_stop = c_FUNLOC(f2c_timer_stop), &
         & f2c_timer_value = c_FUNLOC(f2c_timer_value), &
         & f2c_get_timers_level = c_FUNLOC(f2c_get_timers_level) &
         & ))
    CALL init_ragnarok_timer()
    CALL message(context,'ragnarok_support is initialized.')
    is_initialized = .TRUE.
  END SUBROUTINE init_ragnarok_support

  FUNCTION get_process_info() RESULT(process_info) BIND(c)
    USE mo_mpi, ONLY: my_process_is_mpi_parallel
    TYPE(t_process_info) :: process_info
    process_info = t_process_info(is_mpi_parallel = my_process_is_mpi_parallel())
  END FUNCTION get_process_info

  SUBROUTINE f2c_exchange_data_r3d(pat_descr, lacc, recv, recv_shape) BIND(c)
    ! note: we do not implement the full interface (yet)
    TYPE(t_comm_pattern_cdescr), VALUE :: pat_descr
    LOGICAL(c_bool), VALUE :: lacc ! If true, use openacc
    INTEGER(c_int), INTENT(in) :: recv_shape(3)
    REAL(c_double), INTENT(INOUT) :: recv(recv_shape(1),recv_shape(2),recv_shape(3))

    CALL exchange_data(get_comm_pattern(pat_descr), LOGICAL(lacc), recv)

  END SUBROUTINE f2c_exchange_data_r3d

  SUBROUTINE get_domain_info(dom_info) BIND(c)
    CHARACTER(len=*), PARAMETER :: context = module_name //'::get_domain_info'
    TYPE(t_dom_info), INTENT(out) :: dom_info
    dom_info%nproma = nproma
    IF (ALLOCATED(p_patch)) THEN
      dom_info%id_min = LBOUND(p_patch,1)
      dom_info%id_max = UBOUND(p_patch,1)
    ELSE
      CALL finish(context,'p_patch not allocated')
    ENDIF
  END SUBROUTINE get_domain_info

  SUBROUTINE get_patch_info(icon_patch, patch_info) BIND(c)
    TYPE(t_f2c_patch_descr), VALUE :: icon_patch
    TYPE(t_patch_info), INTENT(out) :: patch_info
    TYPE(t_patch), POINTER :: p
    CALL c_f_POINTER(icon_patch%cptr, p)
    patch_info = t_patch_info( &
         & id = p%id, &
         & nlev = p%nlev, &
         & nblks_c = p%nblks_c, &
         & npromz_c = p%npromz_c, &
         & nblks_e = p%nblks_e, &
         & npromz_e = p%npromz_e, &
         & nblks_v = p%nblks_v, &
         & npromz_v = p%npromz_v &
         & )
  END SUBROUTINE get_patch_info

  SUBROUTINE get_comm_patch(icon_patch, comm_patch) BIND(c)
    TYPE(t_f2c_patch_descr), VALUE, INTENT(in) :: icon_patch
    TYPE(t_comm_patch), INTENT(out) :: comm_patch
    TYPE(t_patch), POINTER :: patch
    CALL c_f_POINTER(icon_patch%cptr, patch)
    comm_patch = t_comm_patch( &
         & comm_pat_c = t_comm_pattern_cdescr(c_LOC(patch%comm_pat_c%descr)), &
         & comm_pat_e = t_comm_pattern_cdescr(c_LOC(patch%comm_pat_e%descr)), &
         & comm_pat_v = t_comm_pattern_cdescr(c_LOC(patch%comm_pat_v%descr)), &
         & comm_pat_c1 = t_comm_pattern_cdescr(c_LOC(patch%comm_pat_c1%descr)) )
  END SUBROUTINE get_comm_patch

  PURE FUNCTION gen_repr5_from_tanv_dp_3d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%tanv_c_double
    TYPE(t_tangent_vectors), TARGET, INTENT(in) :: arr(:,:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0,0], [1,1,1])
    END IF
  END FUNCTION gen_repr5_from_tanv_dp_3d

  PURE FUNCTION gen_repr5_from_geoc_dp_2d(arr, lb) RESULT(r)
    INTEGER(c_int), PARAMETER :: tk = f2c_kind_values%geoc_c_double
    TYPE(t_geographical_coordinates), TARGET, INTENT(in) :: arr(:,:)
    INTEGER, INTENT(IN) :: lb(:)
    TYPE(t_f2c_repr5) :: r
    IF (SIZE(arr) > 0 .AND. IS_CONTIGUOUS(arr)) THEN
      r = make_repr5(C_LOC(arr), tk, SHAPE(arr), lb)
    ELSE
      r = make_repr5(c_null_ptr, tk, [0,0], [1,1])
    END IF
  END FUNCTION gen_repr5_from_geoc_dp_2d

  SUBROUTINE get_f2c_grid_cells(patch_descr, f2c_grid_cells) BIND(c)
    TYPE(t_f2c_patch_descr), VALUE, INTENT(in) :: patch_descr
    TYPE(t_f2c_grid_cells), INTENT(out) :: f2c_grid_cells

    TYPE(t_patch), POINTER :: patch
    CALL c_f_POINTER(patch_descr%cptr, patch)

    f2c_grid_cells = t_f2c_grid_cells( &
         & num_edges = gen_repr5(X_COMMA_LBOUND_X(patch%cells%num_edges)), &
         & area = gen_repr5(X_COMMA_LBOUND_X(patch%cells%area)), &
         & center = gen_repr5(X_COMMA_LBOUND_X(patch%cells%center)), &
         & edge_idx_f = gen_repr5(X_COMMA_LBOUND_X(patch%cells%edge_idx)), &
         & edge_blk_f = gen_repr5(X_COMMA_LBOUND_X(patch%cells%edge_blk)), &
         & start_block_f = gen_repr5(X_COMMA_LBOUND_X(patch%cells%start_block)), &
         & end_block_f = gen_repr5(X_COMMA_LBOUND_X(patch%cells%end_block)), &
         & start_index_f = gen_repr5(X_COMMA_LBOUND_X(patch%cells%start_index)), &
         & end_index_f = gen_repr5(X_COMMA_LBOUND_X(patch%cells%end_index)) &
         & )

  END SUBROUTINE get_f2c_grid_cells

  SUBROUTINE get_f2c_grid_edges(patch_descr, f2c_grid_edges) BIND(c)
    TYPE(t_f2c_patch_descr), VALUE, INTENT(in) :: patch_descr
    TYPE(t_f2c_grid_edges), INTENT(out) :: f2c_grid_edges

    TYPE(t_patch), POINTER :: patch
    CALL c_f_POINTER(patch_descr%cptr, patch)
    f2c_grid_edges = t_f2c_grid_edges( &
         & primal_normal_cell = gen_repr5(X_COMMA_LBOUND_X(patch%edges%primal_normal_cell)), &
         & primal_normal_vert = gen_repr5(X_COMMA_LBOUND_X(patch%edges%primal_normal_vert)), &
         & dual_normal_cell = gen_repr5(X_COMMA_LBOUND_X(patch%edges%dual_normal_cell)), &
         & dual_normal_vert = gen_repr5(X_COMMA_LBOUND_X(patch%edges%dual_normal_vert)), &
         & tangent_orientation = gen_repr5(X_COMMA_LBOUND_X(patch%edges%tangent_orientation)), &
         & inv_primal_edge_length = gen_repr5(X_COMMA_LBOUND_X(patch%edges%inv_primal_edge_length)), &
         & inv_dual_edge_length = gen_repr5(X_COMMA_LBOUND_X(patch%edges%inv_dual_edge_length)), &
         & inv_vert_vert_length = gen_repr5(X_COMMA_LBOUND_X(patch%edges%inv_vert_vert_length)), &
         & edge_cell_length = gen_repr5(X_COMMA_LBOUND_X(patch%edges%edge_cell_length)), &
         & cell_idx_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%cell_idx)), &
         & cell_blk_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%cell_blk)), &
         & vertex_idx_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%vertex_idx)), &
         & vertex_blk_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%vertex_blk)), &
         & start_block_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%start_block)), &
         & end_block_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%end_block)), &
         & start_index_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%start_index)), &
         & end_index_f = gen_repr5(X_COMMA_LBOUND_X(patch%edges%end_index)) &
         & )

  END SUBROUTINE get_f2c_grid_edges

  SUBROUTINE get_f2c_grid_vertices(patch_descr, f2c_grid_vertices) BIND(c)
    TYPE(t_f2c_patch_descr), VALUE, INTENT(in) :: patch_descr
    TYPE(t_f2c_grid_vertices), INTENT(out) :: f2c_grid_vertices

    TYPE(t_patch), POINTER :: patch
    CALL c_f_POINTER(patch_descr%cptr, patch)
    f2c_grid_vertices = t_f2c_grid_vertices( &
         & vertex = gen_repr5(X_COMMA_LBOUND_X(patch%verts%vertex)), &
         & cell_idx_f = gen_repr5(X_COMMA_LBOUND_X(patch%verts%cell_idx)), &
         & cell_blk_f = gen_repr5(X_COMMA_LBOUND_X(patch%verts%cell_blk)), &
         & start_block_f = gen_repr5(X_COMMA_LBOUND_X(patch%verts%start_block)), &
         & end_block_f = gen_repr5(X_COMMA_LBOUND_X(patch%verts%end_block)), &
         & start_index_f = gen_repr5(X_COMMA_LBOUND_X(patch%verts%start_index)), &
         & end_index_f = gen_repr5(X_COMMA_LBOUND_X(patch%verts%end_index)) &
         & )

  END SUBROUTINE get_f2c_grid_vertices


  FUNCTION get_comm_pattern(cdescr) RESULT(pat_ptr)
    CHARACTER(len=*), PARAMETER :: context = module_name //'::get_comm_pattern'
    TYPE(t_comm_pattern_cdescr), INTENT(in) :: cdescr
    CLASS(t_comm_pattern), POINTER :: pat_ptr
    TYPE(t_comm_pattern_descr), pointer :: descr_ptr
    CALL c_f_POINTER(cdescr%cptr, descr_ptr)
    SELECT TYPE (ptr => descr_ptr%ptr)
    CLASS is (t_comm_pattern)
      pat_ptr => ptr
    CLASS DEFAULT
      CALL finish(context,"cannot resolve descr_ptr%ptr")
    END SELECT
  END FUNCTION get_comm_pattern

  FUNCTION get_mo_model_domain_p_patch_descr(id) RESULT(f2c_patch_descr) BIND(c)
    CHARACTER(len=*), PARAMETER :: context = module_name //'::get_mo_model_domain_p_patch_descr'
    INTEGER(c_int), VALUE, INTENT(in) :: id
    TYPE(t_f2c_patch_descr) :: f2c_patch_descr
    f2c_patch_descr%cptr = c_null_ptr
    IF (.NOT. ALLOCATED(p_patch)) CALL finish(context,'p_patch not allocated')
    IF (id >= LBOUND(p_patch,1) .AND. id <= UBOUND(p_patch,1)) THEN
      IF (p_patch(id)%id /= id) CALL finish(context,'unexpected patch%id')
      f2c_patch_descr%cptr = c_LOC(p_patch(id))
    ENDIF
  END FUNCTION get_mo_model_domain_p_patch_descr

  SUBROUTINE f2c_message(cname, cname_len, ctext, ctext_len) BIND(c)
    CHARACTER(c_char), DIMENSION(*), INTENT(in) :: cname, ctext
    INTEGER(c_int), VALUE :: cname_len, ctext_len
    INTEGER :: my_fname_len, my_ftext_len
    my_fname_len = MIN(cname_len, filename_max)
    my_ftext_len = MIN(ctext_len, filename_max)
    CALL my_msg()
  CONTAINS
    SUBROUTINE my_msg
      CHARACTER(len=my_fname_len) :: my_fname
      CHARACTER(len=my_ftext_len) :: my_ftext
      CALL aux_copy_cstring(cname, my_fname, my_fname_len)
      CALL aux_copy_cstring(ctext, my_ftext, my_ftext_len)
      CALL message(my_fname, my_ftext)
    END SUBROUTINE my_msg
  END SUBROUTINE f2c_message

  SUBROUTINE f2c_finish(cname, cname_len, ctext, ctext_len) BIND(c)
    CHARACTER(c_char), DIMENSION(*), INTENT(in) :: cname, ctext
    INTEGER(c_int), VALUE :: cname_len, ctext_len
    INTEGER :: my_fname_len, my_ftext_len
    my_fname_len = MIN(cname_len, filename_max)
    my_ftext_len = MIN(ctext_len, filename_max)
    CALL my_fin()
  CONTAINS
    SUBROUTINE my_fin
      CHARACTER(len=my_fname_len) :: my_fname
      CHARACTER(len=my_ftext_len) :: my_ftext
      CALL aux_copy_cstring(cname, my_fname, my_fname_len)
      IF (my_ftext_len > 0) THEN
        CALL aux_copy_cstring(ctext, my_ftext, my_ftext_len)
        CALL finish(my_fname, my_ftext)
      ELSE
        CALL finish(my_fname)
      ENDIF
    END SUBROUTINE my_fin
  END SUBROUTINE f2c_finish

  SUBROUTINE aux_copy_cstring(cstr, fstr, n)
    INTEGER, INTENT(in) :: n
    CHARACTER(c_char), DIMENSION(n), INTENT(in) :: cstr
    CHARACTER(len=n), INTENT(out) :: fstr
    INTEGER :: i
    DO i = 1, n
      fstr(i:i) = cstr(i)
    ENDDO
  END SUBROUTINE aux_copy_cstring

  FUNCTION f2c_new_timer(cname, cname_len) RESULT(itimer) BIND(c)
    CHARACTER(c_char), DIMENSION(*), INTENT(in) :: cname
    INTEGER(c_int), VALUE, INTENT(in) :: cname_len
    INTEGER(c_int) :: itimer
    INTEGER :: my_fname_len
    my_fname_len = MIN(cname_len, filename_max)
    CALL my_sub()
  CONTAINS
    SUBROUTINE my_sub
      CHARACTER(len=my_fname_len) :: my_fname
      CALL aux_copy_cstring(cname, my_fname, my_fname_len)
      itimer = new_timer(my_fname)
    END SUBROUTINE my_sub
  END FUNCTION f2c_new_timer

  SUBROUTINE f2c_timer_start(it) BIND(c)
    INTEGER(c_int), VALUE, INTENT(in) :: it
    CALL timer_start(it)
  END SUBROUTINE f2c_timer_start

  SUBROUTINE f2c_timer_stop(it) BIND(c)
    INTEGER(c_int), VALUE, INTENT(in) :: it
    CALL timer_stop(it)
  END SUBROUTINE f2c_timer_stop

  REAL(c_double) FUNCTION f2c_timer_value(it) BIND(c)
    INTEGER(c_int), VALUE, INTENT(in) :: it
    f2c_timer_value = timer_val(it)
  END FUNCTION f2c_timer_value

  ! ltimer=.false. means level 0 (mo_run_nml forces that as well); with
  ! ltimer=.true. the plain IF (ltimer) guards measure whatever timers_level
  ! says, so the ragnarok timers mirroring them (threshold 1) need level 1
  ! as the floor
  INTEGER(c_int) FUNCTION f2c_get_timers_level() BIND(c)
    f2c_get_timers_level = INT(MERGE(MAX(timers_level, 1), 0, ltimer), c_int)
  END FUNCTION f2c_get_timers_level

END MODULE mo_ragnarok_support
#endif /* __NO_RAGNAROK__ */
