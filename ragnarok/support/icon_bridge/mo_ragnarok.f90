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

MODULE mo_ragnarok
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_ragnarok, finalize_ragnarok, getKokkosVersion
  PUBLIC :: ragnarok_ensure_empty_memory_pool

  INTERFACE

    SUBROUTINE init_ragnarok_c(is_stdio_proc, device_id, stream_id) BIND(c, name="init_ragnarok_c")
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(c_int), VALUE :: is_stdio_proc
      INTEGER(c_int), VALUE :: device_id
      INTEGER(c_intptr_t), VALUE :: stream_id
    END SUBROUTINE init_ragnarok_c


    SUBROUTINE finalize_ragnarok() BIND(c, name="finalize_ragnarok")
    END SUBROUTINE finalize_ragnarok


    SUBROUTINE ragnarok_ensure_empty_memory_pool() BIND(c)
    END SUBROUTINE ragnarok_ensure_empty_memory_pool


    FUNCTION get_kokkos_version_c(len) RESULT(ret) BIND(c, name="retrieve_kokkos_version_c")
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(c_int), INTENT(INOUT) :: len
      TYPE(c_ptr) :: ret
    END FUNCTION get_kokkos_version_c

  END INTERFACE

  PUBLIC :: get_kokkos_version_c

CONTAINS

  SUBROUTINE getKokkosVersion(str)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    CHARACTER(kind=c_char), POINTER, INTENT(OUT) :: str(:)

    TYPE(c_ptr) :: str_c
    INTEGER(c_int) :: len

    str_c = get_kokkos_version_c(len)
    CALL C_F_POINTER(str_c, str, (/len/))
  END SUBROUTINE getKokkosVersion


  SUBROUTINE init_ragnarok(is_stdio_proc, device_id, stream_id)
    USE, INTRINSIC :: ISO_C_BINDING
    LOGICAL, INTENT(IN) :: is_stdio_proc
    INTEGER, INTENT(IN), OPTIONAL :: device_id
    INTEGER(c_intptr_t), INTENT(IN), OPTIONAL :: stream_id
    INTEGER(c_int) :: c_device_id
    INTEGER(c_intptr_t) :: c_stream_id

    c_device_id = -1_c_int
    IF (PRESENT(device_id) .AND. PRESENT(stream_id)) THEN
      c_device_id = INT(device_id, c_int)
      c_stream_id = stream_id
    ELSE IF (PRESENT(device_id) .OR. PRESENT(stream_id)) THEN
      ERROR STOP "init_ragnarok: device_id and stream_id must be provided together"
    ELSE
      c_device_id = -1_c_int
      c_stream_id = -1_c_intptr_t
    END IF

    CALL init_ragnarok_c(MERGE(1_c_int, 0_c_int, is_stdio_proc), c_device_id, c_stream_id)
  END SUBROUTINE init_ragnarok


END MODULE mo_ragnarok
