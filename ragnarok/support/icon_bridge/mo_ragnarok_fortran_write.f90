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

! Fortran fallback logging support before icon_bridge is ready.
! Uses Fortran runtime buffering for correct message order within existing ICON/Fortran logfile.
MODULE mo_ragnarok_fortran_write
  USE ISO_C_BINDING, ONLY: c_int, c_char
  IMPLICIT NONE
  INTEGER, PARAMETER :: default_funit = 0
  PRIVATE

  ! No Fortran-visible public symbols; BIND(C) subroutines are callable directly from C/C++.

CONTAINS

  ! simplified write to a Fortran unit:
  SUBROUTINE ragnarok_write_to_funit(cname, cname_len, ctext, ctext_len, funit) BIND(c)
    INTEGER, PARAMETER :: max_len = 256
    CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: cname, ctext
    INTEGER(c_int), VALUE :: cname_len, ctext_len, funit
    CHARACTER(len=max_len) :: my_fname, my_ftext
    INTEGER :: i, n_name, n_text
    n_name = MIN(cname_len, max_len)
    n_text = MIN(ctext_len, max_len)
    DO i = 1, n_name
      my_fname(i:i) = cname(i)
    END DO
    DO i = 1, n_text
      my_ftext(i:i) = ctext(i)
    END DO
    IF (n_name > 0) THEN
      WRITE(funit, '(a)') TRIM(my_fname(1:n_name))//': '//TRIM(my_ftext(1:n_text))
    ELSE
      WRITE(funit, '(a)') TRIM(my_ftext(1:n_text))
    END IF
  END SUBROUTINE ragnarok_write_to_funit


  ! wrapper for writing to the default Fortran unit:
  SUBROUTINE ragnarok_write_to_default_funit(cname, cname_len, ctext, ctext_len) BIND(c)
    IMPLICIT NONE
    CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: cname, ctext
    INTEGER(c_int), VALUE :: cname_len, ctext_len
    CALL ragnarok_write_to_funit(cname, cname_len, ctext, ctext_len, default_funit)
  END SUBROUTINE ragnarok_write_to_default_funit

END MODULE mo_ragnarok_fortran_write
