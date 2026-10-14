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

! Copy of ICON's mo_kind to allow for a full encapsulation of granules.
! Including a check to guaranetee consistency between ICON and the granules.

MODULE mo_granules_kind
  USE, INTRINSIC :: iso_c_binding, ONLY: c_float, c_double
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sp, dp, wp, vp, rp
  PUBLIC :: i1, i2, i4, i8
  PUBLIC :: check_granule_precision_consistency
  !--------------------------------------------------------------------
  !
  ! Floating point section
  ! ----------------------
  !
  ! portable floating point kind type parameters: sp, dp
  INTEGER, PARAMETER :: sp = c_float                   !< single precision
  INTEGER, PARAMETER :: dp = c_double                  !< double precision

#if defined(__SINGLE_PRECISION)
  INTEGER, PARAMETER :: wp = sp                        !< selected working precision
  INTEGER, PARAMETER :: vp = sp                        !< selected variable precision
  INTEGER, PARAMETER :: rp = sp                        !< selected nwp radiation precision
#else

#if defined(__MIXED_PRECISION)
  INTEGER, PARAMETER :: wp = dp                        !< selected working precision
  INTEGER, PARAMETER :: vp = sp                        !< selected variable precision
#else
  INTEGER, PARAMETER :: wp = dp                        !< selected working precision
  INTEGER, PARAMETER :: vp = dp                        !< selected variable precision
#endif
#ifdef __SINGLE_PRECISION_ECRAD
  INTEGER, PARAMETER :: rp = sp                        !< selected nwp radiation precision
#else
  INTEGER, PARAMETER :: rp = dp                        !< selected nwp radiation precision
#endif

#endif

  !
  ! Integer section
  ! ---------------
  !
  INTEGER, PARAMETER :: pi1 =  2
  INTEGER, PARAMETER :: pi2 =  4
  INTEGER, PARAMETER :: pi4 =  9
  INTEGER, PARAMETER :: pi8 = 14  ! could be larger, but SX cannot do some operations otherwise
  !
  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(pi1)   !< at least 1 byte integer
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(pi2)   !< at least 2 byte integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)   !< at least 4 byte integer
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)   !< at least 8 byte integer
  !
  !
  !--------------------------------------------------------------------
CONTAINS
  !
  ! check consistency with provided precisions (from ICON)
  !
  ! return true if and only if all checks are passed
  !
  FUNCTION check_granule_precision_consistency(sp_ext, dp_ext, wp_ext, vp_ext, rp_ext) RESULT(is_okay)
    INTEGER, INTENT(IN) :: sp_ext, dp_ext, wp_ext, vp_ext, rp_ext !< external (ICON) precision
    LOGICAL :: is_okay
    is_okay =(sp == sp_ext .AND. dp == dp_ext .AND. wp == wp_ext .AND. vp == vp_ext .AND. rp == rp_ext)
  END FUNCTION check_granule_precision_consistency
END MODULE mo_granules_kind
