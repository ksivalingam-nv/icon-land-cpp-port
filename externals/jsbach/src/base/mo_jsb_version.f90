!> Manages ICON-Land version
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
!>#### ICON-Land version control
!>
!> The module defines the ICON-Land version and contains a function to query it
!>
MODULE mo_jsb_version
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_icon_land_version

  CHARACTER(len=*), PARAMETER :: icon_land_version = '2026.04' !< ICON-Land version string

CONTAINS

  !>
  !> Get ICON-Land version string
  !>
  FUNCTION get_icon_land_version() RESULT(version)

    CHARACTER(len=:), ALLOCATABLE :: version

    version = icon_land_version

  END FUNCTION get_icon_land_version

#endif
END MODULE mo_jsb_version
