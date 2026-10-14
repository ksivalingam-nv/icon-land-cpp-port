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

! Constants used for the computation of lookup tables of the saturation
! mixing ratio over liquid water (*c_les*) or ice(*c_ies*)

MODULE mo_lookup_tables_constants

  USE mo_granules_phys, ONLY:  c1es, c2es, c3les, c3ies, c4les, c4ies,         &
      &                        c5les, c5ies, c5alvcp, c5alscp, alvdcp, alsdcp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies, &
      &       c5alvcp, c5alscp, alvdcp, alsdcp

END MODULE mo_lookup_tables_constants
