!> Contains JSBACH vertical axes structure and methods
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
!>#### Management of vertical axis in ICON-Land
!>
MODULE mo_jsb_vertical_axes
#ifndef __NO_JSBACH__

  USE mo_jsb_vertical_axes_iface, ONLY: t_verticalAxisList, Setup_jsb_vertical_axis, &
    &                                   t_Vgrid, Set_jsb_restart_vgrid
  USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid, t_jsb_vgrid_ptr
  !USE mo_jsb_grid,                ONLY: jsbach_vgrids

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: setup_zaxes_jsbach, set_vertical_grids_jsbach, jsbach_vgrids

  ! register of all vertical grids (global var)
  TYPE(t_jsb_vgrid_ptr), ALLOCATABLE, SAVE :: jsbach_vgrids(:)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_vertical_axes'

CONTAINS

  ! ====================================================================================================== !
  !>
  !>#### Add ICON-Land vertical grids to ICON grids (for output)
  !>
  ! Note: This subroutine is only called from ICON
  !
  SUBROUTINE setup_zaxes_jsbach(vertical_axis_list)

    TYPE(t_verticalAxisList), INTENT(inout) :: vertical_axis_list  !< ICON list of known vertical axis

    TYPE(t_jsb_vgrid), POINTER :: vgrid
    INTEGER :: i

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::setup_zaxes_jsbach"

    IF (.NOT. ALLOCATED(jsbach_vgrids)) RETURN

    ! Append the ICON-Land vgrids to the general ICON registry of vertical axes
    DO i=1,SIZE(jsbach_vgrids)
      vgrid => jsbach_vgrids(i)%p
      CALL Setup_jsb_vertical_axis(vertical_axis_list, &
        &                          vgrid%ZaxisID,      &
        &                          vgrid%name,         &
        &                          vgrid%n_levels,     &
        &                          vgrid%longname,     &
        &                          vgrid%units,        &
        &                          vgrid%levels,       &
        &                          vgrid%lbounds,      &
        &                          vgrid%ubounds       &
        &                         )
    END DO

  END SUBROUTINE setup_zaxes_jsbach

  ! ======================================================================================================= !
  ! Obsolete: not needed anymore after ICON restructuring of restart file reading
  !   TODO: Remove `Set_jsb_restart_vgrid` in ICON mo_util_jsbach.f90
  !
  SUBROUTINE set_vertical_grids_jsbach(vgrid_defs, count)

    TYPE(t_Vgrid), INTENT(inout) :: vgrid_defs(:)
    INTEGER,       INTENT(inout) :: count

    TYPE(t_jsb_vgrid), POINTER :: jsb_vgrid
    INTEGER :: i

    IF (.NOT. ALLOCATED(jsbach_vgrids)) RETURN

    DO i=1,SIZE(jsbach_vgrids)
      jsb_vgrid => jsbach_vgrids(i)%p
      CALL Set_jsb_restart_vgrid(vgrid_defs, count,  &
        &                        jsb_vgrid%ZaxisID,  &
        &                        jsb_vgrid%levels    &
        &                       )
    END DO

  END SUBROUTINE set_vertical_grids_jsbach

#endif
END MODULE mo_jsb_vertical_axes