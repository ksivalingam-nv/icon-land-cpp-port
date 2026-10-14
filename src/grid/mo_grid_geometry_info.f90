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

! Contains the definition of basic structures and geometry parameters
! These are included in the grid/patch info

MODULE mo_grid_geometry_info

  USE mo_kind, ONLY: wp
  USE mo_grid_geometry_info_types
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish
  USE mo_netcdf_errhandler,  ONLY: nf
  ! Access NetCDF file metadata using mo_netcdf_parallel
  ! --> single file (p_pe_work == p_io) accesses file and broadcasts result
  USE mo_netcdf_parallel, ONLY:                            &
    & nf90_open               => p_nf90_open,              &
    & nf90_inquire_attribute  => p_nf90_inquire_attribute, &
    & nf90_inq_varid          => p_nf90_inq_varid,         &
    & nf90_get_var            => p_nf90_get_var,           &
    & nf90_get_att            => p_nf90_get_att,           &
    & nf90_close              => p_nf90_close,             &
    & nf90_nowrite, nf90_global, nf90_noerr

  IMPLICIT NONE

  PRIVATE

  ! public parameters
  PUBLIC :: sphere_geometry, planar_torus_geometry, planar_channel_geometry, planar_geometry
  PUBLIC :: triangular_cell, hexagonal_cell

  ! public structures
  PUBLIC :: t_grid_geometry_info

  ! public methods
  PUBLIC :: set_default_geometry_info, copy_grid_geometry_info, &
    & set_grid_geometry_derived_info, read_geometry_info

  CHARACTER(*), PARAMETER :: modname = "mo_grid_geometry_info"

CONTAINS
  !--------------------------------------------------------------
  !>
  ! The default is sphere triangular geometry
  ! Note the cell characteristice (mean area and lenght) are unknown
  ! and set to 0
  SUBROUTINE set_default_geometry_info(to_geometry_info)
    TYPE(t_grid_geometry_info) :: to_geometry_info

    to_geometry_info%cell_type                  = triangular_cell
    to_geometry_info%geometry_type              = sphere_geometry
    to_geometry_info%grid_creation_process      = undefined
    to_geometry_info%grid_optimization_process  = undefined
    to_geometry_info%center %x(:)               = 0.0_wp
    to_geometry_info%mean_edge_length           = 0.0_wp
    to_geometry_info%mean_dual_edge_length      = 0.0_wp
    to_geometry_info%mean_cell_area             = 0.0_wp
    to_geometry_info%mean_dual_cell_area        = 0.0_wp
    to_geometry_info%domain_length              = 2.0_wp * pi * earth_radius
    to_geometry_info%domain_height              = 2.0_wp * pi * earth_radius
    to_geometry_info%sphere_radius              = earth_radius
    to_geometry_info%grid_mean_resolution       = 0.0_wp
    to_geometry_info%mean_characteristic_length = 0.0_wp

  END SUBROUTINE set_default_geometry_info
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  !>
  SUBROUTINE copy_grid_geometry_info(from_geometry_info, to_geometry_info)
    TYPE(t_grid_geometry_info) :: from_geometry_info, to_geometry_info

    to_geometry_info%cell_type                  = from_geometry_info%cell_type
    to_geometry_info%geometry_type              = from_geometry_info%geometry_type
    to_geometry_info%grid_creation_process      = from_geometry_info%grid_creation_process
    to_geometry_info%grid_optimization_process  = from_geometry_info%grid_optimization_process
    to_geometry_info%center                     = from_geometry_info%center
    to_geometry_info%mean_edge_length           = from_geometry_info%mean_edge_length
    to_geometry_info%mean_dual_edge_length      = from_geometry_info%mean_dual_edge_length
    to_geometry_info%mean_cell_area             = from_geometry_info%mean_cell_area
    to_geometry_info%mean_dual_cell_area        = from_geometry_info%mean_dual_cell_area
    to_geometry_info%domain_length              = from_geometry_info%domain_length
    to_geometry_info%domain_height              = from_geometry_info%domain_height
    to_geometry_info%sphere_radius              = from_geometry_info%sphere_radius
    to_geometry_info%grid_mean_resolution       = from_geometry_info%grid_mean_resolution
    to_geometry_info%mean_characteristic_length = from_geometry_info%grid_mean_resolution

  END SUBROUTINE copy_grid_geometry_info
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE set_grid_geometry_derived_info( to_geometry_info )
    TYPE(t_grid_geometry_info) :: to_geometry_info

    ! derived geometry parameters
    to_geometry_info%grid_mean_resolution = SQRT(to_geometry_info%mean_cell_area)
    to_geometry_info%mean_characteristic_length = to_geometry_info%grid_mean_resolution
!     write(0,*) "------------------------------------------------"
!     write(0,*) "mean_cell_area:", to_geometry_info%mean_cell_area
!     write(0,*) "mean_characteristic_length:", to_geometry_info%mean_characteristic_length
!     write(0,*) "------------------------------------------------"

  END SUBROUTINE set_grid_geometry_derived_info
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION read_geometry_info(ncid, geometry_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_grid_geometry_info) :: geometry_info

    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_geometry_info"

    INTEGER :: geometry_type
    INTEGER :: cell_type
    REAL(wp) :: mean_edge_length
    REAL(wp) :: mean_dual_edge_length
    REAL(wp) :: mean_cell_area
    REAL(wp) :: mean_dual_cell_area
    REAL(wp) :: domain_length
    REAL(wp) :: domain_height
    REAL(wp) :: sphere_radius
    REAL(wp) :: center_x(3)

    CALL set_default_geometry_info(geometry_info)
    read_geometry_info = -1

    geometry_info%grid_file_has_geometry = gridfile_has_geometry(ncid)

    netcd_status = nf90_get_att(ncid, nf90_global,'grid_geometry', geometry_type)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","grid_geometry")
      RETURN
    ENDIF

    geometry_info%geometry_type = geometry_type

    netcd_status = nf90_get_att(ncid, nf90_global,'grid_cell_type', cell_type)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","grid_geometry")
      RETURN
    ENDIF

    geometry_info%cell_type = cell_type

    netcd_status = nf90_get_att(ncid, nf90_global,'mean_edge_length', mean_edge_length)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","mean_edge_length")
      RETURN
    ENDIF

    geometry_info%mean_edge_length = mean_edge_length

    netcd_status = nf90_get_att(ncid, nf90_global,'mean_dual_edge_length', mean_dual_edge_length)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","mean_dual_edge_length")
      RETURN
    ENDIF

    geometry_info%mean_dual_edge_length = mean_dual_edge_length

    netcd_status = nf90_get_att(ncid, nf90_global,'mean_cell_area', mean_cell_area)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","mean_cell_area")
      RETURN
    ENDIF

    geometry_info%mean_cell_area = mean_cell_area

    netcd_status = nf90_get_att(ncid, nf90_global,'mean_dual_cell_area', mean_dual_cell_area)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","mean_dual_cell_area")
      RETURN
    ENDIF

    geometry_info%mean_dual_cell_area = mean_dual_cell_area

    netcd_status = nf90_get_att(ncid, nf90_global,'domain_length', domain_length)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","domain_length")
      RETURN
    ENDIF

    geometry_info%domain_length = domain_length

    netcd_status = nf90_get_att(ncid, nf90_global,'domain_height', domain_height)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","domain_height")
      RETURN
    ENDIF

    geometry_info%domain_height = domain_height

    netcd_status = nf90_get_att(ncid, nf90_global,'sphere_radius', sphere_radius)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","sphere_radius")
      RETURN
    ENDIF

    geometry_info%sphere_radius = sphere_radius

    netcd_status = nf90_get_att(ncid, nf90_global,'domain_cartesian_center', center_x)
    IF (netcd_status /= nf90_noerr) THEN
!       CALL finish("Cannot read","domain_cartesian_center")
      RETURN
    ENDIF

    geometry_info%center%x = center_x

    read_geometry_info = 0

  END FUNCTION read_geometry_info
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION gridfile_has_geometry(ncid)
    INTEGER, INTENT(in)    :: ncid

    INTEGER :: varid

    gridfile_has_geometry = &
      nf90_inq_varid(ncid, 'cell_area', varid) == nf90_noerr

  END FUNCTION gridfile_has_geometry
  !-------------------------------------------------------------------------

END MODULE mo_grid_geometry_info
!----------------------------------------------------------------------------
