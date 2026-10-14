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

MODULE mo_grid_geometry_info_types

  USE mo_kind, ONLY: wp
  USE mo_math_types,  ONLY: t_cartesian_coordinates

  IMPLICIT NONE

  PUBLIC
  !--------------------------------------------------------------
  !> Holds the grid geometry parameters
  TYPE t_grid_geometry_info
    INTEGER :: cell_type          ! triangular_cell, etc
    INTEGER :: geometry_type      ! sphere_geometry, etc
    !> The creation process of the grid (cut_off, refined, etc, see parameters).
    INTEGER :: grid_creation_process
    !> The grid optimization procces
    INTEGER :: grid_optimization_process

    TYPE(t_cartesian_coordinates) :: center
    REAL(wp) :: mean_edge_length       ! (meters)
    REAL(wp) :: mean_dual_edge_length  ! (meters)
    REAL(wp) :: mean_cell_area         ! (meters^2)
    REAL(wp) :: mean_dual_cell_area    ! (meters^2)
    REAL(wp) :: domain_length          ! (meters)
    REAL(wp) :: domain_height          ! (meters)

    !> Sphere/Cylinder parameters
    REAL(wp) :: sphere_radius

    !> derived info, the following can be calclulated from the previous
    REAL(wp) :: grid_mean_resolution ! the sqrt(mean_cell_area)
    REAL(wp) :: mean_characteristic_length ! the same as grid_mean_resolution, deprecated, do NOT use

    ! grid geometry source info
    LOGICAL :: grid_file_has_geometry = .true.
    LOGICAL :: geometry_is_computed_onthefly = .false.
    LOGICAL :: grid_file_has_vertex_cartesian = .false.

  END TYPE t_grid_geometry_info
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! public parameters
  ! -----------------------------
  ! types of grid geometries
  INTEGER, PARAMETER ::  sphere_geometry         = 1
  INTEGER, PARAMETER ::  planar_torus_geometry   = 2
  INTEGER, PARAMETER ::  planar_channel_geometry = 3
  INTEGER, PARAMETER ::  planar_geometry         = 4

  ! -----------------------------
  ! types of grids
  INTEGER, PARAMETER ::  triangular_cell = 3
  INTEGER, PARAMETER ::  hexagonal_cell  = 6

  !--------------------------------------------------------------
  INTEGER, PARAMETER ::  undefined    = -1

END MODULE mo_grid_geometry_info_types
!----------------------------------------------------------------------------
