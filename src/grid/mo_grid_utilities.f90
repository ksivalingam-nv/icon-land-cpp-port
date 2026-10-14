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

!>
!! Contains grid-geometry-aware utilities: coordinate conversions,
!! arc-length computations, and vector-component transformations
!! for sphere and plane-torus geometries.
!!

! ---------------------------------------------------------------
! LL: xlc has trouble optimizing routines with implicit shaped parameters
! #ifdef __xlC__
! @process nohot
! ! @PROCESS NOOPTIMIZE
! #endif
MODULE mo_grid_utilities

  USE mo_kind, ONLY: wp
  USE mo_math_constants, ONLY: pi, pi_2, dbl_eps
  USE mo_gridman_constants, ONLY: TORUS_MAX_LAT
  USE mo_exception, ONLY: finish
  USE mo_grid_geometry_info, ONLY: t_grid_geometry_info, planar_torus_geometry, sphere_geometry
  USE mo_math_types, ONLY: t_cartesian_coordinates, t_geographical_coordinates, &
    &                       t_tangent_vectors
  USE mo_math_utilities, ONLY: cc_dot_product

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cc2gc
  PUBLIC :: gc2cc
  PUBLIC :: cc2tv
  PUBLIC :: gvec2cvec
  PUBLIC :: cvec2gvec
  PUBLIC :: arc_length
  PUBLIC :: arc_length_v
  PUBLIC :: arc_length_on_unitsphere
  PUBLIC :: cos_arc_length
  PUBLIC :: plane_torus_distance, plane_torus_closest_coordinates
  PUBLIC :: sphere_cartesian_midpoint
  PUBLIC :: sphere_tangent_coordinates
  PUBLIC :: cartesian_to_geographical
  PUBLIC :: geographical_to_cartesian
  PUBLIC :: disp_new
  PUBLIC :: disp_new_vect
  PUBLIC :: spherical_intersection2, spherical_intersection2_v

  INTERFACE arc_length
    MODULE PROCEDURE arc_length_sphere
    MODULE PROCEDURE arc_length_generic
  END INTERFACE
  INTERFACE arc_length_v
    MODULE PROCEDURE arc_length_v_sphere
    MODULE PROCEDURE arc_length_v_generic
  END INTERFACE
  INTERFACE gc2cc
    MODULE PROCEDURE gc2cc_sphere
    MODULE PROCEDURE gc2cc_generic
  END INTERFACE
  INTERFACE cc2gc
    MODULE PROCEDURE cc2gc_sphere
    MODULE PROCEDURE cc2gc_generic
  END INTERFACE

!-----------------------------------------------------------------------
! Basic geometry functions definitions
#define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_grid_utilities'

CONTAINS

  !-------------------------------------------------------------------------
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian.
  !!
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian
  !! ones @f$(p\_cu, p\_cv, p\_cw)@f$ using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & -\sin p\_lat \cdot \cos p\_long \\\
  !! \cos p\_long & -\sin p\_lat \cdot \sin p\_long \\\
  !! 0 & -\cos p\_lat
  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix}
  !! @f}
  !!
  ELEMENTAL SUBROUTINE gvec2cvec(p_gu, p_gv, p_long, p_lat, p_cu, p_cv, p_cw, geometry_info)
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(in)  :: p_gu, p_gv ! zonal and meridional vec. component
    REAL(wp), INTENT(in)  :: p_long, p_lat ! geo. coord. of data point
    TYPE(t_grid_geometry_info), INTENT(in), OPTIONAL :: geometry_info

    REAL(wp), INTENT(out) :: p_cu, p_cv, p_cw ! Cart. vector

    REAL(wp)              :: z_cln, z_sln, z_clt, z_slt ! sin and cos of
    INTEGER               :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grid_utilities:gvec2cvec'
    !-------------------------------------------------------------------------

    geometry_type = sphere_geometry
    IF (PRESENT(geometry_info)) &
      geometry_type = geometry_info%geometry_type

    SELECT CASE (geometry_type)

    CASE (planar_torus_geometry)
      p_cu = p_gu
      p_cv = p_gv
      p_cw = 0._wp
    CASE (sphere_geometry)
      z_sln = SIN(p_long)
      z_cln = COS(p_long)
      z_slt = SIN(p_lat)
      z_clt = COS(p_lat)

      p_cu = z_sln*p_gu + z_slt*z_cln*p_gv
      p_cu = -1._wp*p_cu
      p_cv = z_cln*p_gu - z_slt*z_sln*p_gv
      p_cw = z_clt*p_gv
    CASE DEFAULT
      ! Commented by GZ because this destroys vectorization
      !    CALL finish(method_name, "Undefined geometry type")
      p_cu = p_gu
      p_cv = p_gv
      p_cw = 0._wp
    END SELECT

  END SUBROUTINE gvec2cvec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$.
  !!
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$
  !! into zonal @f$p\_gu@f$ and meridional @f$g\_v@f$ velocity components
  !! using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & \cos p\_long & 0 \\\
  !! -\sin p\_lat \cdot \cos p\_long & -\sin p\_lat \cdot \sin p\_long & \cos p\_lat
  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix}
  !! @f}
  !!
  ELEMENTAL SUBROUTINE cvec2gvec(p_cu, p_cv, p_cw, p_long, p_lat, p_gu, p_gv, geometry_info)
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(in)  :: p_cu, p_cv, p_cw ! Cart. vector
    REAL(wp), INTENT(in)  :: p_long, p_lat ! geo. coord. of data point
    TYPE(t_grid_geometry_info), INTENT(in), OPTIONAL :: geometry_info

    REAL(wp), INTENT(out) :: p_gu, p_gv ! zonal and meridional vec. comp.

    REAL(wp)              :: z_cln, z_clt, z_sln, z_slt ! sin and cos of
    INTEGER               :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grid_utilities:cvec2gvec'
    !-------------------------------------------------------------------------

    geometry_type = sphere_geometry
    IF (PRESENT(geometry_info)) &
      geometry_type = geometry_info%geometry_type

    SELECT CASE (geometry_type)

    CASE (planar_torus_geometry)
      p_gu = p_cu
      p_gv = p_cv
    CASE (sphere_geometry)
      ! p_long and p_lat
      z_sln = SIN(p_long)
      z_cln = COS(p_long)
      z_slt = SIN(p_lat)
      z_clt = COS(p_lat)

      p_gu = z_cln*p_cv - z_sln*p_cu
      p_gv = z_cln*p_cu + z_sln*p_cv
      p_gv = z_slt*p_gv
      p_gv = z_clt*p_cw - p_gv
    CASE DEFAULT
      ! Commented by GZ because this destroys vectorization
      !    CALL finish(method_name, "Undefined geometry type")
      p_gu = p_cu
      p_gv = p_cv
    END SELECT

  END SUBROUTINE cvec2gvec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts  vectors (in)  cartesian coordinate representation
  !! to tangent vectors.
  FUNCTION cc2tv(xx, position) RESULT(tt)

    TYPE(t_cartesian_coordinates), INTENT(in) :: xx
    TYPE(t_geographical_coordinates), INTENT(in)  :: position
    TYPE(t_tangent_vectors):: tt

    REAL(wp) :: z_sinlo, z_coslo, z_sinlacoslo, z_sinlasinlo

    z_sinlo = SIN(position%lon)
    z_coslo = COS(position%lon)
    z_sinlacoslo = SIN(position%lat)*z_coslo
    z_sinlasinlo = SIN(position%lat)*z_sinlo
    tt%v1 = -z_sinlo*xx%x(1) + z_coslo*xx%x(2)
    tt%v2 = -z_sinlacoslo*xx%x(1) - z_sinlasinlo*xx%x(2) + COS(position%lat)*xx%x(3)

  END FUNCTION cc2tv
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE cartesian_to_geographical(cartesian, no_of_points, geo_coordinates, &
    & start_point, end_point)
    TYPE(t_cartesian_coordinates), POINTER :: cartesian(:)
    TYPE(t_geographical_coordinates), POINTER :: geo_coordinates(:)
    INTEGER, INTENT(in) :: no_of_points
    INTEGER, INTENT(in), OPTIONAL :: start_point, end_point

    INTEGER :: i, start_p, end_p

    IF (PRESENT(start_point)) THEN
      start_p = start_point
      end_p = end_point
    ELSE
      start_p = 1
      end_p = no_of_points
    END IF
    IF (.NOT. ASSOCIATED(geo_coordinates)) THEN
      ALLOCATE (geo_coordinates(start_p:end_p), stat=i)
      IF (i > 0) THEN
        CALL finish('cartesian_to_geographical', 'Problem in allocating local arrays')
      END IF
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    DO i = start_p, end_p
      geo_coordinates(i) = cc2gc(cartesian(i))
    END DO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE cartesian_to_geographical
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE geographical_to_cartesian(geo_coordinates, no_of_points, cartesian)
    TYPE(t_geographical_coordinates), POINTER :: geo_coordinates(:)
    TYPE(t_cartesian_coordinates), POINTER :: cartesian(:)
    INTEGER, INTENT(in) :: no_of_points

    INTEGER :: i
    REAL(wp) :: cos_lat
    REAL(wp) :: check

    IF (.NOT. ASSOCIATED(cartesian)) THEN
      ALLOCATE (cartesian(no_of_points), stat=i)
      IF (i > 0) THEN
        CALL finish('geographical_to_cartesian', 'Problem in allocating local arrays')
      END IF
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(i, cos_lat, check)
    DO i = 1, no_of_points
      cos_lat = COS(geo_coordinates(i)%lat)
      cartesian(i)%x(1) = COS(geo_coordinates(i)%lon)*cos_lat
      cartesian(i)%x(2) = SIN(geo_coordinates(i)%lon)*cos_lat
      cartesian(i)%x(3) = SIN(geo_coordinates(i)%lat)

      check = d_norma_3d(cartesian(i))
      IF (ABS(check - 1.0_wp) > 1.0e-8_wp) &
        & CALL finish('geographicalToCartesian', 'Error in calculation')

      !    WRITE(*,*) i, " LonLat=", geoCoordinates(i)
      !    WRITE(*,*) i, " Cart=", cartesian(i)

    END DO ! i=1,noOfPoints
!$OMP END DO
!$OMP END PARALLEL

    RETURN

  END SUBROUTINE geographical_to_cartesian
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns a coordinate system (x,y),
  !! for the tangent plane at the sphere on the given point (lon,lat)
  PURE SUBROUTINE sphere_tangent_coordinates(geo_point, x, y)
    TYPE(t_geographical_coordinates), INTENT(in)  :: geo_point
    TYPE(t_cartesian_coordinates), INTENT(out) :: x, y

    REAL(wp) :: lon, lat
    !-------------------------------------------------------------------------
    lon = geo_point%lon
    lat = geo_point%lat
    x%x(1) = -SIN(lon)
    x%x(2) = COS(lon)
    x%x(3) = 0._wp
    y%x(1) = -COS(lon)*SIN(lat)
    y%x(2) = -SIN(lon)*SIN(lat)
    y%x(3) = COS(lat)

  END SUBROUTINE sphere_tangent_coordinates
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  ELEMENTAL FUNCTION arc_length_sphere(p_x, p_y) RESULT(p_arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x, p_y ! endpoints

    REAL(wp)            :: p_arc ! length of geodesic arc

    REAL(wp)            :: z_lx, z_ly ! length of vector p_x and p_y
    REAL(wp)            :: z_cc ! cos of angle between endpoints

    !-----------------------------------------------------------------------

    z_lx = d_norma_3d(p_x)
    z_ly = d_norma_3d(p_y)

    z_cc = DOT_PRODUCT(p_x%x, p_y%x)/(z_lx*z_ly)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp) z_cc = 1._wp
    IF (z_cc < -1._wp) z_cc = -1._wp

    p_arc = ACOS(z_cc)

  END FUNCTION arc_length_sphere

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  FUNCTION arc_length_generic(p_x, p_y, geometry_info) RESULT(p_arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x, p_y ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    REAL(wp)            :: p_arc ! length of geodesic arc

    INTEGER             :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grid_utilities:arc_length'
    !-----------------------------------------------------------------------

    geometry_type = geometry_info%geometry_type

    SELECT CASE (geometry_type)

    CASE (planar_torus_geometry)
      ! arc_length_sphere returns the radius-normalized arc length,
      ! so for consistency plane_torus_distance is normalized too.
      p_arc = plane_torus_distance(p_x%x, p_y%x, geometry_info)/ &
              geometry_info%sphere_radius
    CASE (sphere_geometry)
      !
      p_arc = arc_length_sphere(p_x, p_y)
    CASE DEFAULT
      !
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END FUNCTION arc_length_generic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints x0,x1.
  ELEMENTAL FUNCTION arc_length_on_unitsphere(x0, x1) RESULT(arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: arc
    REAL(wp) :: z_cc

    z_cc = DOT_PRODUCT(x0%x, x1%x)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp) z_cc = 1._wp
    IF (z_cc < -1._wp) z_cc = -1._wp

    arc = ACOS(z_cc)

  END FUNCTION arc_length_on_unitsphere
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !! Vectorizable version
  !!
  PURE FUNCTION arc_length_v_sphere(p_x, p_y) RESULT(p_arc)
    REAL(wp), INTENT(in) :: p_x(3), p_y(3) ! endpoints

    REAL(wp)            :: p_arc ! length of geodesic arc

    REAL(wp)            :: z_lx, z_ly ! length of vector p_x and p_y
    REAL(wp)            :: z_cc ! cos of angle between endpoints

    !-----------------------------------------------------------------------

    z_lx = SQRT(DOT_PRODUCT(p_x, p_x))
    z_ly = SQRT(DOT_PRODUCT(p_y, p_y))

    z_cc = DOT_PRODUCT(p_x, p_y)/(z_lx*z_ly)

    ! in case we get numerically incorrect solutions

    IF (z_cc > 1._wp) z_cc = 1._wp
    IF (z_cc < -1._wp) z_cc = -1._wp

    p_arc = ACOS(z_cc)

  END FUNCTION arc_length_v_sphere

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !! Vectorizable version
  FUNCTION arc_length_v_generic(p_x, p_y, geometry_info) RESULT(p_arc)
    REAL(wp), INTENT(in) :: p_x(3), p_y(3) ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    REAL(wp)            :: p_arc ! length of geodesic arc

    INTEGER             :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grid_utilities:arc_length'
    !-----------------------------------------------------------------------

    geometry_type = geometry_info%geometry_type

    SELECT CASE (geometry_type)

    CASE (planar_torus_geometry)
      ! arc_length_v_sphere returns the radius-normalized arc length,
      ! so for consistency plane_torus_distance is normalized too.
      p_arc = plane_torus_distance(p_x, p_y, geometry_info)/ &
              geometry_info%sphere_radius
    CASE (sphere_geometry)
      !
      p_arc = arc_length_v_sphere(p_x, p_y)
    CASE DEFAULT
      !
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END FUNCTION arc_length_v_generic
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !! Computes the cosine of the length of geodesic arc with endpoints x0,x1.
  ELEMENTAL FUNCTION cos_arc_length(x0, x1) RESULT(carc)
    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: carc
    REAL(wp) :: z_d0, z_d1, z_cc

    z_d0 = d_norma_3d(x0)
    z_d1 = d_norma_3d(x1)

    z_cc = DOT_PRODUCT(x0%x, x1%x)/(z_d0*z_d1)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp) z_cc = 1._wp
    IF (z_cc < -1._wp) z_cc = -1._wp

    carc = z_cc

  END FUNCTION cos_arc_length
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! returns the torus modulo coordinates of v1 that are
  ! closest to v0
  FUNCTION plane_torus_closest_coordinates(v0, v1, &
    & geometry_info) RESULT(new_v1_coord)

    REAL(wp), INTENT(in) :: v0(3), v1(3)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    REAL(wp) :: length_of_torus, height_of_torus
    TYPE(t_cartesian_coordinates) :: new_v1_coord

    length_of_torus = geometry_info%domain_length
    height_of_torus = geometry_info%domain_height
    ! check the x coordinate
    IF (ABS(v0(1) - v1(1)) > length_of_torus*0.5_wp) THEN
      ! we will wrap around + or -  length_of_torus
      IF (v0(1) > v1(1)) THEN
        new_v1_coord%x(1) = v1(1) + length_of_torus
      ELSE
        new_v1_coord%x(1) = v1(1) - length_of_torus
      END IF
    ELSE
      ! keep the same
      new_v1_coord%x(1) = v1(1)
    END IF

    ! check the y coordinate
    IF (ABS(v0(2) - v1(2)) > height_of_torus*0.5_wp) THEN
      ! we will wrap around + or -  length_of_torus
      IF (v0(2) > v1(2)) THEN
        new_v1_coord%x(2) = v1(2) + height_of_torus
      ELSE
        new_v1_coord%x(2) = v1(2) - height_of_torus
      END IF
    ELSE
      ! keep the same
      new_v1_coord%x(2) = v1(2)
    END IF

    ! this should be zero
    new_v1_coord%x(3) = v1(3)

  END FUNCTION plane_torus_closest_coordinates
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  REAL(wp) FUNCTION plane_torus_distance(v0, v1, &
    & geometry_info)
    REAL(wp), INTENT(in) :: v0(3), v1(3)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    TYPE(t_cartesian_coordinates) :: dv

    dv = plane_torus_closest_coordinates(v0, v1, geometry_info)
    dv%x = dv%x - v0
    plane_torus_distance = d_norma_3d(dv)

  END FUNCTION plane_torus_distance
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION sphere_cartesian_midpoint(p1, p2) RESULT(mid_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1, p2
    TYPE(t_cartesian_coordinates) :: mid_point

    mid_point%x = (p1%x + p2%x)
    mid_point%x = mid_point%x/ &
      & d_norma_3d(mid_point)

  END FUNCTION sphere_cartesian_midpoint
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Finds the intersection of two great circles.
  !!
  !! Cannot be made elemental, as long as the
  !! subroutine calls to message and finish are maintained.
  !!
  ELEMENTAL FUNCTION spherical_intersection2(p0, p1, v0, v1) RESULT(p)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0, p1
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1
    TYPE(t_cartesian_coordinates) :: p

    TYPE(t_geographical_coordinates) :: g1, g2
    TYPE(t_cartesian_coordinates) :: en, cn

    !!$    REAL(wp), PARAMETER :: voronoi_tolerance = 1.0e-9_wp

    REAL(wp) :: z_pn

    ! normal to plane where edge lies
    !  en = vector_product (v0, v1)

    g1 = cc2gc(v0)
    g2 = cc2gc(v1)

    en = vecpro(g1, g2)

    ! normal tp plane where triangle centers lies
    ! cn = vector_product (p0, p1)

    g1 = cc2gc(p0)
    g2 = cc2gc(p1)

    cn = vecpro(g1, g2)

    ! check Voronoi/Delaunay properties
    !!$    IF (ABS(DOT_PRODUCT(en%x, cn%x)) > voronoi_tolerance) THEN
    !!$      WRITE (message_text,'(a,e12.6)') &
    !!$           'Orthogonality = ',                &
    !!$           ABS(DOT_PRODUCT(en%x, cn%x))
    !!$      CALL message ('inter_sect', TRIM(message_text))
    !!$      CALL finish ('inter_sect','Voronoi orthogonality violated.')
    !!$    END IF

    g1 = cc2gc(en)
    g2 = cc2gc(cn)

    !  p = vector_product (en, cn)
    p = vecpro(g1, g2)

    ! rescale on sphere surface
    z_pn = d_norma_3d(p)

    p%x = p%x/z_pn

    IF (DOT_PRODUCT(p0%x, p%x) < 0.0_wp) THEN
      p%x = -p%x
    END IF

  END FUNCTION spherical_intersection2

  !>
  !! Vectorized version of spherical_intersection2. Originally developed for the NEC,
  !! but also used for the Cray compiler because the original code (above) generates a segfault.
  !!
  SUBROUTINE spherical_intersection2_v(p0, p1, v0, v1, p, nlen)
    INTEGER, INTENT(in)                       :: nlen
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0(0:nlen - 1), p1(0:nlen - 1)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0(0:nlen - 1), v1(0:nlen - 1)
    TYPE(t_cartesian_coordinates), INTENT(out):: p(0:nlen - 1)

    TYPE(t_geographical_coordinates) :: g1, g2
    TYPE(t_cartesian_coordinates) :: en, cn

    REAL(wp) :: z_pn
    INTEGER  :: ji

!$OMP PARALLEL DO PRIVATE(g1,g2,en,cn,z_pn)
    DO ji = 0, nlen - 1
      g1 = cc2gc(v0(ji))
      g2 = cc2gc(v1(ji))

      en = vecpro(g1, g2)

      g1 = cc2gc(p0(ji))
      g2 = cc2gc(p1(ji))

      cn = vecpro(g1, g2)

      g1 = cc2gc(en)
      g2 = cc2gc(cn)

      !  p = vector_product (en, cn)
      p(ji) = vecpro(g1, g2)

      ! rescale on sphere surface
      z_pn = d_norma_3d(p(ji))

      p(ji)%x = p(ji)%x/z_pn

      IF (DOT_PRODUCT(p0(ji)%x, p(ji)%x) < 0.0_wp) THEN
        p(ji)%x = -p(ji)%x
      END IF
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE spherical_intersection2_v

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the vector product of of 2 unit vectors in lon,lat
  ELEMENTAL FUNCTION vecpro(geog1, geog2) RESULT(xx)

    TYPE(t_geographical_coordinates), INTENT(in) :: geog2, geog1

    TYPE(t_cartesian_coordinates)  :: xx
    REAL(wp)  :: long1, lat1, long2, lat2
    REAL(wp) :: sn1, sn2, sn3, sn4, cs1, cs2

    long1 = geog1%lon
    lat1 = geog1%lat
    long2 = geog2%lon
    lat2 = geog2%lat

    sn1 = SIN(lat1 - lat2)
    sn2 = SIN(lat1 + lat2)
    sn3 = SIN(0.5_wp*(long1 - long2))
    sn4 = SIN(0.5_wp*(long1 + long2))
    cs1 = COS(0.5_wp*(long1 - long2))
    cs2 = COS(0.5_wp*(long1 + long2))

    xx%x(1) = sn1*sn4*cs1 - sn2*cs2*sn3
    xx%x(2) = sn1*cs2*cs1 + sn2*sn4*sn3
    xx%x(3) = COS(lat1)*COS(lat2)*SIN(long1 - long2)

  END FUNCTION vecpro
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  Calculates coordinates of a grid point (longin, latin) in relation to
  !!  a rotated coordinate system with origin in (reflongin, reflatin)
  !!  with rotation matrix and Euler angles. <br>
  !!  The rotation is done by two basic rotations. First around unit vector in
  !!  north pole direction (z-direction), second around new (y-direction).
  !!  The basic rotation matrices are:
  !!  @f{align*}{
  !!  R_z (\lambda) &=
  !!  \begin{pmatrix}
  !!   \cos\lambda & \sin\lambda & 0 \\\
  !!   -\sin\lambda & \cos\lambda & 0 \\\
  !!   0 & 0 & 1
  !!  \end{pmatrix}\\\
  !!  R_y (\varphi) &=
  !!  \begin{pmatrix}
  !!   \cos\varphi & 0 & \sin\varphi \\\
  !!   0 & 1 & 0 \\\
  !!   -\sin\varphi & 0 &\cos\varphi
  !!  \end{pmatrix}
  !!  @f}
  !!  Thus, the rotation matrix @f$R@f$ is
  !!  @f{align*}{
  !!  R(\lambda,\varphi) ~=~  R_y (\varphi) \cdot R_z (\lambda) &=
  !!  \begin{pmatrix}
  !!   \cos\lambda \cos\varphi & \sin\lambda \cos\varphi & \sin\varphi\\\
  !!   -\sin\lambda & \cos\lambda & 0\\\
  !!   -\cos\lambda \sin\varphi & -\sin\lambda \sin\varphi & \cos\varphi
  !!  \end{pmatrix}
  !!  @f}
  !!  Given the Cartesian coordinates of point
  !!  @f$(\lambda, \varphi) \mapsto (x,y,z)^T@f$
  !!  and applying the rotation matrix
  !!  @f{align*}{
  !!  R(\lambda,\varphi) \cdot (x,y,z)^T &= (1,0,0)^T
  !!  @f}
  !!  this point is the origin of the rotated coordinate system. All other points
  !!  are then rotated to this new system in the similar way by applying the
  !!  rotation matrix to their Cartesian coordinates.<br>
  !!  Algorithm:
  !!  <ol>
  !!  <li>  Given the geographical coordinates of a
  !!   point @f$(\lambda_{in}, \varphi_{in})@f$ and the origin of the rotated system
  !!  @f$(\lambda_{in}^{ref}, \varphi_{in}^{ref})@f$
  !!  <li>  Get Cartesian coordinates of @f$(\lambda_{in}, \varphi_{in})@f$ by
  !!  <i>gc2cc</i>
  !!  <li>  Apply rotation matrix @f$R(\lambda_{in}^{ref}, \varphi_{in}^{ref})@f$
  !!  <li>  Get rotated geographical coordinates @f$(\lambda_{dis}, \varphi_{dis})@f$
  !!  of new Cartesian ones by <i>cc2gc</i>
  !!  </ol>
  !!
  PURE SUBROUTINE disp_new(longin, latin, reflongin, reflatin, dislong, dislat)

    REAL(wp), INTENT(in)   :: longin, latin, reflongin, reflatin
    REAL(wp), INTENT(out)  :: dislong, dislat

    TYPE(t_geographical_coordinates) :: gpoint, gpoint_rot
    TYPE(t_cartesian_coordinates)    :: cpoint, row, cpoint_rot

    REAL(wp)               :: z_long, z_lat, z_reflong, z_reflat
    REAL(wp)               :: z_cos_a, z_cos_b, z_sin_a, z_sin_b
    REAL(wp)               :: z_twopi

    z_twopi = 2._wp*pi

    ! check if all angles are in the wright range
    z_long = longin
    IF (z_long < 0._wp) z_long = z_long + z_twopi

    z_reflong = reflongin
    IF (z_reflong < 0._wp) z_reflong = z_reflong + z_twopi

    z_lat = latin
    z_reflat = reflatin

    ! get cartesian coordinates of (z_long, z_lat)

    gpoint%lon = z_long
    gpoint%lat = z_lat

    cpoint = gc2cc(gpoint)

    ! calculate rotation matrix
    z_cos_a = COS(z_reflong)
    z_sin_a = SIN(z_reflong)
    z_cos_b = COS(z_reflat)
    z_sin_b = SIN(z_reflat)

    ! first row stored in a cartesian vector (just for using cc_dot_product)
    row%x(1) = z_cos_a*z_cos_b
    row%x(2) = z_sin_a*z_cos_b
    row%x(3) = z_sin_b

    ! first component of rotated cpoint
    cpoint_rot%x(1) = cc_dot_product(row, cpoint)

    ! second row stored in a cartesian vector (just for using cc_dot_product)
    row%x(1) = -1._wp*z_sin_a
    row%x(2) = z_cos_a
    row%x(3) = 0._wp

    ! second component of rotated cpoint
    cpoint_rot%x(2) = cc_dot_product(row, cpoint)

    ! third row stored in a cartesian vector (just for using cc_dot_product)
    row%x(1) = -1._wp*z_cos_a*z_sin_b
    row%x(2) = -1._wp*z_sin_a*z_sin_b
    row%x(3) = z_cos_b

    ! third component of rotated cpoint
    cpoint_rot%x(3) = cc_dot_product(row, cpoint)

    ! get geographic coordinates of cpoint_rot
    gpoint_rot = cc2gc(cpoint_rot)

    dislong = gpoint_rot%lon
    dislat = gpoint_rot%lat

  END SUBROUTINE disp_new
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Compute the zonal and meridional components of a vector with respect.
  !!
  !! Compute the zonal and meridional components of a vector with respect
  !! to a rotated coordinate system. This function uses the same
  !! conventions and formalism as "disp_new". The algorithm is as
  !! follows:
  !! 1) convert the vector into Cartesian coordinates
  !! 2) apply the rotation to the vector in cartesian components
  !! 3) transform the rotated vector in geographical coordinates, using
  !! the new values of latitude and longitude.
  !!
  SUBROUTINE disp_new_vect(p_gu, p_gv, lon, lat, barlon, barlat, &
    & new_lon, new_lat, new_p_gu, new_p_gv)

    REAL(wp), INTENT(in) :: &
      & p_gu, p_gv,       & ! original lon-lat components
      & lon, lat,         & ! original lon and lat
      & barlon, barlat,   & ! new origin
      & new_lon, new_lat ! point (lon,lat) in the new system (this
    ! value should be computed with disp_new)
    !
    REAL(wp), INTENT(out) :: &
      & new_p_gu, new_p_gv ! new geographical components
    !
    REAL(wp) :: p_c(3), r1(3, 3), r2(3, 3), rot_p_c(3)
    !
    !--------------------------------------------------------------------

    ! 1) convert to Cartesian coordinates
    CALL gvec2cvec(p_gu, p_gv, lon, lat, p_c(1), p_c(2), p_c(3))

    ! 2) apply the rotation
    r1(1, :) = (/COS(barlon), SIN(barlon), 0.0_wp/)
    r1(2, :) = (/-SIN(barlon), COS(barlon), 0.0_wp/)
    r1(3, :) = (/0.0_wp, 0.0_wp, 1.0_wp/)

    r2(1, :) = (/COS(barlat), 0.0_wp, SIN(barlat)/)
    r2(2, :) = (/0.0_wp, 1.0_wp, 0.0_wp/)
    r2(3, :) = (/-SIN(barlat), 0.0_wp, COS(barlat)/)

    rot_p_c = MATMUL(r2, MATMUL(r1, p_c))

    ! 3) back to geographical coordinates
    CALL cvec2gvec(rot_p_c(1), rot_p_c(2), rot_p_c(3), &
      & new_lon, new_lat, new_p_gu, new_p_gv)

  END SUBROUTINE disp_new_vect
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts cartesian coordinates to geographical.
  !!
  ELEMENTAL FUNCTION cc2gc_sphere(p_x) RESULT(p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x ! Cart. coordinates
    TYPE(t_geographical_coordinates)          :: p_pos ! geo. coordinates

    REAL(wp)                                :: z_x, z_y, z_z, z_r
    !-----------------------------------------------------------------------

    z_x = p_x%x(1)
    z_y = p_x%x(2)
    z_z = p_x%x(3)

    z_r = z_x*z_x + z_y*z_y
    z_r = SQRT(z_r)

    IF (ABS(z_r) < dbl_eps) THEN ! one of the poles

      IF (z_z > 0.0_wp) THEN
        p_pos%lat = pi_2
      ELSE
        p_pos%lat = -1._wp*pi_2
      END IF
      p_pos%lon = 0._wp

    ELSE

      p_pos%lat = ATAN2(z_z, z_r)

      IF (ABS(z_x) < dbl_eps) THEN ! z_x == 0 ?

        IF (z_y >= 0.0_wp) THEN
          p_pos%lon = pi_2
        ELSE
          p_pos%lon = -1._wp*pi_2
        END IF

      ELSE
        p_pos%lon = ATAN2(z_y, z_x)
      END IF

    END IF

  END FUNCTION cc2gc_sphere
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !! Implementation for plane torus geometry.
  !!
  !! The lon-lat parameterization of the torus is
  !!    (lon,lat) = [0, 2*pi] x [-max_lat, max_lat]
  !! where max_lat := pi/18 = 10 degrees  (hard-coded in the torus grid generator)
  !!
  !! The Cartesian coordinates of the torus grid are: v = (x,y,0)
  !! where
  !!    (x,y) = [0,domain_length] x [0,domain_height]
  !! and the lengths are stored as global attributes in the grid file.
  !!
  FUNCTION cc2gc_plane_torus(p_x, geometry_info) RESULT(p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x ! Cart. coordinates
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_geographical_coordinates)          :: p_pos ! geo. coordinates

    p_pos%lon = p_x%x(1)*2._wp*pi/geometry_info%domain_length
    p_pos%lat = p_x%x(2)*2._wp*TORUS_MAX_LAT/geometry_info%domain_height - TORUS_MAX_LAT
  END FUNCTION cc2gc_plane_torus
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !! Generic implementation that handles sphere and other geometries.
  !!
  FUNCTION cc2gc_generic(p_x, geometry_info) RESULT(p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x ! Cart. coordinates
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_geographical_coordinates)          :: p_pos ! geo. coordinates

    CHARACTER(LEN=*), PARAMETER :: method_name = modname//':cc2gc_generic'

    SELECT CASE (geometry_info%geometry_type)
    CASE (planar_torus_geometry)
      p_pos = cc2gc_plane_torus(p_x, geometry_info)
      !
    CASE (sphere_geometry)
      p_pos = cc2gc_sphere(p_x)
      !
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT
  END FUNCTION cc2gc_generic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts geographical to cartesian coordinates.
  !!
  ELEMENTAL FUNCTION gc2cc_sphere(p_pos) RESULT(p_x)
    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos ! geo. coordinates
    TYPE(t_cartesian_coordinates)                :: p_x ! Cart. coordinates
    REAL(wp) :: z_cln, z_sln, z_clt, z_slt

    z_sln = SIN(p_pos%lon)
    z_cln = COS(p_pos%lon)
    z_slt = SIN(p_pos%lat)
    z_clt = COS(p_pos%lat)

    p_x%x(1) = z_cln*z_clt
    p_x%x(2) = z_sln*z_clt
    p_x%x(3) = z_slt

  END FUNCTION gc2cc_sphere
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !! Implementation for plane torus geometry.
  !!
  !! The lon-lat parameterization of the torus is
  !!    (lon,lat) = [0, 2*pi] x [-max_lat, max_lat]
  !! where max_lat := pi/18 = 10 degrees  (hard-coded in the torus grid generator)
  !!
  !! The Cartesian coordinates of the torus grid are: v = (x,y,0)
  !! where
  !!    (x,y) = [0,domain_length] x [0,domain_height]
  !! and the lengths are stored as global attributes in the grid file.
  !!
  FUNCTION gc2cc_plane_torus(p_pos, geometry_info) RESULT(p_x)

    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos ! geo. coordinates
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates)                :: p_x ! Cart. coordinates

    p_x%x(:) = (/p_pos%lon*geometry_info%domain_length/(2._wp*pi),                              &
      &           (p_pos%lat + TORUS_MAX_LAT)*geometry_info%domain_height/(2._wp*TORUS_MAX_LAT), &
      &           0._wp/)
  END FUNCTION gc2cc_plane_torus
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !! Generic implementation that handles sphere and other geometries.
  !!
  FUNCTION gc2cc_generic(p_pos, geometry_info) RESULT(p_x)

    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos ! geo. coordinates
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates)                :: p_x ! Cart. coordinates

    CHARACTER(LEN=*), PARAMETER :: method_name = modname//':gc2cc_generic'

    SELECT CASE (geometry_info%geometry_type)
    CASE (planar_torus_geometry)
      p_x = gc2cc_plane_torus(p_pos, geometry_info)
      !
    CASE (sphere_geometry)
      p_x = gc2cc_sphere(p_pos)
      !
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT
  END FUNCTION gc2cc_generic
  !-------------------------------------------------------------------------

END MODULE mo_grid_utilities
