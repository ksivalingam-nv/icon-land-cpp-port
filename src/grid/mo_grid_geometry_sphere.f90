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

! Contains the routines for calculating grid geometry properties derived from the vertex cartesian coordinates
! Currently only the geometry on the shpere for trangle grids is avaliable
! Requires the position of the vertices to be known.
!-------------------------------------------------------------------------------------
MODULE mo_grid_geometry_sphere
#include "icon_definitions.inc"
#define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))

  !-------------------------------------------------------------------------
  USE mo_kind,            ONLY: wp, dp, qp
  USE mo_exception,       ONLY: finish, message

  USE mo_math_types,      ONLY: t_geographical_coordinates, t_cartesian_coordinates
  USE mo_math_constants,  ONLY: pi, pi_2

  USE mo_grid_subset,     ONLY: t_subset_range, get_index_range
  USE mo_model_domain,    ONLY: t_patch, t_grid_cells, t_grid_edges, t_grid_vertices, sphere_geometry, triangular_cell
  USE mo_grid_config,     ONLY: grid_sphere_radius, grid_length_rescale_factor

  USE mo_statistics,      ONLY: global_minmaxmean
  USE mo_sync,            ONLY: sync_c, sync_e, sync_v, sync_patch_array
  USE mo_communication,   ONLY: exchange_data
  USE mo_mpi,             ONLY: my_process_is_mpi_parallel
  USE mo_util_dbg_prnt,   ONLY: dbg_print
  USE iso_fortran_env

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: compute_geometry

#if defined(__QP_GEOMETRY_CALCULATION) && defined(__HAVE_QUAD_PRECISION)
  INTEGER, PARAMETER :: wqp = qp
#else
  INTEGER, PARAMETER :: wqp = dp
#endif

  TYPE t_cartesian_coordinates_qp
    REAL(wqp) :: x(3)
  END TYPE t_cartesian_coordinates_qp

  !---------------------------------------------------------------------------

CONTAINS


  !---------------------------------------------------------------------------
  SUBROUTINE compute_geometry(patch )

    TYPE(t_patch), INTENT(inout), TARGET ::  patch  ! patch data structure

    CHARACTER(LEN=*), PARAMETER :: method_name = 'compute_geometry'

    CALL compute_sphere_geometry_triangles(patch)

  END SUBROUTINE compute_geometry
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! computes the grid geometry porperties based on the vertex cartesian coordinates
  ! and fills the related geo-coordinates
  !
  ! Note 1: if the vertex cartesian coordinates do not exist in the grid file, then they
  !         are calculated from the vertex lon, lat
  !
  ! Note 2: edges%primal_normal_cell, dual_normal_cell, primal_normal_vert, dual_normal_vert are calculated in in mo_intp_coeffs::complete_patchinfo
  !         edges%quad and butterfly variables are calculated in mo_model_domimp_setup
  !
  SUBROUTINE compute_sphere_geometry_triangles(patch)
    TYPE(t_patch), TARGET :: patch

    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    !-------------------------------------------------------------------------

    TYPE(t_cartesian_coordinates_qp) :: cartesian_v(3)
    TYPE(t_cartesian_coordinates_qp) :: cartesian_c(2)
    TYPE(t_cartesian_coordinates_qp) :: edge_cartesian_center,circle_cartesian_center,edge_vector, edge_normal_vector, x,y
    TYPE(t_cartesian_coordinates_qp) :: circumcenters_vector

    REAL(wp) :: sphere_radius_squared,lon,lat
    REAL(wp) :: minmaxmean(3)
    REAL(wqp) :: dual_area

    INTEGER :: jb, je, i
    INTEGER :: start_vertex_index, end_vertex_index
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: edge_idx,edge_blk
    INTEGER :: cell_idx(2), cell_blk(2)

    CHARACTER(LEN=*), PARAMETER :: method_name = 'compute sphere geometry on triangles'
    !----------------------------------------------------------

! to be computed:
!
! cells%edge_orientation
! cells%area
! cells%cartesian_center
!
! edges%cartesian_center
! edges%cartesian_dual_middle
! edges%primal_cart_normal
! edges%dual_cart_normal
! edges%inv_primal_edge_length
! edges%tangent_orientation
! edges%center
! edges%primal_normal
! edges%dual_normal
! edges%primal_edge_length
! edges%dual_edge_length
! edges%inv_dual_edge_length
! edges%area_edge
! edges%edge_vert_length(:,:,1,2)
! edges%edge_cell_length(:,:,1,2)
!
! verts%cartesian
! verts%edge_orientation
! verts%dual_area

    CALL message(method_name, TRIM(patch%grid_filename))

    verts=>patch%verts
    edges=>patch%edges
    cells=>patch%cells
    !----------------------------------------------------------

   ! works only for sphere geometry and triangles
    IF (patch%geometry_info%geometry_type /= sphere_geometry) &
      CALL finish(method_name, "geometry_type /= sphere_geometry")

    IF (cells%max_connectivity /= 3) THEN
      CALL finish(method_name,'Not a triangular grid!')
    ENDIF
    IF (patch%geometry_info%cell_type /= triangular_cell) THEN
      CALL finish(method_name,'Not a triangular grid!')
    ENDIF
    IF (patch%geometry_info%geometry_type /= sphere_geometry) THEN
      CALL finish(method_name,'Grid not on the sphere!')
    ENDIF
    !----------------------------------------------------------
    patch%geometry_info%geometry_is_computed_onthefly = .true.
    !----------------------------------------------------------

!$OMP PARALLEL

! !$  PRINT*,  "OMP_GET_NUM_THREADS=", OMP_GET_NUM_THREADS()

    IF (.not. patch%geometry_info%grid_file_has_vertex_cartesian) THEN
!$OMP DO PRIVATE(jb,je,i,start_vertex_index,end_vertex_index)
      DO jb = verts%all%start_block, verts%all%end_block
        CALL get_index_range(verts%all, jb, start_vertex_index, end_vertex_index)
        DO je = start_vertex_index, end_vertex_index
          verts%cartesian(je,jb) = gc2cc(verts%vertex(je,jb))
          verts%cartesian(je,jb)%x = verts%cartesian(je,jb)%x / &
            & d_norma_3d(verts%cartesian(je,jb))
        ENDDO
      ENDDO
!$OMP BARRIER
    ENDIF

!   compute cell circumecneters and areas
!$OMP DO PRIVATE(jb, je, start_cell_index, end_cell_index, cartesian_v, circle_cartesian_center)
    DO jb = cells%all%start_block, cells%all%end_block
      CALL get_index_range(cells%all, jb, start_cell_index, end_cell_index)
      DO je = start_cell_index, end_cell_index

        cartesian_v(1)%x = verts%cartesian(cells%vertex_idx(je,jb,1), cells%vertex_blk(je,jb,1))%x
        cartesian_v(2)%x = verts%cartesian(cells%vertex_idx(je,jb,2), cells%vertex_blk(je,jb,2))%x
        cartesian_v(3)%x = verts%cartesian(cells%vertex_idx(je,jb,3), cells%vertex_blk(je,jb,3))%x

        circle_cartesian_center = &
          & circum_center(cartesian_v(1), cartesian_v(2), cartesian_v(3))
        cells%cartesian_center(je, jb)%x =  REAL(circle_cartesian_center%x, wp)
        cells%area(je, jb) = REAL(sphere_triangle_area_Huilier(cartesian_v(1),cartesian_v(2),cartesian_v(3)),wp)
        cells%center(je, jb) = cc2gc(circle_cartesian_center) ! this is read from the file
      ENDDO
    ENDDO ! all cells
!$OMP END DO

    ! compute cell-edge orientation
!$OMP DO PRIVATE(jb, je, i, start_cell_index, end_cell_index, edge_idx,edge_blk)
    DO jb = cells%owned%start_block, cells%owned%end_block
      CALL get_index_range(cells%owned, jb, start_cell_index, end_cell_index)
      DO je = start_cell_index, end_cell_index
        ! compute the cell  orientation
        DO i=1,cells%max_connectivity
          edge_idx=cells%edge_idx(je,jb,i)
          edge_blk=cells%edge_blk(je,jb,i)
          ! compute edge orientation according to Gauss formula (positive points outwards)
          ! if this cell is the first of the edge, it points outwards
          IF (edge_idx > 0) THEN
            IF (je == edges%cell_idx(edge_idx,edge_blk,1) .AND. &
              & jb == edges%cell_blk(edge_idx,edge_blk,1) ) THEN
              cells%edge_orientation(je,jb,i) = 1.0_wp
            ELSE
              cells%edge_orientation(je,jb,i) = -1.0_wp
            ENDIF
          ENDIF
        ENDDO

      ENDDO
    ENDDO ! cells owned, cells%edge_orientation need to be synced
!$OMP END DO


    ! compute edge lengths and orientation vectors
!$OMP DO PRIVATE(jb, je, start_edge_index, end_edge_index, cartesian_v, edge_cartesian_center, &
!$OMP cell_idx, cell_blk, cartesian_c,edge_vector, edge_normal_vector, circumcenters_vector, x, y, lon, lat)
    DO jb = edges%all%start_block, edges%all%end_block
      CALL get_index_range(edges%all, jb, start_edge_index, end_edge_index)
      DO je = start_edge_index, end_edge_index

        ! the middle point of the edge
        cartesian_v(1)%x = verts%cartesian(edges%vertex_idx(je, jb, 1), edges%vertex_blk(je, jb, 1) )%x
        cartesian_v(2)%x = verts%cartesian(edges%vertex_idx(je, jb, 2), edges%vertex_blk(je, jb, 2) )%x
        edge_cartesian_center = sphere_cartesian_midpoint(cartesian_v(1), cartesian_v(2))
        edges%cartesian_center(je, jb)%x = REAL(edge_cartesian_center%x,wp)
        edges%center(je, jb) = cc2gc(edge_cartesian_center)
        !------------------------------------------

        edges%primal_edge_length(je, jb) = &
          & REAL(arc_length_on_unitsphere (cartesian_v(1), cartesian_v(2)),wp)

        ! these are calculated in the sphere rescaling at the end
!         edges%inv_primal_edge_length(je,jb) = &
!           & 1._wp/edges%primal_edge_length(je,jb)
!         edges%edge_vert_length(je,jb,1) = 0.5_wp*edges%primal_edge_length(je,jb)
!         edges%edge_vert_length(je,jb,2) = edges%edge_vert_length(je,jb,1)

        !------------------------------------------
        ! cell related edge geometry
        cell_idx(1) = edges%cell_idx(je,jb,1)
        cell_blk(1) = edges%cell_blk(je,jb,1)
        cell_idx(2) = edges%cell_idx(je,jb,2)
        cell_blk(2) = edges%cell_blk(je,jb,2)
        !------------------------------------------
!         IF ( cell_idx(1) < 1 .and.  cell_idx(2) < 1 ) &
!           & CALL finish(method_name, "edge found without any neigboring cells!")

        IF (cell_idx(1) > 0) THEN
          cartesian_c(1)%x = cells%cartesian_center(cell_idx(1), cell_blk(1))%x
          edges%edge_cell_length(je,jb,1) = &
            & REAL(arc_length_on_unitsphere (edge_cartesian_center, cartesian_c(1)),wp)
        ELSE
          cartesian_c(1) = edge_cartesian_center
          edges%edge_cell_length(je,jb,1) = 0.0_wp
        ENDIF

        IF (cell_idx(2) > 0) THEN
          cartesian_c(2)%x = cells%cartesian_center(cell_idx(2), cell_blk(2))%x
          edges%edge_cell_length(je,jb,2) = &
            & arc_length_on_unitsphere (edge_cartesian_center, cartesian_c(2))
        ELSE
          cartesian_c(2) = edge_cartesian_center
          edges%edge_cell_length(je,jb,2) = 0.0_wp
        ENDIF

        x = sphere_cartesian_midpoint(cartesian_c(1), cartesian_c(2)) ! temporarily use x to cast to the cartesian_dual_middle
        edges%cartesian_dual_middle(je,jb)%x = REAL(x%x,wp)
        !------------------------------------------
        ! these are calculated in the sphere rescaling at the end
!         edges%dual_edge_length(je,jb) = &
!             & edges%edge_cell_length(je,jb,1) + edges%edge_cell_length(je,jb,2)
!         edges%inv_dual_edge_length(je,jb) = &
!           & 1._wp/edges%dual_edge_length(je,jb)

        !------------------------------------------
        ! compute normals
        !------------------------------------------
        ! get the vector tangent on the sphere at the edge center, and normal to the edge
        ! this is the primal_cart_normal variable (times 1 or -1, depending the orientation of the cell centers)
        !-------------------------------------------------------------------------
        edge_vector%x = cartesian_v(2)%x - cartesian_v(1)%x
        ! the normal to the edge, as a vector product of the normal to the sphere
        ! times the edge vector, so that the normal-tangent forms a left-handed system
        edge_normal_vector   = vector_product(edge_cartesian_center, edge_vector)

        !------------------------------------------
        ! compute edges vector orientation
        ! tangent_orientation=1 if the tangent vector points from v1 to v2
        ! otherwise = -1
        !------------------------------------------
        edges%tangent_orientation(je,jb) = 1        ! edge_vector at the moment points from v1 to v2
        circumcenters_vector%x = cartesian_c(2)%x - cartesian_c(1)%x
        IF (DOT_PRODUCT(edge_normal_vector%x,circumcenters_vector%x) < 0.0_wp) THEN
          edge_vector%x = - edge_vector%x
          edge_normal_vector%x = - edge_normal_vector%x
          edges%tangent_orientation(je,jb) = -1  ! edge_vector now points from v2 to v1 (this is not needed if we just switch v1 and v2)
        ENDIF
        ! normalize
        edge_vector%x = edge_vector%x /d_norma_3d(edge_vector)
        edge_normal_vector%x = edge_normal_vector%x / &
        & d_norma_3d(edge_normal_vector)

        edges%dual_cart_normal(je,jb)%x = REAL(edge_vector%x,wp)
        edges%primal_cart_normal(je,jb)%x = REAL(edge_normal_vector%x,wp)

        ! the lon,lat directions as cartesian vectors at the edge center
        lon = edges%center(je,jb)%lon
        lat = edges%center(je,jb)%lat
        x%x(1) = -SIN(lon)
        x%x(2) =  COS(lon)
        x%x(3) =  0._wp
        y%x(1) = -COS(lon) * SIN(lat)
        y%x(2) = -SIN(lon) * SIN(lat)
        y%x(3) =  COS(lat)
        ! the normal and tangent vectors in terms of local lon,lat
        ! these define the projections to the lon lat coordiantes
        edges%dual_normal(je,jb)%v1 = REAL(DOT_PRODUCT(edge_vector%x, x%x),wp)
        edges%dual_normal(je,jb)%v2 = REAL(DOT_PRODUCT(edge_vector%x, y%x),wp)
        edges%primal_normal(je,jb)%v2 = edges%dual_normal(je,jb)%v1
        edges%primal_normal(je,jb)%v1 = -edges%dual_normal(je,jb)%v2

       !------------------------------------------
      ENDDO
    ENDDO  ! all edges
!$OMP END DO
!$OMP END PARALLEL

!     CALL dbg_print('all cell center 1:',cells%cartesian_center(:,:)%x(1),"",1, in_subset=cells%all)
!     CALL dbg_print('all cell center 2:',cells%cartesian_center(:,:)%x(2),"",1, in_subset=cells%all)
!     CALL dbg_print('all cell center 3:',cells%cartesian_center(:,:)%x(3),"",1, in_subset=cells%all)
!     CALL dbg_print('c.geo bfsync.:edge-cell1 all',edges%edge_cell_length(:,:,1),"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo bfsync.:edge-cell2 all',edges%edge_cell_length(:,:,2),"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo bfsync:edge-cell1 owned',edges%edge_cell_length(:,:,1),"",1, in_subset=edges%owned)
!     CALL dbg_print('c.geo bfsync:edge-cell2 owned',edges%edge_cell_length(:,:,2),"",1, in_subset=edges%owned)

    !------------------------------------------
    ! here we need to sync the edges%tangent_orientation as it is used in the verts%edge_orientation
    ! so we will just sync all that were computed so sfar
    !------------------------------------------
!     CALL sync_patch_geometry(patch)
    CALL sync_patch_array(SYNC_C, patch, cells%edge_orientation(:,:,1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_C, patch, cells%edge_orientation(:,:,2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_C, patch, cells%edge_orientation(:,:,3), lacc=.FALSE.)

    CALL sync_patch_array(SYNC_E, patch, edges%edge_cell_length(:,:,1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%edge_cell_length(:,:,2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%cartesian_dual_middle(:,:)%x(1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%cartesian_dual_middle(:,:)%x(2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%cartesian_dual_middle(:,:)%x(3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%primal_cart_normal(:,:)%x(1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%primal_cart_normal(:,:)%x(2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%primal_cart_normal(:,:)%x(3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%dual_cart_normal(:,:)%x(1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%dual_cart_normal(:,:)%x(2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%dual_cart_normal(:,:)%x(3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%dual_normal(:,:)%v1, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%dual_normal(:,:)%v2, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%primal_normal(:,:)%v1, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%primal_normal(:,:)%v2, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, edges%tangent_orientation(:,:), lacc=.FALSE.)


!     CALL dbg_print('c.geo afsync.:edge-cell1 all',edges%edge_cell_length(:,:,1),"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo afsync.:edge-cell2 all',edges%edge_cell_length(:,:,2),"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo afsync:edge-cell1 owned',edges%edge_cell_length(:,:,1),"",1, in_subset=edges%owned)
!     CALL dbg_print('c.geo afsync:edge-cell2 owned',edges%edge_cell_length(:,:,2),"",1, in_subset=edges%owned)
    !------------------------------------------
    ! compute verts%dual_area and verts%edge_orientation
!$OMP PARALLEL DO PRIVATE(jb,je,i,start_vertex_index,end_vertex_index,edge_idx,edge_blk,cell_idx, cell_blk, &
!$OMP dual_area, cartesian_c, cartesian_v)
    DO jb = verts%owned%start_block, verts%owned%end_block
      CALL get_index_range(verts%owned, jb, start_vertex_index, end_vertex_index)
      verts%dual_area(:,jb) = 0.0_wqp
      DO je = start_vertex_index, end_vertex_index
        dual_area = 0.0_wqp
!         DO i=1,verts%max_connectivity
        DO i=1,verts%num_edges(je,jb) ! in case of use_duplicated_connectivity and pentagons
          edge_idx = verts%edge_idx(je,jb,i)
          edge_blk = verts%edge_blk(je,jb,i)
          cartesian_v(1)%x = verts%cartesian(je,jb)%x
          IF (edge_idx > 0) THEN
            ! add to the dual area corresponding to the edge
            cell_idx(1) = edges%cell_idx(edge_idx,edge_blk,1)
            cell_blk(1) = edges%cell_blk(edge_idx,edge_blk,1)
            cell_idx(2) = edges%cell_idx(edge_idx,edge_blk,2)
            cell_blk(2) = edges%cell_blk(edge_idx,edge_blk,2)
            !------------------------------------------
            IF (cell_idx(1) > 0) THEN
              cartesian_c(1)%x = cells%cartesian_center(cell_idx(1), cell_blk(1))%x
            ELSE
              cartesian_c(1)%x = edges%cartesian_center(edge_idx,edge_blk)%x
            ENDIF
            IF (cell_idx(2) > 0) THEN
              cartesian_c(2)%x = cells%cartesian_center(cell_idx(2), cell_blk(2))%x
            ELSE
              cartesian_c(2)%x = edges%cartesian_center(edge_idx,edge_blk)%x
            ENDIF

            dual_area = dual_area + &
              sphere_triangle_area_Huilier(cartesian_c(1), cartesian_v(1), cartesian_c(2))

            ! compute edge orientation with respect to this vertex
            ! ie, if the edge tangent points from this vertex to the other.
            ! it requires that the edges%tangent_orientation has been calculated !
            IF (edges%vertex_idx(edge_idx,edge_blk,1) == je .AND. &
              & edges%vertex_blk(edge_idx,edge_blk,1) == jb) THEN
              verts%edge_orientation(je,jb,i) = edges%tangent_orientation(edge_idx,edge_blk)
            ELSE
              verts%edge_orientation(je,jb,i) = -edges%tangent_orientation(edge_idx,edge_blk)
            ENDIF

          ENDIF ! edge_idx > 0

        ENDDO ! i=1,verts%max_connectivity
        verts%dual_area(je,jb) = REAL(dual_area,wp)

      ENDDO
    ENDDO ! verts owned, needs sync...
!$OMP END PARALLEL DO

    ! sync the vertices
    CALL sync_patch_array(SYNC_V, patch, verts%dual_area(:,:), lacc=.FALSE.)
    DO i=1,verts%max_connectivity
      CALL  sync_patch_array(SYNC_V, patch, verts%edge_orientation(:,:,i), lacc=.FALSE.)
    ENDDO
    !----------------------------------------------------------
    ! Finally,rescale distances by the radius of the sphere
!$OMP PARALLEL PRIVATE(sphere_radius_squared)
    sphere_radius_squared = grid_sphere_radius**2
!$OMP DO PRIVATE(jb, je, start_cell_index, end_cell_index)
    DO jb = cells%all%start_block, cells%all%end_block
      CALL get_index_range(cells%all, jb, start_cell_index, end_cell_index)
      DO je = start_cell_index, end_cell_index
        cells%area(je, jb) = cells%area(je, jb) * sphere_radius_squared
      ENDDO
    ENDDO ! all cells

!$OMP DO PRIVATE(jb, je, start_edge_index, end_edge_index)
    DO jb = edges%all%start_block, edges%all%end_block
      CALL get_index_range(edges%all, jb, start_edge_index, end_edge_index)
      DO je = start_edge_index, end_edge_index

        edges%primal_edge_length(je, jb) = edges%primal_edge_length(je, jb) * grid_sphere_radius
        edges%inv_primal_edge_length(je,jb) = &
          & 1._wp/edges%primal_edge_length(je,jb)
        edges%edge_vert_length(je,jb,1) = 0.5_wp*edges%primal_edge_length(je,jb) ! the center of the edges is in the middle
        edges%edge_vert_length(je,jb,2) = edges%edge_vert_length(je,jb,1)

        edges%edge_cell_length(je,jb,1) = edges%edge_cell_length(je,jb,1) * grid_sphere_radius
        edges%edge_cell_length(je,jb,2) = edges%edge_cell_length(je,jb,2) * grid_sphere_radius
        edges%dual_edge_length(je,jb) = &
            & edges%edge_cell_length(je,jb,1) + edges%edge_cell_length(je,jb,2)
!       IF (edges%dual_edge_length(je,jb) > 0.0_wp) &
          edges%inv_dual_edge_length(je,jb) = &
            & 1._wp/edges%dual_edge_length(je,jb)

        edges%area_edge(je,jb) =  &
          &    edges%primal_edge_length(je,jb)  &
          &  * edges%dual_edge_length(je,jb)

      ENDDO
    ENDDO ! all edges

!$OMP DO PRIVATE(jb,je,i,start_vertex_index,end_vertex_index)
    DO jb = verts%all%start_block, verts%all%end_block
      CALL get_index_range(verts%all, jb, start_vertex_index, end_vertex_index)
      DO je = start_vertex_index, end_vertex_index
        verts%dual_area(je,jb) = verts%dual_area(je,jb) * sphere_radius_squared
      ENDDO
    ENDDO ! all vertices

!$OMP END PARALLEL

    !----------------------------------
    ! write sphere geometry properties
    patch%geometry_info%sphere_radius = grid_sphere_radius

    ! -----------------------------------------------------
    ! the following info needs to be read from the file in order
    ! to maintain binary-identical results on different number of procs
    !
!     minmaxmean = global_minmaxmean(values=edges%primal_edge_length, in_subset=edges%owned)
!     patch%geometry_info%mean_edge_length = minmaxmean(3)
!     minmaxmean = global_minmaxmean(values=edges%dual_edge_length, in_subset=edges%owned)
!     patch%geometry_info%mean_dual_edge_length = minmaxmean(3)
!
!     minmaxmean = global_minmaxmean(values=cells%area, in_subset=cells%owned)
!     patch%geometry_info%mean_cell_area = minmaxmean(3)
!     patch%geometry_info%mean_characteristic_length = SQRT(minmaxmean(3))
!
!     minmaxmean = global_minmaxmean(values=verts%dual_area, in_subset=verts%owned)
!     patch%geometry_info%mean_dual_cell_area = minmaxmean(3)

    ! but since we do not calcultae them, we need to rescale
    patch%geometry_info%mean_edge_length           = patch%geometry_info%mean_edge_length           * grid_length_rescale_factor
    patch%geometry_info%mean_cell_area             = patch%geometry_info%mean_cell_area             * grid_length_rescale_factor**2
    patch%geometry_info%mean_dual_edge_length      = patch%geometry_info%mean_dual_edge_length      * grid_length_rescale_factor
    patch%geometry_info%mean_dual_cell_area        = patch%geometry_info%mean_dual_cell_area        * grid_length_rescale_factor**2
    patch%geometry_info%domain_length              = patch%geometry_info%domain_length              * grid_length_rescale_factor
    patch%geometry_info%domain_height              = patch%geometry_info%domain_height              * grid_length_rescale_factor
    patch%geometry_info%sphere_radius              = patch%geometry_info%sphere_radius              * grid_length_rescale_factor
    patch%geometry_info%grid_mean_resolution       = patch%geometry_info%grid_mean_resolution       * grid_length_rescale_factor
    patch%geometry_info%mean_characteristic_length = patch%geometry_info%grid_mean_resolution
    ! syn all (for testing)
!     CALL sync_patch_geometry(patch)
    !----------------------------------
    ! some debug info
!     CALL dbg_print('c.geo.:edge length',edges%primal_edge_length,"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo.:edge dual length',edges%dual_edge_length,"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo.:edges abs orientation', abs(edges%tangent_orientation(:,:)),"",1, in_subset=edges%all)
!     CALL dbg_print('c.geo.:cell area',cells%area,"",1, in_subset=cells%all)
!     CALL dbg_print('c.geo.:dual area',verts%dual_area,"",1, in_subset=verts%all)
!     CALL dbg_print('c.geo.:vert abs orientation',abs(verts%edge_orientation(:,:,1)),"",1, in_subset=verts%all)

  END SUBROUTINE compute_sphere_geometry_triangles
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE sync_patch_geometry(patch)
    TYPE(t_patch), INTENT(inout), TARGET ::  patch  ! patch data structure

    INTEGER :: i

    CALL sync_patch_array(SYNC_C, patch, patch%cells%cartesian_center(:,:)%x(1), lacc=.FALSE.) ! for testing
    CALL sync_patch_array(SYNC_C, patch, patch%cells%cartesian_center(:,:)%x(2), lacc=.FALSE.) ! for testing
    CALL sync_patch_array(SYNC_C, patch, patch%cells%cartesian_center(:,:)%x(3), lacc=.FALSE.) ! for testing
    CALL sync_patch_array(SYNC_C, patch, patch%cells%edge_orientation(:,:,1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_C, patch, patch%cells%edge_orientation(:,:,2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_C, patch, patch%cells%edge_orientation(:,:,3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_C, patch, patch%cells%area(:,:), lacc=.FALSE.) ! for testing

    CALL sync_patch_array(SYNC_E, patch, patch%edges%cartesian_center(:,:)%x(1), lacc=.FALSE.)! for testing
    CALL sync_patch_array(SYNC_E, patch, patch%edges%cartesian_center(:,:)%x(2), lacc=.FALSE.)! for testing
    CALL sync_patch_array(SYNC_E, patch, patch%edges%cartesian_center(:,:)%x(3), lacc=.FALSE.)! for testing
    CALL sync_patch_array(SYNC_E, patch, patch%edges%edge_cell_length(:,:,1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%edge_cell_length(:,:,2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%cartesian_dual_middle(:,:)%x(1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%cartesian_dual_middle(:,:)%x(2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%cartesian_dual_middle(:,:)%x(3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%primal_cart_normal(:,:)%x(1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%primal_cart_normal(:,:)%x(2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%primal_cart_normal(:,:)%x(3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%dual_cart_normal(:,:)%x(1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%dual_cart_normal(:,:)%x(2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%dual_cart_normal(:,:)%x(3), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%dual_normal(:,:)%v1, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%dual_normal(:,:)%v2, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%primal_normal(:,:)%v1, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%primal_normal(:,:)%v2, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%tangent_orientation(:,:), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%center(:,:)%lon, lacc=.FALSE.) ! for testing
    CALL sync_patch_array(SYNC_E, patch, patch%edges%center(:,:)%lat, lacc=.FALSE.) ! for tesing

    ! some syncs for testestin: Remove them after
    CALL sync_patch_array(SYNC_E, patch, patch%edges%primal_edge_length, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%inv_primal_edge_length, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%edge_vert_length(:,:,1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%edge_vert_length(:,:,2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%edge_cell_length(:,:,1), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%edge_cell_length(:,:,2), lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%dual_edge_length, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%inv_dual_edge_length, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E, patch, patch%edges%area_edge, lacc=.FALSE.)

    CALL sync_patch_array(SYNC_V, patch, patch%verts%cartesian(:,:)%x(1), lacc=.FALSE.) ! for tesing
    CALL sync_patch_array(SYNC_V, patch, patch%verts%cartesian(:,:)%x(2), lacc=.FALSE.) ! for tesing
    CALL sync_patch_array(SYNC_V, patch, patch%verts%cartesian(:,:)%x(3), lacc=.FALSE.) ! for tesing
    CALL sync_patch_array(SYNC_V, patch, patch%verts%dual_area, lacc=.FALSE.)
    DO i=1,patch%verts%max_connectivity
      CALL sync_patch_array(SYNC_V, patch, patch%verts%edge_orientation(:,:,i), lacc=.FALSE.)
    ENDDO
    !----------------------------------
  END SUBROUTINE sync_patch_geometry
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints x0,x1 on a unit sphere
  REAL(wqp) ELEMENTAL FUNCTION arc_length_on_unitsphere (x0, x1)

    TYPE(t_cartesian_coordinates_qp), INTENT(in) :: x0, x1
    REAL(wqp) :: z_cc


    z_cc = DOT_PRODUCT(x0%x, x1%x)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    arc_length_on_unitsphere = ACOS(z_cc)

  END FUNCTION arc_length_on_unitsphere
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL TYPE(t_cartesian_coordinates_qp) FUNCTION circum_center (v0, v1, v2)
   !> the coordinates of the three triangle vertices (unit vectors)
    !! in counter clockwise order.
    TYPE(t_cartesian_coordinates_qp), INTENT(in) :: v0, v1, v2

    !  Local variables:
    TYPE(t_cartesian_coordinates_qp) :: e1, e2 ! edges of the underlying planar triangle:
    ! v1-v0 ands v2-v0, respectively
    TYPE(t_cartesian_coordinates_qp) ::  cu0, cu1, cu2    ! vector product of center:  e1 x e2
    REAL(wqp) :: z_cnorm, cu0_norm, cu1_norm,  cu2_norm
    !-----------------------------------------------------------------------

    e1%x = v1%x - v0%x
    e2%x = v2%x - v0%x
    cu0 = vector_product (e1, e2)
    IF (DOT_PRODUCT(cu0%x, v0%x) < 0.0_wp) THEN
      cu0%x = -cu0%x
    END IF

    e1%x = v2%x - v1%x
    e2%x = v0%x - v1%x
    cu1 = vector_product (e1, e2)
    IF (DOT_PRODUCT(cu1%x, v1%x) < 0.0_wp) THEN
      cu1%x = -cu1%x
    END IF

    e1%x = v0%x - v2%x
    e2%x = v1%x - v2%x
    cu2 = vector_product (e1, e2)
    IF (DOT_PRODUCT(cu2%x, v2%x) < 0.0_wp) THEN
      cu2%x = -cu2%x
    END IF

    cu0_norm = d_norma_3d(cu0)
    cu1_norm = d_norma_3d(cu1)
    cu2_norm = d_norma_3d(cu2)

    IF (cu0_norm > cu1_norm) THEN
      IF (cu0_norm > cu2_norm) THEN
        circum_center%x = cu0%x / cu0_norm
      ELSE
        circum_center%x = cu2%x / cu2_norm
      ENDIF
    ELSE
      IF (cu1_norm > cu2_norm) THEN
        circum_center%x = cu1%x / cu1_norm
      ELSE
        circum_center%x = cu2%x / cu2_norm
      ENDIF
    ENDIF

  END FUNCTION circum_center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION sphere_cartesian_midpoint (p1, p2) result(mid_point)
    TYPE(t_cartesian_coordinates_qp), INTENT(in) :: p1, p2
    TYPE(t_cartesian_coordinates_qp) :: mid_point

      mid_point%x = (p1%x + p2%x)
      mid_point%x = mid_point%x / &
        & d_norma_3d(mid_point)

  END FUNCTION sphere_cartesian_midpoint
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!  Computes area of a spherical triangle based on L'Huilier's Theorem
  !\
  !  /*
  !   * source code is taken from code by Robert Oehmke of Earth System Modeling
  !   * Framework (www.earthsystemmodeling.org) and adjusted to CDO data structures
  !   *
  !   * the license statement for this routine is as follows:
  !   * Earth System Modeling Framework
  !   * Copyright 2002-2013, University Corporation for Atmospheric Research,
  !   * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
  !   * Laboratory, University of Michigan, National Centers for Environmental
  !   * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
  !   * NASA Goddard Space Flight Center.
  !   * Licensed under the University of Illinois-NCSA License.
  !   */
  ELEMENTAL FUNCTION sphere_triangle_area_Huilier (x0, x1, x2) result(area)

    TYPE(t_cartesian_coordinates_qp), INTENT(in) :: x0, x1, x2

    REAL(wqp) :: area
    REAL(wqp) :: a, b, c, s, t

    TYPE(t_cartesian_coordinates_qp) :: productVector

    productVector = vector_product (x0, x1)
    a = ASIN(d_norma_3d(productVector))

    productVector = vector_product (x1, x2)
    b = ASIN(d_norma_3d(productVector))

    productVector = vector_product (x2, x0)
    c = ASIN(d_norma_3d(productVector))

    s = 0.5_wp * (a+b+c)

    t = TAN(s*0.5_wp) * TAN((s - a)*0.5_wp) * TAN((s - b)*0.5_wp) * TAN((s - c)*0.5_wp)
    area = ABS(4.0_wp * ATAN(SQRT(ABS(t))))

  END FUNCTION sphere_triangle_area_Huilier
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION vector_product (v0, v1) result(v2)

    TYPE(t_cartesian_coordinates_qp), INTENT(in) :: v0, v1

    TYPE(t_cartesian_coordinates_qp) :: v2

    !-----------------------------------------------------------------------
    v2%x(1) = v0%x(2)*v1%x(3) - v0%x(3)*v1%x(2)
    v2%x(2) = v0%x(3)*v1%x(1) - v0%x(1)*v1%x(3)
    v2%x(3) = v0%x(1)*v1%x(2) - v0%x(2)*v1%x(1)

  END FUNCTION vector_product
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION cc2gc(p_x) result (p_pos)

    TYPE(t_cartesian_coordinates_qp), INTENT(in) :: p_x          ! Cart. coordinates
    TYPE(t_geographical_coordinates)          :: p_pos        ! geo. coordinates

    REAL (wqp), PARAMETER ::  dbl_eps   = EPSILON(1._wp)
    REAL(wqp)                                :: z_x, z_y, z_z, z_r
    !-----------------------------------------------------------------------

    z_x = p_x%x(1)
    z_y = p_x%x(2)
    z_z = p_x%x(3)

    z_r = z_x * z_x + z_y * z_y
    z_r = SQRT(z_r)

    IF (ABS(z_r) < dbl_eps) THEN    ! one of the poles

      IF (z_z > 0.0_wp) THEN
        p_pos%lat = pi_2
      ELSE
        p_pos%lat = -1._wp * pi_2
      END IF
      p_pos%lon = 0._wp

    ELSE

      p_pos%lat = REAL(ATAN2 ( ABS(z_z), z_r),wp)
      IF (z_z < 0.0_wp) p_pos%lat = -p_pos%lat

      IF (ABS(z_x) < dbl_eps) THEN    ! z_x == 0 ?

        IF (z_y >= 0.0_wp) THEN
          p_pos%lon = pi_2
        ELSE
          p_pos%lon = -1._wp * pi_2
        END IF

      ELSE
        p_pos%lon = REAL(ATAN2( z_y, z_x),wp)
      END IF

    END IF

  END FUNCTION cc2gc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION gc2cc (p_pos)  result(p_x)

    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos     ! geo. coordinates

    TYPE(t_cartesian_coordinates)                :: p_x       ! Cart. coordinates

    REAL (wp)                                  :: z_cln, z_sln, z_clt, z_slt

    !-----------------------------------------------------------------------
    z_sln = SIN(p_pos%lon)
    z_cln = COS(p_pos%lon)
    z_slt = SIN(p_pos%lat)
    z_clt = COS(p_pos%lat)

    p_x%x(1) = z_cln*z_clt
    p_x%x(2) = z_sln*z_clt
    p_x%x(3) = z_slt

  END FUNCTION gc2cc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
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
  !! Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
  !! Currenlty noyt in use, but kept it in comments for future use
!   ELEMENTAL SUBROUTINE gvec2cvec_wqp (p_gu, p_gv, p_long, p_lat, p_cu, p_cv, p_cw)
!     !
!     REAL(wp), INTENT(in)  :: p_gu, p_gv     ! zonal, meridional vec. component
!     REAL(wp), INTENT(in)  :: p_long, p_lat  ! geo. coord. of data point
!     REAL(wp), INTENT(out) :: p_cu, p_cv, p_cw            ! Cart. vector
!
!     REAL(wqp)              :: p_gu_qp, p_gv_qp, p_long_qp, p_lat_qp
!     REAL(wqp)              :: z_cln, z_sln, z_clt, z_slt  ! sin and cos of
!     !-------------------------------------------------------------------------
!
!     p_gu_qp = p_gu
!     p_gv_qp = p_gv
!     p_long_qp = p_long
!     p_lat_qp  = p_lat
!
!
!     z_sln = SIN(p_long_qp)
!     z_cln = COS(p_long_qp)
!     z_slt = SIN(p_lat_qp)
!     z_clt = COS(p_lat_qp)
!
!     p_cu = -REAL(z_sln * p_gu_qp + z_slt * z_cln * p_gv_qp, wp)
!     ! p_cu = -1._wqp * p_cu
!     p_cv = REAL(z_cln * p_gu_qp - z_slt * z_sln * p_gv_qp, wp)
!     p_cw = REAL(z_clt * p_gv_qp, wp)
!
!   END SUBROUTINE gvec2cvec_wqp
  !-------------------------------------------------------------------------


END MODULE mo_grid_geometry_sphere
