!> Contains the memory class for phenology.
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
MODULE mo_pheno_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of processes in this module
  dsl4jsb_Use_processes HYDRO_, RAD_
  ! Use process configurations
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(RAD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  ! ======================================================================================================= !
  !>Type definition for pheno memory
  !>
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_pheno_memory

    TYPE(t_jsb_var_real2d) :: &
      & lai                               !< actual leaf area index [m2 m-2]
    TYPE(t_jsb_var_real3d) :: &
      & lai_cl                            !< leaf area index per canopy layer [m2 m-2]

    ! variables NOT added to varlist
    REAL(wp), POINTER :: lai_mon(:,:,:)       => NULL()   !< Monthly leaf area index (per year or climatology) <nproma x nblks x 14> !
    REAL(wp), POINTER :: fract_fpc_mon(:,:,:) => NULL()   !< Monthly fpc (per year or climatology) <nproma x nblks x 14> !

    ! Diagnostic land mean phenology e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) ::           &
      & lai_ta_gmean,                   & !< Global Land mean LAI (box tile)
      & fract_fpc_gmean                   !< Global Land mean foliage cover (box tile)

    TYPE(t_jsb_var_real2d)           :: &
      & lai_ta,                         & !< Leaf area index relative to the tile area
      & fract_fpc,                      & !< Foliage projected cover of the pft
                                          !  (fraction of ground area occupied by the vertical projection of foliage)
      & fract_fpc_max,                  & !< Maximum fraction of tile that can be covered by leaves (assuming infinite lai)
      & fract_forest,                   & !< Forest fraction
      & maxLAI_allom,                   & !< dynamical maximum LAI [m^2(leaf)/m^2(ground)]:
                                          !< Current max LAI according to allometric relationships --> only changes ones a year!

      & veg_fract_correction              !< Factor representing clumpyness of vegetation (for each tile type).
                                          !< To write it in the output file we define it as 2d even when there is no change
                                          !< between gridcells as it is in fact a 1d constant for each tile.
                                          !   JN: with l_forestRegrowth = true veg_fract_correction can change over time
                                          !       and is spatially variable
                                          !       if l_forestRegrowth = false there is still no change between gridcells
                                          !       and it is a 1d constant for each tile in this case.
                                          ! Gives the fraction of canopy area from the tile area. Between 0 and 1.
                                          ! canopy area = veg_fract_correction * tile area.
                                          ! var (tile area) = veg_fract_correction * var(canopy area).
                                          ! ACHTUNG: da wir vorerst fract_fpc_max zur Abgrenzung zwischen "Wueste" und Vegetation
                                          !          in den JS4 Berechnungen haben, gilt eigentlich dies hier:
                                          ! var (tile area) = veg_fract_correction * var(canopy area)  * fract_fpc_max.

  CONTAINS
    PROCEDURE :: Init => Init_pheno_memory
  END TYPE t_pheno_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: phenology
  !>
  SUBROUTINE Init_pheno_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_varlist,       ONLY: MINIMUM, BASIC, FULL  ! MEDIUM
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_pheno_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: lib_id     !< id of primary lct in lctlib
    INTEGER,               INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(RAD_)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    INTEGER                     :: table                        ! ...
    TYPE(t_grib2)               :: grib2_desc                   ! used for lai settings
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pheno_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model        => Get_model(model_id)
    table        = tables(1)
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(RAD_)

    IF ( (     One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &  ) ) THEN

      ! IF (TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'climatology') THEN   ! TODO
        ALLOCATE(mem%lai_mon      (hgrid%nproma, hgrid%nblks, 0:13))
        ALLOCATE(mem%fract_fpc_mon(hgrid%nproma, hgrid%nblks, 0:13))
        !$ACC ENTER DATA CREATE(mem%lai_mon, mem%fract_fpc_mon)
        !$ACC ENTER DATA ATTACH(mem%lai_mon, mem%fract_fpc_mon)
      ! END IF

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(2,0,28, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'lai', mem%lai,                                        &
        & hgrid, surface,                                                      &
        & t_cf('leaf_area_index', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), grib2_desc,                          &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lai_ta', mem%lai_ta,                                  &
        & hgrid, surface,                                                      &
        & t_cf('leaf_area_index', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=MINIMUM,                                                &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (TRIM(mem%owner_tile_name) == 'veg') THEN
        grib2_desc = t_grib2(2,0,228, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'fract_fpc', mem%fract_fpc,                            &
        & hgrid, surface,                                                      &
        & t_cf('foliage_projected_cover', '', ''),                             &
        & t_grib1(table, 255, grib_bits), grib2_desc,                          &
        & prefix, suffix,                                                      &
        & output_level=MINIMUM,                                                &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_fpc_max', mem%fract_fpc_max,                    &
        & hgrid, surface,                                                      &
        & t_cf('maximum_fpc', '', ''),                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      IF (((dsl4jsb_Config(HYDRO_)%l_organic .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_socmap)) .OR. &
        &   dsl4jsb_Config(RAD_)%use_alb_veg_simple) THEN
        CALL mem%Add_var( 'fract_forest', mem%fract_forest,                      &
          & hgrid, surface,                                                      &
          & t_cf('forest_fraction', '', ''),                                     &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & initval_r=0.0_wp )
      END IF

      CALL mem%Add_var( 'veg_fract_correction', mem%veg_fract_correction,      &
        & hgrid, surface,                                                      &
        & t_cf('Correction factor to convert tile to canopy area', '', '') ,   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=FULL,                                                   &
        & lrestart=.TRUE.,                                                     &
        & initval_r=0.0_wp )

      ! Diagnostic global carbon sums for experiment monitoring
      !   (1d stream variables are not supported with echam)
      IF (TRIM(suffix) == 'box') THEN

        CALL mem%Add_var( 'lai_ta_gmean', mem%lai_ta_gmean,                           &
          & hgrid, surface,                                                           &
          & t_cf('lai_ta_gmean', '1', 'Global land mean leaf_area_index'),            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
          & prefix, suffix,                                                           &
          & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )

        CALL mem%Add_var( 'fract_fpc_gmean', mem%fract_fpc_gmean,                     &
          & hgrid, surface,                                                           &
          & t_cf('fract_fpc_gmean', '1', 'Global land mean foilage projected cover'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
          & prefix, suffix,                                                           &
          & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )

      END IF
    END IF !  VEG / LAND

  END SUBROUTINE Init_pheno_memory

#endif
END MODULE mo_pheno_memory_class
