!> Contains the memory class for the surface turbulence process.
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
MODULE mo_turb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_lct_class,    ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE, GLACIER_TYPE !, LAKE_TYPE
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_turb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 15

  TYPE, EXTENDS(t_jsb_memory) :: t_turb_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      & rough_m,              & !< Surface roughness length for momentum
      & rough_h,              & !< Surface roughness length for heat
      ! The following two variables are transformations of the corresponding roughness length variables
      ! and are used in the aggregation of roughness length:
      !   z0_star = 1 / (ln(blending_height/z0))^2
      ! For the average over surface types, we have:
      !   z0_star_avg = sum(f_i * z0_star_i)
      ! where _i indicates the different surface types
      & rough_m_star,         &
      & rough_h_star,         &
      ! For everything but lakes
      & fact_q_air,           & !< Factor for lowest atm. level in formulation of surface fluxes
                                !! of specific humidity and dry static energy to account for
                                !! evapotranspiration                                              []
      & fact_qsat_srf,        & !< Similar factor for surface values of specific humidity
                                !! and dry static energy
      & fact_qsat_trans_srf     !< Similar factor for transpiration

CONTAINS
    PROCEDURE :: Init => Init_turb_memory
  END TYPE t_turb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_memory_class'

CONTAINS

  SUBROUTINE Init_turb_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,         ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid

    CLASS(t_turb_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: lib_id     !< id of primary lct in lctlib
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_turb_memory'

    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used
    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL mem%Add_var( 'rough_m', mem%rough_m,                                            &
      & hgrid, surface,                                                                  &
      & t_cf('roughness_length_momentum', 'm', 'surface roughness length for momentum'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & output_level=FULL,                                                               &
      & initval_r=1.0_wp )

    CALL mem%Add_var( 'rough_h', mem%rough_h,                                    &
      & hgrid, surface,                                                          &
      & t_cf('roughness_length_heat', 'm', 'surface roughness length for heat'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
      & prefix, suffix,                                                          &
      & lrestart=.TRUE.,                                                         &
      & output_level=FULL,                                                       &
      & initval_r=1.0_wp )

    CALL mem%Add_var( 'rough_m_star', mem%rough_m_star,                          &
      & hgrid, surface,                                                          &
      & t_cf('rough_m_star', '', 'rough_m_star'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
      & prefix, suffix,                                                          &
      & lrestart=.FALSE.,                                                        &
      & loutput=.FALSE.,                                                         &
      & initval_r=1.0_wp )

    CALL mem%Add_var( 'rough_h_star', mem%rough_h_star,                          &
      & hgrid, surface,                                                          &
      & t_cf('rough_h_star', '', 'rough_h_star'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
      & prefix, suffix,                                                          &
      & lrestart=.FALSE.,                                                        &
      & loutput=.FALSE.,                                                         &
      & initval_r=1.0_wp )

    CALL mem%Add_var( 'fact_q_air', mem%fact_q_air,                        &
      & hgrid, surface,                                                    &
      & t_cf('fact_q_air', '', ''),                                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & output_level=FULL,                                                 &
      & initval_r=0.5_wp )

    CALL mem%Add_var( 'fact_qsat_srf', mem%fact_qsat_srf,                  &
      & hgrid, surface,                                                    &
      & t_cf('fact_qsat_srf', '', ''),                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & output_level=FULL,                                                 &
      & initval_r=0.5_wp )

    CALL mem%Add_var( 'fact_qsat_trans_srf', mem%fact_qsat_trans_srf,      &
      & hgrid, surface,                                                    &
      & t_cf('fact_qsat_trans_srf', '', ''),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & output_level=FULL,                                                 &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_turb_memory

#endif
END MODULE mo_turb_memory_class
