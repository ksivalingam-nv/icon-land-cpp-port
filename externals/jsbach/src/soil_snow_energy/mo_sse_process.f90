!> Contains the process routines of the soil and snow energy process.
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
!> The module provides the process routines for the soil and snow energy scheme. It
!> contains subroutines for computing thermal properties of soil and snow, as well as
!> updating snow and soil temperature including the computation of soil moisture
!> phase change and the active layer depth.
!>
MODULE mo_sse_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control,            ONLY: acc_stream
  USE mo_jsb_physical_constants, ONLY: rhoh2o
  USE mo_exception, ONLY: message, finish
  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: relative_humidity_soil, init_soil_temperature,                       &
    & calc_vol_heat_capacity, calc_thermal_conductivity, add_heat_capacity_pond, &
    & calc_snow_abcoeff, calc_soil_temperature, calc_snow_temperature, Get_liquid_max

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_process' !< Module name

CONTAINS

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Calculate the volumetric heat capacity of a soil layer.
  !>
  !> The formulation follows de Vries (1963) and considers the heat storage of
  !> mineral and organic solids as well as liquid water, ice, and the gas phase in
  !> the pore space. If the layer extends to bedrock, the result is mixed with the
  !> bedrock heat capacity to represent the soil-bedrock composite.
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION calc_vol_heat_capacity( &
    & vol_heat_cap_soil_dry_mineral,              &
    & vol_heat_cap_soil_dry_organic,              &
    & porosity,                                   &
    & fract_organic,                              &
    & fract_water,                                &
    & fract_ice,                                  &
    & layer_depth,                                &
    & dz                                          &
    & ) RESULT(vol_heat_cap)

    !$ACC ROUTINE SEQ

    USE mo_sse_constants, ONLY: vol_hcap_water, vol_hcap_ice, vol_hcap_air, vol_hcap_bedrock

    REAL(wp), INTENT(in) :: porosity                      !< Soil porosity []
    REAL(wp), INTENT(in) :: vol_heat_cap_soil_dry_mineral !< Volumetric heat capacity of dry mineral soil [J m-3 K-1]
    REAL(wp), INTENT(in) :: vol_heat_cap_soil_dry_organic !< Volumetric heat capacity of dry organic soil [J m-3 K-1]
    REAL(wp), INTENT(in) :: fract_organic                 !< Organic fraction of soil solids []
    REAL(wp), INTENT(in) :: fract_water                   !< Volumetric water fraction []
    REAL(wp), INTENT(in) :: fract_ice                     !< Volumetric ice fraction []
    REAL(wp), INTENT(in) :: layer_depth                   !< Soil depth above bedrock [m]
    REAL(wp), INTENT(in) :: dz                            !< Soil layer thickness [m]
    REAL(wp)             :: vol_heat_cap                  !< Volumetric heat capacity of the soil layer [J m-3 K-1]

    IF (layer_depth > 0._wp) THEN
      ! above bedrock compute the fraction weighted volumetric heat capacity based on the individual contributions
      vol_heat_cap = (1._wp - porosity)                                            & ! vol. fraction of solid (mineral + organic)
        &              * ( (1._wp - fract_organic) * vol_heat_cap_soil_dry_mineral & ! heat cap of mineral soil
        &                  +        fract_organic  * vol_heat_cap_soil_dry_organic & ! heat cap of organic soil
        &                )                                                         &
        &            + fract_water * vol_hcap_water                                & ! heat cap of liquid water
        &            + fract_ice   * vol_hcap_ice                                  & ! heat cap of frozen water
        ! @todo : porosity - fract_water - fract_ice should not become lower than zero.
        &            + MAX((porosity - fract_water - fract_ice),0._wp) * vol_hcap_air  ! heat cap of remaining free pores
      vol_heat_cap = ( vol_heat_cap * layer_depth + vol_hcap_bedrock * (dz - layer_depth) ) / dz ! soil and bedrock mixture
    ELSE
      ! Below bedrock, the volumetric heat capacity is set to the bedrock value as no water, ice or organic matter is present
      vol_heat_cap = vol_hcap_bedrock
    END IF

  END FUNCTION calc_vol_heat_capacity

! --------------------------------------------------------------------------------------------------------- !
!>
!> Include surface water heat capacity in top soil layer heat capacity.
!>
!> The routine calculates an effective volumetric heat capacity of the top soil layer
!> that includes the additional heat capacity of the surface water and ice in ponds.
!> Since ponds do not have a dedicated energy content in this model, their effect on
!> surface and soil temperatures is otherwise not accounted for.
!>
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION  add_heat_capacity_pond( &
    & dz,                                 &
    & wtr_pond,                           &
    & ice_pond,                           &
    & vol_heat_cap                        &
    & ) RESULT(vol_heat_cap_eff)

    !$ACC ROUTINE SEQ

    USE mo_sse_constants, ONLY: vol_hcap_water, vol_hcap_ice
    USE mo_jsb_physical_constants, ONLY: rhoh2o, rhoi

    REAL(wp), INTENT(in)  :: dz,          &   !< Layer thickness [m]
      &                      wtr_pond,    &   !< Water content of surface reservoir [m]
      &                      ice_pond,    &   !< Ice content of surface reservoir [m]
      &                      vol_heat_cap     !< Volumetric heat capacity of top soil layer [J m-3 K-1]
    REAL(wp)              :: vol_heat_cap_eff !< Effective volumetric heat capacity of top soil layer [J m-3 K-1]

    ! We approximate the potential energy storage by calculating an effective heat
    ! capacity of the top soil layer that includes the additional heat capacity of the surface reservoir.

    vol_heat_cap_eff = vol_heat_cap &
      &              + (wtr_pond * vol_hcap_water + ice_pond * rhoh2o / rhoi * vol_hcap_ice) / dz

  END FUNCTION add_heat_capacity_pond

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Calculate the soil thermal conductivity following Johansen (1977).
  !>
  !> The routine accounts for the dry conductivity of the soil matrix, the effect of
  !> liquid water and ice, and the soil-bedrock mixture when the layer reaches the
  !> bedrock boundary.
  !>
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION calc_thermal_conductivity( &
    & heat_cond_soil_mineral,                        &
    & water, ice,                                    &
    & porosity,                                      &
    & fract_organic,                                 &
    & vol_porosity_org, hcond_org, hcond_dry_org,    &
    & layer_depth,                                   &
    & dz                                             &
    & ) RESULT(heat_cond)

    USE mo_jsb_physical_constants, ONLY: dens_soil
    USE mo_sse_constants,          ONLY: hcond_water, hcond_ice, hcond_bedrock


    REAL(wp), INTENT(in) :: heat_cond_soil_mineral !< Thermal conductivity of dry mineral soil [J m-1 s-1 K-1]
    REAL(wp), INTENT(in) :: water                  !< Liquid water content [m]
    REAL(wp), INTENT(in) :: ice                    !< Ice content [m]
    REAL(wp), INTENT(in) :: porosity               !< Soil porosity []
    REAL(wp), INTENT(in) :: fract_organic          !< Organic fraction of soil solids []
    REAL(wp), INTENT(in) :: vol_porosity_org       !< Porosity of organic soil []
    REAL(wp), INTENT(in) :: hcond_org              !< Thermal conductivity of organic soil [J m-1 s-1 K-1]
    REAL(wp), INTENT(in) :: hcond_dry_org          !< Dry thermal conductivity of organic soil [J m-1 s-1 K-1]
    REAL(wp), INTENT(in) :: layer_depth            !< Soil depth above bedrock [m]
    REAL(wp), INTENT(in) :: dz                     !< Soil layer thickness [m]
    REAL(wp)             :: heat_cond              !< Effective thermal conductivity of the soil layer [J m-1 s-1 K-1]

    REAL(wp) ::           &
      & dens_bulk,        & !< Bulk density of dry mineral soil (including pores) [kg m-3]
      & porosity_mineral, & !< Volumetric porosity of mineral soil part [ ]
      & hcond_dry_min,    & !< Dry thermal conductivity soil (including mineral soil and pores) [J m-1 s-1 K-1]
      & hcond_dry,        & !< Dry thermal conductivity of soil (including mineral and organic soil and pores) [J m-1 s-1 K-1]
      & sat_fact,         & !< Scaling factor to account for empty pore space [ ]
      & saturation,       & !< Degree of saturation [ ]
      & kersten,          & !< Kersten number
      & hcond_soil,       & !< Thermal conductivity of soil solids (mineral and organic) [J m-1 s-1 K-1]
      & hcond_sat,        & !< Thermal conductivity for saturated soil [J m-1 s-1 K-1]
      & fract_water,      & !< Effective fractional water volume for saturated soil computation [ ]
      & fract_ice           !< Effective fractional ice volume for saturated soil computation [ ]

    !$ACC ROUTINE SEQ

    IF (layer_depth > 0._wp) THEN
      ! Above bedrock, compute thermal conductivity of dry mineral and organic soil first
      porosity_mineral = (porosity - fract_organic * vol_porosity_org) / (1._wp - fract_organic)
      dens_bulk        = dens_soil * (1._wp - porosity_mineral)
      hcond_dry_min    = (0.135_wp * dens_bulk + 64.7_wp) / (dens_soil - 0.947_wp * dens_bulk)
      hcond_dry        = (1._wp - fract_organic) * hcond_dry_min + fract_organic * hcond_dry_org

      IF (water + ice <= EPSILON(1._wp)) THEN
        ! No water or ice present, the thermal conductivity is set to the dry value
        heat_cond = hcond_dry
      ELSE
        ! Thermal conductivity of soil solids
        hcond_soil = (1._wp - fract_organic) * heat_cond_soil_mineral + fract_organic * hcond_org

        ! Thermal conductivity of saturated soil
        saturation = (water + ice) / layer_depth / porosity  ! fractional volume of water/ice rel. to pore volume

        IF (ice > 0._wp) THEN  ! Frozen case
          kersten = saturation
        ELSE                   ! Unfrozen case
          ! This is the case for fine unfrozen soil
          IF (saturation > 0.1_wp) THEN
            kersten = LOG10(saturation) + 1._wp
          ELSE
            kersten = 0._wp
          END IF
        END IF
        ! The scaling is changed compared to earlier versions. The current implementation interpretes the fractions
        ! of the individual components such that fract_water + fract_ice equals porosity, which is also supported by
        ! Nikoosokhan et al (2016), equation 4
        sat_fact    = MAX(0._wp, porosity / ((water + ice) / layer_depth)) ! Scaling factor to upscale water/ice to saturation
        fract_water = (water / layer_depth) * sat_fact ! Water fraction rel. to completely saturated soil
        fract_ice   = (  ice / layer_depth) * sat_fact ! Ice fraction rel. to completely saturated soil
        hcond_sat   = hcond_soil      ** (1._wp-fract_water-fract_ice) &
          &         * hcond_water     ** fract_water                   &
          &         * hcond_ice       ** fract_ice

        heat_cond = (hcond_sat - hcond_dry) * kersten + hcond_dry

        ! Finally, consider soil and bedrock mixture
        heat_cond = ( heat_cond * layer_depth + hcond_bedrock * (dz - layer_depth) ) / dz

        heat_cond = MAX(heat_cond, hcond_dry)
      END IF
    ELSE
      ! Below bedrock, the thermal conductivity is set to the bedrock value as no water, ice or organic matter is present
      heat_cond = hcond_bedrock
    END IF

  END FUNCTION calc_thermal_conductivity

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Calculate a smooth soil-moisture humidity factor between dry and saturated
  !> conditions.
  !>
  !> The function is used to scale soil-water-related processes smoothly as the moisture
  !> state moves between the dry and saturated limits.
  !>
  REAL(wp) FUNCTION relative_humidity_soil(ws, ws_min, ws_max)

    !$ACC ROUTINE SEQ

  USE mo_jsb_math_constants, ONLY: pi

    REAL(wp), INTENT(in) :: &
      & ws,                 & !< Soil (layer) moisture content [m]
      & ws_min,             & !< Minimum soil (layer) moisture content [m]
      & ws_max                !< Maximum soil (layer) moisture content [m]

    REAL(wp) :: moisture      !< Moisture content limited to the configured maximum [m]

    moisture = MIN(ws, ws_max)
    IF (moisture > 0._wp) THEN
      ! Soil moisture values below ws_min result in a rel. humidity of 0, for values above ws_max the
      ! rel. humidity is 1. Between ws_min and ws_max rel. humidity increases using the Cosine function
      ! resulting in a smooth transitions at ws_min and was_max.
      relative_humidity_soil = 0.5_wp * (1._wp - COS(MAX(moisture - ws_min, 0._wp) / (ws_max - ws_min) * pi))
    ELSE
      relative_humidity_soil = 0._wp
    END IF

  END FUNCTION relative_humidity_soil

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Initialize soil temperatures from climatological surface temperatures.
  !>
  !> The routine derives an annual-cycle profile for the soil column from the
  !> prescribed surface temperature climatology and the soil layer depths.
  !!
  SUBROUTINE init_soil_temperature(cmid, tclim, tsoil)
    ! Method:
    !
    ! Starting from the tslclim field temperatures are set in relation
    ! to the depth of the soil layer and position of the initial
    ! day in the annual cycle.
    !
    ! tsl is at 0.07 m
    ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m (s. soiltemp)
    !
    USE mo_jsb_time,       ONLY: t_datetime, get_time_start, get_year, deallocateDatetime, &
                                 & get_year_length, get_month_length, get_day_length, get_year_day
    USE mo_jsb_math_constants, ONLY: pi

    REAL(wp), INTENT(in)  :: cmid(:)     !< Depths of soil layer mid points [m]
    REAL(wp), INTENT(in)  :: tclim(:,:)  !< Climatological monthly surface temperature [K]
    REAL(wp), INTENT(out) :: tsoil(:,:)  !< Soil temperatures [K]

    INTEGER                   :: nsoil, ndim !< Number of soil layers and horizontal grid cells
    REAL(wp)                  :: zkap, zsqrt !< Soil thermal diffusivity and characteristic damping depth
    REAL(wp), ALLOCATABLE     :: zrange(:), zdmax(:) !< Annual surface-temperature range and day of maximum temperature
    INTEGER, ALLOCATABLE      :: jmax(:) !< Month of maximum surface temperature
    TYPE(t_datetime), POINTER :: inidate !< Start date of the experiment
    INTEGER                   :: iniyear !< Start year of the experiment
    INTEGER                   :: im, i, ic !< Month, soil-layer, and horizontal-cell indices
    REAL(wp)                  :: year_len, day_in_year !< Length of the year and current day of year
    REAL(wp)                  :: month_len(12), month_mid(12) !< Month lengths and month midpoints [days]

    REAL(wp), PARAMETER :: cmid_offset = 0.07_wp !< Reference depth of the climatological surface temperature [m]

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_temperature'

    ndim  = SIZE(tsoil,1)
    nsoil = SIZE(tsoil,2)

    ALLOCATE(zrange(ndim), jmax(ndim), zdmax(ndim))

    zkap = 7.5e-7_wp

    ! Compute dampening depth based on year and day length
    inidate => get_time_start()
    iniyear = get_year(inidate)

    year_len    = REAL(get_year_length(iniyear), wp)
    day_in_year = INT(get_year_day(inidate))
    CALL deallocateDatetime(inidate)

    DO im=1,12
      month_len(im) = REAL(get_month_length(iniyear, im), wp)
      month_mid(im) = month_len(im) * 0.5_wp
    END DO

    zsqrt = SQRT(zkap * year_len * get_day_length() / pi)

    ! Calendar month of maximum surface temperature
    jmax(:) = MAXLOC(tclim(:,:), DIM=2)

    ! Difference between months of maximum and minimum surface temperature
    zrange(:) = MAXVAL(tclim(:,:), DIM=2) - MINVAL(tclim(:,:), DIM=2)

    ! Calendar day of maximum surface temperature
    zdmax(:) = 0._wp
    DO ic=1,ndim
      DO im=1,jmax(ic)
        zdmax(ic) = zdmax(ic) + month_len(im)
      END DO
      zdmax(ic) = zdmax(ic) - month_mid(jmax(ic))
    END DO

    ! Initialize soil temperatures
    DO i=1,nsoil
      tsoil(:,i) = SUM(tclim(:,:), DIM=2) / 12._wp                            &
                 & + 0.5_wp * zrange(:) * EXP(-(cmid(i)-cmid_offset) / zsqrt) &
                 &   * COS(2._wp * pi * (day_in_year - zdmax(:)) / year_len - (cmid(i)-cmid_offset) / zsqrt)
    END DO

    DEALLOCATE(zrange, jmax, zdmax)

  END SUBROUTINE init_soil_temperature

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Calculate the soil-layer temperatures and the associated implicit coefficients.
  !>
  !> The routine updates the soil temperature profile, the Richtmyer and Morton A/B
  !> coefficients, the ground heat flux, and the ground heat capacity. It also tracks
  !> thaw depth and the freezing or melting of soil water when the corresponding
  !> information is available.
  !!
  SUBROUTINE calc_soil_temperature(       &
    & nc, nsoil, dz, mids, bots,          &
    & delta_time, lstart,                 &
    & l_freeze, l_supercool,              &
    & enforce_water_budget,               &
    & t_soil_top,                         &
    & vol_heat_cap, heat_cond,            &
    & t_soil_sl,                          &
    & t_soil_acoef, t_soil_bcoef,         &
    & hcap_grnd, grnd_hflx,               &
    & thaw_depth, is_snow_covered,        &
    ! optional
    & dsoil, ws_max, matric_pot, bclapp,  &
    & wtr, ice, le_pc_remain,             &
    & liquid_max, frozen, melt            &
    & )

    USE mo_jsb_physical_constants, ONLY: tmelt, alf, rhoh2o, rhoi
    USE mo_jsb_impl_constants,     ONLY: WB_IGNORE, WB_LOGGING, WB_ERROR


    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    ::  nc                          !< length of the vector
    INTEGER,  INTENT(in)    ::  nsoil                       !< number of soil layers
    REAL(wp), INTENT(in)    ::  delta_time                  !< Time step length [s]
    LOGICAL,  INTENT(in)    ::  lstart                      !< if T, start of experiment
    LOGICAL,  INTENT(in)    ::  l_freeze, l_supercool       !< if T, freeze or supercooling is allowed
    INTEGER,  INTENT(IN)    ::  enforce_water_budget        !< Water balance check setting
    REAL(wp), INTENT(in)    ::  dz(:)                       !< soil layer thickness [m]
    REAL(wp), INTENT(in)    ::  mids(:)                     !< depth of mids of soil layers [m]
    REAL(wp), INTENT(in)    ::  bots(:)                     !< lower boundary of soil layers [m]
    REAL(wp), INTENT(in)    ::  t_soil_top(:)               !< surface temperature at top of soil [K]
    REAL(wp), INTENT(in)    ::  vol_heat_cap(:,:)           !< soil heat capacity [J m-3 K-1]
    REAL(wp), INTENT(in)    ::  heat_cond(:,:)              !< soil thermal conductivity [J m-1 s-1 K-1]
    LOGICAL,  INTENT(in)    ::  is_snow_covered(:)          !< T: soil layer is snow covered

    REAL(wp), INTENT(inout) ::  t_soil_sl(:,:)              !< soil temperature [K]
    REAL(wp), INTENT(inout) ::  t_soil_acoef(:,:)           !< soil A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(inout) ::  t_soil_bcoef(:,:)           !< soil B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(out)   ::  hcap_grnd(:)                !< heat capacity of the ground [J m-2 K-1]
    REAL(wp), INTENT(out)   ::  grnd_hflx(:)                !< ground heat flux [J m-2 s-1]
    REAL(wp), INTENT(out)   ::  thaw_depth(:)               !< Thawing depth [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  dsoil(:,:)        !< soil depth until bedrock [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  ws_max(:,:)       !< Maximum water storage [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  matric_pot(:,:)   !< Soil matric potential [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  bclapp(:,:)       !< Clapp and Hornberger parameter [ ]
    REAL(wp), INTENT(inout), OPTIONAL ::  wtr(:,:)          !< Water storage [m]
    REAL(wp), INTENT(inout), OPTIONAL ::  ice(:,:)          !< Ice storage [m]
    REAL(wp), INTENT(inout), OPTIONAL ::  le_pc_remain(:)   !< Energy available for or required from phase change
                                                            !< in top soil layer [J m-2]
    REAL(wp), INTENT(inout), OPTIONAL ::  liquid_max(:,:)   !< Water that remains liquid at sub-zero temps [m]
    REAL(wp), INTENT(out), OPTIONAL   ::  frozen(:,:)       !< Water flux from freezing soil water [kg m-2 s-1]
    REAL(wp), INTENT(out), OPTIONAL   ::  melt(:,:)         !< Water flux from melting soil ice    [kg m-2 s-1]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: ic, is                                     !< Horizontal-cell and soil-layer indices
    REAL(wp) :: le_phase(nc)                               !< Residual phase-change energy in the top soil layer [J m-2]
    REAL(wp) :: z1(nc,nsoil)                               !< Inverse coefficients for the implicit soil-temperature solver
    REAL(wp) :: zd1(nsoil-1)                               !< Inverse distances between soil-layer midpoints [m-1]
    REAL(wp) :: zdz1(nc,nsoil-1), zdz2(nc,nsoil)           !< Conductive and storage coefficients of soil layers
    REAL(wp) :: heat_cap(nc,nsoil)                         !< Heat capacity of soil layers per horizontal cell [J m-2 K-1]
    REAL(wp) :: wtr_old(nc,nsoil), ice_old(nc,nsoil)       !< Water and ice storage before phase change [m]
    REAL(wp) :: freeze_mass(nc,nsoil), melt_mass(nc,nsoil) !< Mass of frozen and melted water by soil layer [kg m-2]
    CHARACTER(len=4096) :: message_text_long               !< Long string for debug output

    ! Temporary index lists for error detection
    INTEGER :: ie             !< Error cell index
    INTEGER :: n_errors       !< Number of error cells
    INTEGER :: error_idx(nc)  !< Cell indices of error cells
    LOGICAL :: error_flag(nc) !< Flag that given cell index is already in list

    !-----------------------------------------------------------------------------------------------
    !  Precompute soil layer spacing and soil heat capacities for the implicit solver.

    !$ACC DATA CREATE(z1, zd1, zdz1, zdz2, heat_cap, wtr_old, ice_old, freeze_mass, melt_mass, le_phase) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO is = 1,nsoil-1
       zd1(is) = 1._wp / (mids(is+1) - mids(is))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO ic = 1, nc
      DO is = 1, nsoil
        heat_cap(ic,is) = dz(is) * vol_heat_cap(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !-----------------------------------------------------------------------------------------------
    !  Initialize variables related to phase change and thaw depth calculations.
    !-----------------------------------------------------------------------------------------------

    IF (PRESENT(wtr)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO is = 1, nsoil
        DO ic = 1, nc
          frozen    (ic,is)  = 0._wp
          melt      (ic,is)  = 0._wp
          freeze_mass(ic,is) = 0._wp
          melt_mass(ic,is)   = 0._wp
          wtr_old(ic,is)     = wtr(ic,is)
          ice_old(ic,is)     = ice(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      thaw_depth(ic) = 0._wp
      le_phase(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    IF (PRESENT(le_pc_remain)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF (.NOT. is_snow_covered(ic)) THEN
          le_phase(ic) = le_pc_remain(ic)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !-----------------------------------------------------------------------------------------------
    !  Update soil temperatures from the current surface temperature and the A/B
    !  coefficients from the previous time step.
    !-----------------------------------------------------------------------------------------------

    IF (.NOT. lstart) THEN

      ! set uppermost layer to surface temperature
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        t_soil_sl(ic,1) = t_soil_top(ic)
      END DO
      !$ACC END PARALLEL LOOP

      ! compute deeper layers for upper layer temperature and A/B coefficients from previous time step
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP SEQ
      DO is = 1, nsoil-1
        !$ACC LOOP GANG VECTOR
        DO ic=1,nc
          t_soil_sl(ic,is+1) = t_soil_acoef(ic,is) + t_soil_bcoef(ic,is) * t_soil_sl(ic,is)
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

      !-----------------------------------------------------------------------------------------------
      !  Compute the amount of freezing or melting water/ice in the soil layers based on its
      !  current temperature and the available energy for phase change (top layer soil).
      !-----------------------------------------------------------------------------------------------

      IF (l_freeze .AND. PRESENT(wtr)) THEN

        ! Note that there is an unavoidable inconsistency between the phase change energy computed
        ! in SEB and its use here, because in the mean time the soil heat capacity might have changed
        ! due to hydrological processes.

        ! Temporarily adding residual energy from phase change (from SEB) to top layer soil temperature.
        ! This allows using the soil temperature difference to tmelt to compute the actual
        ! freezing/melting water/ice mass.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          IF (.NOT. is_snow_covered(ic)) THEN
            t_soil_sl(ic,1) = t_soil_sl(ic,1) + le_phase(ic) / heat_cap(ic,1)
          END IF
        END DO
        !$ACC END PARALLEL LOOP

        ! Compute mass of freezing water or melting ice limited by required/available
        ! energy (only for snow-free top layer soil) and by the temperature difference
        ! to tmelt, and available water/ice.
        ! Note: only either freeze_mass or melt_mass is non-zero at any given point.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is=1,nsoil
          DO ic=1,nc
            IF (.NOT. is_snow_covered(ic) .AND. is == 1) THEN
              ! Without snow cover, the top soil layer took part in the surface energy balance and the
              ! phase change latent heat needs to be accounted for here.
              IF (le_phase(ic) < 0._wp .AND. t_soil_sl(ic,is) < tmelt .AND. wtr(ic,is) > liquid_max(ic,is)) THEN
                freeze_mass(ic,is) = MIN(-le_phase(ic) / alf, heat_cap(ic,is) * (tmelt - t_soil_sl(ic,is)) / alf, &
                  &                     (wtr(ic,is) - liquid_max(ic,is)) * rhoh2o)
              ELSE IF (le_phase(ic) > 0._wp .AND. t_soil_sl(ic,is) > tmelt .AND. ice(ic,is) > 0._wp) THEN
                melt_mass(ic,is)   = MIN(le_phase(ic) / alf,  heat_cap(ic,is) * (t_soil_sl(ic,is) - tmelt) / alf, &
                  &                      ice(ic,is) * rhoi)
              END IF
            ELSE
              ! For snow covered soil, the phase change energy is accounted for in snow layers. Top layer and deeper
              ! soil change phase based only on their temperature difference to tmelt and available water/ice.
              IF (t_soil_sl(ic,is) < tmelt .AND. wtr(ic,is) > liquid_max(ic,is)) THEN
                freeze_mass(ic,is) = MIN(heat_cap(ic,is) * (tmelt - t_soil_sl(ic,is)) / alf, &
                  &                     (wtr(ic,is) - liquid_max(ic,is)) * rhoh2o)
              ELSE IF (t_soil_sl(ic,is) > tmelt .AND. ice(ic,is) > 0._wp) THEN
                melt_mass(ic,is)   = MIN(heat_cap(ic,is) * (t_soil_sl(ic,is) - tmelt) / alf, &
                  &                      ice(ic,is) * rhoi)
              END IF
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        ! Update water/ice storage and soil layer temperature based on the actually
        ! frozen/melted water/ice mass.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is=1,nsoil
          DO ic=1,nc
            wtr(ic,is) = wtr(ic,is) + (melt_mass(ic,is) - freeze_mass(ic,is)) / rhoh2o
            ice(ic,is) = ice(ic,is) + (freeze_mass(ic,is) - melt_mass(ic,is)) / rhoi
            t_soil_sl(ic,is) = t_soil_sl(ic,is) + &
              & (freeze_mass(ic,is) - melt_mass(ic,is)) * alf / heat_cap(ic,is)
            frozen(ic,is) = freeze_mass(ic,is) / (rhoh2o * delta_time)
            melt  (ic,is) = melt_mass(ic,is) / (rhoi * delta_time)
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        ! Diagnose the remaining phase change energy that was not used for freezing/melting, and
        ! and remove its contribution to the top layer soil temperature change.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          IF (.NOT. is_snow_covered(ic)) THEN
            le_phase(ic)    = le_phase(ic) + (freeze_mass(ic,1) - melt_mass(ic,1)) * alf
            t_soil_sl(ic,1) = t_soil_sl(ic,1) - le_phase(ic) / heat_cap(ic,1)
          END IF
        END DO
        !$ACC END PARALLEL LOOP

        ! For unknown reasons, snow layers can have temperatures which are slightly above tmelt. Thus,
        ! also soil layers with freezing water can start into this routine with values > tmelt,
        ! and would be reset to this values after the energy release for freezing is accounted for.
        ! Therefore an additional correction is done, to ensure that the top soil layer cannot cross
        ! the tmelt threshold as long as freezing/melting is ongoing and will use the residual phase
        ! change heat term for this. This correction can potentially be removed again, after
        ! the positive snow temperature issue is resolved.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          IF (.NOT. is_snow_covered(ic)) THEN
            IF ((freeze_mass(ic,1) > 0._wp .AND. t_soil_sl(ic,1) > tmelt) .OR. &
              & (melt_mass(ic,1)   > 0._wp .AND. t_soil_sl(ic,1) < tmelt)) THEN
              le_phase(ic)    = le_phase(ic) + (t_soil_sl(ic,1) - tmelt) * heat_cap(ic,1)
              t_soil_sl(ic,1) = tmelt
            END IF
          END IF
        END DO
        !$ACC END PARALLEL LOOP

#ifndef _OPENACC
        IF (enforce_water_budget /= WB_IGNORE) THEN
          ! Check if phase change causes negative water or ice values
          n_errors = 0
          error_flag(:) = .FALSE.
          DO is=1,nsoil
            DO ic=1,nc
              IF (((wtr(ic,is) < -EPSILON(1._wp) .AND. frozen(ic,is) > EPSILON(1._wp)) .OR. &
                & (ice(ic,is) < -EPSILON(1._wp) .AND. melt(ic,is) > EPSILON(1._wp))) .AND.  &
                & .NOT. error_flag(ic)) THEN
                n_errors = n_errors + 1    ! count the soil moisture errors
                error_idx(n_errors) = ic   ! keep the respective cell indices
                error_flag(ic) = .TRUE.    ! mark gridcell as containing errors
              END IF
            END DO
          END DO

          ! Write error messages in case errors have been found.
          DO ie = 1, n_errors
            ic = error_idx(ie)
            WRITE (message_text_long,*) 'negative soil water or ice content due to phasechange', &
              & ' at index ', ic,'.:',                                 NEW_LINE('a'), &
              & 'Water content old [m]: ', wtr_old(ic,:),              NEW_LINE('a'), &
              & 'Ice content old [m]:   ', ice_old(ic,:),              NEW_LINE('a'), &
              & 'Frozen water [m]:      ', freeze_mass(ic,:) / rhoh2o, NEW_LINE('a'), &
              & 'Melted ice [m]:        ', melt_mass(ic,:)   / rhoi,   NEW_LINE('a'), &
              & 'Water content new [m]: ', wtr(ic,:),                  NEW_LINE('a'), &
              & 'Ice content new [m]:   ', ice(ic,:)
            IF (enforce_water_budget == WB_ERROR) THEN
              CALL finish ('calc_soil_temperature', message_text_long)
            ELSE IF (enforce_water_budget == WB_LOGGING) THEN
              CALL message ('calc_soil_temperature', message_text_long, all_print=.TRUE.)
            END IF
          END DO
        END IF
#endif
        ! Correct for numerical errors
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is=1,nsoil
          DO ic=1,nc
            wtr(ic,is) = MAX(0._wp, wtr(ic,is))
            ice(ic,is) = MAX(0._wp, ice(ic,is))
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        ! Update the le_phase
        IF (PRESENT(le_pc_remain)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1, nc
            IF (.NOT. is_snow_covered(ic)) THEN
              le_pc_remain(ic) = le_phase(ic)
            END IF
          END DO
          !$ACC END PARALLEL LOOP
        END IF

        ! Update supercooled water based on actual soil temperature.
        IF (l_supercool) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
          DO is=1,nsoil
            DO ic=1,nc
              liquid_max(ic,is) = Get_liquid_max( &
                & t_soil_sl(ic,is), ws_max(ic,is), matric_pot(ic,is), bclapp(ic,is))
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END IF

      END IF

      ! Thaw depth
      ! The general idea of the thaw depth calculation is to not simply determine the first layer with sub-zero temperatures,
      ! but to assume a linear temperatur gradient between two layers (either layer centers or layer bottoms) and determine
      ! the depth at which the temperature is zero between an unfrozen and a frozen layer. The underlying assumption is thus,
      ! that a layer with an average temperature below the melting point may contain an unfrozen part if the temperature of
      ! the overlying layer is above the melting point. Here, the unfrozen fraction of an, on average, frozen layer is
      ! determined by scaling the temperature excess of the overlying layer (t_lm1-tmelt) by the temperature differences
      ! between the layer and the overlying layer (t_lm1-t_l). For layers above the bedrock boundary the ice content needs to
      ! be accounted for since this represents an additionional energy sink. This is done including the temperature change
      ! corresponding to the latent energy of thaw at the soils current heat capacity. For layers at the thaw front
      ! (temperature at 0C) it is simply assumed that the unfrozen fraction of the layer is proportional to the liquid
      ! water content, relative to the total water content.
      ! Calculation for full layers (bots), rather than layer centers (mids) appears more consistant with the actual energy
      ! balance in the soil since temperatures can only increase above the melting point when the entirety of the ice in the
      ! respective layer has completly melted. However, thaw depths are quite large, so it may become necessary to switch
      ! back to the assumption that temepratures are valid for the center rather than the whole layer.
      IF (PRESENT(wtr)) THEN              ! Not under glacier
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          thaw_depth(ic) = bots(nsoil)      ! Initialize with maximum thaw depth
          IF (t_soil_sl(ic,1) < tmelt) THEN ! Top soil layer is frozen
            thaw_depth(ic) = 0._wp
          ELSE IF (t_soil_sl(ic,1) == tmelt) THEN ! Top soil layer is thawing
            ! in case water is present we diagnose thaw depth using ration of liquid / total water
            IF (ice(ic,1) > 0._wp) THEN
              thaw_depth(ic) = dsoil(ic,1) * wtr(ic,1) / (wtr(ic,1) + ice(ic,1))
            ELSE
              thaw_depth(ic) = dz(1)
            END IF
          END IF
        END DO
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
        !$ACC LOOP SEQ
        DO is = 2, nsoil
          !$ACC LOOP GANG VECTOR
          DO ic = 1, nc
            IF (thaw_depth(ic) == bots(nsoil)) THEN
              IF  (t_soil_sl(ic,is-1) > tmelt .AND. t_soil_sl(ic,is) == tmelt) THEN ! Lower layer is at zero degrees
                IF (ice(ic,is) > 0._wp) THEN
                ! in case water is present we diagnose thaw depth using ration of liquid / total water
                  thaw_depth(ic) = bots(is-1) + dsoil(ic,is) * wtr(ic,is) / (wtr(ic,is) + ice(ic,is))
                ELSE
                  thaw_depth(ic) = bots(is)
                END IF
              ELSE IF (t_soil_sl(ic,is-1) > tmelt .AND. t_soil_sl(ic,is  ) < tmelt) THEN
                  ! in case water is present this would increase the energy required to bring temperatures above the melting point
                  thaw_depth(ic) = bots(is-1) + dsoil(ic,is) * (t_soil_sl(ic,is-1) - tmelt) &
                    & / (t_soil_sl(ic,is-1) - t_soil_sl(ic,is) + ice(ic,is) * alf * rhoi / heat_cap(ic,is))
              END IF
            END IF
          END DO
          !$ACC END LOOP
        END DO
        !$ACC END LOOP
        !$ACC END PARALLEL
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          thaw_depth(ic) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    ELSE ! lstart
      IF (l_freeze .AND. PRESENT(wtr)) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO is = 1, nsoil
          DO ic = 1, nc
            liquid_max(ic,is) = 0._wp
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF

    !------------------------------------------------------------------------------------------------
    !  Computation of the Richtmyer and Morton A and B coefficients for the next time step
    !------------------------------------------------------------------------------------------------
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil
      DO ic=1,nc
        zdz2(ic,is) = heat_cap(ic,is) / delta_time
        z1(ic,is)   = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1, nsoil-1
      DO ic=1,nc
        zdz1(ic,is) = zd1(is) * heat_cond(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! lowest soil layer uses a zero flux bottom condition
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      z1(ic,nsoil-1) = zdz2(ic,nsoil) + zdz1(ic,nsoil-1)
      t_soil_acoef(ic,nsoil-1) = zdz2(ic,nsoil) * t_soil_sl(ic,nsoil) / z1(ic,nsoil-1)
      t_soil_bcoef(ic,nsoil-1) = zdz1(ic,nsoil-1) / z1(ic,nsoil-1)
    END DO
    !$ACC END PARALLEL LOOP

    ! soil layers above consider the A/B coefficients of the layer below
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR
    DO ic=1,nc
      !$ACC LOOP SEQ
      DO is = nsoil-1, 2, -1
        z1(ic,is-1) = 1._wp / (zdz2(ic,is) + zdz1(ic,is-1) + zdz1(ic,is) * (1._wp - t_soil_bcoef(ic,is)))
        t_soil_acoef(ic,is-1) = (t_soil_sl(ic,is) * zdz2(ic,is) + zdz1(ic,is) * t_soil_acoef(ic,is)) * z1(ic,is-1)
        t_soil_bcoef(ic,is-1) = zdz1(ic,is-1) * z1(ic,is-1)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      grnd_hflx(ic) = zdz1(ic,1) * (t_soil_acoef(ic,1) + (t_soil_bcoef(ic,1) - 1._wp) * t_soil_sl(ic,1))
      hcap_grnd(ic) = (zdz2(ic,1) * delta_time + delta_time * (1._wp - t_soil_bcoef(ic,1)) * zdz1(ic,1))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE calc_soil_temperature

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Calculate the temperature evolution of the snow layers.
  !>
  !> The routine propagates snow temperatures through the snow column and derives the
  !> temperature of the top soil layer from the snow-layer coefficients.
  !>
  SUBROUTINE calc_snow_temperature(         &
    & nc, nsnow,                            &
    & lstart,                               &
    & itop_old, itop,                       &
    & t_srf, t_soil_acoef, t_soil_bcoef,    &
    & t_snow, t_snow_acoef, t_snow_bcoef,   &
    & t_soil_top)

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    :: nc                         !< length of the vector
    INTEGER,  INTENT(in)    :: nsnow                      !< number of (fixed) snow layers
    LOGICAL,  INTENT(in)    :: lstart                     !< if T, start of experiment
    INTEGER,  INTENT(in)    :: itop(:), itop_old(:)       !< Top snow layer for previous and new time step
    REAL(wp), INTENT(in)    :: t_srf(:)                   !< Surface temperature at top of snow [K]
    REAL(wp), INTENT(in)    :: t_soil_acoef(:)            !< Top layer soil A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_soil_bcoef(:)            !< Top layer soil B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(inout) :: t_snow(:,:)                !< Snow temperature [K]
    REAL(wp), INTENT(inout) :: t_snow_acoef(:,:)          !< Snow A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(inout) :: t_snow_bcoef(:,:)          !< Snow B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(out)   :: t_soil_top(:)              !< New temperature of top soil layer [K]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: ic, is !< Horizontal-cell and snow-layer indices

    ! Initialize snow layer temperatures with surface temperature
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO ic=1,nc
      DO is=1,nsnow
        t_snow(ic,is) = t_srf(ic)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !-----------------------------------------------------------------------------------------------
    !  Computation of new snow temperatures
    !-----------------------------------------------------------------------------------------------
    !
    IF (.NOT. lstart) THEN
      ! Initialize a/b coefficients for snow layers which did not exists in the last time step
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO ic=1,nc
        DO is=1,nsnow
          IF (itop(ic) <= is .AND. itop_old(ic) > nsnow) THEN         ! New snow layers on formerly snowfree soil
            t_snow_acoef(ic,is) = t_soil_acoef(ic)                    ! --> use old top soil coefficients
            t_snow_bcoef(ic,is) = t_soil_bcoef(ic)
          ELSE IF (itop(ic) <= is .AND. is <= itop_old(ic)) THEN      ! zero+ new snow layers: from new top layer to top prev layer
            t_snow_acoef(ic,is) = t_snow_acoef(ic,itop_old(ic))       ! --> use old top snow coefficients
            t_snow_bcoef(ic,is) = t_snow_bcoef(ic,itop_old(ic))
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ! Compute snow layer temperatures from a/b coefficients, except for the top snow layer which
      ! retains the surface temperature
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        !$ACC LOOP SEQ
        DO is=2,nsnow
          IF (is > itop(ic)) THEN                                     ! Snow layers below new top layer
            t_snow(ic,is) = t_snow_acoef(ic,is-1) + t_snow_bcoef(ic,is-1) * t_snow(ic,is-1)
          END IF
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

      ! Compute new temperature for top soil layer ( or use the surface temperature where there was or is no snow)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        IF (MAX(itop(ic), itop_old(ic)) <= nsnow) THEN
          t_soil_top(ic) = t_snow_acoef(ic,nsnow) + t_snow_bcoef(ic,nsnow) * t_snow(ic,nsnow)
        ELSE
          t_soil_top(ic) = t_srf(ic)
        END IF
      END DO
      !$ACC END PARALLEL LOOP

    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        t_soil_top(ic) = t_srf(ic)
      END DO
      !$ACC END PARALLEL LOOP

    END IF

  END SUBROUTINE calc_snow_temperature

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Calculate snow-layer A/B coefficients and the associated heat fluxes.
  !>
  !> The routine assembles the implicit coefficients for the snow layers and the
  !> coupling to the upper soil layers, while also providing the ground heat flux and
  !> the heat capacity of the surface layer.
  !>
  SUBROUTINE calc_snow_abcoeff(           &
    & nc, nsnow,                          &
    & dz_soil,                            &
    & delta_time,                         &
    & dz, itop,                           &
    & vol_heat_cap, heat_cond,            &
    & vol_heat_cap_soil, heat_cond_soil,  &
    & t_soil_sl,                          &
    & t_soil_acoef, t_soil_bcoef,         &
    & t_snow, t_snow_acoef, t_snow_bcoef, &
    & hcap_grnd, grnd_hflx)

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    :: nc                         !< length of the vector
    INTEGER,  INTENT(in)    :: nsnow                      !< number of (fixed) snow layers
    REAL(wp), INTENT(in)    :: delta_time                 !< Time step
    REAL(wp), INTENT(in)    :: dz(:,:)                    !< (dynamic) snow layer thickness [m]
    REAL(wp), INTENT(in)    :: dz_soil(2)                 !< Thickness of two uppermost soil layers [m]
    INTEGER,  INTENT(in)    :: itop(:)                    !< Top snow layer for previous and new time step
    REAL(wp), INTENT(in)    :: vol_heat_cap(:)            !< snow volumetric heat capacity [J m-3 K-1]
    REAL(wp), INTENT(in)    :: heat_cond(:)               !< snow thermal conductivity [J m-1 s-1 K-1]
    REAL(wp), INTENT(in)    :: vol_heat_cap_soil(:)       !< top soil layer volumetric heat capacity [J m-3 K-1]
    REAL(wp), INTENT(in)    :: heat_cond_soil(:)          !< top soil layer heat conductivity [J m-1 s-1 K-1]
    REAL(wp), INTENT(in)    :: t_soil_sl(:,:)             !< Temperature of soil layers [K]
    REAL(wp), INTENT(in)    :: t_soil_acoef(:,:)          !< top soil layer A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_soil_bcoef(:,:)          !< top soil layer B coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_snow(:,:)                !< snow temperature [K]

    REAL(wp), INTENT(out)   :: t_snow_acoef(:,:)          !< snow A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(out)   :: t_snow_bcoef(:,:)          !< snow B coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(out)   :: hcap_grnd(:)               !< heat capacity of the top snow layer [J m-2 K-1]
    REAL(wp), INTENT(out)   :: grnd_hflx(:)               !< ground heat flux [J m-2 s-1]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: ic, is !< Horizontal-cell and snow-layer indices
    REAL(wp) :: zmid(nc,nsnow+2) !< Midpoint depths of snow and upper soil layers [m]
    REAL(wp) :: z1 !< Inverse coefficient for the implicit snow-temperature solver
    REAL(wp) :: zd1(nc,nsnow+1) !< Inverse distances between snow and soil-layer midpoints [m-1]
    REAL(wp) :: zdz1(nc,nsnow+1), zdz2(nc,nsnow+1) !< Conductive and storage coefficients of snow layers

    !$ACC DATA CREATE(zmid, z1, zd1, zdz1, zdz2) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      grnd_hflx(ic) = 0._wp
      hcap_grnd(ic) = 0._wp

      !-----------------------------------------------------------------------------------------------
      !  Computation of useful constants
      !
      ! dz and zmid are zero for unused snow layers
      !
      ! mid of snow layers
      zmid(ic,1) = dz(ic,1) / 2._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is=2,nsnow
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        zmid(ic,is) = zmid(ic,is-1) + 0.5_wp * (dz(ic,is-1) + dz(ic,is))
        IF (zmid(ic,is) > zmid(ic,is-1)) THEN
          zd1(ic,is-1) = 1._wp / (zmid(ic,is) - zmid(ic,is-1))
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    ! mid of two uppermost soil layers
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      zmid(ic,nsnow+1) = zmid(ic,nsnow) + 0.5_wp * (dz(ic,nsnow) + dz_soil(1))
      zmid(ic,nsnow+2) = zmid(ic,nsnow+1) + 0.5_wp * (dz_soil(1) + dz_soil(2))
      zd1(ic,nsnow) = 1._wp / (zmid(ic,nsnow+1) - zmid(ic,nsnow))
      zd1(ic,nsnow+1) = 1._wp / (zmid(ic,nsnow+2) - zmid(ic,nsnow+1))
    END DO
    !$ACC END PARALLEL LOOP

    !
    !------------------------------------------------------------------------------------------------
    !  Computation of the Richtmyer and Morton A and B coefficients for the next time step
    !------------------------------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = 1, nsnow
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        IF (is <= itop(ic) - 1) THEN
          t_snow_acoef(ic,is) = 0._wp
          t_snow_bcoef(ic,is) = 0._wp
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    ! Compute thermal storage and conductivity coefficients for snow layers
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO is = 1,nsnow
      DO ic=1, nc
        IF (is >= itop(ic)) THEN
          zdz2(ic,is) = vol_heat_cap(ic) * dz(ic,is) / delta_time
          zdz1(ic,is) = zd1(ic,is) * heat_cond(ic)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Compute thermal storage and conductivity coefficients as well as Richtmyer and Morton A/B coefficients
    ! for the lowest snow layer, which is coupled to the top soil layer
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) PRIVATE(z1)
    DO ic=1,nc
      zdz2(ic,nsnow+1) = vol_heat_cap_soil(ic) * dz_soil(1) / delta_time
      zdz1(ic,nsnow+1) = zd1(ic,nsnow+1) * heat_cond_soil(ic)

      IF (itop(ic) <= nsnow) THEN
        ! Snow is present --> Compute for lowest layer
        z1 = 1._wp / (zdz2(ic,nsnow+1) + zdz1(ic,nsnow) + zdz1(ic,nsnow+1) * (1._wp - t_soil_bcoef(ic,1)))
        t_snow_acoef(ic,nsnow) = (t_soil_sl(ic,1) * zdz2(ic,nsnow+1) + zdz1(ic,nsnow+1) * t_soil_acoef(ic,1)) * z1
        t_snow_bcoef(ic,nsnow) = zdz1(ic,nsnow) * z1
      ELSE
        ! Without snow coefficients are set to zero
        t_snow_acoef(ic,nsnow) = 0._wp
        t_snow_bcoef(ic,nsnow) = 0._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! snow layers above up to top snow layer
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO is = nsnow, 2, -1
      !$ACC LOOP GANG VECTOR PRIVATE(z1)
      DO ic=1,nc
        IF (is > itop(ic)) THEN
          z1 = 1._wp / (zdz2(ic,is) + zdz1(ic,is-1) + zdz1(ic,is) * (1._wp - t_snow_bcoef(ic,is)))
          t_snow_acoef(ic,is-1) = (t_snow(ic,is) * zdz2(ic,is) + zdz1(ic,is) * t_snow_acoef(ic,is)) * z1
          t_snow_bcoef(ic,is-1) = zdz1(ic,is-1) * z1
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) PRIVATE(is)
    DO ic=1,nc
      IF (itop(ic) <= nsnow) THEN
        is = itop(ic)
        grnd_hflx(ic) = zdz1(ic,is) * (t_snow_acoef(ic,is) + (t_snow_bcoef(ic,is) - 1._wp) * t_snow(ic,is))
        hcap_grnd(ic) = (zdz2(ic,is) * delta_time + delta_time * (1._wp - t_snow_bcoef(ic,is)) * zdz1(ic,is))
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE calc_snow_abcoeff

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !> Estimate the amount of water that remains liquid in frozen soil layers.
  !>
  !> The expression follows a supercooling formulation and provides the maximum
  !> liquid water content allowed for a given temperature and soil-water state.
  !> It is based on (Niu & Yang, 2006)
  !>
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION Get_liquid_max( &
    & temp, w_max, matric_pot, bclapp     &
    & ) RESULT(liquid_max)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, alf, grav

    REAL(wp), INTENT(in) :: &
      & temp,               & !< Soil temperature [K]
      & w_max,              & !< Maximum water storage [m]
      & matric_pot,         & !< Soil matric potential [m]
      & bclapp                !< Clapp-Hornberger parameter []
    REAL(wp) :: liquid_max    !< Maximum liquid water content at the given temperature [m]

    IF (temp < tmelt - 1e-10_wp) THEN
      ! supercooled water equation (Niu&Yang,2006)
      liquid_max =                                                               &
        & w_max                                                                  & ! Maximum water storage
        & * ( alf * (tmelt - temp) / MAX(1.e-10_wp,  grav * temp * (-matric_pot)) ) &
        &   ** ( -1._wp / MAX(1._wp, bclapp))
    ELSE
      liquid_max = 0._wp
    END IF

  END FUNCTION Get_liquid_max

#endif
END MODULE mo_sse_process
