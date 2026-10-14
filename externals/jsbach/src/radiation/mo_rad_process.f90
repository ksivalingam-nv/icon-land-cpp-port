!> Contains the routines for the radiation processes
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
MODULE mo_rad_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control, ONLY: acc_stream
  USE mo_kind,        ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC ::                        &
    & calc_radiation_surface_net,  &
    & calc_alb_lwtr, calc_alb_lice, calc_glacier_albedo, Has_minimal_radiation, calc_snow_age,                 &
    & calc_pond_albedo, calc_snow_albedo, &
    & get_surface_albedo_simple

  INTERFACE Has_minimal_radiation
    MODULE PROCEDURE Has_minimal_radiation_1
    MODULE PROCEDURE Has_minimal_radiation_2
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_process'

CONTAINS

  SUBROUTINE calc_radiation_surface_net( &
    & swvis_down,    & ! in
    & swnir_down,    & ! in
    & alb_vis,       & ! in
    & alb_nir,       & ! in
    & lw_down,       & ! in
    & t,             & ! in
    & rad_net,       & ! output
    & swvis_net,     & ! output, optional
    & swnir_net,     & ! output, optional
    & sw_net,        & ! output, optional
    & lw_net         & ! output, optional
    & )

    USE mo_phy_schemes,   ONLY: lwnet_from_lwdown

    REAL(wp),  INTENT(in), DIMENSION(:) :: &
      & swvis_down,          &
      & swnir_down,          &
      & alb_vis,             &
      & alb_nir,             &
      & lw_down,             &
      & t

    REAL(wp),  INTENT(out) :: &
      & rad_net(:)

    REAL(wp),  INTENT(out), OPTIONAL, DIMENSION(:) :: &
      & swvis_net,                      & ! net shortwave radiation in visible range
      & swnir_net,                      & ! net shortwave radiation in near infrared range
      & sw_net,                         &
      & lw_net

    REAL(wp) :: zswvis_net, zswnir_net, zsw_net, zlw_net

    INTEGER :: nc, ic

    nc = SIZE(swvis_down)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(zswvis_net, zswnir_net, zsw_net, zlw_net)
    DO ic=1,nc

      ! Compute net SW radiation from downward SW and albedo
      zswvis_net = swvis_down(ic) * (1._wp - alb_vis(ic))
      zswnir_net = swnir_down(ic) * (1._wp - alb_nir(ic))
      zsw_net    = zswvis_net + zswnir_net

      ! Compute LW net radiation from incoming and the thermal radiation from the surface
      zlw_net = lwnet_from_lwdown(lw_down(ic), t(ic))

      ! Compute net radiation
      rad_net(ic) = zsw_net + zlw_net

      IF (PRESENT(swvis_net)) swvis_net(ic) = zswvis_net
      IF (PRESENT(swnir_net)) swnir_net(ic) = zswnir_net
      IF (PRESENT(sw_net))    sw_net(ic)    = zsw_net
      IF (PRESENT(lw_net))    lw_net(ic)    = zlw_net

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_radiation_surface_net

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Compute non-dimensional snow age
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_snow_age( &
    & dtime,                & ! in
    & temperature_surface,  & ! in
    & snow,                 & ! in
    & fract_snow_soil,      & ! in
    & snow_age              & ! inout
    & )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants,  ONLY: tmelt     ! tmelt =273.15_wp

    REAL(wp), INTENT(in)    :: dtime,               & !< model time step [s]
      &                        temperature_surface, & !< surface temperature [K]
      &                        snow,                & !< snowfall [kg m-2 s-1]
      &                        fract_snow_soil        !< Fraction of snow covered ground [-]

    REAL(wp), INTENT(inout) :: snow_age               !< non-dimensional age of snow

    REAL(wp)                :: hlp1                   !< temporary variable

    ! ---------------------------
    ! Update snow age where ever there is snow on the ground
    ! Snow age increases with time and temperatures closer to the melting point
    ! It decreases with new snowfall or resets to zero in case all snow is melted

    IF (fract_snow_soil > EPSILON(1.0_wp)) THEN
      hlp1     = MIN(1._wp, EXP(5000._wp * ((1._wp / tmelt) - (1._wp / temperature_surface))))
      snow_age = MAX(0.0_wp, (snow_age +                          &
        & ((hlp1**10._wp + hlp1 + 0.3_wp) * dtime * 1.e-06_wp)) * &
        & (1._wp - (snow * dtime * 0.5_wp)))
    ELSE
      snow_age = 0._wp
    END IF

  END SUBROUTINE calc_snow_age

  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_alb_lwtr( &
    & l_day,                               &
    & alb_lwtr                             &
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants, ONLY: AlbedoLakeWater

    ! Arguments
    LOGICAL,  INTENT(in)    :: &
      & l_day
    REAL(wp), INTENT(inout) :: &
      & alb_lwtr

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! @todo: use method from mo_surface_ocean:update_albedo_ocean from ECHAM6
    alb_lwtr = AlbedoLakeWater

  END SUBROUTINE calc_alb_lwtr
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_alb_lice( &
    & l_day,         &
    & t,             &
    & weq_snow_lice, &
    & alb_lice       &
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants,          ONLY: AlbedoLakeIceMin, AlbedoLakeIceMax, AlbedoLakeSnowMin, AlbedoLakeSnowMax
    USE mo_jsb_physical_constants, ONLY: tmelt

    ! Arguments
    LOGICAL,  INTENT(in)    :: &
      & l_day
    REAL(wp), INTENT(in)    :: &
      & t,                     &
      & weq_snow_lice
    REAL(wp), INTENT(inout) :: &
      & alb_lice

    ! Local variables
    REAL(wp) :: t_upper_limit    ! Upper temperature limit for cold snow albedo
    REAL(wp) :: alb_min, alb_max ! Minimum and maximum snow albedo

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    t_upper_limit = tmelt - 1._wp

    IF (weq_snow_lice > 0.01_wp) THEN
      alb_min = AlbedoLakeSnowMin
      alb_max = AlbedoLakeSnowMax
    ELSE
      alb_min = AlbedoLakeIceMin
      alb_max = AlbedoLakeIceMax
    END IF

    ! Temperature-dependent snow albedo
    IF (t >= tmelt) THEN
      alb_lice = alb_min
    ELSE IF (t < t_upper_limit) THEN
      alb_lice = alb_max
    ELSE
      alb_lice = alb_min + &
        &      ( ((alb_max - alb_min) / (tmelt - t_upper_limit)) &  ! Change of snow albedo per deg C
        &        * (tmelt - t))
    END IF

  END SUBROUTINE calc_alb_lice
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_glacier_albedo( &
    & l_day,   &
    & t,       &
    & alb_vis, &
    & alb_nir  &
    & )

    !$ACC ROUTINE SEQ

    ! Use declarations
    USE mo_rad_constants, ONLY: AlbedoGlacierVisMin,  & ! Albedo of glacier in the visible range at the melting point
                                AlbedoGlacierVisMax,  & ! Albedo of glacier in the visible range at hard frost
                                TempAlbedoGlacierMax, & ! Maximum glacier albedo at this temperature below melting point of H2O
                                AlbedoGlacierNirMin,  & ! Albedo of glacier in the NIR range at at the melting point
                                AlbedoGlacierNirMax     ! Albedo of glacier in the NIR range at hard frost
    USE mo_jsb_physical_constants, ONLY: tmelt

    ! Arguments
    LOGICAL,  INTENT(in)    :: &
      & l_day
    REAL(wp), INTENT(in)    :: &
      & t
    REAL(wp), INTENT(inout) :: &
      & alb_vis,               & ! Albedo of canopy in the VIS range without snow cover.
      & alb_nir                  ! Albedo of canopy in the NIR range without snow cover.

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! calculate new glacier albedo (visible and NIR) only at grid points with solar radiation (i.e. during daytime)
    IF (t >= tmelt) THEN
       alb_vis = AlbedoGlacierVisMin
       alb_nir = AlbedoGlacierNirMin
    ELSE IF (t < tmelt - TempAlbedoGlacierMax) THEN
       alb_vis = AlbedoGlacierVisMax
       alb_nir = AlbedoGlacierNirMax
    ELSE
       alb_vis = AlbedoGlacierVisMin +                                  &
            (tmelt - t) * (AlbedoGlacierVisMax - AlbedoGlacierVisMin) / &
            TempAlbedoGlacierMax
       alb_nir = AlbedoGlacierNirMin +                                  &
            (tmelt - t) * (AlbedoGlacierNirMax - AlbedoGlacierNirMin) / &
             TempAlbedoGlacierMax
    END IF

  END SUBROUTINE calc_glacier_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_snow_albedo( &
    & l_day,             & ! in
    & fract_snow_soil,   & ! in
    & snow_age,          & ! in
    & cos_zenith_angle,  & ! in
    & albedo_age_weight, & ! in
    & t_srf,             & ! in
    & alb_vis_snow,      & ! inout
    & alb_nir_snow       & ! inout
    & )

    !$ACC ROUTINE SEQ

    ! Declarations
    USE mo_rad_constants, ONLY: &
      & ZenithAngleFactor,      & ! Factor in solar zenith angle dependence of snow albedo
                                  ! (the increase of snow albedo is the higher this factor)
      & AlbedoSnowVisMax_age,   &
      & AlbedoSnowNirMax_age,   &
      & AlbedoSnowVisAge,       & ! Maximal rel. reduction of snow albedo by aging in the visible range
      & AlbedoSnowNirAge,       & ! Maximal rel. reduction of snow albedo by aging in the NIR range
      & AlbedoSnowAngle,        & ! Maximal rel. reduction of snow absorption by large solar zenith angle
      & TempAlbedoSnowMax,      & ! Maximum snow albedo at this temperature below melting point of H2O
      & AlbedoSnowVisMin_temp,  &
      & AlbedoSnowNirMin_temp,  &
      & AlbedoSnowVisMax_temp,  &
      & AlbedoSnowNirMax_temp

    USE mo_jsb_physical_constants,  ONLY: tmelt     ! tmelt =273.15_wp

    ! Arguments
    LOGICAL,  INTENT(in)  :: l_day
    REAL(wp), INTENT(in)  :: fract_snow_soil,   & !< Fraction of snow covered ground
                             snow_age,          & !< non-dimensional age of snow
                             cos_zenith_angle,  & !< cosinus of zenith angle of incoming radiation
                             t_srf

    REAL(wp), INTENT(in)  ::                    & !
                             albedo_age_weight    ! 0 < albedo_age_weight < 1: snow albedo is calcualted by linearly weighting
                                                  !                            the snow albedo resulting from both age and
                                                  !                            temperature schemes
    REAL(wp), INTENT(inout) ::              &
                             alb_vis_snow,  & ! Albedo of canopy in the VIS range without snow cover.
                             alb_nir_snow     ! Albedo of canopy in the NIR range without snow cover.

    ! Local variables
    REAL(wp)              ::                    &
                             alb_vis_snow_age,  &
                             alb_nir_snow_age,  &
                             alb_vis_snow_temp, &
                             alb_nir_snow_temp, &
                             snow_age_factor,   &
                             snow_albedo_angle_factor

    ! Initialize local variables
    alb_vis_snow_age  = 0.0_wp
    alb_nir_snow_age  = 0.0_wp
    alb_vis_snow_temp = 0.0_wp
    alb_nir_snow_temp = 0.0_wp
    alb_vis_snow      = 0.0_wp
    alb_nir_snow      = 0.0_wp

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    ! calculate albedo of snow dependent on age (snow age scheme, from BATS model of NCAR)
    IF (fract_snow_soil > EPSILON(1._wp)) THEN
       snow_age_factor = snow_age / (1._wp + snow_age)
       IF (cos_zenith_angle < 0.5_wp) THEN
         snow_albedo_angle_factor = ((1.0_wp + ZenithAngleFactor)/ &
                                    (1.0_wp + 2.0_wp * ZenithAngleFactor * cos_zenith_angle) - 1.0_wp) / &
                                  ZenithAngleFactor
       ELSE
          snow_albedo_angle_factor = 0._wp
       END IF
       alb_vis_snow_age = AlbedoSnowVisMax_age * (1._wp - AlbedoSnowVisAge * snow_age_factor)
       alb_vis_snow_age = alb_vis_snow_age + AlbedoSnowAngle * snow_albedo_angle_factor * &
                                                   (1._wp - alb_vis_snow_age)
       alb_nir_snow_age = AlbedoSnowNirMax_age * (1._wp - AlbedoSnowNirAge * snow_age_factor)
       alb_nir_snow_age = alb_nir_snow_age + AlbedoSnowAngle * snow_albedo_angle_factor * &
                                                   (1._wp - alb_nir_snow_age)
    ELSE
       alb_vis_snow_age = 0._wp
       alb_nir_snow_age = 0._wp
    END IF


    ! calculate albedo of snow dependent on temperature (snow temperature scheme, from echam5 model of MPI-M)
    IF (t_srf >= tmelt) THEN
          alb_vis_snow_temp = AlbedoSnowVisMin_temp
          alb_nir_snow_temp = AlbedoSnowNirMin_temp
    ELSE IF (t_srf < tmelt - TempAlbedoSnowMax) THEN
          alb_vis_snow_temp = AlbedoSnowVisMax_temp
          alb_nir_snow_temp = AlbedoSnowNirMax_temp
    ELSE
          alb_vis_snow_temp = AlbedoSnowVisMin_temp +     &
               (tmelt - t_srf) * (AlbedoSnowVisMax_temp - &
               AlbedoSnowVisMin_temp) / TempAlbedoSnowMax
          alb_nir_snow_temp = AlbedoSnowNirMin_temp +     &
               (tmelt - t_srf) * (AlbedoSnowNirMax_temp - &
               AlbedoSnowNirMin_temp) / TempAlbedoSnowMax
    END IF

    ! weight age vs. temperature of snow albedo to receive alb_vis_snow and alb_nir_snow
    alb_vis_snow = albedo_age_weight * alb_vis_snow_age + &
                              (1._wp - albedo_age_weight) * alb_vis_snow_temp
    alb_nir_snow = albedo_age_weight * alb_nir_snow_age + &
                              (1._wp - albedo_age_weight) * alb_nir_snow_temp

  END SUBROUTINE calc_snow_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_pond_albedo( &
    & l_day,             & ! in
    & cos_zenith_angle,  & ! in
    & t,                 & ! in
    & wtr_pond,          & ! in
    & ice_pond,          & ! in
    & alb_pond_vis,      & ! inout
    & alb_pond_nir       & ! inout
    & )

    !$ACC ROUTINE SEQ

    ! Declarations
    USE mo_rad_constants,           ONLY: AlbedoLakeIceMin, AlbedoLakeIceMax

    USE mo_jsb_physical_constants,  ONLY: tmelt     ! tmelt =273.15_wp

    ! Arguments
    LOGICAL,  INTENT(in)  :: l_day
    REAL(wp), INTENT(in)  :: cos_zenith_angle,  & ! cosinus of zenith angle of incoming radiation
                             t,                 &
                             wtr_pond,          &
                             ice_pond

    REAL(wp), INTENT(inout) :: alb_pond_vis,    & ! Albedo of ponds in the VIS range without snow cover
                               alb_pond_nir       ! Albedo of ponds in the NIR range without snow cover

    ! Local variables

    REAL(wp) :: t_upper_limit    ! Upper temperature limit for cold snow albedo

    REAL(wp)              :: alb_min,           &
                             alb_max,           &
                             alb_wtr_pond_vis,  &
                             alb_wtr_pond_nir,  &
                             alb_ice_pond_vis,  &
                             alb_ice_pond_nir,  &
                             fact_ice

    ! ---------------------------
    ! Go

    IF (.NOT. l_day) RETURN

    !!!! Albedo formulation for unfrozen ponds is taken from CLM
    alb_wtr_pond_vis = 0.05_wp /(MAX(0.001_wp,cos_zenith_angle) + 0.15_wp)
    alb_wtr_pond_nir = alb_pond_vis


    ! Temperature-dependent ice albedo (taken from lake ice fomulation)
    t_upper_limit     = tmelt - 1._wp
    IF (t >= tmelt) THEN
      alb_ice_pond_vis = AlbedoLakeIceMin
    ELSE IF (t < t_upper_limit) THEN
      alb_ice_pond_vis = AlbedoLakeIceMax
    ELSE
      alb_ice_pond_vis = AlbedoLakeIceMin + &
        &      ( ((AlbedoLakeIceMax - AlbedoLakeIceMin) / (tmelt - t_upper_limit)) &  ! Change of snow albedo per deg C
        &        * (tmelt - t))
    END IF
    ! Albedo nir from CLM
    alb_ice_pond_nir = 0.4_wp

    ! Agregation of albedo of ice and wtr under the assumption that ice is predominantly
    ! near the surfce hence has a stronger influence on albedo
    IF (ice_pond > 0._wp) THEN
      fact_ice     = MAX(0._wp, MIN(1._wp, 1._wp - (wtr_pond/(wtr_pond+ice_pond))**2._wp))
      alb_pond_vis = fact_ice * alb_ice_pond_vis + (1.0_wp - fact_ice) * alb_wtr_pond_vis
      alb_pond_nir = fact_ice * alb_ice_pond_nir + (1.0_wp - fact_ice) * alb_wtr_pond_nir
    ELSE
      alb_pond_vis = alb_wtr_pond_vis
      alb_pond_nir = alb_wtr_pond_nir
    END IF

  END SUBROUTINE calc_pond_albedo
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
  ! This is the echam5 scheme
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  REAL(wp) FUNCTION get_surface_albedo_simple( &
    & t_srf,             &
    & fract_forest,      &
    & grnd_snow_fract,   & ! R: fract_snow_soil umbenannt, da grnd_snow_fract nicht unbedingt dasselbe sein muss...
    & fract_snow_can,    &
    & fract_pond,        &
    & lai,               &
    & alb_background     &
    & )

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    ! Use declarations
    USE mo_jsb_physical_constants,  ONLY: tmelt
    USE mo_rad_constants,           ONLY: AlbedoCanopySnow_simple, AlbedoSnowMin_simple, AlbedoSnowMax_simple, &
      &                                   AlbedoLakeWater, AlbedoLakeIceMin

    ! Arguments
    REAL(wp), INTENT(in) ::                  &
                            t_srf,           &
                            fract_forest,    &
                            grnd_snow_fract, &
                            fract_snow_can,  &
                            fract_pond,      &
                            lai,             &
                            alb_background

    ! Local variables
    REAL(wp) ::                                   &
                            sky_view_fract,       &  ! R: how much of light reaches the ground through the canopy
                            forest_view_fract,    &  ! R: how much of light is catched only by forest
                            albedo_snow_grnd,     &  ! Albedo of snow on ground
                            min_temp_snow_albedo, &  ! Temperature below which albedo of snow is at maximum
                            albedo_grnd,          &  ! Albedo on ground
                            albedo_can,           &  ! Albedo on canopy
                            albedo_pond              ! Albedo of ponds


    ! ---------------------------
    ! Go

    ! Calculate albedo_grnd incl. snowcover
    min_temp_snow_albedo = tmelt - 5._wp

    IF (t_srf >= tmelt) THEN
      albedo_snow_grnd = AlbedoSnowMin_simple
      albedo_pond      = AlbedoLakeWater
    ELSE IF (t_srf < min_temp_snow_albedo) THEN
      albedo_snow_grnd = AlbedoSnowMax_simple
      albedo_pond      = AlbedoLakeIceMin
    ELSE
      albedo_snow_grnd = AlbedoSnowMin_simple &
                         + (tmelt - t_srf) * (AlbedoSnowMax_simple - AlbedoSnowMin_simple) / (tmelt - min_temp_snow_albedo)
      albedo_pond      = AlbedoLakeIceMin &
                         - (tmelt - t_srf) * (AlbedoLakeIceMin - AlbedoLakeWater) / (tmelt - min_temp_snow_albedo)
    END IF

    ! Calclulate albedo of ground incl. snowcover
    albedo_grnd = grnd_snow_fract                                   * albedo_snow_grnd &
                + (1._wp - grnd_snow_fract) * fract_pond            * albedo_pond      &
                + (1._wp - grnd_snow_fract) * (1.0_wp - fract_pond) * alb_background

    ! Calculate albedo_can incl. snowcover
    albedo_can = fract_snow_can * AlbedoCanopySnow_simple + (1._wp - fract_snow_can) * alb_background

    ! Calculate albedo of surface
    sky_view_fract = EXP(-1.0_wp * MAX(lai, 2._wp)) ! The value 1.0 is the SkyViewFactor.
                                                    ! The value 2.0 represents the stem_area.
                                                    ! Keep these values constant:
                                                    ! The function get_surface_albedo_simple_echam5 is not thought as an actual
                                                    ! alternative albedo scheme but to provide comparability to old ECHAM versions

    forest_view_fract = fract_forest * (1._wp - sky_view_fract)


!   R: falls man Glaciermaske verwenden will, siehe auch radiation interface...
!     If (glaciermaske) THEN
!        get_surface_albedo_simple = albedo_snow_grnd
!     ELSE
    get_surface_albedo_simple = MAX(alb_background, &
                                    (1._wp - forest_view_fract) * albedo_grnd &
                                    + forest_view_fract * albedo_can)

  END FUNCTION get_surface_albedo_simple
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !
  FUNCTION Has_minimal_radiation_2(swvis, swnir) RESULT(l_day)

    !$ACC ROUTINE SEQ

    USE mo_rad_constants, ONLY: MinRadiation ! Minimum amount of radiation for albedo calculations

    REAL(wp), INTENT(in) :: &
      & swvis, &
      & swnir
    LOGICAL :: l_day

    l_day = swvis + swnir > MinRadiation

  END FUNCTION Has_minimal_radiation_2
  !
  FUNCTION Has_minimal_radiation_1(sw) RESULT(l_day)

    !$ACC ROUTINE SEQ

    USE mo_rad_constants, ONLY: MinRadiation ! Minimum amount of radiation for albedo calculations

    REAL(wp), INTENT(in) :: &
      & sw
    LOGICAL :: l_day

    l_day = sw > MinRadiation

  END FUNCTION Has_minimal_radiation_1

#endif
END MODULE mo_rad_process
