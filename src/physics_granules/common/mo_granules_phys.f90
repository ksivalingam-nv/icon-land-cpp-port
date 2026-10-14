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

! Module to allow for a full encapsulation of granules.
! Collection of general physics constants and function for consistency between different granules.
! These are copied from mo_physical_constants and mo_lookup_tables_constants.

MODULE mo_granules_phys

  USE mo_granules_kind,        ONLY: wp

  IMPLICIT NONE

  REAL(wp), PARAMETER :: grav  = 9.80665_wp  !< [m/s2]   av. gravitational acceleration
  REAL(wp), PARAMETER :: rd    = 287.04_wp   !< [J/K/kg] gas constant
  REAL(wp), PARAMETER :: cpd   = 1004.64_wp  !< [J/K/kg] specific heat at constant pressure
  REAL(wp), PARAMETER :: cvd   = cpd - rd    !< [J/K/kg] specific heat at constant volume
  REAL(wp), PARAMETER :: rv    = 461.51_wp   !< [J/K/kg] gas constant for water vapor
  REAL(wp), PARAMETER :: alv   = 2.5008e6_wp !< [J/kg]   latent heat for vaporisation
  REAL(wp), PARAMETER :: als   = 2.8345e6_wp !< [J/kg]   latent heat for sublimation
  REAL(wp), PARAMETER :: tmelt = 273.15_wp   !< [K]      melting temperature of ice/snow
  REAL(wp), PARAMETER :: rdv   = rd / rv     !< [ ]
  !$ACC DECLARE COPYIN(grav, rd, cpd, cvd, alv, als)

  !ACCWA (Cray Fortran => 17.0.1) : non-parameter scalars have zero values in GPU when used via renames CAST-37645
#if defined(_CRAYFTN)
  REAL(wp), PARAMETER :: c1es  = 610.78_wp
  REAL(wp), PARAMETER :: c2es  = c1es * rd / rv
  REAL(wp), PARAMETER :: c3les = 17.269_wp
  REAL(wp), PARAMETER :: c3ies = 21.875_wp
  REAL(wp), PARAMETER :: c4les = 35.86_wp
  REAL(wp), PARAMETER :: c4ies = 7.66_wp
  REAL(wp), PARAMETER :: c5les = c3les * (tmelt - c4les)
  REAL(wp), PARAMETER :: c5ies = c3ies * (tmelt - c4ies)
  REAL(wp), PARAMETER :: c5alvcp = c5les * alv / cpd
  REAL(wp), PARAMETER :: c5alscp = c5ies * als / cpd
  REAL(wp), PARAMETER :: alvdcp  = alv / cpd
  REAL(wp), PARAMETER :: alsdcp  = als / cpd
  !$ACC DECLARE COPYIN(c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies) &
  !$ACC   COPYIN(c5alvcp, c5alscp, alvdcp, alsdcp)
#else
  REAL(wp) :: c1es
  REAL(wp) :: c2es
  REAL(wp) :: c3les
  REAL(wp) :: c3ies
  REAL(wp) :: c4les
  REAL(wp) :: c4ies
  REAL(wp) :: c5les
  REAL(wp) :: c5ies
  REAL(wp) :: c5alvcp
  REAL(wp) :: c5alscp
  REAL(wp) :: alvdcp
  REAL(wp) :: alsdcp
  !$ACC DECLARE CREATE(c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies, c5alvcp, c5alscp, alvdcp, alsdcp)
#endif

  PRIVATE

  ! Public parameters
  PUBLIC :: grav, alv, als, rd, rv, rdv, cpd, cvd
  PUBLIC :: tmelt
  PUBLIC :: c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies, &
      &       c5alvcp, c5alscp, alvdcp, alsdcp

  ! Public functions
  PUBLIC :: init_satpres_coeffs
  PUBLIC :: fgew
  PUBLIC :: fgee
  PUBLIC :: fgqs
  PUBLIC :: dqsdt

CONTAINS

  !---------------------------------
  ! init lookup table constants
  !---------------------------------

  INTEGER FUNCTION init_satpres_coeffs(itype_satpres)
    INTEGER, INTENT(IN), OPTIONAL :: itype_satpres
    INTEGER :: itype

    IF (PRESENT(itype_satpres)) THEN
      itype = itype_satpres
    ELSE
      itype = 1
    END IF

    init_satpres_coeffs = 0 ! default return code

    !ACCWA (Cray Fortran => 17.0.1) : non-parameter scalars have zero values in GPU when used via renames CAST-37645
#if defined(_CRAYFTN)
    IF (itype == 2) THEN
      init_satpres_coeffs = 1 ! itype=2 not supported with Cray compilers, set return code to 1
      RETURN
    END IF
#else
    ! Set coefficients for saturation pressure (originally set as fortran parameters in the declaration above, but this
    ! needed to be moved in order to switch between the old COSMO model coefficients and the more accurate IFS coefficients)

    c1es  = 610.78_wp
    c2es  = c1es * rd / rv
    IF (itype == 1) THEN ! DWD coefficients (namelist default)
      c3les = 17.269_wp
      c3ies = 21.875_wp
      c4les = 35.86_wp
      c4ies = 7.66_wp
    ELSE IF (itype == 2) THEN ! IFS coefficients
      c3les = 17.502_wp
      c3ies = 22.587_wp
      c4les = 32.19_wp
      c4ies = -0.7_wp
    END IF
    c5les = c3les * (tmelt - c4les)
    c5ies = c3ies * (tmelt - c4ies)
    c5alvcp = c5les * alv / cpd
    c5alscp = c5ies * als / cpd
    alvdcp  = alv / cpd
    alsdcp  = als / cpd

    !$ACC UPDATE DEVICE(c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies, c5alvcp, c5alscp, alvdcp, alsdcp)
#endif

  END FUNCTION init_satpres_coeffs


  !---------------------------------
  ! saturation vapour pressure
  !---------------------------------

  ELEMENTAL REAL(wp) FUNCTION fgew(ztt)
    ! function to calculate saturation vapour pressure over water
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(IN) :: ztt ! temperature
    fgew = c1es * EXP(c3les * (ztt - tmelt) / (ztt - c4les))
  END FUNCTION fgew


  ELEMENTAL REAL(wp) FUNCTION fgee(ztt)
    ! function to calculate saturation vapour pressure over ice
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(IN) :: ztt ! temperature
    fgee = c1es * EXP(c3ies * (ztt - tmelt) / (ztt - c4ies))
  END FUNCTION fgee


  PURE REAL(wp) FUNCTION fgqs(zzps,zzpv,zzpa)
    ! function to calculate saturation specific humidities from RH=esat/e (proper and safe)
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(IN) :: zzps ! saturation vapour pressure
    REAL(wp), INTENT(IN) :: zzpv !
    REAL(wp), INTENT(IN) :: zzpa
    fgqs = rdv * zzps / (zzpa - (1._wp - rdv) * zzpv)
  END FUNCTION fgqs


  PURE REAL(wp) FUNCTION dqsdt(ztt,zqs)
    ! function for dq_sat_dT
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(IN) :: ztt ! temperature
    REAL(wp), INTENT(IN) :: zqs !
    dqsdt = c5les * (1._wp - zqs) * zqs / (ztt - c4les)**2
  END FUNCTION dqsdt

END MODULE mo_granules_phys
