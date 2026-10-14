!> Contains parameters for the phenology processes
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
MODULE mo_pheno_parameters
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  ! List of phenology types
  ! evergreen   (EG)
  ! summergreen (SG)
  ! raingreen   (RG)
  ! grasses     (GRS)
  ! crops       (CRP)
  !
  ! public types used for parameters
  !
  ! R: nicht mehr notwendig laut Christian:
  !   TYPE t_bounds  !! For saving parameter boundary values
  !      LOGICAL   ::  exists = .FALSE.  ! indicates whether a boundary value exists
  !      REAL(dp)  ::  value  =  NaN     ! the boundary value if it exists
  !   END TYPE t_bounds
  !  PUBLIC :: t_bounds
  !
  !   TYPE rdata_type  ! To have a real value of a variable, its name and possible upper and lower range together. Should be set by
  !                    ! routine set_rdata().
  !      REAL(dp)          :: v =  NaN   ! the value of the variable
  !      CHARACTER(len=32) :: n =  ""    ! the name of the variable
  ! ! R: nicht mehr notwendig laut Christian    TYPE(t_bounds) :: uBound        ! the upper bound
  ! ! R: nicht mehr notwendig laut Christian    TYPE(t_bounds) :: lBound        ! the lower bound
  !      CHARACTER(len=32) :: units =""  ! the physical units of the stored values
  !   END TYPE rdata_type
  ! !  PUBLIC :: rdata_type

  !> Parameters that are common to the evergreen and raingreen phenology type
  !>
  TYPE t_EG_SG_common_param
    REAL(wp) :: &
      & tau_pseudo_soil = def_parameters      ! [days] Characteristic time for the memory loss for computing the pseudo
                                              ! soil temperature from air temperature [days] (see subroutine
                                              ! "update_pseudo_soil_temp")
  END TYPE t_EG_SG_common_param

  ! ======================================================================================================= !
  !> 2.1 TYPE that collects all JSBACH specific parameters in a single structure
  !>
  TYPE t_PhenologyParameters
    TYPE(t_EG_SG_common_param)   :: EG_SG       ! Parameters common to evergreen (EG) and summergreen (SG) phenology type
    ! TYPE(rdata_type), POINTER     :: memory(:)  => NULL() ! pointer to serial array containing the parameter values physically
    ! logical                       :: initialized =.FALSE. ! .true. indicates that phenology parameters have been initialized.
  END TYPE t_PhenologyParameters

  ! create object
  TYPE(t_PhenologyParameters), SAVE :: pheno_param_jsbach

  !$ACC DECLARE CREATE(pheno_param_jsbach%EG_SG) &
  !$ACC   CREATE(pheno_param_jsbach)

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_pheno_parameters'

CONTAINS

  ! ======================================================================================================= !
  !> initialize JSBACH specific parameters for the process: phenology
  !>
  SUBROUTINE init_pheno_parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_pheno_parameters'

    pheno_param_jsbach%EG_SG%tau_pseudo_soil  = 10.0_wp             !< citation or source ....  [days]

    !$ACC UPDATE DEVICE(pheno_param_jsbach) ASYNC(1)
    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(pheno_param_jsbach%EG_SG)

  END SUBROUTINE init_pheno_parameters

#endif
END MODULE mo_pheno_parameters
