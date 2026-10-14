!> Contains methods for initialization of JSBACH usecases.
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
!! This module defines various model usecases that can be enabled individually via namelist.
!! A model usecase defines:
!!    * the tiles that are enabled in the model and their nested structure
!!    * the processes that are active on each tile or its leaves
!!    * the fraction of each tile
!!    * the land-cover-type of each tile

MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, PHENO_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_ !, ON_SUBTREE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE !, BARE_TYPE   ! define the type of the tile

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  ! ======================================================================================================= !
  !>initialize the usecase that determines the tiles and what processes are running where
  !>
  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase'

    IF (model%config%use_lakes) THEN
      CALL message(TRIM(routine), 'Initializing model usecase "'//TRIM(model%config%usecase)//'" with lakes')
    ELSE
      CALL message(TRIM(routine), 'Initializing model usecase "'//TRIM(model%config%usecase)//'"')
    END IF

    ! select usecase based on namelist option
    SELECT CASE (TRIM(model%config%usecase))

    ! jsbach physics
    CASE ('jsbach_lite')
      CALL init_usecase_lite(model)

    ! catch wrong usecase setting
    CASE DEFAULT
      CALL finish(TRIM(routine), 'No such usecase defined.')
    END SELECT

  END SUBROUTINE init_usecase

  ! ======================================================================================================= !
  !>init routine for usecase lite
  !>
  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |
    !  vegetation  glacier

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
        !$ACC ENTER DATA COPYIN(lake_tile, land_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
        !$ACC ENTER DATA COPYIN(land_tile)
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')
          !$ACC ENTER DATA COPYIN(glacier_tile)
        END IF

          veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
          !$ACC ENTER DATA COPYIN(veg_tile)

    NULLIFY(box_tile, land_tile, veg_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_lite

#endif
END MODULE mo_jsb_model_usecases
