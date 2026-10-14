!> Contains the initialization routines for the soil and snow energy process.
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
!>#### Initialize soil properties and initial conditions of the SSE process
!>
!> The module provides the subroutines used to initialize the soil and snow energy process.
!> The main subroutine is [[sse_init]], which calls the other initialization routines.
!>
!> Depending on the SSE configuration, initial and boundary-condition data is read from input
!> files in subroutine [[sse_read_init_vars]]. The data are then used to set the soil properties
!> in subroutine [[sse_init_bc]] and the initial conditions in subroutine [[sse_init_ic]].
!> After the last tile has been initialized, the temporary initialization data are released in
!> subroutine [[sse_finalize_init_vars]].
!>
MODULE mo_sse_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, soil_init_from_texture, read_bc_data

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sse_init, init_soil_properties

  REAL(wp), PARAMETER :: &
    & fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/), & !< FAO thermal diffusivity by soil class [m2 s-1]
    & fao_vol_hcap    (0:6) = (/2.25e+6_wp,1.93e+6_wp,2.10e+6_wp,2.25e+6_wp,2.36e+6_wp,2.48e+6_wp,-1._wp/)    !< FAO volumetric heat capacity by soil class [J m-3 K-1]

  ! Heat capacity [J m-3 K-1] (see table 2.2 in Reick et al. 2020,
  ! Reports on Earth System Science 240)
  REAL(wp), PARAMETER ::  heat_cap_sand = 1.93e+06_wp !< Volumetric heat capacity of sand [J m-3 K-1]
  REAL(wp), PARAMETER ::  heat_cap_silt = 2.25e+06_wp !< Volumetric heat capacity of silt [J m-3 K-1]
  REAL(wp), PARAMETER ::  heat_cap_clay = 2.48e+06_wp !< Volumetric heat capacity of clay [J m-3 K-1]
  REAL(wp), PARAMETER ::  heat_cap_oc = 2.5e+06_wp    !< Volumetric heat capacity of organic carbon [J m-3 K-1]

  ! Heat conductivity (see table 2.2 in Reick et al. 2020,
  ! Reports on Earth System Science 240)
  REAL(wp), PARAMETER ::  heat_cond_sand = 1.6791_wp !< Thermal conductivity of sand [J m-1 s-1 K-1]
  REAL(wp), PARAMETER ::  heat_cond_silt = 1.665_wp  !< Thermal conductivity of silt [J m-1 s-1 K-1]
  REAL(wp), PARAMETER ::  heat_cond_clay = 1.6616_wp !< Thermal conductivity of clay [J m-1 s-1 K-1]
  REAL(wp), PARAMETER ::  heat_cond_oc = 0.25_wp     !< Thermal conductivity of organic carbon [J m-1 s-1 K-1]

  !> Temporary container for auxiliary variables that are read from input files and used to
  !> initialize the SSE state.
  TYPE t_sse_init_vars
    REAL(wp), POINTER ::     &
      & fao         (:,:  ) => NULL(), & !< FAO soil classification index
      & vol_heat_cap(:,:  ) => NULL(), & !< Volumetric heat capacity of dry soil [J m-3 K-1]
      & heat_cond   (:,:  ) => NULL(), & !< Thermal conductivity of dry soil [J m-1 s-1 K-1]
      & fr_sand     (:,:  ) => Null(), & !< Sand fraction in upper soil []
      & fr_silt     (:,:  ) => Null(), & !< Silt fraction in upper soil []
      & fr_clay     (:,:  ) => Null(), & !< Clay fraction in upper soil []
      & fr_sand_deep(:,:  ) => Null(), & !< Sand fraction in deep soil []
      & fr_silt_deep(:,:  ) => Null(), & !< Silt fraction in deep soil []
      & fr_clay_deep(:,:  ) => Null(), & !< Clay fraction in deep soil []
      & t_srf_mon   (:,:,:) => NULL(), & !< Monthly climatological surface temperature [K]
      & ifs_tsoil   (:,:,:) => NULL(), & !< Soil temperatures from IFS [K]
      & ifs_tskin   (:,:  ) => NULL()    !< Skin temperatures from IFS [K]
  END TYPE t_sse_init_vars

  TYPE(t_sse_init_vars) :: sse_init_vars !< Temporary initialization data

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_init'

CONTAINS

  ! ====================================================================================================== !
  !>
  !>#### Initialize the soil and snow energy process
  !>
  !> The routine is called once the process memory has already been set up. It delegates to the
  !> routines that read input data, populate the boundary-condition fields, and initialize the
  !> soil and snow state according to the current configuration and tile structure.
  !>
  SUBROUTINE sse_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Current tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init'

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL sse_read_init_vars(tile)
    END IF

    CALL sse_init_bc(tile)
    CALL sse_init_ic(tile)

    IF (tile%Is_last_process_tile(SSE_)) THEN
      CALL sse_finalize_init_vars()
    END IF
  END SUBROUTINE sse_init

  ! ====================================================================================================== !
  !>
  !>#### Finalize the temporary SSE initialization state
  !>
  !> Release the temporary arrays that were created while reading and processing the SSE
  !> initialization data.
  SUBROUTINE sse_finalize_init_vars

    DEALLOCATE(                     &
      & sse_init_vars%vol_heat_cap, &
      & sse_init_vars%heat_cond     &
      & )

    IF (ASSOCIATED(sse_init_vars%fao))          DEALLOCATE(sse_init_vars%fao)
    IF (ASSOCIATED(sse_init_vars%fr_sand))      DEALLOCATE(sse_init_vars%fr_sand)
    IF (ASSOCIATED(sse_init_vars%fr_silt))      DEALLOCATE(sse_init_vars%fr_silt)
    IF (ASSOCIATED(sse_init_vars%fr_clay))      DEALLOCATE(sse_init_vars%fr_clay)
    IF (ASSOCIATED(sse_init_vars%fr_sand_deep)) DEALLOCATE(sse_init_vars%fr_sand_deep)
    IF (ASSOCIATED(sse_init_vars%fr_silt_deep)) DEALLOCATE(sse_init_vars%fr_silt_deep)
    IF (ASSOCIATED(sse_init_vars%fr_clay_deep)) DEALLOCATE(sse_init_vars%fr_clay_deep)

    IF (ASSOCIATED(sse_init_vars%t_srf_mon))    DEALLOCATE(sse_init_vars%t_srf_mon)

    IF (ASSOCIATED(sse_init_vars%ifs_tsoil))    DEALLOCATE(sse_init_vars%ifs_tsoil)
    IF (ASSOCIATED(sse_init_vars%ifs_tskin))    DEALLOCATE(sse_init_vars%ifs_tskin)

  END SUBROUTINE sse_finalize_init_vars

  ! ====================================================================================================== !
  !>
  !>#### Read initial and boundary-condition data for SSE
  !>
  !> Read the soil-energy input fields and climatological data needed to initialize the SSE
  !> state for the current tile.
  SUBROUTINE sse_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Current tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< Temporary input-data pointers

    TYPE(t_input_file) :: input_file !< Current boundary-condition input file

    TYPE(t_jsb_model), POINTER :: model  !< Current model instance
    TYPE(t_jsb_grid),  POINTER :: grid   !< Model horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e !< Soil energy vertical grid
    INTEGER :: nproma, nblks, isoil !< Horizontal block dimensions and IFS soil-layer index
    LOGICAL :: l_organic !< Whether organic soil properties are enabled

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_e => Get_vgrid('soil_depth_energy')

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Reading soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get config parameters
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    l_organic = dsl4jsb_Config(HYDRO_)%l_organic

    ALLOCATE(                                          &
      & sse_init_vars%vol_heat_cap(nproma, nblks    ), &
      & sse_init_vars%heat_cond   (nproma, nblks    )  &
      & )

    sse_init_vars%vol_heat_cap(:,:)   = missval
    sse_init_vars%heat_cond   (:,:)   = missval

    IF (      .NOT. (dsl4jsb_Config(SSE_)%l_heat_cap_map .AND. dsl4jsb_Config(SSE_)%l_heat_cond_map) &
      & .AND. .NOT. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      ALLOCATE(sse_init_vars%fao(nproma, nblks))
      sse_init_vars%fao(:,:)   = missval
    END IF

    IF (dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      ALLOCATE(                                          &
        & sse_init_vars%fr_sand     (nproma, nblks    ), &
        & sse_init_vars%fr_silt     (nproma, nblks    ), &
        & sse_init_vars%fr_clay     (nproma, nblks    ), &
        & sse_init_vars%fr_sand_deep(nproma, nblks    ), &
        & sse_init_vars%fr_silt_deep(nproma, nblks    ), &
        & sse_init_vars%fr_clay_deep(nproma, nblks    )  &
        & )

      sse_init_vars%fr_sand     (:,:)   = missval
      sse_init_vars%fr_silt     (:,:)   = missval
      sse_init_vars%fr_clay     (:,:)   = missval
      sse_init_vars%fr_sand_deep(:,:)   = missval
      sse_init_vars%fr_silt_deep(:,:)   = missval
      sse_init_vars%fr_clay_deep(:,:)   = missval
    END IF ! l_soil_texture

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SSE_)%bc_filename), model%grid_id)

    ! soil texture information
    !
    IF (dsl4jsb_Config(SSE_)%l_soil_texture) THEN

      !> Fraction of sand in the upper soil
      CALL read_bc_data(input_file, 'FR_SAND', tile%name, sse_init_vars%fr_sand, llimit=0._wp)

      !> Fraction of silt in the upper soil
      CALL read_bc_data(input_file, 'FR_SILT', tile%name, sse_init_vars%fr_silt, llimit=0._wp)

      !> Fraction of clay in the upper soil
      CALL read_bc_data(input_file, 'FR_CLAY', tile%name, sse_init_vars%fr_clay, llimit=0._wp)

      !> Fraction of sand in the deep soil
      CALL read_bc_data(input_file, 'SUB_FR_SAND', tile%name, sse_init_vars%fr_sand_deep, llimit=0._wp)

      !> Fraction of silt in the deep soil
      CALL read_bc_data(input_file, 'SUB_FR_SILT', tile%name, sse_init_vars%fr_silt_deep, llimit=0._wp)

      !> Fraction of clay in the deep soil
      CALL read_bc_data(input_file, 'SUB_FR_CLAY', tile%name, sse_init_vars%fr_clay_deep, llimit=0._wp)

      ! read 'bulk density' values from quincy specific input file
      ! calculate soil heat capacity [J m-3 K-1] incl. organic fraction
      CALL soil_init_from_texture( &
         heat_cap_sand, heat_cap_silt, heat_cap_clay, heat_cap_oc,                                          &
         sse_init_vars%fr_sand(:,:), sse_init_vars%fr_silt(:,:), sse_init_vars%fr_clay(:,:),                &
         sse_init_vars%fr_sand_deep(:,:), sse_init_vars%fr_silt_deep(:,:), sse_init_vars%fr_clay_deep(:,:), &
         sse_init_vars%vol_heat_cap(:,:))

      ! calculate soil heat conductivity incl. organic fraction
      CALL soil_init_from_texture( &
         heat_cond_sand, heat_cond_silt, heat_cond_clay, heat_cond_oc,                                      &
         sse_init_vars%fr_sand(:,:), sse_init_vars%fr_silt(:,:), sse_init_vars%fr_clay(:,:),                &
         sse_init_vars%fr_sand_deep(:,:), sse_init_vars%fr_silt_deep(:,:), sse_init_vars%fr_clay_deep(:,:), &
         sse_init_vars%heat_cond(:,:))
    ELSE
      IF (.NOT. dsl4jsb_Config(SSE_)%l_heat_cap_map .OR. .NOT. dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
        ! FAO
        CALL read_bc_data(input_file, 'fao', tile%name, sse_init_vars%fao, default=0._wp)
        sse_init_vars%fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), &
                                     MIN(REAL(UBOUND(fao_thermal_diff,1),wp), sse_init_vars%fao(:,:)))
      END IF

      IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
        ! Volumetric heat capacity of dry soil [J m-3 K-1]
        CALL read_bc_data(input_file, 'heat_capacity', tile%name, sse_init_vars%vol_heat_cap, &
          &  l_organic=l_organic, default=2.e6_wp)
      END IF

      IF (dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
        ! Heat conductivity of mineral soil
        CALL read_bc_data(input_file, 'heat_conductivity', tile%name, sse_init_vars%heat_cond, &
          &  l_organic=l_organic, default=0.2_wp)
      END IF

    END IF ! l_soil_texture

    CALL input_file%Close()

    IF (model%config%init_from_ifs) THEN

      ALLOCATE(sse_init_vars%ifs_tsoil(nproma,ifs_nsoil,nblks))
      sse_init_vars%ifs_tsoil(:,:,:) = missval
      ALLOCATE(sse_init_vars%ifs_tskin(nproma,nblks))
      sse_init_vars%ifs_tskin(:,:) = missval
      DO isoil=1,ifs_nsoil
        ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
          & variable_name='STL'//int2string(isoil))
        sse_init_vars%ifs_tsoil(:,isoil,:) = ptr_2D(:,:)
        DEALLOCATE(ptr_2D)
      END DO
      ptr_3D => model%config%ifs_input_file%Read_2d_time( &
        & variable_name='SKT', start_time_step=1, end_time_step=1)
      sse_init_vars%ifs_tskin(:,:) = ptr_3D(:,:,1)
      DEALLOCATE(ptr_3D)

    ELSE

      ALLOCATE(sse_init_vars%t_srf_mon(nproma, nblks, 12))
      sse_init_vars%t_srf_mon(:,:,:) = missval

      input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%ic_filename, model%grid_id)

      ! Read climatological surface temperature
      ptr_3D => input_file%Read_2d_time( &
        & variable_name='surf_temp')
      sse_init_vars%t_srf_mon(:,:,:) = MERGE(ptr_3D, 280._wp, ptr_3D > 0._wp)
      DEALLOCATE(ptr_3D)

      CALL input_file%Close()

    END IF

  END SUBROUTINE sse_read_init_vars

  ! ====================================================================================================== !
  !>
  !>#### Set SSE boundary conditions for the tiles
  !>
  !> Transfer the initialized soil-property fields into the tile memory so that the SSE state
  !> can be advanced during the subsequent model time steps.
  SUBROUTINE sse_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Current tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    INTEGER :: iblk !< Horizontal block index
    INTEGER :: nsoil !< Number of soil energy layers

    TYPE(t_jsb_model), POINTER :: model !< Current model instance
    TYPE(t_jsb_vgrid), POINTER :: soil_e !< Soil energy vertical grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_bc'

    model => get_model(tile%owner_model_id)
    soil_e => Get_vgrid('soil_depth_energy')

    nsoil = soil_e%n_levels

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Setting soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get process config
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)
    ! Get process memory
    dsl4jsb_Get_memory(SSE_)

    ! Volumetric heat capacity of dry soil [J m-3 K-1]
    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map .OR. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      dsl4jsb_var_ptr(SSE_, vol_heat_cap) = sse_init_vars%vol_heat_cap
    ELSE
      DO iblk=1,SIZE(sse_init_vars%fao,2)
        dsl4jsb_var_ptr(SSE_,vol_heat_cap)(:,iblk) = fao_vol_hcap(NINT(sse_init_vars%fao(:,iblk)))
      END DO
    END IF

    ! Heat conductivity of mineral soil
    IF (dsl4jsb_Config(SSE_)%l_heat_cond_map .OR. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      dsl4jsb_var_ptr(SSE_, heat_cond) = sse_init_vars%heat_cond
    ELSE
      ! Thermal diffusivity of soil from FAO
      DO iblk=1,SIZE(sse_init_vars%fao,2)
        dsl4jsb_var_ptr(SSE_,thermal_diffusivity)(:,iblk) = fao_thermal_diff(NINT(sse_init_vars%fao(:,iblk)))
      END DO
      dsl4jsb_var2D_onDomain(SSE_, heat_cond) = &
        & dsl4jsb_var2D_onDomain(SSE_, vol_heat_cap) * dsl4jsb_var2D_onDomain(SSE_, thermal_diffusivity)
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(SSE_,thermal_diffusivity)) ASYNC(1)
    END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SSE_,vol_heat_cap)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SSE_,heat_cond))
    !$ACC WAIT(1)

  END SUBROUTINE sse_init_bc

  ! ====================================================================================================== !
  !>
  !>#### Initialize the SSE initial conditions
  !>
  !> Populate the initial soil and surface temperatures in the SSE memory from either the
  !> climatological forcing or the IFS input data, depending on the model configuration.
  SUBROUTINE sse_init_ic(tile)

    USE mo_sse_process, ONLY: init_soil_temperature
    USE mo_util,        ONLY: ifs2soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Current tile

    TYPE(t_jsb_model), POINTER :: model  !< Current model instance
    TYPE(t_jsb_vgrid), POINTER :: soil_e !< Soil energy vertical grid

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)

    INTEGER :: iblk !< Horizontal block index

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_ic'

    model => get_model(tile%owner_model_id)
    soil_e => Get_vgrid('soil_depth_energy')

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

    ! Soil and surface temperatures
    !
    IF (model%config%init_from_ifs) THEN

      ! CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), &
      !   & dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:), &
      !   & ifs_sfc=sse_init_vars%ifs_tskin(:,:))
      CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:))

    ELSE
      DO iblk=1,SIZE(sse_init_vars%t_srf_mon, DIM=2)
        CALL init_soil_temperature(                   &
          & soil_e%levels(:),                         &
          & sse_init_vars%t_srf_mon(:,iblk,1:12),     &
          & dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,iblk))
      END DO
    END IF

    !$ACC UPDATE DEVICE(dsl4jsb_var3D_onDomain(SSE_,t_soil_sl)) ASYNC(1)

    ! IF (TRIM(tile%name) /= 'box') THEN
    IF (TRIM(tile%name) == 'land') THEN
      dsl4jsb_var2D_onDomain(SEB_,t)        = dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,1,:)
      dsl4jsb_var2D_onDomain(SEB_,t_unfilt) = dsl4jsb_var2D_onDomain(SEB_,t)
      dsl4jsb_var2D_onDomain(SEB_,t_srf)    = dsl4jsb_var2D_onDomain(SEB_,t)

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SEB_,t_srf)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SEB_,t_unfilt)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SEB_,t))
    END IF
    !$ACC WAIT(1)

    ! Note: called from model_init now
    ! CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  ! ====================================================================================================== !
  !>
  !>#### Initialize the soil thermal and hydrological properties
  !>
  !> Derive the initial soil thermal properties, snow properties, and supercooled-water limits
  !> from the current soil state and the selected SSE configuration.
  SUBROUTINE init_soil_properties(tile)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min, tmelt, alf, rhoh2o
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, &   ! vol_hcap_snow_min,
      &                                  vol_hcap_org_top, vol_hcap_org_below, &
      &                                  hcond_ice, hcond_snow, &         !, hcond_snow_min, hcond_org
      &                                  hcond_org_top, hcond_dry_org_top, hcond_org_below, hcond_dry_org_below
    USE mo_hydro_constants,        ONLY: vol_porosity_org_top, vol_porosity_org_below
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity, Get_liquid_max
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & vol_heat_cap,          & !< Volumetric heat capacity of dry soil [J m-3 K-1]
      & heat_cond,             & !< Thermal conductivity of dry soil [J m-1 s-1 K-1]
      & vol_heat_cap_snow,     & !< Volumetric heat capacity of snow [J m-3 K-1]
      & heat_cond_snow           !< Thermal conductivity of snow [J m-1 s-1 K-1]
    dsl4jsb_Real3D_onDomain :: &
      & soil_depth_sl,         & !< Soil layer depths [m]
      & wtr_soil_sl,           & !< Liquid water content by soil layer [m]
      & ice_soil_sl,           & !< Ice content by soil layer [m]
      & t_soil_sl,             & !< Soil layer temperature [K]
      & matric_pot_sl,         & !< Soil matric potential [m]
      & bclapp_sl,             & !< Clapp-Hornberger parameter []
      & fract_org_sl,          & !< Organic matter fraction by soil layer []
      & vol_heat_cap_sl,       & !< Dynamic volumetric heat capacity [J m-3 K-1]
      & heat_cond_sl,          & !< Dynamic thermal conductivity [J m-1 s-1 K-1]
      & vol_porosity_sl          !< Volumetric soil porosity []

    ! Locally allocated vectors
    !
    dsl4jsb_Real3D_onDomain :: &
      & liquid_max,            & !< Maximum liquid water content by soil layer [m]
      & fract_water,           & !< Volumetric water fraction by soil layer []
      & fract_ice,             & !< Volumetric ice fraction by soil layer []
      & vol_heat_cap_tmp,      & !< Temporary soil heat capacity [J m-3 K-1]
      & vol_heat_cap_org_tmp,  & !< Temporary organic heat capacity [J m-3 K-1]
      & heat_cond_tmp,         & !< Temporary soil thermal conductivity [J m-1 s-1 K-1]
      & vol_porosity_org_tmp,  & !< Temporary organic soil porosity []
      & hcond_org_tmp,         & !< Temporary organic thermal conductivity [J m-1 s-1 K-1]
      & hcond_dry_org_tmp        !< Temporary dry organic thermal conductivity [J m-1 s-1 K-1]

    TYPE(t_jsb_model), POINTER :: model !< Current model instance
    TYPE(t_jsb_vgrid), POINTER :: soil_e !< Soil energy vertical grid

    INTEGER  :: nsoil, nproma, nblks, ic, is, ib !< Soil layers, horizontal dimensions, and loop indices
    REAL(wp), POINTER :: soil_e_dz(:) !< Soil energy-layer thicknesses [m]
    REAL(wp) :: heat_cap_sl, ws_max_sl !< Soil-layer heat capacity and maximum water storage [J m-2 K-1, m]

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    soil_e_dz => soil_e%dz

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var3D_onDomain(SSE_, vol_heat_cap_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_, heat_cond_sl)
    dsl4jsb_Get_var2D_onDomain(SSE_, vol_heat_cap_snow)
    dsl4jsb_Get_var2D_onDomain(SSE_, heat_cond_snow)

    nproma = SIZE(heat_cond_snow,1)
    nblks  = SIZE(heat_cond_snow,2)

    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(3)
      DO is=1,nsoil
        DO ib=1,nblks
          DO ic=1,nproma
            vol_heat_cap_sl(ic,is,ib) = vol_hcap_ice
            heat_cond_sl   (ic,is,ib) = hcond_ice
          END DO
        END DO
      END DO
      !$ACC END LOOP
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO ib=1,nblks
        DO ic=1,nproma
          vol_heat_cap_snow(ic  ,ib) = vol_hcap_snow
          heat_cond_snow   (ic  ,ib) = hcond_snow
        END DO
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL
      RETURN
    END IF

    dsl4jsb_Get_var2D_onDomain(SSE_,   vol_heat_cap)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)

    IF (dsl4jsb_Config(SSE_)%l_freeze) THEN

      dsl4jsb_Get_var3D_onDomain(SSE_,   t_soil_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_, matric_pot_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)

      ALLOCATE(liquid_max (nproma, nsoil, nblks))
      !$ACC ENTER DATA CREATE(liquid_max)

      IF (dsl4jsb_Config(SSE_)%l_supercool) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1) PRIVATE(ws_max_sl)
        DO is = 1, nsoil
          DO ib = 1, nblks
            DO ic = 1, nproma
              ws_max_sl = soil_depth_sl(ic,is,ib) * vol_porosity_sl(ic,is,ib)
              liquid_max(ic,is,ib) = Get_liquid_max( &
                & t_soil_sl    (ic,is,ib), &
                & ws_max_sl              , &
                & matric_pot_sl(ic,is,ib), &
                & bclapp_sl    (ic,is,ib)  &
                & )
            END DO
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
        DO is = 1, nsoil
          DO ib = 1, nblks
            DO ic = 1, nproma
              liquid_max(ic,is,ib) = 0._wp
            END DO
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            heat_cap_sl = soil_e_dz(is) * vol_heat_cap(ic,ib)
            ice_soil_sl(ic,is,ib) = 0._wp
            IF (t_soil_sl(ic,is,ib) < tmelt .AND. wtr_soil_sl(ic,is,ib) > liquid_max(ic,is,ib)) THEN
              ice_soil_sl(ic,is,ib) = MIN(wtr_soil_sl(ic,is,ib) - liquid_max(ic,is,ib),                &
                &                         heat_cap_sl * (tmelt - t_soil_sl(ic,is,ib)) / (alf * rhoh2o) &
                &                        )
            END IF
            wtr_soil_sl(ic,is,ib) = wtr_soil_sl(ic,is,ib) - ice_soil_sl(ic,is,ib)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ! t_soil_sl(:,:,:) = t_soil_sl(:,:,:) + ice_soil_sl(:,:,:) * alf * rhoh2o / heat_cap_sl(:,:,:)
      ! IF (TRIM(tile%name) == 'land') THEN
      !   dsl4jsb_var_ptr(SEB_,t)(:,:)        = t_soil_sl(:,1,:)
      !   dsl4jsb_var_ptr(SEB_,t_unfilt)(:,:) = t_soil_sl(:,1,:)
      ! END IF

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(liquid_max)
      DEALLOCATE(liquid_max)

      ! dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = &
      !   & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl), dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl))
      ! !$acc update device(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))

    END IF

    dsl4jsb_Get_var2D_onDomain(SSE_,   heat_cond)

    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nproma, nsoil, nblks))
      !$ACC ENTER DATA CREATE(fract_org_sl)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            fract_org_sl(ic,is,ib) = 0._wp
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ALLOCATE(vol_heat_cap_tmp(nproma,nsoil,nblks), heat_cond_tmp(nproma,nsoil,nblks))
    !$ACC ENTER DATA CREATE(vol_heat_cap_tmp, heat_cond_tmp)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
    DO is = 1, nsoil
      DO ib = 1, nblks
        DO ic = 1, nproma
          vol_heat_cap_tmp(ic,is,ib) = vol_heat_cap(ic,ib)
          heat_cond_tmp   (ic,is,ib) = heat_cond(ic,ib)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nproma, nsoil, nblks), fract_ice(nproma, nsoil, nblks), &
        &      vol_heat_cap_org_tmp(nproma, nsoil, nblks))
      !$ACC ENTER DATA CREATE(fract_water, fract_ice, vol_heat_cap_org_tmp)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO ic = 1, nproma
        DO is = 1, nsoil
          DO ib = 1, nblks
            IF (soil_depth_sl(ic,is,ib) > 0._wp) THEN
              fract_water(ic,is,ib) = wtr_soil_sl(ic,is,ib) / soil_depth_sl(ic,is,ib)
              fract_ice  (ic,is,ib) = ice_soil_sl(ic,is,ib) / soil_depth_sl(ic,is,ib)
            ELSE
              fract_water(ic,is,ib) = 0._wp
              fract_ice  (ic,is,ib) = 0._wp
            END IF
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC UPDATE HOST(fract_water, fract_ice)
      !$ACC WAIT(1)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ib = 1, nblks
        DO ic = 1, nproma
           vol_heat_cap_org_tmp(ic,1,ib) = vol_hcap_org_top
        END DO
      END DO
      !$ACC END PARALLEL LOOP
     !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 2, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            vol_heat_cap_org_tmp(ic,is,ib) = vol_hcap_org_below
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            vol_heat_cap_sl (ic,is,ib) = calc_vol_heat_capacity( &
              &                                 vol_heat_cap_tmp (ic,is,ib)      , &
              &                                 vol_heat_cap_org_tmp(ic,is,ib)   , &
              &                                 vol_porosity_sl  (ic,is,ib)      , &
              &                                 fract_org_sl     (ic,is,ib)      , &
              &                                 fract_water      (ic,is,ib)      , &
              &                                 fract_ice        (ic,is,ib)      , &
              &                                 soil_depth_sl    (ic,is,ib)      , &
              &                                 soil_e_dz        (   is   )        &
              & )
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(fract_water, fract_ice, vol_heat_cap_org_tmp)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_org_tmp)
    ELSE

      !$ACC UPDATE DEVICE(vol_heat_cap) ASYNC(1)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib=1,nblks
          DO ic = 1, nproma
            vol_heat_cap_sl(ic,is,ib) = vol_heat_cap(ic,ib)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN
      ALLOCATE( &
        & vol_porosity_org_tmp(nproma,nsoil,nblks), &
        & hcond_org_tmp(nproma,nsoil,nblks),        &
        & hcond_dry_org_tmp(nproma,nsoil,nblks))
      !$ACC ENTER DATA CREATE(vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ib = 1, nblks
        DO ic = 1, nproma
          vol_porosity_org_tmp(ic,1,ib) = vol_porosity_org_top
          hcond_org_tmp(ic,1,ib)        = hcond_org_top
          hcond_dry_org_tmp(ic,1,ib)    = hcond_dry_org_top
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 2, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            vol_porosity_org_tmp(ic,is,ib) = vol_porosity_org_below
            hcond_org_tmp(ic,is,ib)        = hcond_org_below
            hcond_dry_org_tmp(ic,is,ib)    = hcond_dry_org_below
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 2, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            heat_cond_sl(ic,is,ib) = calc_thermal_conductivity( &
              &                                  heat_cond_tmp       (ic,is,ib)   , &
              &                                  wtr_soil_sl         (ic,is,ib)   , &
              &                                  ice_soil_sl         (ic,is,ib)   , &
              &                                  vol_porosity_sl     (ic,is,ib)   , &
              &                                  fract_org_sl        (ic,is,ib)   , &
              &                                  vol_porosity_org_tmp(ic,is,ib)   , &
              &                                  hcond_org_tmp       (ic,is,ib)   , &
              &                                  hcond_dry_org_tmp   (ic,is,ib)   , &
              &                                  soil_depth_sl       (ic,is,ib)   , &
              &                                  soil_e_dz           (   is   )     &
              & )
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp)
      DEALLOCATE(vol_porosity_org_tmp,hcond_org_tmp,hcond_dry_org_tmp)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib=1,nblks
          DO ic = 1, nproma
            heat_cond_sl(ic,is,ib) = heat_cond(ic,ib)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
      ! Use snow density values for fresh snow for initialization
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ib=1,nblks
        DO ic = 1, nproma
          vol_heat_cap_snow(ic,ib) = ci * dens_snow_min              ! Use specific heat of ice (Verseghy 1991)
          heat_cond_snow   (ic,ib) =   2.5e-6_wp * dens_snow_min**2._wp &
            &                        - 0.123E-3_wp *dens_snow_min + 0.024_wp ! Calonne et al. (2011)
          ! vol_heat_cap_snow(ic,ib) = vol_hcap_snow_min
          ! heat_cond_snow   (ic,ib) = hcond_snow_min
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ib=1,nblks
        DO ic = 1, nproma
          vol_heat_cap_snow(ic,ib) = vol_hcap_snow
          heat_cond_snow   (ic,ib) = hcond_snow
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC WAIT(1)

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_sse_init
