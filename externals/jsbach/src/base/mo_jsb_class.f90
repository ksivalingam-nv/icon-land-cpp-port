!> Contains basic model structure
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
!> Contains basic model structure
!>
!> In most usecases one object of this class is used.
!> However, ICON-land can use more than one independent model at the same time having different configurations.
!> This is used for nested domains.
MODULE mo_jsb_class
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class, ONLY: t_jsb_model, t_jsb_model_m

  IMPLICIT NONE
  PRIVATE

  ! t_jsb needs to be public for the NEC compiler
  PUBLIC :: t_jsb, jsbach, get_model

  !> Type to contain a registry of models
  TYPE t_jsb
     INTEGER :: no_of_models                    !< Number of models - in most cases this would be 1
     TYPE(t_jsb_model_m), POINTER :: models(:)  !< List of models in use
  END TYPE t_jsb

  TYPE(t_jsb) :: jsbach !< The model registry. This is outdated, does not need to be called jsbach, it could be called "model registry" or the like

  INTERFACE get_model
     MODULE PROCEDURE get_model_by_id
  END INTERFACE get_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_class' !< Name of the module

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Get the pointer to a certain model from the registry given the model ID
  !>
  FUNCTION get_model_by_id(id) RESULT(model)

    INTEGER, INTENT(in)        :: id     !< ID of the queried model
    TYPE(t_jsb_model), POINTER :: model  !< Pointer to the queried model

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_model_by_id' !< Name of the function

    model => jsbach%models(id)%m

  END FUNCTION get_model_by_id

#endif
END MODULE mo_jsb_class
