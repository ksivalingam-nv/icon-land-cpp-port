!> Definition of the abstract node type of the tile tree.
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
!>### Definition of abstract class and methods for nodes of the tile tree
!>
!> The module defines the abstract type [[t_jsb_tile_tree_node]] holding the structural information of a
!> node in the tile tree (parent, children, siblings, level and path), and the methods to query and
!> traverse this structure. It is the base type for tiles ([[mo_jsb_tile_class:t_jsb_tile_abstract]]).
!>

MODULE mo_jsb_tile_tree_class
#ifndef __NO_JSBACH__

  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_exception,          ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile_tree_node, t_jsb_tile_tree_node_p, MAX_TILE_TREE_DEPTH

  INTEGER, PARAMETER :: MAX_TILE_TREE_DEPTH = 10 !< Maximum nesting depth of the tile tree

  !> Abstract type for a node of the tile tree
  !>
  !> The type holds the structure of the tile tree (links to parent, children and siblings, level and
  !> path) and the methods to query it; tiles ([[mo_jsb_tile_class:t_jsb_tile_abstract]]) extend it.
  !> The link pointers defined here are of the node type. Tiles additionally keep the links to parent,
  !> first child and next sibling as pointers of the tile type (`parent_tile`, `first_child_tile`,
  !> `next_sibling_tile`): these point to the same tiles, but allow to access the tile components and
  !> methods without "SELECT TYPE".
  TYPE, ABSTRACT :: t_jsb_tile_tree_node
    CHARACTER(len=SHORT_NAME_LEN) :: name
    CHARACTER(len=50)             :: description
    INTEGER                       :: level = 0                       !< Level in hierarchy (top=1)
    INTEGER                       :: path(MAX_TILE_TREE_DEPTH) = -1  !< Position of the node in the tile tree
                                                                     !< (child index on each level, -1 below node)
    INTEGER, PRIVATE              :: no_of_children = -1
    !
    CLASS(t_jsb_tile_tree_node), POINTER :: parent       => NULL()
    CLASS(t_jsb_tile_tree_node), POINTER :: first_child  => NULL()
    CLASS(t_jsb_tile_tree_node), POINTER :: last_child   => NULL()
    CLASS(t_jsb_tile_tree_node), POINTER :: next_sibling => NULL()
    CLASS(t_jsb_tile_tree_node), POINTER :: prev_sibling => NULL()
  CONTAINS
    PROCEDURE :: Init_node
    PROCEDURE :: Get_parent
    PROCEDURE :: Get_parent_name
    PROCEDURE :: Get_no_of_children
    PROCEDURE :: Get_pos_of_child
    PROCEDURE :: Has_children
    PROCEDURE :: Has_siblings
    PROCEDURE :: Is_first_child
    PROCEDURE :: Is_last_child
    PROCEDURE :: Is_last_leaf
    PROCEDURE :: Print_node
    PROCEDURE :: Get_subtree_nodes
    PROCEDURE :: Is_equal     => Is_equal_node
    GENERIC   :: OPERATOR(==) => Is_equal
    PROCEDURE :: Is_not_equal => Is_not_equal_node
    GENERIC   :: OPERATOR(/=) => Is_not_equal
  END TYPE t_jsb_tile_tree_node

  !> A type to hold vectors of pointers to nodes of the tile tree
  TYPE t_jsb_tile_tree_node_p
    CLASS(t_jsb_tile_tree_node), POINTER :: p => NULL()
  END TYPE t_jsb_tile_tree_node_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_tree_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initialize a node and link it into the tile tree below `parent` (if present) as last child
  !>
  SUBROUTINE Init_node(this, name, description, parent)

    CLASS(t_jsb_tile_tree_node),          INTENT(inout), TARGET :: this
    CHARACTER(len=*),                INTENT(in)            :: name
    CHARACTER(len=*),                INTENT(in)            :: description
    CLASS(t_jsb_tile_tree_node), POINTER, OPTIONAL              :: parent

    this%name = TRIM(name)
    this%description = TRIM(description)

    IF (.NOT. PRESENT(parent)) THEN
      this%path(:) = -1
      this%path(1) = 1
      this%level = 1
      RETURN
    END IF

    this%path(:) = parent%path(:)
    this%level = parent%level + 1

    IF (.NOT. ASSOCIATED(parent%first_child)) THEN
      parent%first_child => this
      this%path(this%level) = 1
    END IF

    IF (ASSOCIATED(parent%last_child)) THEN
      parent%last_child%next_sibling => this
      this%prev_sibling              => parent%last_child
      this%path(this%level) = this%prev_sibling%path(this%level) + 1
    END IF

    parent%last_child => this
    this%parent       => parent

  END SUBROUTINE Init_node

  FUNCTION Get_parent(this) RESULT(return_ptr)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this
    CLASS(t_jsb_tile_tree_node), POINTER    :: return_ptr

    return_ptr => this%parent

  END FUNCTION Get_parent

  FUNCTION Get_parent_name(this) RESULT(return_value)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this
    CHARACTER(len=SHORT_NAME_LEN)      :: return_value

    IF (ASSOCIATED(this%parent)) THEN
      return_value = this%parent%name
    ELSE
      return_value = ''
    END IF

  END FUNCTION Get_parent_name

  LOGICAL FUNCTION Has_children(this)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    Has_children = ASSOCIATED(this%first_child)

  END FUNCTION Has_children

  LOGICAL FUNCTION Has_siblings(this)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    Has_siblings = ASSOCIATED(this%next_sibling) .OR. ASSOCIATED(this%prev_sibling)

  END FUNCTION Has_siblings

  LOGICAL FUNCTION Is_first_child(this)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    Is_first_child = .NOT. ASSOCIATED(this%prev_sibling)

  END FUNCTION Is_first_child

  LOGICAL FUNCTION Is_last_child(this)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    Is_last_child = .NOT. ASSOCIATED(this%next_sibling)

  END FUNCTION Is_last_child

  !> True if node is a leaf and the last child of a parent that is itself a last child.
  !> Note: This only identifies the last leaf of the tile tree if the tree has at most three levels.
  LOGICAL FUNCTION Is_last_leaf(this)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    IF (ASSOCIATED(this%parent)) THEN
      Is_last_leaf = .NOT. this%Has_children() .AND. this%Is_last_child() .AND. this%parent%Is_last_child()
    ELSE
      Is_last_leaf = .NOT. this%Has_children() .AND. this%Is_last_child()
    END IF

  END FUNCTION Is_last_leaf

  INTEGER FUNCTION Get_no_of_children(this)

    CLASS(t_jsb_tile_tree_node), INTENT(inout) :: this

    CLASS(t_jsb_tile_tree_node), POINTER :: current

    IF (this%no_of_children < 0) THEN
      this%no_of_children = 0
      current => this%first_child
      DO WHILE (ASSOCIATED(current))
        this%no_of_children = this%no_of_children + 1
        current => current%next_sibling
      END DO
    END IF

    Get_no_of_children = this%no_of_children

  END FUNCTION Get_no_of_children

  !> Position of this node among the children of its parent
  INTEGER FUNCTION Get_pos_of_child(this)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    Get_pos_of_child = this%path(this%level)

  END FUNCTION Get_pos_of_child

  SUBROUTINE Print_node(this)

    USE mo_util, ONLY: int2string

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this

    CHARACTER(len=4*MAX_TILE_TREE_DEPTH) :: c_path
    INTEGER :: i

    c_path = int2string(this%path(1))
    DO i=2,this%level
      c_path = TRIM(c_path)//'-'//int2string(this%path(i))
    END DO

    CALL message('Tile', this%name)
    CALL message('  Description', this%description)
    CALL message('  Tree path', TRIM(c_path))

  END SUBROUTINE Print_node

  ! ======================================================================================================= !
  !>
  !> Collect all nodes of the sub-tree rooted at this node (including this node) in pre-order (each node
  !> before its children) and in post-order (each node after its children). Siblings are always
  !> visited in the order they were added to the tile tree.
  !>
  SUBROUTINE Get_subtree_nodes(this, preorder, postorder)

    CLASS(t_jsb_tile_tree_node), TARGET,       INTENT(in)  :: this
    TYPE(t_jsb_tile_tree_node_p), ALLOCATABLE, INTENT(out) :: preorder(:)  !< Nodes in pre-order
    TYPE(t_jsb_tile_tree_node_p), ALLOCATABLE, INTENT(out) :: postorder(:) !< Nodes in post-order

    INTEGER :: nnodes, npre, npost

    ! First pass counts the nodes, second pass fills the lists
    CALL Walk(.FALSE.)
    nnodes = npre
    ALLOCATE(preorder(nnodes), postorder(nnodes))
    CALL Walk(.TRUE.)

  CONTAINS

    !> Depth-first walk through the sub-tree, adding each node when first reached (pre-order) and
    !> when left for good (post-order)
    SUBROUTINE Walk(l_fill)

      LOGICAL, INTENT(in) :: l_fill !< .FALSE.: only count the nodes (`npre`, `npost`), since the lists are
                                    !< not allocated yet; .TRUE.: also store the nodes in the lists

      CLASS(t_jsb_tile_tree_node), POINTER :: node

      npre  = 0
      npost = 0
      node => this
      DO WHILE (ASSOCIATED(node))
        npre = npre + 1
        IF (l_fill) preorder(npre)%p => node
        IF (ASSOCIATED(node%first_child)) THEN
          node => node%first_child
        ELSE
          ! Leaf: go up until a node with a next sibling is found; all nodes passed on the way are done
          DO WHILE (ASSOCIATED(node))
            npost = npost + 1
            IF (l_fill) postorder(npost)%p => node
            IF (ASSOCIATED(node, this)) THEN
              node => NULL()
            ELSE IF (ASSOCIATED(node%next_sibling)) THEN
              node => node%next_sibling
              EXIT
            ELSE
              node => node%parent
            END IF
          END DO
        END IF
      END DO

    END SUBROUTINE Walk

  END SUBROUTINE Get_subtree_nodes

  !> Only valid for two nodes in the same tile tree
  LOGICAL FUNCTION Is_equal_node(this, that)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this, that

    Is_equal_node = ALL(this%path == that%path)

  END FUNCTION Is_equal_node

  LOGICAL FUNCTION Is_not_equal_node(this, that)

    CLASS(t_jsb_tile_tree_node), INTENT(in) :: this, that

    Is_not_equal_node = ANY(this%path /= that%path)

  END FUNCTION Is_not_equal_node

#endif
END MODULE mo_jsb_tile_tree_class
