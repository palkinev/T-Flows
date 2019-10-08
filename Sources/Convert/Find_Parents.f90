!==============================================================================!
  subroutine Find_Parents(grid)
!------------------------------------------------------------------------------!
!   Looks boundary cells' parents for meshes in which they are not given       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include 'Cell_Numbering_Neu.f90'
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fn(6,4), j, n1, n2, c1, c2, cb, run
  integer              :: n_match, n_found_parents, n_cells_fraction
  integer              :: n_bnd_cells  ! near boundary cells
  integer              :: n_bnd_nodes  ! near boundary nodes
  integer, allocatable :: bnd_cells(:) ! near boundary cells
  logical, allocatable :: is_node_bnd(:)
!==============================================================================!

  print *, '#================================================='
  print *, '# Parent information not given in CGNS file!'
  print *, '# Looking for parents. This may take a few minutes'
  print *, '#-------------------------------------------------'

  allocate(is_node_bnd(grid % n_nodes))
  is_node_bnd(:) = .false.

  ! Mark all boundary nodes
  do c2 = -grid % n_bnd_cells, -1
    do n2 = 1, grid % cells_n_nodes(c2)  ! 3 or 4
      is_node_bnd( grid % cells_n(n2, c2) ) = .true.
    end do
  end do

  ! Count boundary nodes
  n_bnd_nodes = 0
  do n1 = 1, grid % n_nodes
    if( is_node_bnd(n1) ) n_bnd_nodes = n_bnd_nodes + 1
  end do
  print *, '# Total    nodes: ', grid % n_nodes
  print *, '# Boundary nodes: ', n_bnd_nodes

  ! Store near boundary cells
  do run = 1, 2
    n_bnd_cells = 0
    do c1 = 1, grid % n_cells
      do n1 = 1, grid % cells_n_nodes(c1)  ! 4 to 8
        if( is_node_bnd( grid % cells_n(n1, c1) ) ) then
          n_bnd_cells = n_bnd_cells + 1
          if(run .eq. 2) bnd_cells(n_bnd_cells) = c1
          exit
        end if
      end do
    end do
    if(run .eq. 1) allocate(bnd_cells(n_bnd_cells))
  end do
  print *, '# Total number of cells:         ', grid % n_cells
  print *, '# Number of near boundary cells: ', n_bnd_cells

  !----------------------------------------!
  !   Double loop through boundary cells   !
  !     and near boundary inside cells     !
  !----------------------------------------!
  n_found_parents  = 0
  n_cells_fraction = n_bnd_cells / 20
  do cb = 1, n_bnd_cells

    ! Print some statistics on the screen
    if(mod( cb, n_cells_fraction ) .eq. 0) then
      print '(a2, f5.0, a14)',              &
        ' #',                               &
        (100. * cb / (1.0*(n_bnd_cells))),  &
        ' % complete...'
    end if ! each 5%

    ! Get near boundary cell value
    c1 = bnd_cells(cb)

    if(grid % cells_n_nodes(c1) .eq. 4) fn = neu_tet
    if(grid % cells_n_nodes(c1) .eq. 5) fn = neu_pyr
    if(grid % cells_n_nodes(c1) .eq. 6) fn = neu_wed
    if(grid % cells_n_nodes(c1) .eq. 8) fn = neu_hex

    ! Brows through bundary cells
    do c2 = -grid % n_bnd_cells, -1

      !------------------------------!
      !   Number of matching nodes   !
      !------------------------------!
      do j = 1, 6  ! from 1 to 6th face (6 = hexahedron)
        n_match = 0
        do n1 = 1, 4  ! four is maximum number of cell nodes (quad)
          do n2 = 1, grid % cells_n_nodes(c2)  ! from 1 to 3/4th node in face
            if(fn(j, n1) > 0) then  ! if this node exists
              if(grid % cells_n(fn(j, n1), c1) .eq.  &
                 grid % cells_n(n2,        c2)) then
                n_match = n_match + 1
              end if
            end if
          end do ! n2
        end do ! n1
        if(n_match >= 3) then
          grid % cells_bnd_color(j, c1) = grid % bnd_cond % color(c2)
          n_found_parents = n_found_parents + 1
        end if ! n_match >= 3
      end do ! j

    end do
  end do

  end subroutine
