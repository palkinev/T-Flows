!==============================================================================!
  subroutine Cgns_Mod_Merge_Blocks(grid)
!------------------------------------------------------------------------------!
!  Merges blocks by merging common surface without changing structure of node  !
!  and cells connection  tables                                                !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, i, cnt_node
  real,    allocatable :: criterion(:) ! sorting criterion
  integer, allocatable :: old_seq(:)
  real,    allocatable :: x_new(:), y_new(:), z_new(:)
  logical, allocatable :: nodes_to_remove(:) ! marked duplicated nodes to remove
  real,    parameter   :: BIG   = 104651.    ! prime number
  real,    parameter   :: SMALL = 1.0e-6     ! precision for sorted criterion
!==============================================================================!

  print *, '# Merging blocks since they have duplicating nodes '
  print *, '# Hint: Join blocks in mesh builder to avoid any problems'
  print *, '# Old number of nodes: ', grid % n_nodes

  ! Allocate memory
  allocate(criterion      (grid % n_nodes));  criterion       = 0.
  allocate(old_seq        (grid % n_nodes));  old_seq         = 0
  allocate(nodes_to_remove(grid % n_nodes));  nodes_to_remove = .false.

  !--------------------------------------!
  !   Prescribe some sorting criterion   !
  !--------------------------------------!
  do n = 1, grid % n_nodes
    criterion(n) = grid % xn(n) + grid % yn(n) * BIG + grid % zn(n) * BIG ** 2
    old_seq(n) = n
  end do

  !----------------------------------------------------------------------------!
  !   Original block structure with duplicate nodes:                           !
  !        x  y  z                                                             !
  !    1:  1  1  1         13--------9    <->      5--------1                  !
  !    2:  1  1  0        /|        /|    <->     /|       /|                  !
  !    3:  1 -1  1       / |       / |    <->    / |      / |                  !
  !    4:  1 -1  0      14-|------10 |    <->   6--|-----2--|                  !
  !    5:  0  1  1      |  15     |  11   <->   |  7     |  3                  !
  !    6:  0  1  0      |  /      | /     <->   | /      | /                   !
  !    7:  0 -1  1      | /       |/      <->   |/       |/                    !
  !    8:  0 -1  0      16--------12      <->   8--------4                     !
  !    9:  0  1  1                                                             !
  !   10:  0  1  0        cell 2                    cell 1                     !
  !   11:  0 -1  1    y                                                        !
  !   12:  0 -1  0    ^                                                        !
  !   13: -1  1  1    |   ^ z                                                  !
  !   14: -1  1  0    |  /                                                     !
  !   15: -1 -1  1    | /                                                      !
  !   16: -1 -1  0    --------> x                                              !
  !                                                                            !
  !   cell 1: 1   2   3   4   5   6   7   8                                    !
  !   cell 2: 9  10  11  12  13  14  15  16                                    !
  !                                                                            !
  !   after function:                                                          !
  !        x  y  z                                                             !
  !    1:  1  1  1       9----------5   <->      5--------1                    !
  !    2:  1  1  0      / |        /|   <->     /|       /|                    !
  !    3:  1 -1  1     /  |       / |   <->    / |      / |                    !
  !    4:  1 -1  0    10--|------6--|   <->   6--|-----2--|                    !
  !    5:  0  1  1    |   11     |  7   <->   |  7     |  3                    !
  !    6:  0  1  0    |  /       | /    <->   | /      | /                     !
  !    7:  0 -1  1    | /        |/     <->   |/       |/                      !
  !    8:  0 -1  0    12---------8      <->   8--------4                       !
  !    9: -1  1  1                                                             !
  !   10: -1  1  0    cell 2                    cell 1                         !
  !   11: -1 -1  1                                                             !
  !   12: -1 -1  0                                                             !
  !                                                                            !
  !   cell 1: 1  2  3  4   5   6   7   8                                       !
  !   cell 2: 5  6  7  8   9  10  11  12                                       !
  !                                                                            !
  !   old_seq:                                                                 !
  !    1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16                   !
  !----------------------------------------------------------------------------!

  !----------------------------------!
  !   Sort nodes by this criterion   !
  !----------------------------------!
  print *, '# Sorting the nodes ...'
  call Sort_Real_Carry_Int_Heapsort(criterion(1), old_seq(1), grid % n_nodes)
  print *, '# ... done!'

  print *, '# cells before function'
  do c = 1, grid % n_cells
    print *, (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
  end do

  !----------------------------------------------------------------------------!
  !   old_seq now became:                                                      !
  !   16 12--8  4 14  6--10 2 15   7--11   3  13   5--9   1                    !
  !   "--" means that criterion has same value for these elements              !
  !----------------------------------------------------------------------------!
  cnt_node = 1
  if (verbose) print *,'Uniting nodes:'

  do n = 2, grid % n_nodes
    ! if node is unique
    if( .not. Approx(criterion(n), criterion(n-1), SMALL) ) then
      cnt_node = cnt_node + 1
    else ! if node is duplicated

      ! substitute dup. nodes to the lowest id
      do c = 1, grid % n_cells
        do i = 1, grid % cells_n_nodes(c)
          if (grid % cells_n(i,c) .eq.  old_seq(n)   .or. &
              grid % cells_n(i,c) .eq.  old_seq(n-1)      ) then

            grid % cells_n(i,c) = min(old_seq(n), old_seq(n-1))

          end if
        end do
      end do

      ! mark node to remove
      nodes_to_remove(max(old_seq(n-1), old_seq(n))) = .true.

      if (verbose) then
        print *,'----------------------'
        print '(a,i14,a,i14)','n: ', min(old_seq(n-1), old_seq(n)), '<-', &
                                     max(old_seq(n-1), old_seq(n))
        print '(a,es14.7,a,es14.7)','c: ', criterion(n), ' -', criterion(n-1)
        print '(a,es14.7,a,es14.7)','x: ', grid % xn(old_seq(n)), ' -', &
                                           grid % xn(old_seq(n-1))
        print '(a,es14.7,a,es14.7)','y: ', grid % yn(old_seq(n)), ' -', &
                                           grid % yn(old_seq(n-1))
        print '(a,es14.7,a,es14.7)','z: ', grid % zn(old_seq(n)), ' -', &
                                           grid % zn(old_seq(n-1))

      end if
    end if
  end do

  print *, '# New number of nodes: ', cnt_node

  !--------------------------------------------!
  !   Reconstruct new nodes and cells arrays   !
  !--------------------------------------------!

  allocate(x_new(1:cnt_node))
  allocate(y_new(1:cnt_node))
  allocate(z_new(1:cnt_node))

  cnt_node = 1
  do n = 1, grid % n_nodes
    if (.not. nodes_to_remove(n)) then ! if node is unique
      x_new(cnt_node) = grid % xn(n)
      y_new(cnt_node) = grid % yn(n)
      z_new(cnt_node) = grid % zn(n)

      cnt_node = cnt_node + 1

    else ! if node is duplicated

      ! shift nodes is cell
      do c = 1, grid % n_cells
        do i = 1, grid % cells_n_nodes(c)
          if (grid % cells_n(i,c) > cnt_node) then ! >= ???
            grid % cells_n(i,c) = grid % cells_n(i,c) - 1
          end if
        end do ! i
      end do ! c

    end if
  end do ! n
  cnt_node = cnt_node - 1

  !print *, '# nodes after function'
  !do n = 1, cnt_node
  !  print '(a,i14)' ,'n: ', n
  !  print '(a,f5.3,a,f5.3,a,f5.3)', &
  !    'x: ', x_new(n), ' y: ', y_new(n), ' z: ', z_new(n)
  !end do

  if (verbose) then
    print *, '# cells after function'
    do c = 1, grid % n_cells
      print *, (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  !------------------------------------------!
  !   Construct new nodes and cells arrays   !
  !------------------------------------------!

  !-------------------------------!
  !   Final touch: reinitialize   !
  !-------------------------------!
  grid % n_nodes = cnt_node

  deallocate(grid % xn)
  deallocate(grid % yn)
  deallocate(grid % zn)

  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)

  grid % xn(1: grid % n_nodes) = x_new(1: grid % n_nodes)
  grid % yn(1: grid % n_nodes) = y_new(1: grid % n_nodes)
  grid % zn(1: grid % n_nodes) = z_new(1: grid % n_nodes)

  deallocate(x_new)
  deallocate(y_new)
  deallocate(z_new)

  deallocate(criterion)
  deallocate(old_seq)
  deallocate(nodes_to_remove)

  end subroutine