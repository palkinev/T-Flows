!==============================================================================!
  subroutine Grid_Mod_Create_Buffers(grid)
!------------------------------------------------------------------------------!
!   Reads: name.buf                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: sub, subo
  character(len=80) :: name_in
  integer           :: c1, c2, s, buf_cnt
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.      !
!   A receive buffer will be stored as aditional cells for each subdomain.     !
!   So each subdomain will have grid % n_cells cells, which entials physical   !
!   and buffer cells. It is handy to do it that way, because most of the       !
!   algorythms can remain the same as in sequential run. On the other hand,    !
!   a sending buffer has to be allocated in a new separate array called        !
!   simply buffer(). An additional array is needed to keep track of all the    !
!   indexes. That one is called buffind().                                     !
!------------------------------------------------------------------------------!

  if(n_proc < 2) return

  allocate (grid % comm % buff_s_cell(0:n_proc))  ! why from zero?
  allocate (grid % comm % buff_e_cell(0:n_proc))

  ! Initialize
  do sub = 0, n_proc
    grid % comm % buff_s_cell(sub) = grid % n_cells - grid % comm % n_buff_cells
    grid % comm % buff_e_cell(sub) = grid % n_cells - grid % comm % n_buff_cells
  end do

  ! Browse through other sub-domains
  do subo = 1, n_proc

    ! Connection with new domain will start ...
    ! ... where the domain with previous ended
    grid % comm % buff_s_cell(subo) = grid % comm % buff_e_cell(subo-1) + 1
    grid % comm % buff_e_cell(subo) = grid % comm % buff_e_cell(subo-1)

    ! Initialize buffer connections with this subdomain
    buf_cnt = 0

    if(subo .ne. this_proc) then

      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 > 0) then

          ! If c2 is a buffer face
          if( (grid % comm % cell_proc(c1) .eq. this_proc) .and.  &
              (grid % comm % cell_proc(c2) .eq. subo) ) then
            buf_cnt = buf_cnt + 1                ! increase buffer cell count
            grid % comm % buff_index(c2) = c1    ! buffer send index
          end if
          if( (grid % comm % cell_proc(c2) .eq. this_proc) .and.  &
              (grid % comm % cell_proc(c1) .eq. subo) ) then
            print *, '# FATAL ERROR! Shouldn''t be here at all!'
            print *, '# Exiting now!'
            stop
          end if
        end if  ! c2 > 0
      end do    ! through faces

    end if

    ! Set the end of the current buffer
    grid % comm % buff_e_cell(subo) =  &
    grid % comm % buff_e_cell(subo) + buf_cnt

  end do

  end subroutine
