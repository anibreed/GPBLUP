module link_heap
  implicit none
  public :: Heap, heap_init, heap_clear, insertL, getL, deleteL, heap_empty

  type :: Heap
    integer, allocatable :: data(:)
    integer :: size = 0
  end type Heap

contains

  subroutine heap_init(h, cap)
    type(Heap), intent(inout) :: h
    integer, intent(in), optional :: cap
    integer :: init_cap

    init_cap = 128
    if (present(cap)) init_cap = max(16, cap)

    if (allocated(h%data)) deallocate(h%data)
    allocate(h%data(init_cap))
    h%size = 0
  end subroutine heap_init

  subroutine heap_clear(h)
    type(Heap), intent(inout) :: h
    if (allocated(h%data)) deallocate(h%data)
    h%size = 0
  end subroutine heap_clear

  logical function heap_empty(h)
    type(Heap), intent(in) :: h
    heap_empty = (h%size <= 0)
  end function heap_empty

  subroutine heap_grow(h)
    type(Heap), intent(inout) :: h
    integer, allocatable :: new_data(:)
    integer :: new_cap

    if (.not. allocated(h%data)) then
      call heap_init(h, 128)
      return
    end if

    new_cap = max(2 * size(h%data), size(h%data) + 128)
    allocate(new_data(new_cap))
    new_data(1:h%size) = h%data(1:h%size)
    call move_alloc(new_data, h%data)
  end subroutine heap_grow

  subroutine insertL(h, id)
    type(Heap), intent(inout) :: h
    integer, intent(in) :: id
    integer :: i, parent

    if (id <= 0) return
    if (.not. allocated(h%data)) call heap_init(h, 128)
    if (h%size + 1 > size(h%data)) call heap_grow(h)

    h%size = h%size + 1
    i = h%size
    h%data(i) = id

    do while (i > 1)
      parent = i / 2
      if (h%data(parent) >= h%data(i)) exit
      call swap_int(h%data(parent), h%data(i))
      i = parent
    end do
  end subroutine insertL

  integer function getL(h) result(id)
    type(Heap), intent(in) :: h
    if (h%size <= 0) then
      id = 0
    else
      id = h%data(1)
    end if
  end function getL

  subroutine deleteL(h)
    type(Heap), intent(inout) :: h
    integer :: i, left, right, largest

    if (h%size <= 0) return

    h%data(1) = h%data(h%size)
    h%size = h%size - 1

    i = 1
    do
      left = 2 * i
      right = left + 1
      largest = i

      if (left <= h%size) then
        if (h%data(left) > h%data(largest)) largest = left
      end if
      if (right <= h%size) then
        if (h%data(right) > h%data(largest)) largest = right
      end if

      if (largest == i) exit
      call swap_int(h%data(i), h%data(largest))
      i = largest
    end do
  end subroutine deleteL

  subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: t
    t = a
    a = b
    b = t
  end subroutine swap_int

end module link_heap
