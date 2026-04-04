module M_ped
use M_param
implicit none
 integer, public, parameter:: Male=2, Female=1

type, public :: CPED
     character(len=LEN_STR):: ASD(3)
     character(len=LEN_STR):: BREED, ID, ARN, SIRE, DAM, LOC
     integer::  NASD(3)
     integer:: SEX, BYR, BDate, GEN, CG, OBS, NMat
     real(kind=r8):: F
     real(kind=r8):: EBV
end type CPED

type, public :: VPED
     character(len=LEN_STR):: ID
     integer::  NID
end type VPED

type, public :: str_hash_t
     integer :: n_buckets = 0
     integer :: n_items = 0
     integer :: n_deleted = 0
     real(kind=r8) :: max_load_factor = 0.80_r8
     integer, allocatable :: state(:)  ! 0=empty, 1=used, 2=deleted
     character(len=LEN_STR), allocatable :: keys(:)
     integer, allocatable :: vals(:)
contains
     procedure :: reset => str_hash_reset
     procedure :: get => str_hash_get
     procedure :: put => str_hash_put
     procedure :: remove => str_hash_remove
     procedure, private :: ensure_capacity => str_hash_ensure_capacity
     procedure, private :: rehash_to => str_hash_rehash_to
end type str_hash_t

type, public :: int_hash_t
     integer :: n_buckets = 0
     integer :: n_items = 0
     integer :: n_deleted = 0
     real(kind=r8) :: max_load_factor = 0.80_r8
     integer, allocatable :: state(:)  ! 0=empty, 1=used, 2=deleted
     integer, allocatable :: keys(:)
     integer, allocatable :: vals(:)
contains
     procedure :: reset => int_hash_reset
     procedure :: get => int_hash_get
     procedure :: put => int_hash_put
     procedure :: remove => int_hash_remove
     procedure, private :: ensure_capacity => int_hash_ensure_capacity
     procedure, private :: rehash_to => int_hash_rehash_to
end type int_hash_t

type, public :: ped_hash_t
     integer :: n_buckets = 0
     integer :: n_keys_stored = 0
     integer :: n_occupied = 0
     integer :: hash_mask = 0
     real(kind=r8) :: max_load_factor = 0.70_r8
     integer, allocatable :: flags(:)   ! 0=empty, 1=used, 2=deleted
     character(len=LEN_STR), allocatable :: keys(:)
     type(CPED), allocatable :: vals(:)
     type(str_hash_t) :: key_to_slot
contains
     procedure :: get_index
     procedure :: valid_index
     procedure :: store_value
     procedure :: ustore_value
     procedure :: get_value
     procedure :: delete_key
     procedure :: udelete_key
     procedure :: reset
     procedure, private :: ensure_capacity
     procedure, private :: rehash_to
end type ped_hash_t

! Integer-view PED record to avoid duplicating full payload fields.
! CPED remains the source-of-truth payload; IPED stores only key links.
type, public :: IPED
     integer :: KEY_ID = missI
     integer :: SIRE_ID = missI
     integer :: DAM_ID = missI
     integer :: SLOT = -1
end type IPED

integer, save :: cache_n = 0
character(len=LEN_STR), allocatable, save :: cache_key_by_id(:)
integer, allocatable, save :: cache_slot_by_id(:)
integer, allocatable, save :: cache_id_by_ffh_slot(:)
type(str_hash_t), save :: cache_key_to_slot
type(int_hash_t), save :: cache_id_to_slot
type(int_hash_t), save :: cache_slot_to_id

public :: ped_build_index_cache, ped_clear_index_cache
public :: ped_slot_from_key, ped_id_from_key, ped_slot_from_id
public :: cped_to_iped, iped_to_key
public :: cped_sync_links

contains

subroutine init_cped_record(x)
     type(CPED), intent(inout) :: x

     x%ASD = missA
     x%BREED = ''
     x%ID = ''
     x%ARN = ''
     x%SIRE = ''
     x%DAM = ''
     x%LOC = ''
     x%NASD = missI
     x%SEX = missI
     x%BYR = missI
     x%BDate = missI
     x%GEN = -1
     x%CG = missI
     x%OBS = missI
     x%NMat = missI
     x%F = missR
     x%EBV = missR
end subroutine init_cped_record

subroutine cped_sync_links(x)
     type(CPED), intent(inout) :: x

     if (len_trim(x%ARN) > 0 .and. trim(x%ARN) /= trim(missA)) then
          x%ASD(1) = trim(x%ARN)
     else if (x%ASD(1) /= missA) then
          x%ARN = trim(x%ASD(1))
     end if

     if (len_trim(x%ID) == 0 .or. trim(x%ID) == trim(missA)) then
          if (len_trim(x%ARN) > 0 .and. trim(x%ARN) /= trim(missA)) then
               x%ID = trim(x%ARN)
          else if (x%ASD(1) /= missA) then
               x%ID = trim(x%ASD(1))
          end if
     end if

     if (len_trim(x%SIRE) > 0 .and. trim(x%SIRE) /= trim(missA)) then
          x%ASD(2) = trim(x%SIRE)
     else if (x%ASD(2) /= missA) then
          x%SIRE = trim(x%ASD(2))
     end if

     if (len_trim(x%DAM) > 0 .and. trim(x%DAM) /= trim(missA)) then
          x%ASD(3) = trim(x%DAM)
     else if (x%ASD(3) /= missA) then
          x%DAM = trim(x%ASD(3))
     end if
end subroutine cped_sync_links

integer function hash_index_from_key(key, n_buckets) result(ix)
     character(len=*), intent(in) :: key
     integer, intent(in) :: n_buckets
     integer(kind=8) :: h
     integer :: i

     ix = 0
     if (n_buckets <= 0) return

     h = 1469598103934665603_8
     do i = 1, len_trim(key)
          h = ieor(h, int(iachar(key(i:i)), kind=8))
          h = h * 1099511628211_8
     end do
     ix = modulo(int(modulo(h, int(huge(1), kind=8))), n_buckets)
end function hash_index_from_key

integer function hash_index_from_int(key, n_buckets) result(ix)
     integer, intent(in) :: key
     integer, intent(in) :: n_buckets
     integer(kind=8) :: h

     ix = 0
     if (n_buckets <= 0) return

     h = int(key, kind=8)
     h = h * 6364136223846793005_8 + 1442695040888963407_8
     ix = modulo(int(modulo(h, int(huge(1), kind=8))), n_buckets)
end function hash_index_from_int

subroutine str_hash_reset(h)
     class(str_hash_t), intent(inout) :: h

     if (allocated(h%state)) deallocate(h%state)
     if (allocated(h%keys)) deallocate(h%keys)
     if (allocated(h%vals)) deallocate(h%vals)
     h%n_buckets = 0
     h%n_items = 0
     h%n_deleted = 0
end subroutine str_hash_reset

subroutine str_hash_rehash_to(h, requested_size)
     class(str_hash_t), intent(inout) :: h
     integer, intent(in) :: requested_size
     integer :: new_size, i
     integer, allocatable :: old_state(:), old_vals(:)
     character(len=LEN_STR), allocatable :: old_keys(:)

     new_size = max(64, requested_size)

     if (allocated(h%state)) then
          allocate(old_state(0:ubound(h%state, 1)))
          allocate(old_keys(0:ubound(h%keys, 1)))
          allocate(old_vals(0:ubound(h%vals, 1)))
          old_state = h%state
          old_keys = h%keys
          old_vals = h%vals
          deallocate(h%state, h%keys, h%vals)
     else
          allocate(old_state(0:-1))
          allocate(old_keys(0:-1))
          allocate(old_vals(0:-1))
     end if

     allocate(h%state(0:new_size-1))
     allocate(h%keys(0:new_size-1))
     allocate(h%vals(0:new_size-1))
     h%state = 0
     h%keys = missC
     h%vals = missI
     h%n_buckets = new_size
     h%n_items = 0
     h%n_deleted = 0

     do i = lbound(old_state, 1), ubound(old_state, 1)
          if (old_state(i) == 1) call h%put(old_keys(i), old_vals(i))
     end do

     if (allocated(old_state)) deallocate(old_state)
     if (allocated(old_keys)) deallocate(old_keys)
     if (allocated(old_vals)) deallocate(old_vals)
end subroutine str_hash_rehash_to

subroutine str_hash_ensure_capacity(h, extra)
     class(str_hash_t), intent(inout) :: h
     integer, intent(in), optional :: extra
     integer :: need

     need = 1
     if (present(extra)) need = max(1, extra)

     if (.not. allocated(h%state)) then
          call h%rehash_to(max(256, need * 4))
          return
     end if

     if (real(h%n_items + h%n_deleted + need, r8) >= real(h%n_buckets, r8) * h%max_load_factor) then
          call h%rehash_to(max(2 * h%n_buckets, h%n_buckets + need * 4))
          return
     end if

     if (h%n_deleted > h%n_items) call h%rehash_to(h%n_buckets)
end subroutine str_hash_ensure_capacity

subroutine str_hash_get(h, key, val, found)
     class(str_hash_t), intent(in) :: h
     character(len=*), intent(in) :: key
     integer, intent(out) :: val
     logical, intent(out) :: found
     integer :: slot, step

     val = missI
     found = .false.
     if (.not. allocated(h%state)) return
     if (len_trim(key) == 0) return

     slot = hash_index_from_key(trim(key), h%n_buckets)
     do step = 0, h%n_buckets - 1
          if (h%state(slot) == 0) return
          if (h%state(slot) == 1 .and. trim(h%keys(slot)) == trim(key)) then
               val = h%vals(slot)
               found = .true.
               return
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do
end subroutine str_hash_get

subroutine str_hash_put(h, key, val)
     class(str_hash_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     integer, intent(in) :: val
     integer :: slot, step, first_deleted

     if (len_trim(key) == 0) return
     call h%ensure_capacity(1)

     slot = hash_index_from_key(trim(key), h%n_buckets)
     first_deleted = -1
     do step = 0, h%n_buckets - 1
          if (h%state(slot) == 1) then
               if (trim(h%keys(slot)) == trim(key)) then
                    h%vals(slot) = val
                    return
               end if
          else if (h%state(slot) == 2) then
               if (first_deleted == -1) first_deleted = slot
          else
               if (first_deleted /= -1) slot = first_deleted
               h%state(slot) = 1
               h%keys(slot) = trim(key)
               h%vals(slot) = val
               h%n_items = h%n_items + 1
               if (first_deleted /= -1) h%n_deleted = max(0, h%n_deleted - 1)
               return
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do

     call h%rehash_to(2 * h%n_buckets)
     call h%put(key, val)
end subroutine str_hash_put

subroutine str_hash_remove(h, key, found)
     class(str_hash_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     logical, intent(out), optional :: found
     integer :: slot, step
     logical :: ok

     ok = .false.
     if (.not. allocated(h%state)) then
          if (present(found)) found = ok
          return
     end if

     slot = hash_index_from_key(trim(key), h%n_buckets)
     do step = 0, h%n_buckets - 1
          if (h%state(slot) == 0) exit
          if (h%state(slot) == 1 .and. trim(h%keys(slot)) == trim(key)) then
               h%state(slot) = 2
               h%keys(slot) = missC
               h%vals(slot) = missI
               h%n_items = h%n_items - 1
               h%n_deleted = h%n_deleted + 1
               ok = .true.
               exit
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do

     if (present(found)) found = ok
end subroutine str_hash_remove

subroutine int_hash_reset(h)
     class(int_hash_t), intent(inout) :: h

     if (allocated(h%state)) deallocate(h%state)
     if (allocated(h%keys)) deallocate(h%keys)
     if (allocated(h%vals)) deallocate(h%vals)
     h%n_buckets = 0
     h%n_items = 0
     h%n_deleted = 0
end subroutine int_hash_reset

subroutine int_hash_rehash_to(h, requested_size)
     class(int_hash_t), intent(inout) :: h
     integer, intent(in) :: requested_size
     integer :: new_size, i
     integer, allocatable :: old_state(:), old_keys(:), old_vals(:)

     new_size = max(64, requested_size)

     if (allocated(h%state)) then
          allocate(old_state(0:ubound(h%state, 1)))
          allocate(old_keys(0:ubound(h%keys, 1)))
          allocate(old_vals(0:ubound(h%vals, 1)))
          old_state = h%state
          old_keys = h%keys
          old_vals = h%vals
          deallocate(h%state, h%keys, h%vals)
     else
          allocate(old_state(0:-1))
          allocate(old_keys(0:-1))
          allocate(old_vals(0:-1))
     end if

     allocate(h%state(0:new_size-1))
     allocate(h%keys(0:new_size-1))
     allocate(h%vals(0:new_size-1))
     h%state = 0
     h%keys = missI
     h%vals = missI
     h%n_buckets = new_size
     h%n_items = 0
     h%n_deleted = 0

     do i = lbound(old_state, 1), ubound(old_state, 1)
          if (old_state(i) == 1) call h%put(old_keys(i), old_vals(i))
     end do

     if (allocated(old_state)) deallocate(old_state)
     if (allocated(old_keys)) deallocate(old_keys)
     if (allocated(old_vals)) deallocate(old_vals)
end subroutine int_hash_rehash_to

subroutine int_hash_ensure_capacity(h, extra)
     class(int_hash_t), intent(inout) :: h
     integer, intent(in), optional :: extra
     integer :: need

     need = 1
     if (present(extra)) need = max(1, extra)

     if (.not. allocated(h%state)) then
          call h%rehash_to(max(256, need * 4))
          return
     end if

     if (real(h%n_items + h%n_deleted + need, r8) >= real(h%n_buckets, r8) * h%max_load_factor) then
          call h%rehash_to(max(2 * h%n_buckets, h%n_buckets + need * 4))
          return
     end if

     if (h%n_deleted > h%n_items) call h%rehash_to(h%n_buckets)
end subroutine int_hash_ensure_capacity

subroutine int_hash_get(h, key, val, found)
     class(int_hash_t), intent(in) :: h
     integer, intent(in) :: key
     integer, intent(out) :: val
     logical, intent(out) :: found
     integer :: slot, step

     val = missI
     found = .false.
     if (.not. allocated(h%state)) return

     slot = hash_index_from_int(key, h%n_buckets)
     do step = 0, h%n_buckets - 1
          if (h%state(slot) == 0) return
          if (h%state(slot) == 1 .and. h%keys(slot) == key) then
               val = h%vals(slot)
               found = .true.
               return
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do
end subroutine int_hash_get

subroutine int_hash_put(h, key, val)
     class(int_hash_t), intent(inout) :: h
     integer, intent(in) :: key
     integer, intent(in) :: val
     integer :: slot, step, first_deleted

     call h%ensure_capacity(1)

     slot = hash_index_from_int(key, h%n_buckets)
     first_deleted = -1
     do step = 0, h%n_buckets - 1
          if (h%state(slot) == 1) then
               if (h%keys(slot) == key) then
                    h%vals(slot) = val
                    return
               end if
          else if (h%state(slot) == 2) then
               if (first_deleted == -1) first_deleted = slot
          else
               if (first_deleted /= -1) slot = first_deleted
               h%state(slot) = 1
               h%keys(slot) = key
               h%vals(slot) = val
               h%n_items = h%n_items + 1
               if (first_deleted /= -1) h%n_deleted = max(0, h%n_deleted - 1)
               return
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do

     call h%rehash_to(2 * h%n_buckets)
     call h%put(key, val)
end subroutine int_hash_put

subroutine int_hash_remove(h, key, found)
     class(int_hash_t), intent(inout) :: h
     integer, intent(in) :: key
     logical, intent(out), optional :: found
     integer :: slot, step
     logical :: ok

     ok = .false.
     if (.not. allocated(h%state)) then
          if (present(found)) found = ok
          return
     end if

     slot = hash_index_from_int(key, h%n_buckets)
     do step = 0, h%n_buckets - 1
          if (h%state(slot) == 0) exit
          if (h%state(slot) == 1 .and. h%keys(slot) == key) then
               h%state(slot) = 2
               h%keys(slot) = missI
               h%vals(slot) = missI
               h%n_items = h%n_items - 1
               h%n_deleted = h%n_deleted + 1
               ok = .true.
               exit
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do

     if (present(found)) found = ok
end subroutine int_hash_remove

subroutine reset(h)
     class(ped_hash_t), intent(inout) :: h

     if (allocated(h%flags)) deallocate(h%flags)
     if (allocated(h%keys)) deallocate(h%keys)
     if (allocated(h%vals)) deallocate(h%vals)
     call h%key_to_slot%reset()
     h%n_buckets = 0
     h%n_keys_stored = 0
     h%n_occupied = 0
     h%hash_mask = 0
end subroutine reset

logical function valid_index(h, i) result(ok)
     class(ped_hash_t), intent(in) :: h
     integer, intent(in) :: i

     ok = .false.
     if (.not. allocated(h%flags)) return
     if (i < 0 .or. i > ubound(h%flags, 1)) return
     ok = (h%flags(i) == 1)
end function valid_index

integer function get_index(h, key) result(ix)
     class(ped_hash_t), intent(in) :: h
     character(len=*), intent(in) :: key
     logical :: found

     ix = -1
     if (.not. allocated(h%flags)) return
     if (len_trim(key) == 0) return

     call h%key_to_slot%get(trim(key), ix, found)
     if (.not. found) ix = -1
     if (ix < 0) return
     if (.not. h%valid_index(ix)) ix = -1
end function get_index

subroutine ensure_capacity(h, extra)
     class(ped_hash_t), intent(inout) :: h
     integer, intent(in), optional :: extra
     integer :: need

     need = 1
     if (present(extra)) need = max(1, extra)

     if (.not. allocated(h%flags)) then
          call h%rehash_to(max(1024, need * 2))
          return
     end if

     if (real(h%n_occupied + need, r8) >= real(h%n_buckets, r8) * h%max_load_factor) then
          call h%rehash_to(max(2 * h%n_buckets, h%n_buckets + need * 2))
     end if

     call h%key_to_slot%ensure_capacity(need)
end subroutine ensure_capacity

subroutine rehash_to(h, requested_size)
     class(ped_hash_t), intent(inout) :: h
     integer, intent(in) :: requested_size
     integer :: new_size, i
     integer, allocatable :: old_flags(:)
     character(len=LEN_STR), allocatable :: old_keys(:)
     type(CPED), allocatable :: old_vals(:)

     if (allocated(h%flags)) then
          allocate(old_flags(0:ubound(h%flags, 1)))
          allocate(old_keys(0:ubound(h%keys, 1)))
          allocate(old_vals(0:ubound(h%vals, 1)))
          old_flags = h%flags
          old_keys = h%keys
          old_vals = h%vals
          deallocate(h%flags, h%keys, h%vals)
     else
          allocate(old_flags(0:-1))
          allocate(old_keys(0:-1))
          allocate(old_vals(0:-1))
     end if

     new_size = max(101, requested_size)
     allocate(h%flags(0:new_size-1))
     allocate(h%keys(0:new_size-1))
     allocate(h%vals(0:new_size-1))
     h%flags = 0
     h%keys = missC
     do i = 0, new_size - 1
          call init_cped_record(h%vals(i))
     end do

     h%n_buckets = new_size
     h%hash_mask = new_size - 1
     h%n_keys_stored = 0
     h%n_occupied = 0

     call h%key_to_slot%reset()
     call h%key_to_slot%ensure_capacity(max(256, new_size / 2))

     do i = lbound(old_flags, 1), ubound(old_flags, 1)
          if (old_flags(i) == 1) call h%ustore_value(old_keys(i), old_vals(i))
     end do

     if (allocated(old_flags)) deallocate(old_flags)
     if (allocated(old_keys)) deallocate(old_keys)
     if (allocated(old_vals)) deallocate(old_vals)
end subroutine rehash_to

subroutine store_value(h, key, val, status)
     class(ped_hash_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     type(CPED), intent(in) :: val
     integer, intent(out), optional :: status
     integer :: slot
     logical :: found

     if (present(status)) status = -1
     if (len_trim(key) == 0) return

     call h%ensure_capacity(1)

     call h%key_to_slot%get(trim(key), slot, found)
     if (found .and. h%valid_index(slot)) then
          h%vals(slot) = val
          h%keys(slot) = trim(key)
          if (present(status)) status = 0
          return
     end if

     do slot = 0, h%n_buckets - 1
          if (h%flags(slot) /= 1) then
               if (h%flags(slot) == 0) h%n_occupied = h%n_occupied + 1
               h%flags(slot) = 1
               h%keys(slot) = trim(key)
               h%vals(slot) = val
               h%n_keys_stored = h%n_keys_stored + 1
               call h%key_to_slot%put(trim(key), slot)
               if (present(status)) status = 0
               return
          end if
     end do

     call h%rehash_to(2 * h%n_buckets)
     call h%store_value(key, val, status)
end subroutine store_value

subroutine ustore_value(h, key, val)
     class(ped_hash_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     type(CPED), intent(in) :: val
     integer :: status

     call h%store_value(key, val, status)
end subroutine ustore_value

subroutine get_value(h, key, val, status)
     class(ped_hash_t), intent(in) :: h
     character(len=*), intent(in) :: key
     type(CPED), intent(out) :: val
     integer, intent(out) :: status
     integer :: ix

     call init_cped_record(val)
     ix = h%get_index(key)
     if (ix == -1) then
          status = -1
     else
          val = h%vals(ix)
          status = 0
     end if
end subroutine get_value

subroutine delete_key(h, key, status)
     class(ped_hash_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     integer, intent(out) :: status
     integer :: ix
     logical :: found

     status = -1
     ix = h%get_index(key)
     if (ix == -1) return
     h%flags(ix) = 2
     h%keys(ix) = missC
     call init_cped_record(h%vals(ix))
     h%n_keys_stored = h%n_keys_stored - 1
     call h%key_to_slot%remove(trim(key), found)
     status = 0
end subroutine delete_key

subroutine udelete_key(h, key)
     class(ped_hash_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     integer :: status

     call h%delete_key(key, status)
end subroutine udelete_key

subroutine ped_clear_index_cache()
     if (allocated(cache_key_by_id)) deallocate(cache_key_by_id)
     if (allocated(cache_slot_by_id)) deallocate(cache_slot_by_id)
     if (allocated(cache_id_by_ffh_slot)) deallocate(cache_id_by_ffh_slot)
     call cache_key_to_slot%reset()
     call cache_id_to_slot%reset()
     call cache_slot_to_id%reset()
     cache_n = 0
end subroutine ped_clear_index_cache

subroutine ped_build_index_cache(h)
     type(ped_hash_t), intent(in) :: h
     integer :: i, id

     call ped_clear_index_cache()
     if (h%n_keys_stored <= 0) return

     cache_n = h%n_keys_stored
     allocate(cache_slot_by_id(1:h%n_keys_stored))
     allocate(cache_key_by_id(1:h%n_keys_stored))
     allocate(cache_id_by_ffh_slot(0:h%n_buckets-1))

     cache_slot_by_id = -1
     cache_key_by_id = missA
     cache_id_by_ffh_slot = missI

     call cache_key_to_slot%ensure_capacity(h%n_keys_stored)
     call cache_id_to_slot%ensure_capacity(h%n_keys_stored)
     call cache_slot_to_id%ensure_capacity(h%n_keys_stored)

     cache_n = 0
     do i = 0, h%n_buckets-1
          if (.not. h%valid_index(i)) cycle
          cache_n = cache_n + 1
          id = cache_n
          cache_slot_by_id(id) = i
          cache_key_by_id(id) = h%keys(i)
          cache_id_by_ffh_slot(i) = id
          call cache_key_to_slot%put(h%keys(i), i)
          call cache_id_to_slot%put(id, i)
          call cache_slot_to_id%put(i, id)
     end do
end subroutine ped_build_index_cache

integer function ped_slot_from_key(key) result(slot_ffh)
     character(len=*), intent(in) :: key
     logical :: found

     slot_ffh = -1
     if (cache_n <= 0) return
     if (len_trim(key) == 0) return

     call cache_key_to_slot%get(trim(key), slot_ffh, found)
     if (.not. found) slot_ffh = -1
end function ped_slot_from_key

integer function ped_id_from_key(key) result(key_id)
     character(len=*), intent(in) :: key
     integer :: slot_ffh
     logical :: found

     key_id = missI
     slot_ffh = ped_slot_from_key(key)
     if (slot_ffh < 0) return

     call cache_slot_to_id%get(slot_ffh, key_id, found)
     if (found) return

     if (.not. allocated(cache_id_by_ffh_slot)) return
     if (slot_ffh > ubound(cache_id_by_ffh_slot, 1)) return
     key_id = cache_id_by_ffh_slot(slot_ffh)
end function ped_id_from_key

integer function ped_slot_from_id(key_id) result(slot_ffh)
     integer, intent(in) :: key_id
     logical :: found

     slot_ffh = -1
     call cache_id_to_slot%get(key_id, slot_ffh, found)
     if (found) return

     if (.not. allocated(cache_slot_by_id)) return
     if (key_id < 1 .or. key_id > size(cache_slot_by_id)) return
     slot_ffh = cache_slot_by_id(key_id)
end function ped_slot_from_id

subroutine cped_to_iped(src, dst)
     type(CPED), intent(in) :: src
     type(IPED), intent(out) :: dst
     type(CPED) :: tmp
     character(len=LEN_STR) :: key_main

     tmp = src
     call cped_sync_links(tmp)

     key_main = trim(tmp%ARN)
     if (len_trim(key_main) == 0 .or. trim(key_main) == trim(missA)) key_main = trim(tmp%ID)
     if (len_trim(key_main) == 0 .or. trim(key_main) == trim(missA)) key_main = trim(tmp%ASD(1))

     dst%KEY_ID = ped_id_from_key(key_main)
     dst%SIRE_ID = src%NASD(2)
     if (dst%SIRE_ID == missI) then
          if (len_trim(tmp%SIRE) > 0 .and. trim(tmp%SIRE) /= trim(missA)) then
               dst%SIRE_ID = ped_id_from_key(tmp%SIRE)
          else if (tmp%ASD(2) /= missA) then
               dst%SIRE_ID = ped_id_from_key(tmp%ASD(2))
          end if
     end if
     dst%DAM_ID = src%NASD(3)
     if (dst%DAM_ID == missI) then
          if (len_trim(tmp%DAM) > 0 .and. trim(tmp%DAM) /= trim(missA)) then
               dst%DAM_ID = ped_id_from_key(tmp%DAM)
          else if (tmp%ASD(3) /= missA) then
               dst%DAM_ID = ped_id_from_key(tmp%ASD(3))
          end if
     end if
     dst%SLOT = ped_slot_from_id(dst%KEY_ID)
end subroutine cped_to_iped

subroutine iped_to_key(src, key, status)
     type(IPED), intent(in) :: src
     character(len=*), intent(out) :: key
     integer, intent(out) :: status
     integer :: slot_ffh

     key = missA
     status = -1
     slot_ffh = ped_slot_from_id(src%KEY_ID)
     if (slot_ffh < 0) return
     if (.not. allocated(cache_key_by_id)) return
     if (src%KEY_ID < 1 .or. src%KEY_ID > size(cache_key_by_id)) return
     key = cache_key_by_id(src%KEY_ID)
     status = 0
end subroutine iped_to_key

end module M_ped

