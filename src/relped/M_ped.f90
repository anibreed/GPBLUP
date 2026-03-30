module M_ped
use M_param
use M_Hashtable, only: htab_get_prime_size, htab_hash_string
implicit none
 integer, public, parameter:: Male=2, Female=1

type, public :: CPED
     character(len=LEN_STR):: ASD(3)
     integer::  NASD(3)
     integer:: SEX, BYR, BDate, GEN, CG, OBS, NMat
     real(kind=r8):: F
     real(kind=r8):: EBV
end type CPED

type, public :: VPED
     character(len=LEN_STR):: ID
     integer::  NID
end type VPED

type, public :: ffh_t
     integer :: n_buckets = 0
     integer :: n_keys_stored = 0
     integer :: n_occupied = 0
     integer :: hash_mask = 0
     real(kind=r8) :: max_load_factor = 0.70_r8
     integer, allocatable :: flags(:)   ! 0=empty, 1=used, 2=deleted
     character(len=LEN_STR), allocatable :: keys(:)
     type(CPED), allocatable :: vals(:)
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
end type ffh_t

! Integer-view PED record to avoid duplicating full payload fields.
! CPED remains the source-of-truth payload; IPED stores only key links.
type, public :: IPED
     integer :: KEY_ID = missI
     integer :: SIRE_ID = missI
     integer :: DAM_ID = missI
     integer :: SLOT = -1
end type IPED

integer, save :: cache_size = 0
integer, save :: cache_n = 0
character(len=LEN_STR), allocatable, save :: cache_keys(:)
character(len=LEN_STR), allocatable, save :: cache_key_by_id(:)
integer, allocatable, save :: cache_used(:)
integer, allocatable, save :: cache_ffh_slot(:)
integer, allocatable, save :: cache_slot_by_id(:)
integer, allocatable, save :: cache_id_by_ffh_slot(:)

public :: ped_build_index_cache, ped_clear_index_cache
public :: ped_slot_from_key, ped_id_from_key, ped_slot_from_id
public :: cped_to_iped, iped_to_key

contains

subroutine init_cped_record(x)
     type(CPED), intent(inout) :: x

     x%ASD = missA
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

subroutine reset(h)
     class(ffh_t), intent(inout) :: h

     if (allocated(h%flags)) deallocate(h%flags)
     if (allocated(h%keys)) deallocate(h%keys)
     if (allocated(h%vals)) deallocate(h%vals)
     h%n_buckets = 0
     h%n_keys_stored = 0
     h%n_occupied = 0
     h%hash_mask = 0
end subroutine reset

logical function valid_index(h, i) result(ok)
     class(ffh_t), intent(in) :: h
     integer, intent(in) :: i

     ok = .false.
     if (.not. allocated(h%flags)) return
     if (i < 0 .or. i > ubound(h%flags, 1)) return
     ok = (h%flags(i) == 1)
end function valid_index

integer function get_index(h, key) result(ix)
     class(ffh_t), intent(in) :: h
     character(len=*), intent(in) :: key
     integer :: slot, step

     ix = -1
     if (.not. allocated(h%flags)) return
     if (len_trim(key) == 0) return

     slot = htab_hash_string(trim(key), h%n_buckets)
     do step = 0, h%n_buckets - 1
          if (h%flags(slot) == 0) return
          if (h%flags(slot) == 1 .and. trim(h%keys(slot)) == trim(key)) then
               ix = slot
               return
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do
end function get_index

subroutine ensure_capacity(h, extra)
     class(ffh_t), intent(inout) :: h
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
end subroutine ensure_capacity

subroutine rehash_to(h, requested_size)
     class(ffh_t), intent(inout) :: h
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

     new_size = htab_get_prime_size(max(101, requested_size))
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

     do i = lbound(old_flags, 1), ubound(old_flags, 1)
          if (old_flags(i) == 1) call h%ustore_value(old_keys(i), old_vals(i))
     end do

     if (allocated(old_flags)) deallocate(old_flags)
     if (allocated(old_keys)) deallocate(old_keys)
     if (allocated(old_vals)) deallocate(old_vals)
end subroutine rehash_to

subroutine store_value(h, key, val, status)
     class(ffh_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     type(CPED), intent(in) :: val
     integer, intent(out), optional :: status
     integer :: slot, step, first_deleted

     if (present(status)) status = -1
     if (len_trim(key) == 0) return

     call h%ensure_capacity(1)

     slot = htab_hash_string(trim(key), h%n_buckets)
     first_deleted = -1
     do step = 0, h%n_buckets - 1
          if (h%flags(slot) == 1) then
               if (trim(h%keys(slot)) == trim(key)) then
                    h%vals(slot) = val
                    if (present(status)) status = 0
                    return
               end if
          else if (h%flags(slot) == 2) then
               if (first_deleted == -1) first_deleted = slot
          else
               if (first_deleted /= -1) slot = first_deleted
               h%flags(slot) = 1
               h%keys(slot) = trim(key)
               h%vals(slot) = val
               h%n_keys_stored = h%n_keys_stored + 1
               if (first_deleted == -1) h%n_occupied = h%n_occupied + 1
               if (present(status)) status = 0
               return
          end if
          slot = modulo(slot + 1, h%n_buckets)
     end do

     call h%rehash_to(2 * h%n_buckets)
     call h%store_value(key, val, status)
end subroutine store_value

subroutine ustore_value(h, key, val)
     class(ffh_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     type(CPED), intent(in) :: val
     integer :: status

     call h%store_value(key, val, status)
end subroutine ustore_value

subroutine get_value(h, key, val, status)
     class(ffh_t), intent(in) :: h
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
     class(ffh_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     integer, intent(out) :: status
     integer :: ix

     status = -1
     ix = h%get_index(key)
     if (ix == -1) return
     h%flags(ix) = 2
     h%n_keys_stored = h%n_keys_stored - 1
     status = 0
end subroutine delete_key

subroutine udelete_key(h, key)
     class(ffh_t), intent(inout) :: h
     character(len=*), intent(in) :: key
     integer :: status

     call h%delete_key(key, status)
end subroutine udelete_key

subroutine ped_clear_index_cache()
     if (allocated(cache_keys)) deallocate(cache_keys)
     if (allocated(cache_key_by_id)) deallocate(cache_key_by_id)
     if (allocated(cache_used)) deallocate(cache_used)
     if (allocated(cache_ffh_slot)) deallocate(cache_ffh_slot)
     if (allocated(cache_slot_by_id)) deallocate(cache_slot_by_id)
     if (allocated(cache_id_by_ffh_slot)) deallocate(cache_id_by_ffh_slot)
     cache_size = 0
     cache_n = 0
end subroutine ped_clear_index_cache

subroutine ped_build_index_cache(h)
     type(ffh_t), intent(in) :: h
     integer :: i, slot, id

     call ped_clear_index_cache()
     if (h%n_keys_stored <= 0) return

     cache_size = htab_get_prime_size(max(101, int(real(h%n_keys_stored) * 1.3)))
     allocate(cache_keys(0:cache_size-1))
     allocate(cache_used(0:cache_size-1))
     allocate(cache_ffh_slot(0:cache_size-1))
     allocate(cache_slot_by_id(1:h%n_keys_stored))
     allocate(cache_key_by_id(1:h%n_keys_stored))
     allocate(cache_id_by_ffh_slot(0:h%n_buckets-1))

     cache_keys = missC
     cache_used = 0
     cache_ffh_slot = -1
     cache_slot_by_id = -1
     cache_key_by_id = missA
     cache_id_by_ffh_slot = missI
     cache_n = 0

     do i = 0, h%n_buckets-1
          if (.not. h%valid_index(i)) cycle
          slot = htab_hash_string(trim(h%keys(i)), cache_size)
          do
               if (cache_used(slot) == 0) then
                    cache_used(slot) = 1
                    cache_keys(slot) = h%keys(i)
                    cache_ffh_slot(slot) = i
                    cache_n = cache_n + 1
                    id = cache_n
                    cache_slot_by_id(id) = i
                    cache_key_by_id(id) = h%keys(i)
                    cache_id_by_ffh_slot(i) = id
                    exit
               else if (trim(cache_keys(slot)) == trim(h%keys(i))) then
                    cache_ffh_slot(slot) = i
                    if (cache_id_by_ffh_slot(i) == missI) then
                         cache_n = cache_n + 1
                         id = cache_n
                         cache_slot_by_id(id) = i
                         cache_key_by_id(id) = h%keys(i)
                         cache_id_by_ffh_slot(i) = id
                    end if
                    exit
               end if
               slot = modulo(slot + 1, cache_size)
          end do
     end do
end subroutine ped_build_index_cache

integer function ped_slot_from_key(key) result(slot_ffh)
     character(len=*), intent(in) :: key
     integer :: slot, probe

     slot_ffh = -1
     if (cache_size <= 0) return
     if (len_trim(key) == 0) return

     slot = htab_hash_string(trim(key), cache_size)
     do probe = 1, cache_size
          if (cache_used(slot) == 0) return
          if (trim(cache_keys(slot)) == trim(key)) then
               slot_ffh = cache_ffh_slot(slot)
               return
          end if
          slot = modulo(slot + 1, cache_size)
     end do
end function ped_slot_from_key

integer function ped_id_from_key(key) result(key_id)
     character(len=*), intent(in) :: key
     integer :: slot_ffh

     key_id = missI
     slot_ffh = ped_slot_from_key(key)
     if (slot_ffh < 0) return
     if (.not. allocated(cache_id_by_ffh_slot)) return
     if (slot_ffh > ubound(cache_id_by_ffh_slot, 1)) return
     key_id = cache_id_by_ffh_slot(slot_ffh)
end function ped_id_from_key

integer function ped_slot_from_id(key_id) result(slot_ffh)
     integer, intent(in) :: key_id

     slot_ffh = -1
     if (.not. allocated(cache_slot_by_id)) return
     if (key_id < 1 .or. key_id > size(cache_slot_by_id)) return
     slot_ffh = cache_slot_by_id(key_id)
end function ped_slot_from_id

subroutine cped_to_iped(src, dst)
     type(CPED), intent(in) :: src
     type(IPED), intent(out) :: dst

     dst%KEY_ID = ped_id_from_key(src%ASD(1))
     dst%SIRE_ID = src%NASD(2)
     if (dst%SIRE_ID == missI .and. src%ASD(2) /= missA) dst%SIRE_ID = ped_id_from_key(src%ASD(2))
     dst%DAM_ID = src%NASD(3)
     if (dst%DAM_ID == missI .and. src%ASD(3) /= missA) dst%DAM_ID = ped_id_from_key(src%ASD(3))
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

