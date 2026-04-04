module M_Hashtable
  use M_Kinds
  use M_Variables, only: LEN_STR
  implicit none
  
  ! 공통 해시 함수 및 소수 계산 기능을 제공하는 기본 모듈
  
  public :: htab_get_prime_size
  public :: htab_hash_string
  public :: htab_hash_string2
  public :: htab_hash_int
  
contains

  !**********************************************************************
  ! 소수 계산 - 해시 테이블 크기로 적합
  !**********************************************************************
  integer function htab_get_prime_size(min_size) result(prime)
    integer, intent(in) :: min_size
    integer :: candidate, i, is_prime
    
    candidate = max(101, min_size)
    if (mod(candidate, 2) == 0) candidate = candidate + 1
    
    do while (.true.)
      is_prime = 1
      do i = 3, int(sqrt(real(candidate))) + 1, 2
        if (mod(candidate, i) == 0) then
          is_prime = 0
          exit
        end if
      end do
      
      if (is_prime == 1) then
        prime = candidate
        return
      end if
      candidate = candidate + 2
    end do
  end function htab_get_prime_size

  !**********************************************************************
  ! 공통 문자열 해시 함수 (DJB2 알고리즘 변형)
  !**********************************************************************
  integer function htab_hash_string(key, table_size) result(hash_val)
    character(len=*), intent(in) :: key
    integer, intent(in) :: table_size
    integer :: i
    integer(kind=ki8) :: hash_i8
    
    if (len_trim(key) == 0) then
      hash_val = 0
      return
    end if

    hash_i8 = 5381_ki8  ! Magic number for DJB2
    do i = 1, len_trim(key)
      hash_i8 = mod(hash_i8 * 33_ki8 + int(iachar(key(i:i)), ki8), int(table_size, ki8))
    end do
    hash_val = int(abs(hash_i8))
  end function htab_hash_string

  !**********************************************************************
  ! Secondary string hash for double hashing probe step
  ! Returns value in [1, table_size-1] when table_size > 1
  !**********************************************************************
  integer function htab_hash_string2(key, table_size) result(step_val)
    character(len=*), intent(in) :: key
    integer, intent(in) :: table_size
    integer :: i
    integer(kind=ki8) :: h2_i8

    if (table_size <= 1) then
      step_val = 1
      return
    end if

    if (len_trim(key) == 0) then
      step_val = 1
      return
    end if

    h2_i8 = 0_ki8
    do i = 1, len_trim(key)
      h2_i8 = mod(h2_i8 * 131_ki8 + int(iachar(key(i:i)), ki8) + int(i, ki8), int(table_size - 1, ki8))
    end do
    step_val = 1 + int(h2_i8)
  end function htab_hash_string2

  !**********************************************************************
  ! 공통 정수 해시 함수
  !**********************************************************************
  integer function htab_hash_int(key, table_size) result(hash_val)
    integer, intent(in) :: key
    integer, intent(in) :: table_size
    integer(kind=ki8) :: x

    if (table_size <= 0) then
      hash_val = 0
      return
    end if

    x = int(key, ki8)
    x = ieor(x, ishft(x, -16))
    x = x * 2246822519_ki8
    x = ieor(x, ishft(x, -13))
    x = x * 3266489917_ki8
    x = ieor(x, ishft(x, -16))

    hash_val = int(mod(abs(x), int(table_size, ki8)))
  end function htab_hash_int

end module M_Hashtable
