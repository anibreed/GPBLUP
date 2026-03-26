module M_Hashtable
  use M_Kinds
  use M_Variables, only: LEN_STR
  implicit none
  
  ! 공통 해시 함수 및 소수 계산 기능을 제공하는 기본 모듈
  
  public :: htab_get_prime_size
  public :: htab_hash_string
  
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

end module M_Hashtable
