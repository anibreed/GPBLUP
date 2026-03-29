module H_PedRead
  use M_Kinds
  use M_Variables
  use M_Hashtable
  implicit none
  
  ! PEDInfo를 저장하는 해시 테이블 노드
  type, public :: PEDHashNode
    character(len=LEN_STR) :: ANIMAL_key   ! ARN을 키로 사용 (문자열)
    type(PEDInfo) :: ped_data           ! 실제 PED 정보
    type(PEDHashNode), pointer :: next => null()
  end type PEDHashNode
  
  ! PEDInfo 해시 테이블
  type, public :: PEDHashTable
    type(PEDHashNode), pointer :: table(:) => null()
    integer :: size = 0
    integer :: count = 0
  end type PEDHashTable
  
  public :: pht_create, pht_insert_with_key, pht_search, pht_delete, pht_free
  public :: pht_print_stats
  public :: pht_load_ped_file
  
contains

  !**********************************************************************
  ! 가계 파일(PED)을 읽어 해시 테이블에 로드 (범용 서브루틴)
  !**********************************************************************
  subroutine pht_load_ped_file(pht, n_rec, unit_fd, key_name, quiet)
    use M_ReadFile, only: fopen, readline, N_recf
    type(PEDHashTable), intent(out) :: pht
    integer, intent(out) :: n_rec, unit_fd
    character(len=*), intent(in) :: key_name
    logical, intent(in), optional :: quiet
    
    integer :: i, n, j
    integer(kind=ki4) :: XI(MAX_VAR)
    character(len=LEN_STR) :: XC(MAX_VAR)
    type(PEDInfo) :: ped_rec
    logical :: do_log

    do_log = .true.
    if (present(quiet)) do_log = .not. quiet
    
    n_rec = N_recf(PEDFile%FileName) - PEDFile%Header
    
    if (do_log) then
      print *, " M_PEDHashTable: PED 레코드를 읽는 중 (", trim(PEDFile%FileName), ")"
    end if
    call pht_create(pht, int(n_rec * 1.3), quiet=.not.do_log)
    
    unit_fd = fopen(trim(PEDFile%FileName))
    if (unit_fd < 0) return
    
    do i = 1, PEDFile%Header
      call readline(unit_fd, n, XC, XI, dlm_str=trim(PEDFile%Delim_char))
    end do
    
    do i = 1, n_rec
      call readline(unit_fd, n, XC, XI, dlm_str=trim(PEDFile%Delim_char))
      if (n < 0) exit
      call init_Ped(ped_rec)
      do j = 1, PEDFile%NVAR
         select case(trim(PEDFile%FieldName(j)))
         case("BREED"); ped_rec%BREED = trim(XC(PEDFile%FieldLoc(j)))
         case("ID","ANIMAL_ID"); ped_rec%ID = trim(XC(PEDFile%FieldLoc(j)))
         case("ARN","ANIMAL_ARN"); ped_rec%ARN = trim(XC(PEDFile%FieldLoc(j)))
         case("SIRE","SRN","SIRE_ID"); ped_rec%SIRE = trim(XC(PEDFile%FieldLoc(j)))
         case("DAM","DRN","DAM_ID"); ped_rec%DAM = trim(XC(PEDFile%FieldLoc(j)))
         case("BIRTH_DATE", "BIRHT_DATE", "BIRTH", "BDATE", "BIRTHDATE"); ped_rec%BDate = XI(PEDFile%FieldLoc(j))
         case("SEX"); ped_rec%SEX = XI(PEDFile%FieldLoc(j))
         case("LOCATION"); ped_rec%LOC = adjustl(trim(XC(PEDFile%FieldLoc(j))))
         end select
      end do
      call pht_insert_with_key(pht, ped_rec, key_name)
    end do
    close(unit=unit_fd)
    n_rec = pht%count
  end subroutine pht_load_ped_file

  !**********************************************************************
  ! PED 해시 테이블 생성
  !**********************************************************************
  subroutine pht_create(pht, table_size, quiet)
    type(PEDHashTable), intent(out) :: pht
    integer, intent(in) :: table_size
    logical, intent(in), optional :: quiet
    integer :: i, actual_size
    logical :: do_log

    do_log = .true.
    if (present(quiet)) do_log = .not. quiet
    actual_size = htab_get_prime_size(table_size)
    allocate(pht%table(0:actual_size-1))
    pht%size = actual_size
    pht%count = 0
    do i = 0, actual_size - 1
      nullify(pht%table(i)%next)
    end do
    if (do_log) then
      print '(A,I10,A,I10)', "  PED 해시 테이블 생성: 요청 크기=", table_size, &
                             "  실제 크기(소수)=", actual_size
    end if
  end subroutine pht_create

  !**********************************************************************
  ! PED 정보를 해시 테이블에 삽입 (키 지정)
  !**********************************************************************
  subroutine pht_insert_with_key(pht, ped, key_type)
    type(PEDHashTable), intent(inout) :: pht
    type(PEDInfo), intent(in) :: ped
    character(len=*), intent(in) :: key_type
    integer :: hash_idx
    type(PEDHashNode), pointer :: node, current
    character(len=LEN_STR) :: key_value
    if (index(trim(key_type), 'ID') > 0 .and. index(trim(key_type), 'ARN') == 0) then
      key_value = trim(ped%ID)
    else
      key_value = trim(ped%ARN)
    end if
    if (len_trim(key_value) == 0 .or. trim(key_value) == "0") return
    hash_idx = htab_hash_string(trim(key_value), pht%size)
    current => pht%table(hash_idx)%next
    do while (associated(current))
      if (trim(current%ANIMAL_key) == trim(key_value)) then
        current%ped_data = ped
        return
      end if
      current => current%next
    end do
    allocate(node)
    node%ANIMAL_key = trim(key_value)
    node%ped_data = ped
    node%next => pht%table(hash_idx)%next
    pht%table(hash_idx)%next => node
    pht%count = pht%count + 1
  end subroutine pht_insert_with_key

  !**********************************************************************
  ! PED 정보 검색
  !**********************************************************************
  logical function pht_search(pht, ANIMAL_key, ped) result(found)
    type(PEDHashTable), intent(in) :: pht
    character(len=*), intent(in) :: ANIMAL_key
    type(PEDInfo), intent(out) :: ped
    integer :: hash_idx
    type(PEDHashNode), pointer :: current
    found = .false.
    if (.not. associated(pht%table)) return
    hash_idx = htab_hash_string(trim(ANIMAL_key), pht%size)
    current => pht%table(hash_idx)%next
    do while (associated(current))
      if (trim(current%ANIMAL_key) == trim(ANIMAL_key)) then
        ped = current%ped_data
        found = .true.
        return
      end if
      current => current%next
    end do
  end function pht_search

  !**********************************************************************
  ! PED 정보 삭제
  !**********************************************************************
  logical function pht_delete(pht, ANIMAL_key) result(deleted)
    type(PEDHashTable), intent(inout) :: pht
    character(len=*), intent(in) :: ANIMAL_key
    integer :: hash_idx
    type(PEDHashNode), pointer :: current, prev
    deleted = .false.
    if (.not. associated(pht%table)) return
    hash_idx = htab_hash_string(trim(ANIMAL_key), pht%size)
    prev => null()
    current => pht%table(hash_idx)%next
    do while (associated(current))
      if (trim(current%ANIMAL_key) == trim(ANIMAL_key)) then
        if (associated(prev)) then
          prev%next => current%next
        else
          pht%table(hash_idx)%next => current%next
        end if
        deallocate(current)
        pht%count = pht%count - 1
        deleted = .true.
        return
      end if
      prev => current
      current => current%next
    end do
  end function pht_delete

  !**********************************************************************
  ! 메모리 해제
  !**********************************************************************
  subroutine pht_free(pht)
    type(PEDHashTable), intent(inout) :: pht
    integer :: i
    type(PEDHashNode), pointer :: current, next_node
    if (.not. associated(pht%table)) return
    do i = 0, pht%size - 1
      current => pht%table(i)%next
      do while (associated(current))
        next_node => current%next
        deallocate(current)
        current => next_node
      end do
    end do
    deallocate(pht%table)
    pht%size = 0
    pht%count = 0
  end subroutine pht_free

  !**********************************************************************
  ! 통계 출력
  !**********************************************************************
  subroutine pht_print_stats(pht)
    type(PEDHashTable), intent(in) :: pht
    integer :: i, chain_len, max_chain
    type(PEDHashNode), pointer :: current
    max_chain = 0
    do i = 0, pht%size - 1
      chain_len = 0
      current => pht%table(i)%next
      do while (associated(current))
        chain_len = chain_len + 1
        current => current%next
      end do
      max_chain = max(max_chain, chain_len)
    end do
    print *, "Count: ", pht%count, " / Size: ", pht%size, " / Max Chain: ", max_chain
  end subroutine pht_print_stats

end module H_PedRead
