module M_HashSort
  !---------------------------------------------------------------------------
  ! 범용 희소행렬 유틸리티 모듈
  !
  !   해시 테이블 (단일 정수 키, open addressing)
  !   해시 → CSR 변환
  !   다중 키 quicksort
  !   원본 호환 루틴 (hashvr, hashia1, clusterr)
  !---------------------------------------------------------------------------
  use M_Kinds, only: ki4, ki8, r8
  implicit none
  private

  ! === 타입 ===
  public :: hash_table_t, sparse_csr_t

  ! === 해시 테이블 API ===
  public :: hash_init, hash_free
  public :: hash_add            ! (row,col)에 값 누적
  public :: hash_get            ! (row,col) 값 조회
  public :: hash_find           ! 내부 슬롯 탐색
  public :: hash_to_csr         ! 해시 → CSR 변환
  public :: encode_key, decode_key

  ! === CSR API ===
  public :: csr_free

  ! === 정렬/군집 ===
  public :: sortqr

  ! === 원본 호환 ===
  public :: hashvr, clusterr, hashia1

  integer, parameter, public :: HASH_NOT_FOUND = 0
  integer, parameter, public :: HASH_FULL      = -1

  integer, parameter :: PROBE_STEP = 17
  integer, parameter :: MAX_PROBE  = 5000
  integer, parameter :: NP = 10
  integer, parameter :: prime(NP) = [433,1039,53,761,197,31,887,389,277,97]

  !--- 해시 테이블 타입 (단일 정수 키) ---
  type :: hash_table_t
    integer              :: tblsz = 0
    integer              :: nmax  = 0
    integer              :: count = 0
    integer(ki8), allocatable :: keys(:)
    real(r8),     allocatable :: vals(:)
  end type hash_table_t

  !--- CSR 희소행렬 타입 ---
  type :: sparse_csr_t
    integer              :: n   = 0
    integer              :: nnz = 0
    integer, allocatable :: ia(:)
    integer, allocatable :: ja(:)
    real(r8), allocatable :: a(:)
  end type sparse_csr_t

contains

  !==========================================================================
  !  키 인코딩/디코딩
  !==========================================================================
  pure integer(ki8) function encode_key(row, col, nmax)
    integer, intent(in) :: row, col, nmax
    encode_key = int(row - 1, ki8) * int(nmax, ki8) + int(col, ki8)
  end function encode_key

  pure subroutine decode_key(key, nmax, row, col)
    integer(ki8), intent(in)  :: key
    integer,      intent(in)  :: nmax
    integer,      intent(out) :: row, col
    row = int((key - 1_ki8) / int(nmax, ki8)) + 1
    col = int(mod(key - 1_ki8, int(nmax, ki8))) + 1
  end subroutine decode_key

  !==========================================================================
  !  hash_init / hash_free
  !==========================================================================
  subroutine hash_init(ht, tblsz, nmax)
    type(hash_table_t), intent(out) :: ht
    integer, intent(in) :: tblsz, nmax
    ht%tblsz = tblsz;  ht%nmax = nmax;  ht%count = 0
    allocate(ht%keys(tblsz), ht%vals(tblsz))
    ht%keys = 0_ki8;  ht%vals = 0.0_r8
  end subroutine hash_init

  subroutine hash_free(ht)
    type(hash_table_t), intent(inout) :: ht
    if (allocated(ht%keys)) deallocate(ht%keys)
    if (allocated(ht%vals)) deallocate(ht%vals)
    ht%tblsz = 0;  ht%count = 0
  end subroutine hash_free

  !==========================================================================
  !  hash_find — 슬롯 탐색 (mode=1: 없으면 삽입, mode=0: 조회만)
  !==========================================================================
  function hash_find(ht, key, mode) result(slot)
    type(hash_table_t), intent(inout) :: ht
    integer(ki8), intent(in) :: key
    integer, intent(in)      :: mode
    integer :: slot, addr, k

    addr = int(mod(abs(key * 2654435761_ki8), int(ht%tblsz, ki8))) + 1
    do k = 1, MAX_PROBE
      if (ht%keys(addr) == 0_ki8) then
        if (mode == 1) then
          ht%keys(addr) = key;  ht%vals(addr) = 0.0_r8
          ht%count = ht%count + 1;  slot = addr
        else
          slot = HASH_NOT_FOUND
        end if
        return
      end if
      if (ht%keys(addr) == key) then;  slot = addr;  return;  end if
      addr = mod(addr + PROBE_STEP - 1, ht%tblsz) + 1
    end do
    slot = HASH_FULL
  end function hash_find

  !==========================================================================
  !  hash_add — (row,col)에 값 누적 (대칭: 상삼각 자동)
  !==========================================================================
  subroutine hash_add(ht, row, col, val)
    type(hash_table_t), intent(inout) :: ht
    integer, intent(in)  :: row, col
    real(r8), intent(in) :: val
    integer(ki8) :: key
    integer :: r, c, slot
    if (row <= col) then;  r = row;  c = col
    else;                  r = col;  c = row;  end if
    key = encode_key(r, c, ht%nmax)
    slot = hash_find(ht, key, mode=1)
    if (slot > 0) ht%vals(slot) = ht%vals(slot) + val
  end subroutine hash_add

  !==========================================================================
  !  hash_get — (row,col) 값 조회. 없으면 0 반환
  !==========================================================================
  function hash_get(ht, row, col) result(val)
    type(hash_table_t), intent(inout) :: ht
    integer, intent(in) :: row, col
    real(r8) :: val
    integer(ki8) :: key
    integer :: r, c, slot
    if (row <= col) then;  r = row;  c = col
    else;                  r = col;  c = row;  end if
    key = encode_key(r, c, ht%nmax)
    slot = hash_find(ht, key, mode=0)
    if (slot > 0) then;  val = ht%vals(slot)
    else;                 val = 0.0_r8;  end if
  end function hash_get

  !==========================================================================
  !  hash_to_csr — 해시 → CSR 변환
  !==========================================================================
  subroutine hash_to_csr(ht, n, csr)
    type(hash_table_t), intent(inout) :: ht
    integer, intent(in) :: n
    type(sparse_csr_t), intent(out) :: csr
    integer :: rowcnt(n), i, row, col, k
    integer(ki8) :: key

    csr%n = n;  rowcnt = 0
    do i = 1, ht%tblsz
      key = ht%keys(i)
      if (key == 0_ki8) cycle
      call decode_key(key, ht%nmax, row, col)
      if (row < 1 .or. row > n .or. col < 1 .or. col > n) then
        ht%keys(i) = 0_ki8;  cycle
      end if
      rowcnt(row) = rowcnt(row) + 1
    end do

    allocate(csr%ia(n + 1))
    csr%ia(1) = 1
    do i = 1, n;  csr%ia(i+1) = csr%ia(i) + rowcnt(i);  end do
    csr%nnz = csr%ia(n + 1) - 1
    allocate(csr%ja(csr%nnz), csr%a(csr%nnz))

    rowcnt = 0
    do i = 1, ht%tblsz
      key = ht%keys(i)
      if (key == 0_ki8) cycle
      call decode_key(key, ht%nmax, row, col)
      k = csr%ia(row) + rowcnt(row)
      csr%ja(k) = col;  csr%a(k) = ht%vals(i)
      rowcnt(row) = rowcnt(row) + 1
    end do
  end subroutine hash_to_csr

  subroutine csr_free(csr)
    type(sparse_csr_t), intent(inout) :: csr
    if (allocated(csr%ia)) deallocate(csr%ia)
    if (allocated(csr%ja)) deallocate(csr%ja)
    if (allocated(csr%a))  deallocate(csr%a)
    csr%n = 0;  csr%nnz = 0
  end subroutine csr_free

  !==========================================================================
  !  sortqr — 원본 인터페이스 (a, n, m, ix, iy, col, icol)
  !==========================================================================
  subroutine sortqr(a, n, m, ix, iy, col, icol)
    integer,  intent(in)    :: n, m, ix, iy, icol
    integer,  intent(in)    :: col(:)
    real(r8), intent(inout) :: a(:,:)
    integer, parameter :: STKMAX = 64
    integer  :: stkL(STKMAX), stkR(STKMAX), sp
    integer  :: lo, hi, i, j, mid, k, j2, a1, a2
    real(r8) :: x(icol), x1
    logical  :: adv

    if (ix <= 1) return
    sp = 1;  stkL(1) = 1;  stkR(1) = ix
    do while (sp > 0)
      lo = stkL(sp);  hi = stkR(sp);  sp = sp - 1
      do while (lo < hi)
        i = lo;  j = hi;  mid = lo + ishft(hi - lo, -1)
        do k = 1, icol;  x(k) = a(col(k), mid);  end do
        do
          scan_l: do
            adv = .false.
            do k = 1, icol
              a1 = a(col(k), i);  a2 = x(k)
              if (a1 < a2) then;  adv = .true.;  exit
              else if (a1 > a2) then;  exit;  end if
            end do
            if (.not. adv) exit scan_l;  i = i + 1
          end do scan_l
          scan_r: do
            adv = .false.
            do k = 1, icol
              a1 = a(col(k), j);  a2 = x(k)
              if (a1 > a2) then;  adv = .true.;  exit
              else if (a1 < a2) then;  exit;  end if
            end do
            if (.not. adv) exit scan_r;  j = j - 1
          end do scan_r
          if (i <= j) then
            do j2 = 1, iy
              x1 = a(j2,i);  a(j2,i) = a(j2,j);  a(j2,j) = x1
            end do
            i = i + 1;  j = j - 1
          end if
          if (i > j) exit
        end do
        if (j-lo < hi-i) then
          if (i<hi) then; sp=sp+1; stkL(sp)=i; stkR(sp)=hi; end if; hi=j
        else
          if (lo<j) then; sp=sp+1; stkL(sp)=lo; stkR(sp)=j; end if; lo=i
        end if
      end do
    end do
  end subroutine sortqr

  !==========================================================================
  !  hashvr — 원본 인터페이스 (dat, r, ind, m, n, mode, nr)
  !==========================================================================
  integer function hashvr(dat, r, ind, m, n, mode, nr)
    integer,  intent(in)    :: r, m, n, mode
    integer,  intent(inout) :: nr
    real(r8), intent(inout) :: ind(:,:)
    real(r8), intent(in)    :: dat(r)
    integer(ki8) :: iaddr
    integer :: iaddress, k, i, i1
    logical :: is_empty, is_match
    iaddr = 0_ki8
    do i = 1, r
      iaddr = iaddr + int(prime(mod(i-1,NP)+1), ki8) * int(dat(i), ki8)
    end do
    iaddress = int(abs(mod(iaddr, int(m, ki8)))) + 1
    do k = 1, MAX_PROBE
      is_empty = .true.;  is_match = .true.
      do i = 1, r
        i1 = ind(i, iaddress)
        if (i1 /= dat(i))          is_match = .false.
        if (ind(i, iaddress) /= 0) is_empty = .false.
        if (.not. is_empty .and. .not. is_match) exit
      end do
      if (is_empty .or. is_match) then
        if (is_empty .and. mode == 1) then
          ind(1:r, iaddress) = dat(1:r);  nr = nr + 1
        end if
        if (mode == 0 .and. is_empty) then; hashvr = 0
        else; hashvr = iaddress; end if
        return
      end if
      iaddress = mod(iaddress + PROBE_STEP - 1, m) + 1
    end do
    hashvr = -1
  end function hashvr

  !==========================================================================
  !  clusterr — 원본 인터페이스 (x, m1, m2, ix)
  !==========================================================================
  subroutine clusterr(x, m1, m2, ix)
    integer,  intent(in)    :: m1, m2, ix
    real(r8), intent(inout) :: x(:,:)
    integer :: i, j1, k1, k2
    j1 = 1;  k1 = 1
    do i = 1, ix
      do while (j1 <= m1 .and. x(1,j1) /= 0);  j1 = j1 + 1;  end do
      do while (k1 <= m1 .and. x(1,k1) == 0);  k1 = k1 + 1;  end do
      if (k1 > j1) then
        do k2 = 1, m2;  x(k2,j1) = x(k2,k1);  x(k2,k1) = 0;  end do
      else;  k1 = k1 + 1;  end if
    end do
  end subroutine clusterr

  !==========================================================================
  !  hashia1 — 원본 인터페이스 (x, nhash, n, ia, ja, a, m)
  !==========================================================================
  subroutine hashia1(x, nhash, n, ia, ja, a, m)
    integer,  intent(in)    :: nhash, n, m
    real(r8), intent(inout) :: x(:,:)
    integer,  intent(out)   :: ia(:), ja(:)
    real(r8), intent(out)   :: a(:)
    integer :: tmp(n), i, j, k, cut
    cut = 0;  tmp = 0
    do i = 1, nhash
      j = x(1,i)
      if (j /= 0) then
        if (j > n .or. x(2,i) > n) then;  x(1,i) = 0;  cut = cut + 1
        else;  tmp(j) = tmp(j) + 1;  end if
      end if
    end do
    ia(1) = 1
    do i = 1, n;  ia(i+1) = ia(i) + tmp(i);  end do
    if (ia(n+1)-1 > m) then
      print *, 'Too small parameter m in hashia: should be > ', ia(n+1);  stop
    end if
    tmp = 0
    do i = 1, nhash
      j = x(1,i)
      if (j /= 0) then
        k = ia(j)+tmp(j);  ja(k) = x(2,i);  a(k) = x(3,i);  tmp(j) = tmp(j)+1
      end if
    end do
    if (cut > 0) print *, cut, ' elements in HASH->IJA conversion had index(es) over ', n
  end subroutine hashia1

end module M_HashSort
