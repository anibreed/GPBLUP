module M_IO_Map
  use M_Kinds
  use M_Variables
  use M_Hashtable
  implicit none

  type, public :: MAPHashNode
    character(len=LEN_STR) :: snp_id_key
    integer :: map_index
    type(MAPHashNode), pointer :: next => null()
  end type MAPHashNode

  type, public :: MAPHashTable
    type(MAPHashNode), pointer :: table(:) => null()
    integer :: size = 0
    integer :: count = 0
  end type MAPHashTable

  public :: mht_create, mht_insert, mht_search, mht_free
  public :: mht_load_map_file

contains

  subroutine mht_create(mht, table_size)
    type(MAPHashTable), intent(out) :: mht
    integer, intent(in) :: table_size
    integer :: i, actual_size

    actual_size = htab_get_prime_size(table_size)
    allocate(mht%table(0:actual_size-1))
    mht%size = actual_size
    mht%count = 0

    do i = 0, actual_size - 1
      nullify(mht%table(i)%next)
    end do
  end subroutine mht_create

  subroutine mht_insert(mht, snp_id, map_idx)
    type(MAPHashTable), intent(inout) :: mht
    character(len=*), intent(in) :: snp_id
    integer, intent(in) :: map_idx
    integer :: h_idx
    type(MAPHashNode), pointer :: node

    h_idx = htab_hash_string(trim(snp_id), mht%size)

    allocate(node)
    node%snp_id_key = trim(snp_id)
    node%map_index = map_idx
    node%next => mht%table(h_idx)%next
    mht%table(h_idx)%next => node
    mht%count = mht%count + 1
  end subroutine mht_insert

  function mht_search(mht, snp_id) result(map_idx)
    type(MAPHashTable), intent(in) :: mht
    character(len=*), intent(in) :: snp_id
    integer :: map_idx
    integer :: h_idx
    type(MAPHashNode), pointer :: current

    map_idx = 0
    if (.not. associated(mht%table)) return

    h_idx = htab_hash_string(trim(snp_id), mht%size)
    current => mht%table(h_idx)%next

    do while (associated(current))
      if (trim(current%snp_id_key) == trim(snp_id)) then
        map_idx = current%map_index
        return
      end if
      current => current%next
    end do
  end function mht_search

  subroutine mht_load_map_file(mht, map_array, n_snp, quiet)
    use M_ReadFile, only: fopen, readline, N_recf
    type(MAPHashTable), intent(out) :: mht
    type(SNPInfo), allocatable, intent(out) :: map_array(:)
    integer, intent(out) :: n_snp
    logical, intent(in), optional :: quiet

    integer :: i, n, u_fd, j
    integer(kind=ki4) :: XI(MAX_VAR)
    character(len=MAX_STR) :: XC(MAX_VAR)
    logical :: do_log

    do_log = .true.
    if (present(quiet)) do_log = .not. quiet

    n_snp = N_recf(trim(MAPFile%FileName)) - MAPFile%Header
    allocate(map_array(n_snp))

    call mht_create(mht, int(n_snp * 1.5))

    u_fd = fopen(trim(MAPFile%FileName))
    do i = 1, MAPFile%Header
      call readline(u_fd, n, XC, XI, dlm_str=trim(MAPFile%Delim_char))
    end do

    if (do_log) then
      print *, "  M_MAPHashTable: MAP loading (", trim(MAPFile%FileName), ")"
    end if

    do i = 1, n_snp
      call readline(u_fd, n, XC, XI, dlm_str=trim(MAPFile%Delim_char))
      if (n < 0) exit

      do j = 1, MAPFile%NVAR
         select case(trim(MAPFile%FieldName(j)))
         case("SNP_ID")    ; map_array(i)%SNP_ID    = trim(XC(MAPFile%FieldLoc(j)))
         case("CHR")       ; map_array(i)%Chr       = XI(MAPFile%FieldLoc(j))
         case("POS")       ; map_array(i)%Pos       = XI(MAPFile%FieldLoc(j))
         case("ARRAY_ALL") ; map_array(i)%Array_All = XI(MAPFile%FieldLoc(j))
         case("ARRAY_CHR") ; map_array(i)%Array_Chr = XI(MAPFile%FieldLoc(j))
         end select
      end do

      call mht_insert(mht, map_array(i)%SNP_ID, i)
    end do

    close(u_fd)
    if (do_log) then
      print *, "  MAP hash built: ", mht%count, " SNPs"
    end if
  end subroutine mht_load_map_file

  subroutine mht_free(mht)
    type(MAPHashTable), intent(inout) :: mht
    integer :: i
    type(MAPHashNode), pointer :: cur, nxt

    if (.not. associated(mht%table)) return
    do i = 0, mht%size - 1
      cur => mht%table(i)%next
      do while (associated(cur))
        nxt => cur%next
        deallocate(cur)
        cur => nxt
      end do
    end do
    deallocate(mht%table)
    mht%size = 0
    mht%count = 0
  end subroutine mht_free

end module M_IO_Map
