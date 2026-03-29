program diag_colleau
  use colleau_mod
  use M_param
  implicit none
  integer :: n, i
  integer, allocatable :: sire(:), dam(:)
  real(kind=r8), allocatable :: F(:)
  character(len=256) :: logfile
  integer :: unit

  ! We'll build small arrays by reading `Geno_animal.txt_ren` and mapping IDs
  call build_from_ref_file('Geno_animal.txt_ren', n, sire, dam)
  allocate(F(n))
  write(*,*) 'Sample mapping (first 20):'
  do i = 1, min(20,n)
    write(*,'(I6,2X,I6,2X,I6)') i, sire(i), dam(i)
  end do
  call set_colleau_verbose(.true.)
  logfile = 'diag_colleau.log'
  open(newunit=unit, file=logfile, status='replace', action='write')
  write(unit,*) 'Running compute_inbreeding diagnostic, n=', n

  call compute_inbreeding(n, sire, dam, F)

  write(unit,*) 'Finished compute_inbreeding'
  do i = 1, n
    write(unit,'(A,I6,2X,A,ES18.12)') 'F[', i, ']=', F(i)
  end do
  close(unit)

contains
  subroutine build_from_ref_file(fname, n_out, sire_out, dam_out)
    character(len=*), intent(in) :: fname
    integer, intent(out) :: n_out
    integer, allocatable, intent(out) :: sire_out(:), dam_out(:)
    character(len=200) :: line
    integer :: iostat, unit, i
    character(len=LEN_STR) :: id, s, d
    integer, allocatable :: ids_idx(:)
    integer :: count
    character(len=LEN_STR), allocatable :: id_list(:), sire_list(:), dam_list(:)
    integer :: j, j2

    open(newunit=unit, file=fname, status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
      write(*,*) 'Cannot open', trim(fname)
      stop
    end if
    count = 0
    do
      read(unit,'(A)', iostat=iostat) line
      if (iostat /= 0) exit
      count = count + 1
    end do
    rewind(unit)
    n_out = count
    allocate(id_list(n_out), sire_list(n_out), dam_list(n_out))
    do i = 1, n_out
      read(unit,'(A)', iostat=iostat) line
      if (iostat /= 0) then
        id_list(i) = '0'
        sire_list(i) = '0'
        dam_list(i) = '0'
      else
        call parse_line(line, id, s, d)
        id_list(i) = id
        sire_list(i) = s
        dam_list(i) = d
      end if
    end do
    ! debug: print first few parsed tokens
    write(*,*) 'Parsed tokens sample (first 10):'
    do i = 1, min(10,n_out)
      write(*,'(I6,2X,A,2X,A,2X,A)') i, trim(id_list(i)), trim(sire_list(i)), trim(dam_list(i))
    end do
    close(unit)

    ! build mapping from ID -> index
    allocate(ids_idx(n_out))
    ids_idx = 0
    ! naive map: for small n, use O(n^2) search; acceptable here
    allocate(sire_out(n_out))
    allocate(dam_out(n_out))
    do i = 1, n_out
      sire_out(i) = 0
      dam_out(i) = 0
    end do
    do i = 1, n_out
      if (len_trim(sire_list(i)) > 0 .and. sire_list(i) /= '0') then
        do j = 1, n_out
          if (trim(id_list(j)) == trim(sire_list(i))) then
            sire_out(i) = j
            exit
          end if
        end do
      end if
      if (len_trim(dam_list(i)) > 0 .and. dam_list(i) /= '0') then
        do j2 = 1, n_out
          if (trim(id_list(j2)) == trim(dam_list(i))) then
            dam_out(i) = j2
            exit
          end if
        end do
      end if
    end do
  end subroutine

  subroutine parse_line(line, id, s, d)
    character(len=*), intent(in) :: line
    character(len=LEN_STR), intent(out) :: id, s, d
    integer :: n, pos, start, tok
    character(len=1) :: ch
    n = len_trim(line)
    id = '0'; s = '0'; d = '0'
    pos = 1
    tok = 0
    do while (pos <= n)
      ! skip whitespace
      do while (pos <= n)
        ch = line(pos:pos)
        if (ch /= ' ' .and. ch /= char(9)) exit
        pos = pos + 1
      end do
      if (pos > n) exit
      start = pos
      do while (pos <= n)
        ch = line(pos:pos)
        if (ch == ' ' .or. ch == char(9)) exit
        pos = pos + 1
      end do
      tok = tok + 1
      select case (tok)
      case (9)
        id = adjustl(line(start:pos-1))
      case (10)
        s = adjustl(line(start:pos-1))
      case (11)
        d = adjustl(line(start:pos-1))
      end select
      if (tok >= 11) exit
    end do
  end subroutine


end program diag_colleau
