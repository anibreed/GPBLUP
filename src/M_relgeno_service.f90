module M_relgeno_service
  use M_Kinds, only: r4, r8
  use M_ReadFile, only: fopen, N_recf
  use M_readpar_relgeno, only: relgeno_params_t, FMT_BINARY, FMT_TEXT, FMT_BOTH
  use M_StrEdit, only: to_upper, split_string, get_field_pos, remove_duplicate_spaces
  use LAPACKinv, only: LAPACK_inv
  implicit none
  private

  integer, parameter :: LEN_STR = 20
  integer, parameter :: MAX_STR = 128
  integer, parameter :: MAX_COL = 256
  character(len=1), parameter :: missC = achar(32)
  real(kind=r8), parameter :: missR = 0.0_r8

  public :: run_relgeno
  public :: build_g_matrices
  public :: write_relgeno_outputs

contains

  subroutine run_relgeno(params)
    type(relgeno_params_t), intent(in) :: params
    character(len=LEN_STR), allocatable :: id(:)
    real(kind=r8), allocatable :: g(:,:), ginv(:,:)
    logical :: has_ginv

    call validate_params(params)
    call print_param_summary(params)
    call build_g_matrices(params, id, g, ginv, has_ginv)
    if (has_ginv) then
      call write_relgeno_outputs(params, id, g, has_ginv, ginv)
    else
      call write_relgeno_outputs(params, id, g, has_ginv)
    end if
  end subroutine run_relgeno

  subroutine build_g_matrices(params, id, g, ginv, has_ginv)
    use, intrinsic :: iso_fortran_env, only: int64
    type(relgeno_params_t), intent(in) :: params
    character(len=LEN_STR), allocatable, intent(out) :: id(:)
    real(kind=r8), allocatable, intent(out) :: g(:,:), ginv(:,:)
    logical, intent(out) :: has_ginv

    integer, parameter :: dp = kind(1d0)
    integer, allocatable :: geno(:,:)
    character(len=MAX_STR) :: genofile, informat
    real(kind=dp), allocatable :: freq(:)
    real(kind=dp), allocatable :: gz(:,:), gw(:,:), gk(:)
    real(kind=r4), allocatable :: diag(:), offdiag(:)
    double precision :: t1, t2
    character(len=:), allocatable :: line
    character(len=:), allocatable :: line_ws
    integer :: recsize, unt, io
    integer :: i, j, k, l, nrec, nsnp, nv, ip_snp, cv
    integer :: missing_code
    real(kind=r8) :: weight, diag_weight
    integer(int64) :: missing_total
    integer(int64) :: total_calls
    character(len=LEN_STR) :: tokens(MAX_COL)
    character(len=1) :: dlm
    integer :: ntok, geno_col, id_col, data_nrec, geno_len, header_skip
    integer :: geno_s, geno_e, id_s, id_e, line_len
    logical :: has_header, is_popqc, found_field, first_line_is_data
    external dgemm

    genofile = trim(params%genofile)
    missing_code = params%missing_code
    weight = params%weight
    diag_weight = params%diag_weight

    recsize = chkRECL(genofile, params%header_lines)
    print *, 'records length=', recsize
    allocate(character(len=recsize+1) :: line)
    allocate(character(len=recsize+1) :: line_ws)

    nrec = N_recf(genofile, recsize)
    if (nrec <= 0) then
      write(*,*) 'No records found in GENOFILE: ', trim(genofile)
      stop
    end if

    unt = fopen(genofile, recsize)
    read(unt, '(a)', iostat=io) line
    if (io /= 0) then
      write(*,*) 'Reading error on Genotype file: ', trim(genofile)
      stop
    end if

    call detect_delim(line, dlm)
    if (dlm == ' ') then
      call remove_duplicate_spaces(line, line_ws)
    else
      line_ws = line
    end if
    call split_string(line_ws, tokens, ntok, dlm)
    geno_col = 0
    id_col = 0
    do i = 1, min(ntok, MAX_COL)
      if (trim(to_upper(tokens(i))) == 'GENO') geno_col = i
      if (trim(to_upper(tokens(i))) == 'ANIMAL_ID' .or. trim(to_upper(tokens(i))) == 'ID') id_col = i
    end do
    has_header = (geno_col > 0)
    is_popqc = has_header
    first_line_is_data = .false.

    if (.not. has_header) then
      if (ntok >= 2) then
        geno_col = 0
        id_col = 1
        line_len = len_trim(line)
        call get_last_field(line, line_len, geno_s, geno_e, found_field)
        if (.not. found_field .or. geno_e < geno_s) then
          write(*,*) 'Invalid GENO length in GENOFILE data: ', trim(genofile)
          stop
        end if
        nsnp = geno_e - geno_s + 1
        if (nsnp <= 0) then
          write(*,*) 'Invalid GENO length in GENOFILE data: ', trim(genofile)
          stop
        end if
        is_popqc = .true.
        first_line_is_data = .true.
        data_nrec = nrec
      else
        is_popqc = .false.
      end if
    end if

    if (is_popqc) then
      if (nrec <= 0) then
        write(*,*) 'No data records found in GENOFILE: ', trim(genofile)
        stop
      end if

      if (.not. first_line_is_data) then
        read(unt, '(a)', iostat=io) line
        if (io /= 0) then
          write(*,*) 'Reading error on Genotype file: ', trim(genofile)
          stop
        end if

        line_len = len_trim(line)
        call get_field_pos(line, line_len, dlm, geno_col, geno_s, geno_e, found_field)
        if (.not. found_field .or. geno_e < geno_s) then
          call get_last_field(line, line_len, geno_s, geno_e, found_field)
          if (.not. found_field .or. geno_e < geno_s) then
            write(*,*) 'GENO column not found in GENOFILE data: ', trim(genofile)
            stop
          end if
          nsnp = geno_e - geno_s + 1
          if (nsnp <= 0) then
            write(*,*) 'Invalid GENO length in GENOFILE data: ', trim(genofile)
            stop
          end if
          first_line_is_data = .true.
        else
          nsnp = geno_e - geno_s + 1
          if (nsnp <= 0) then
            write(*,*) 'Invalid GENO length in GENOFILE data: ', trim(genofile)
            stop
          end if
        end if

        header_skip = 0
        if (params%header_lines > 0) then
          header_skip = params%header_lines
        else
          header_skip = 1
        end if

        data_nrec = nrec - header_skip
        nrec = data_nrec

        print '(10x,a,2i7,2x,a)', 'No. SNP, Pos. of 1st Marker, FORMAT: ', nsnp, 1, 'DELIMITED'
        rewind(unt)
        if (header_skip > 0) then
          do j = 1, header_skip
            read(unt, '(a)', iostat=io) line
            if (io /= 0) exit
          end do
        end if
        print *, 'NREC= ', nrec, 'nSNP= ', nsnp
      else
        header_skip = 0
        if (params%header_lines > 0) header_skip = params%header_lines
        if (header_skip > 0) then
          data_nrec = nrec - header_skip
          nrec = data_nrec
          rewind(unt)
          do j = 1, header_skip
            read(unt, '(a)', iostat=io) line
            if (io /= 0) exit
          end do
          first_line_is_data = .false.
        end if
        print '(10x,a,2i7,2x,a)', 'No. SNP, Pos. of 1st Marker, FORMAT: ', nsnp, 1, 'DELIMITED'
        print *, 'NREC= ', nrec, 'nSNP= ', nsnp
      end if

      allocate(id(nrec), geno(nrec, nsnp), freq(nsnp))
    else
      ip_snp = index(line(1:30), ' ', back=.true.) + 1
      nsnp = index(line(ip_snp:), ' ') - 1
      if (nsnp <= 0) then
        write(*,*) 'Invalid SNP count in GENOFILE header: ', trim(genofile)
        stop
      end if

      write(informat, '(a,i0,a,i0,a)') '(a', ip_snp-2, ',1x,', nsnp, 'i1)'
      print '(10x,a,2i7,2x,a)', 'No. SNP, Pos. of 1st Marker, FORMAT: ', nsnp, ip_snp, informat
      rewind(unt)
      print *, 'NREC= ', nrec, 'nSNP= ', nsnp

      allocate(id(nrec), geno(nrec, nsnp), freq(nsnp))
    end if
    allocate(gw(3, nsnp), gk(nrec), gz(nrec, nsnp))
    allocate(g(nrec, nrec))

    freq = 0d0
    gk = 0d0
    gz = 0d0
    g = 0d0
    gw(1,:) = 0d0
    gw(2,:) = 1d0
    gw(3,:) = 2d0

    missing_total = 0_int64
    if (is_popqc) then
      do i = 1, nrec
        if (first_line_is_data .and. i == 1) then
          ! use the already-read first data line
        else
          read(unt, '(a)', iostat=io) line
          if (io /= 0) exit
        end if
        line_len = len_trim(line)
        call get_field_pos(line, line_len, dlm, geno_col, geno_s, geno_e, found_field)
        if (.not. found_field .or. geno_e < geno_s) then
          call get_last_field(line, line_len, geno_s, geno_e, found_field)
          if (.not. found_field .or. geno_e < geno_s) then
            write(*,*) 'GENO column missing at record ', i
            exit
          end if
        end if

        if (id_col > 0) then
          call get_field_pos(line, line_len, dlm, id_col, id_s, id_e, found_field)
          if (found_field .and. id_e >= id_s) then
            id(i) = trim(adjustl(line(id_s:id_e)))
          else
            id(i) = ''
          end if
        else
          call get_field_pos(line, line_len, dlm, 1, id_s, id_e, found_field)
          if (found_field .and. id_e >= id_s) then
            id(i) = trim(adjustl(line(id_s:id_e)))
          else
            id(i) = ''
          end if
        end if

        geno_len = geno_e - geno_s + 1
        if (geno_len <= 0) then
          geno(i,:) = missing_code
        else
          do k = 1, min(geno_len, nsnp)
            cv = ichar(line(geno_s + k - 1:geno_s + k - 1)) - 48
            if (cv >= 0 .and. cv <= 9) then
              geno(i,k) = cv
            else
              geno(i,k) = missing_code
            end if
          end do
          if (geno_len < nsnp) geno(i, geno_len+1:nsnp) = missing_code
        end if

        missing_total = missing_total + int(count(geno(i,:) == missing_code), int64)
        if (mod(i, 100) == 0) print *, i, id(i)
      end do
    else
      do i = 1, nrec
        read(unt, informat, iostat=io) id(i), geno(i,:)
        if (io /= 0) exit
        id(i) = trim(adjustl(id(i)))
        missing_total = missing_total + int(count(geno(i,:) == missing_code), int64)
        if (mod(i, 100) == 0) print *, i, id(i)
      end do
    end if

    call cpu_time(t1)

    do j = 1, nsnp
      nv = 0
      do i = 1, nrec
        if (geno(i,j) /= missing_code) then
          freq(j) = freq(j) + geno(i,j)
          nv = nv + 1
        end if
      end do

      if (nv <= 0) then
        write(*,*) 'All missing for SNP index ', j
        stop
      end if
      freq(j) = freq(j) / (2d0 * nv)

      do i = 1, 3
        if (geno(i,j) /= missing_code) then
          gw(i,j) = gw(i,j) - 2d0 * freq(j)
        else
          gw(i,j) = 0d0
        end if
      end do
      if (mod(j, 10000) == 0) print *, j, 'ith SNP'
    end do

    do i = 1, nrec
      do j = 1, nsnp
        if (geno(i,j) /= missing_code) gz(i,j) = gw(geno(i,j) + 1, j)
      end do
    end do
    print *, 'Finished calculating Z matrix (M-2pj)'

    call dgemm('N', 'T', nrec, nrec, nsnp, 1.0d0, gz, nrec, gz, nrec, 0.0d0, g, nrec)
    print *, 'Finished calculating ZZ'''

    do i = 1, nrec
      do j = 1, nsnp
        if (geno(i,j) /= missing_code) gk(i) = gk(i) + freq(j) * (1 - freq(j))
      end do
      gk(i) = gk(i) * 2d0
      if (gk(i) <= 0d0) then
        write(*,*) 'Non-positive sigma(2pq) for individual ', i
        stop
      end if
    end do
    print *, 'Finished calculating sigma(2pq)'

    do i = 1, nrec
      do j = i, nrec
        g(i,j) = g(i,j) / sqrt(gk(i) * gk(j))
        g(j,i) = g(i,j)
      end do
    end do
    print *, 'Finished calculating G-matrix divided by sigma(2pq)'

    call cpu_time(t2)
    print *, 'Time G-Matrix ', t2 - t1
    call cpu_time(t1)

    g = g * weight
    do i = 1, nrec
      g(i,i) = g(i,i) + diag_weight
    end do

    total_calls = int(nrec, int64) * int(nsnp, int64)
    if (total_calls > 0_int64) then
      print '(a,f8.4)', 'Missing rate = ', real(missing_total, r8) / real(total_calls, r8)
    end if

    if (params%compute_ginv) then
      ginv = LAPACK_inv(g)
      print *, 'Finished inverse G'
      has_ginv = .true.
    else
      has_ginv = .false.
    end if

    call cpu_time(t2)
    if (has_ginv) print *, 'Time Ginv-Matrix ', t2 - t1

    allocate(diag(nrec), offdiag(nrec * (nrec - 1) / 2))
    diag = missR
    offdiag = missR
    k = 0
    l = 0
    do i = 1, nrec
      do j = i, nrec
        l = l + 1
        if (j == i) then
          diag(i) = g(j,j)
        else
          k = k + 1
          offdiag(k) = g(i,j)
        end if
      end do
    end do

    print '(24x,a32)', 'Mean   SD     Min.   Max.'
    print '(a27,4f8.4)', 'Diagonal on G-matrix', mean(diag), sqrt(var(diag)), minval(diag), maxval(diag)
    print '(a27,4f8.4)', 'Offdiag  on G-matrix', mean(offdiag), sqrt(var(offdiag)), minval(offdiag), maxval(offdiag)

    if (has_ginv) then
      diag = missR
      offdiag = missR
      k = 0
      l = 0
      do i = 1, nrec
        do j = i, nrec
          l = l + 1
          if (j == i) then
            diag(i) = ginv(j,j)
          else
            k = k + 1
            offdiag(k) = ginv(i,j)
          end if
        end do
      end do

      print '(24x,a32)', 'Mean   SD     Min.   Max.'
      print '(a27,4f8.4)', 'Diagonal on Ginv-matrix', mean(diag), sqrt(var(diag)), minval(diag), maxval(diag)
      print '(a27,4f8.4)', 'Offdiag  on Ginv-matrix', mean(offdiag), sqrt(var(offdiag)), minval(offdiag), maxval(offdiag)
    end if

  end subroutine build_g_matrices

  subroutine write_relgeno_outputs(params, id, g, has_ginv, ginv)
    type(relgeno_params_t), intent(in) :: params
    character(len=LEN_STR), intent(in) :: id(:)
    real(kind=r8), intent(in) :: g(:,:)
    logical, intent(in) :: has_ginv
    real(kind=r8), intent(in), optional :: ginv(:,:)

    character(len=MAX_STR) :: outfile1, outfile2, outfile_gib
    integer :: i, j, l, unt_gdiag, unt_gib, unt_g_txt, unt_ginv_txt
    logical :: write_bin, write_txt
    real(kind=r8) :: fib

    write_bin = (params%output_format == FMT_BINARY .or. params%output_format == FMT_BOTH)
    write_txt = (params%output_format == FMT_TEXT .or. params%output_format == FMT_BOTH)

    outfile1 = trim(adjustl(params%outfile_g))
    if (write_bin) then
      open(2, file=outfile1, form='unformatted', access='direct', recl=2 * LEN_STR + 8, status='replace')
    end if
    if (write_txt) then
      open(newunit=unt_g_txt, file=trim(adjustl(params%outfile_g)) // '.txt', status='replace', action='write')
    end if

    unt_gdiag = fopen(trim(adjustl(params%outfile_gdiag)))
    outfile_gib = trim(adjustl(params%outfile_gib))
    open(newunit=unt_gib, file=outfile_gib, status='replace', action='write')
    l = 0
    do i = 1, size(id)
      do j = i, size(id)
        l = l + 1
        if (write_bin) write(2, rec=l) id(i), id(j), g(i,j)
        if (write_txt) write(unt_g_txt, '(a,1x,a,1x,es24.16)') trim(id(i)), trim(id(j)), g(i,j)
        if (j == i) then
          write(unt_gdiag, fmt='(i8,2x,a16,12f12.4)') j, id(j), g(j,j)
          fib = g(j,j) - 1.0_r8
          write(unt_gib, '(a,1x,es24.16)') trim(id(j)), fib
        end if
      end do
    end do
    if (write_bin) close(unit=2)
    if (write_txt) close(unit=unt_g_txt)
    close(unit=unt_gdiag)
    close(unit=unt_gib)

    if (.not. has_ginv) then
      print *, 'INFO: COMPUTE_GINV=FALSE. Skipping G inverse output.'
      return
    end if

    if (.not. present(ginv)) then
      write(*,*) 'GINV output requested but not provided'
      stop
    end if

    outfile2 = trim(adjustl(params%outfile_ginv))
    if (write_bin) then
      open(3, file=outfile2, form='unformatted', access='direct', recl=2 * LEN_STR + 8, status='replace')
    end if
    if (write_txt) then
      open(newunit=unt_ginv_txt, file=trim(adjustl(params%outfile_ginv)) // '.txt', status='replace', action='write')
    end if

    l = 0
    do i = 1, size(id)
      do j = i, size(id)
        l = l + 1
        if (write_bin) write(3, rec=l) id(i), id(j), ginv(i,j)
        if (write_txt) write(unt_ginv_txt, '(a,1x,a,1x,es24.16)') trim(id(i)), trim(id(j)), ginv(i,j)
      end do
    end do
    if (write_bin) close(unit=3)
    if (write_txt) close(unit=unt_ginv_txt)
  end subroutine write_relgeno_outputs

  subroutine validate_params(params)
    type(relgeno_params_t), intent(in) :: params
    logical :: exists

    inquire(file=trim(params%genofile), exist=exists)
    if (.not. exists) then
      write(*,*) 'GENOFILE does not exist: ', trim(params%genofile)
      stop
    end if

    if (params%compute_ginv .and. .not. params%compute_g) then
      write(*,*) 'COMPUTE_GINV requires COMPUTE_G=TRUE'
      stop
    end if

    if (.not. params%compute_g .and. len_trim(params%outfile_gib) > 0) then
      write(*,*) 'OUT_GIB requires COMPUTE_G=TRUE'
      stop
    end if

    if (params%weight <= 0.0_r8) then
      write(*,*) 'WEIGHT must be positive'
      stop
    end if

    if (params%diag_weight < 0.0_r8) then
      write(*,*) 'DIAG_WEIGHT must be non-negative'
      stop
    end if

    if (params%missing_code < 0) then
      write(*,*) 'MISSING code must be non-negative'
      stop
    end if

    if (len_trim(params%outfile_g) == 0) then
      write(*,*) 'OUT_G is empty'
      stop
    end if

    if (len_trim(params%outfile_gdiag) == 0) then
      write(*,*) 'OUT_GDIAG is empty'
      stop
    end if

    if (len_trim(params%outfile_gib) == 0) then
      write(*,*) 'OUT_GIB is empty'
      stop
    end if

    if (params%output_format < FMT_BINARY .or. params%output_format > FMT_BOTH) then
      write(*,*) 'OUTPUT FORMAT must be BINARY, TEXT, or BOTH'
      stop
    end if

    if (params%compute_ginv) then
      if (len_trim(params%outfile_ginv) == 0) then
        write(*,*) 'OUT_GINV is empty'
        stop
      end if
    end if
  end subroutine validate_params

  subroutine print_param_summary(params)
    type(relgeno_params_t), intent(in) :: params

    print *, 'RELGENO parameters:'
    print '(a,1x,a)', '  GENOFILE:', trim(params%genofile)
    print '(a,1x,i0)', '  MISSING:', params%missing_code
    print '(a,1x,f8.4)', '  WEIGHT:', params%weight
    print '(a,1x,f8.4)', '  DIAG_WEIGHT:', params%diag_weight
    print '(a,1x,l1)', '  COMPUTE_G:', params%compute_g
    print '(a,1x,l1)', '  COMPUTE_GINV:', params%compute_ginv
    print '(a,1x,a)', '  OUT_G:', trim(params%outfile_g)
    print '(a,1x,a)', '  OUT_GDIAG:', trim(params%outfile_gdiag)
    print '(a,1x,a)', '  OUT_GIB:', trim(params%outfile_gib)
    print '(a,1x,i0)', '  OUTPUT_FORMAT:', params%output_format
    print '(a,1x,a)', '  OUT_GINV:', trim(params%outfile_ginv)
  end subroutine print_param_summary

  subroutine detect_delim(str, dlm)
    character(len=*), intent(in) :: str
    character(len=1), intent(out) :: dlm

    if (index(str, achar(9)) > 0) then
      dlm = achar(9)
    else if (index(str, ',') > 0) then
      dlm = ','
    else if (index(str, '|') > 0) then
      dlm = '|'
    else if (index(str, ';') > 0) then
      dlm = ';'
    else
      dlm = ' '
    end if
  end subroutine detect_delim

  subroutine get_last_field(str, str_len, out_start, out_end, out_found)
    character(len=*), intent(in) :: str
    integer, intent(in) :: str_len
    integer, intent(out) :: out_start, out_end
    logical, intent(out) :: out_found

    integer :: i

    out_found = .false.
    out_start = 0
    out_end = -1
    if (str_len <= 0) return

    i = str_len
    do while (i >= 1 .and. (str(i:i) == ' ' .or. str(i:i) == achar(9)))
      i = i - 1
    end do
    if (i < 1) return
    out_end = i
    do while (i >= 1 .and. str(i:i) /= ' ' .and. str(i:i) /= achar(9))
      i = i - 1
    end do
    out_start = i + 1
    if (out_start <= out_end) out_found = .true.
  end subroutine get_last_field

  function var(x) result(v)
    implicit none
    real(kind=r4), intent(in) :: x(:)
    real(kind=r4) :: v
    v = (dot_product(x, x) - sum(x) * sum(x) / size(x)) / (size(x) - 1.d0)
  end function var

  function mean(x) result(m)
    implicit none
    real(kind=r4), intent(in) :: x(:)
    real(kind=r4) :: m
    m = sum(x) / size(x)
  end function mean

  integer function chkRECL(filename, header_lines) result(rec_len)
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: header_lines
    integer :: unitf, ios
    character(len=200000) :: linebuf
    integer, parameter :: max_scan = 1000
    integer, parameter :: target_len = 50000
    integer :: scan_count, line_len, skip_lines, i

    rec_len = 0
    unitf = fopen(filename, 200000)
    if (unitf < 0) return
    skip_lines = 0
    if (present(header_lines)) skip_lines = max(0, header_lines)
    do i = 1, skip_lines
      read(unitf, '(A)', iostat=ios) linebuf
      if (ios /= 0) exit
    end do

    scan_count = 0
    do
      read(unitf, '(A)', iostat=ios) linebuf
      if (ios /= 0) exit
      line_len = len_trim(linebuf)
      if (line_len > rec_len) rec_len = line_len
      scan_count = scan_count + 1
      if (rec_len >= target_len) exit
      if (scan_count >= max_scan) exit
    end do
    close(unit=unitf)
    if (rec_len <= 0) rec_len = 1024
  end function chkRECL

end module M_relgeno_service
