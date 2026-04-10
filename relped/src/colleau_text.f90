subroutine write_triplets_text_mod(n, sire, dam, F, orig_id, filename)
  use M_param
  use omp_lib
  use colleau_mod, only: colleau_apply
  implicit none
  real(kind=r8), parameter :: TEXT_TRIPLET_EPS = 1.0d-6
  integer, intent(in) :: n
  integer, intent(in) :: sire(n), dam(n)
  real(kind=r8), intent(in) :: F(n)
  character(len=*), intent(in) :: orig_id(n)
  character(len=*), intent(in) :: filename
  integer(kind=ki4) :: j, i, ios
  integer(kind=ki4) :: row_out, col_out
  integer :: unit, merge_u, nthreads, tid
  character(len=MAX_STR), allocatable :: part_files(:)
  character(len=MAX_STR) :: part_name, line_buf
  real(kind=r8), allocatable :: x(:), y(:)

  nthreads = max(1, omp_get_max_threads())
  allocate(part_files(0:nthreads-1))

!$omp parallel default(shared) private(tid,unit,part_name,j,i,x,y,row_out,col_out,ios,line_buf)
  allocate(x(n), y(n))
  tid = omp_get_thread_num()
  write(part_name,'(A,".part.",I0)') trim(filename), tid
  part_files(tid) = trim(part_name)

  open(newunit=unit, file=trim(part_name), status='replace', action='write', &
    form='formatted', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open temporary triplet file:', trim(part_name)
    stop 7
  end if

!$omp do schedule(static)
  do j = 1, n
    x = 0.0_r8
    x(j) = 1.0_r8
    call colleau_apply(n, sire, dam, F, x, y)
    do i = j, n
      if (abs(y(i)) > TEXT_TRIPLET_EPS) then
        row_out = int(i, kind=ki4)
        col_out = int(j, kind=ki4)
        write(unit, '(I0,1X,I0,1X,F12.6,1X,A,1X,A)') row_out, col_out, y(i), &
          trim(orig_id(row_out)), trim(orig_id(col_out))
      end if
    end do
  end do
!$omp end do

  close(unit)
  deallocate(x, y)
!$omp end parallel

  open(newunit=merge_u, file=filename, status='replace', action='write', &
    form='formatted', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open triplet text file:', trim(filename)
    stop 7
  end if

  do tid = 0, nthreads - 1
    open(newunit=unit, file=trim(part_files(tid)), status='old', action='read', &
      form='formatted', iostat=ios)
    if (ios /= 0) cycle
    do
      read(unit, '(A)', iostat=ios) line_buf
      if (ios /= 0) exit
      write(merge_u, '(A)') trim(line_buf)
    end do
    close(unit, status='delete')
  end do

  close(merge_u)
  deallocate(part_files)
end subroutine write_triplets_text_mod

subroutine write_ainv_triplets_text_mod(n, sire, dam, F, orig_id, filename)
  use M_param
  use colleau_mod, only: add_lower_term, sort_key_values
  implicit none
  real(kind=r8), parameter :: TEXT_TRIPLET_EPS = 1.0d-6
  integer(kind=ki4), intent(in) :: n
  integer(kind=ki4), intent(in) :: sire(n), dam(n)
  real(kind=r8), intent(in) :: F(n)
  character(len=*), intent(in) :: orig_id(n)
  character(len=*), intent(in) :: filename

  integer(kind=ki4) :: i, s, d, unit, ios
  integer(kind=ki4) :: nt, max_terms
  integer(kind=ki4), allocatable :: rows(:), cols(:)
  integer(kind=ki8), allocatable :: keys(:)
  real(kind=r8), allocatable :: vals(:)
  real(kind=r8) :: di, bi
  real(kind=r8) :: v_acc
  integer(kind=ki8) :: key_acc
  integer(kind=ki4) :: row_acc, col_acc

  max_terms = 8 * n + 16
  allocate(rows(max_terms), cols(max_terms), vals(max_terms))
  nt = 0

  do i = 1, n
    s = sire(i)
    d = dam(i)

    if (s > 0 .and. d > 0) then
      di = 0.5_r8 - 0.25_r8 * (F(s) + F(d))
    else if (s > 0 .or. d > 0) then
      if (s > 0) then
        di = 0.75_r8 - 0.25_r8 * F(s)
      else
        di = 0.75_r8 - 0.25_r8 * F(d)
      end if
    else
      di = 1.0_r8
    end if

    if (di <= 1.0d-14) cycle
    bi = 1.0_r8 / di

    call add_lower_term(i, i, 1.0_r8 * bi, rows, cols, vals, nt, max_terms)

    if (s > 0) then
      call add_lower_term(i, s, -0.5_r8 * bi, rows, cols, vals, nt, max_terms)
      call add_lower_term(s, s,  0.25_r8 * bi, rows, cols, vals, nt, max_terms)
    end if
    if (d > 0) then
      call add_lower_term(i, d, -0.5_r8 * bi, rows, cols, vals, nt, max_terms)
      call add_lower_term(d, d,  0.25_r8 * bi, rows, cols, vals, nt, max_terms)
    end if
    if (s > 0 .and. d > 0) then
      call add_lower_term(s, d, 0.25_r8 * bi, rows, cols, vals, nt, max_terms)
    end if
  end do

  if (nt <= 0) then
    open(newunit=unit, file=filename, status='replace', action='write', &
      form='formatted', iostat=ios)
    if (ios == 0) close(unit)
    deallocate(rows, cols, vals)
    return
  end if

  allocate(keys(nt))
  do i = 1, nt
    keys(i) = int(rows(i), kind=ki8) * int(n + 1, kind=ki8) + int(cols(i), kind=ki8)
  end do
  call sort_key_values(keys, vals, nt)

  open(newunit=unit, file=filename, status='replace', action='write', &
    form='formatted', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open inverse triplet text file:', trim(filename)
    deallocate(rows, cols, vals, keys)
    return
  end if

  key_acc = keys(1)
  v_acc = vals(1)
  do i = 2, nt
    if (keys(i) == key_acc) then
      v_acc = v_acc + vals(i)
    else
      if (abs(v_acc) > TEXT_TRIPLET_EPS) then
        row_acc = int(key_acc / int(n + 1, kind=ki8), kind=ki4)
        col_acc = int(mod(key_acc, int(n + 1, kind=ki8)), kind=ki4)
        if (row_acc >= 1 .and. col_acc >= 1) then
          write(unit, '(I0,1X,I0,1X,F12.6,1X,A,1X,A)') row_acc, col_acc, v_acc, &
            trim(orig_id(row_acc)), trim(orig_id(col_acc))
        end if
      end if
      key_acc = keys(i)
      v_acc = vals(i)
    end if
  end do

  if (abs(v_acc) > TEXT_TRIPLET_EPS) then
    row_acc = int(key_acc / int(n + 1, kind=ki8), kind=ki4)
    col_acc = int(mod(key_acc, int(n + 1, kind=ki8)), kind=ki4)
    if (row_acc >= 1 .and. col_acc >= 1) then
      write(unit, '(I0,1X,I0,1X,F12.6,1X,A,1X,A)') row_acc, col_acc, v_acc, &
        trim(orig_id(row_acc)), trim(orig_id(col_acc))
    end if
  end if

  close(unit)
  deallocate(rows, cols, vals, keys)
end subroutine write_ainv_triplets_text_mod

subroutine triplets_bin_to_text(bin_file, txt_file)
  use M_param
  implicit none
  real(kind=r8), parameter :: TEXT_TRIPLET_EPS = 1.0d-6
  character(len=*), intent(in) :: bin_file, txt_file
  integer :: u_bin, u_txt, ios
  integer(kind=ki4) :: ver, nfile, row, col
  integer(kind=ki8) :: nnz, k
  real(kind=r8) :: val
  character(len=4) :: magic

  open(newunit=u_bin, file=trim(bin_file), status='old', action='read', &
    form='unformatted', access='stream', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open triplet binary file:', trim(bin_file)
    return
  end if

  read(u_bin, iostat=ios) magic, ver, nfile, nnz
  if (ios /= 0) then
    write(*,*) 'Invalid triplet binary header:', trim(bin_file)
    close(u_bin)
    return
  end if

  open(newunit=u_txt, file=trim(txt_file), status='replace', action='write', &
    form='formatted', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open triplet text file:', trim(txt_file)
    close(u_bin)
    return
  end if

  do k = 1, nnz
    read(u_bin, iostat=ios) row, col, val
    if (ios /= 0) exit
    if (abs(val) > TEXT_TRIPLET_EPS) then
      write(u_txt, '(I0,1X,I0,1X,F12.4)') row, col, val
    end if
  end do

  close(u_bin)
  close(u_txt)
end subroutine triplets_bin_to_text
