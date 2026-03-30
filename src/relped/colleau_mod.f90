module colleau_mod
  use M_param
  use omp_lib
  implicit none
  integer(kind=ki4), save :: call_count = 0
  integer(kind=ki4), parameter :: PROG_STEP = 100
  integer(kind=ki4), save :: max_ns = 0
  integer(kind=ki4), save :: max_nd = 0
  logical, save :: verbose = .false.
  
  
  ! Colleau 관련 서브루틴 모듈
  ! use real(kind=r8) for double precision
contains

  subroutine set_colleau_verbose(flag)
    logical, intent(in) :: flag
    verbose = flag
  end subroutine set_colleau_verbose



  subroutine build_contrib_sparse_local(start, sire, dam, anc, val, nanc, maxanc)
    implicit none
    integer(kind=ki4), intent(in) :: start
    integer(kind=ki4), intent(in) :: sire(:), dam(:)
    integer(kind=ki4), intent(out) :: anc(:), nanc
    real(kind=r8), intent(out) :: val(:)
    integer(kind=ki4), intent(in) :: maxanc

    integer(kind=ki4), allocatable :: stack(:)
    real(kind=r8), allocatable :: wstack(:)
    integer(kind=ki4) :: top, k, s, d, i
    real(kind=r8) :: curr_w

    allocate(stack(maxanc))
    allocate(wstack(maxanc))

    ! initialize output buffers
    anc = 0
    val = 0.0_r8

    ! We'll traverse ancestors starting at `start`, but only record
    ! contributions for founders (animals with no parents). This builds
    ! k-vectors over founders which matches the tabular method.
    nanc = 0

    top = 1
    stack(1) = start
    wstack(1) = 1.0_r8

    do while (top > 0)
      k = stack(top)
      ! current weight at top
      curr_w = wstack(top)
      ! pop
      top = top - 1
      s = 0; d = 0
      if (k >= 1 .and. k <= size(sire)) then
        s = sire(k)
        d = dam(k)
      end if

      ! Accumulate contribution for founders always; for non-founders,
      ! record contributions for all ancestors except the starting node.
      ! This ensures founders (including when start is a founder) are
      ! captured while avoiding duplicating the start's self-contribution
      ! when not appropriate.
      ! Record this node's contribution (including the start node).
      call accum_local(k, curr_w, anc, val, nanc)
      ! If node has parents, push them for further expansion with half weight
      if (s > 0) then
        top = top + 1
        if (top > maxanc) stop 'stack overflow in build_contrib_sparse_local'
        stack(top) = s
        wstack(top) = 0.5_r8*curr_w
      end if
      if (d > 0) then
        top = top + 1
        if (top > maxanc) stop 'stack overflow in build_contrib_sparse_local'
        stack(top) = d
        wstack(top) = 0.5_r8*curr_w
      end if

    end do

  end subroutine build_contrib_sparse_local

  subroutine accum_local(id, w, anc, val, nanc)
    implicit none
    integer(kind=ki4), intent(in) :: id
    real(kind=r8), intent(in) :: w
    integer(kind=ki4), intent(inout) :: anc(:), nanc
    real(kind=r8), intent(inout) :: val(:)
    integer(kind=ki4) :: ii
    do ii = 1, nanc
      if (anc(ii) == id) then
        val(ii) = val(ii) + w
        return
      end if
    end do

    nanc = nanc + 1
    if (nanc > size(anc)) stop 'anc buffer overflow in accum_local'
    anc(nanc) = id
    val(nanc) = w

  end subroutine accum_local

  real(kind=r8) function dot_sparse_local(a1, v1, n1, a2, v2, n2)
    implicit none
    integer(kind=ki4), intent(in) :: a1(:), a2(:)
    real(kind=r8), intent(in) :: v1(:), v2(:)
    integer(kind=ki4), intent(in) :: n1, n2
    integer(kind=ki4) :: i, j

    dot_sparse_local = 0.0_r8
    do i = 1, n1
      do j = 1, n2
        if (a1(i) == a2(j)) then
          dot_sparse_local = dot_sparse_local + v1(i)*v2(j)
          exit
        end if
      end do
    end do

  end function dot_sparse_local

  subroutine compute_inbreeding(n, sire, dam, F)
    ! simple inbreeding calculation using ancestor sparse expansion
    implicit none
    integer(kind=ki4), intent(in) :: n
    integer(kind=ki4), intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(out) :: F(n)

    integer(kind=ki4) :: i, n_s, n_d, maxanc
    integer(kind=ki4) :: k
    integer(kind=ki4) :: ttt, jj
    integer(kind=ki4) :: anc_id
    real(kind=r8) :: dval
    integer :: outunit, ios
    character(len=128) :: fname
    integer(kind=ki4), allocatable :: anc_s(:), anc_d(:)
    real(kind=r8), allocatable :: val_s(:), val_d(:)
    real(kind=r8) :: a_sd
    real(kind=r8) :: t_start, t_end, t_elapsed, t_compute_total
    t_compute_total = 0.0_r8

    F = 0.0_r8

    maxanc = min(n, MAXANC_LIMIT)

!$omp parallel do default(shared) private(i,n_s,n_d,anc_s,anc_d,val_s,val_d,a_sd)
    do i = 1, n
      if (sire(i) == 0 .or. dam(i) == 0) then
        F(i) = 0.0_r8
      else
        ! time the per-target compute (high-resolution)
        t_start = omp_get_wtime()
        n_s = 0; n_d = 0
        allocate(anc_s(maxanc)); allocate(anc_d(maxanc))
        allocate(val_s(maxanc)); allocate(val_d(maxanc))
        anc_s = 0; anc_d = 0
        val_s = 0.0_r8; val_d = 0.0_r8
        ! normal execution (no per-target debug)
        call build_contrib_sparse_local(sire(i), sire, dam, anc_s, val_s, n_s, maxanc)
        call build_contrib_sparse_local(dam(i),  sire, dam, anc_d, val_d, n_d, maxanc)
        ! (no temporary debug prints here)
        ! (diagnostic sums removed)
        ! Compute a_sd = sum_{anc} val_s(anc)*val_d(anc)*D(anc)
        ! where D(anc) depends on the inbreeding of anc's parents
        a_sd = 0.0_r8
        if (n_s > 0 .and. n_d > 0) then
          do jj = 1, n_s
            do k = 1, n_d
              if (anc_s(jj) == anc_d(k)) then
                anc_id = anc_s(jj)
                if (anc_id >= 1 .and. anc_id <= size(sire)) then
                  if (sire(anc_id) /= 0 .and. dam(anc_id) /= 0) then
                    dval = 0.5_r8 - 0.25_r8*(F(sire(anc_id)) + F(dam(anc_id)))
                  else if (sire(anc_id) /= 0 .or. dam(anc_id) /= 0) then
                    if (sire(anc_id) /= 0) then
                      dval = 0.75_r8 - 0.25_r8*F(sire(anc_id))
                    else
                      dval = 0.75_r8 - 0.25_r8*F(dam(anc_id))
                    end if
                  else
                    dval = 1.0_r8
                  end if
                else
                  dval = 1.0_r8
                end if
                a_sd = a_sd + val_s(jj) * val_d(k) * dval
                exit
              end if
            end do
          end do
        end if
        ! Now a_sd represents A(s,d); inbreeding F = 0.5 * A(s,d)
        F(i) = 0.5_r8 * a_sd
        t_end = omp_get_wtime()
        t_elapsed = t_end - t_start
      !$omp atomic
        t_compute_total = t_compute_total + t_elapsed
      ! (DBG_ASSIGN removed)
        ! (per-target debug logic removed)
        ! update global maxima (thread-safe)
      !$omp critical
        if (n_s > max_ns) max_ns = n_s
        if (n_d > max_nd) max_nd = n_d
      !$omp end critical

        ! periodic progress and basic error checks (synchronized)
        if (mod(i, PROG_STEP) == 0) then
      !$omp critical
          if (verbose) then
            write(*,'(A,I8,2X,A,I8,2X,A,I8,2X,A,I8,A,F10.4)') 'compute_inbreeding: processed', i, 'of', n, 'n_s=', n_s, 'n_d=', n_d,' F=', F(i)
            if (n_s >= maxanc .or. n_d >= maxanc) then
              write(*,'(A)') 'WARNING: ancestor expansion reached maxanc limit; results may be truncated'
            end if
          end if
      !$omp end critical
        end if
        ! (debug prints removed)
        deallocate(anc_s, anc_d, val_s, val_d)
      end if
    end do
!$omp end parallel do
  ! After parallel region, print maximum ancestor counts observed and timing
  if (verbose) then
    write(*,'(A,F8.4)') 'compute_inbreeding: total compute time (sum threads) =', t_compute_total
  end if
  ! After parallel region, print maximum ancestor counts observed
  if (verbose) then
    call timestamp()
    write(*,'(A,I8,2X,A,I8)') 'compute_inbreeding: maximum n_s observed =', max_ns, 'maximum n_d observed =', max_nd
  end if

  end subroutine compute_inbreeding

  subroutine colleau_apply(n, sire, dam, F, x, y)
    ! Compute y = A * x using Colleau indirect method (T' D T) with precomputed F
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(in) :: F(n)
    real(kind=r8), intent(in) :: x(n)
    real(kind=r8), intent(out) :: y(n)

    real(kind=r8), allocatable :: w(:)
    real(kind=r8), allocatable :: d_vec(:)
    integer :: i, s, d, ii

    allocate(w(n), d_vec(n))

    ! 1) w = T' * x  (reverse order)
    w = x
    do i = n, 1, -1
      s = sire(i)
      d = dam(i)
      if (s /= 0) w(s) = w(s) + 0.5_r8 * w(i)
      if (d /= 0) w(d) = w(d) + 0.5_r8 * w(i)
    end do

    ! 2) build D from F: diagonal elements d_i
    do i = 1, n
      s = sire(i)
      d = dam(i)
      if (s /= 0 .and. d /= 0) then
        d_vec(i) = 0.5_r8 - 0.25_r8*(F(s) + F(d))
      else if (s /= 0 .or. d /= 0) then
        if (s /= 0) then
          d_vec(i) = 0.75_r8 - 0.25_r8*F(s)
        else
          d_vec(i) = 0.75_r8 - 0.25_r8*F(d)
        end if
      else
        d_vec(i) = 1.0_r8
      end if
    end do

      ! progress/counting for calls (use saved counter)
      call_count = call_count + 1
      if (mod(call_count, PROG_STEP) == 0 .and. verbose) then
        call timestamp()
        write(*,'(A,I8,2X,A,I8)') 'colleau_apply: call', call_count, ' n=', n
        write(*,'(A)') 'sample d_vec[1:min(10,n)]:'
        do ii = 1, min(10,n)
          write(*,'(F10.6,1x)',advance='no') d_vec(ii)
        end do
        write(*,*)
        ! basic check: any non-positive diagonal entry
        do ii = 1, n
          if (d_vec(ii) <= 0.0_r8) then
            write(*,'(A,I8)') 'WARNING: non-positive D element at index', ii
            exit
          end if
        end do
      end if

    ! 3) z = D * w (store in w)
    do i = 1, n
      w(i) = d_vec(i) * w(i)
    end do

    ! 4) y = T * z (forward order)
    y = w
    do i = 1, n
      s = sire(i)
      d = dam(i)
      if (s /= 0) y(i) = y(i) + 0.5_r8 * y(s)
      if (d /= 0) y(i) = y(i) + 0.5_r8 * y(d)
    end do

    deallocate(w, d_vec)
  end subroutine colleau_apply

  subroutine write_triplets(n, sire, dam, F, filename)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(in) :: F(n)
    character(len=*), intent(in) :: filename
    integer(kind=ki4) :: j, i
    real(kind=r8), allocatable :: x(:), y(:)
    integer :: unit

    allocate(x(n), y(n))
    open(newunit=unit, file=filename, status='replace', action='write', iostat=i)
    if (i /= 0) then
      write(*,*) 'Cannot open triplet file:', trim(filename)
      return
    end if

      do j = 1, n
      x = 0.0_r8
      x(j) = 1.0_r8
      call colleau_apply(n, sire, dam, F, x, y)
      do i = j, n
        if (abs(y(i)) > 0.0_r8) write(unit,'(I8,1x,I8,1x,ES16.8)') i, j, y(i)
      end do
      if (verbose) then
        write(*,*) 'Write relationships: Processed column', j, 'of', n
      end if
    end do

    close(unit)
    deallocate(x, y)
  end subroutine write_triplets

  subroutine write_ainv_triplets(n, sire, dam, F, filename)
    implicit none
    integer(kind=ki4), intent(in) :: n
    integer(kind=ki4), intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(in) :: F(n)
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
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios == 0) close(unit)
      deallocate(rows, cols, vals)
      return
    end if

    allocate(keys(nt))
    do i = 1, nt
      keys(i) = int(rows(i), kind=ki8) * int(n + 1, kind=ki8) + int(cols(i), kind=ki8)
    end do
    call sort_key_values(keys, vals, nt)

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open inverse triplet file:', trim(filename)
      deallocate(rows, cols, vals, keys)
      return
    end if

    key_acc = keys(1)
    v_acc = vals(1)
    do i = 2, nt
      if (keys(i) == key_acc) then
        v_acc = v_acc + vals(i)
      else
        if (abs(v_acc) > 1.0d-16) then
          row_acc = int(key_acc / int(n + 1, kind=ki8), kind=ki4)
          col_acc = int(mod(key_acc, int(n + 1, kind=ki8)), kind=ki4)
          if (row_acc >= 1 .and. col_acc >= 1) then
            write(unit,'(I8,1x,I8,1x,ES16.8)') row_acc, col_acc, v_acc
          end if
        end if
        key_acc = keys(i)
        v_acc = vals(i)
      end if
    end do

    if (abs(v_acc) > 1.0d-16) then
      row_acc = int(key_acc / int(n + 1, kind=ki8), kind=ki4)
      col_acc = int(mod(key_acc, int(n + 1, kind=ki8)), kind=ki4)
      if (row_acc >= 1 .and. col_acc >= 1) then
        write(unit,'(I8,1x,I8,1x,ES16.8)') row_acc, col_acc, v_acc
      end if
    end if

    close(unit)
    deallocate(rows, cols, vals, keys)
  end subroutine write_ainv_triplets

  subroutine add_lower_term(i, j, v, rows, cols, vals, nt, max_terms)
    implicit none
    integer(kind=ki4), intent(in) :: i, j, max_terms
    real(kind=r8), intent(in) :: v
    integer(kind=ki4), intent(inout) :: nt
    integer(kind=ki4), intent(inout) :: rows(:), cols(:)
    real(kind=r8), intent(inout) :: vals(:)
    integer(kind=ki4) :: r, c

    if (i <= 0 .or. j <= 0) return
    if (abs(v) <= 1.0d-20) return

    r = i
    c = j
    if (r < c) then
      r = j
      c = i
    end if

    nt = nt + 1
    if (nt > max_terms) stop 'write_ainv_triplets: term buffer overflow'
    rows(nt) = r
    cols(nt) = c
    vals(nt) = v
  end subroutine add_lower_term

  subroutine sort_key_values(keys, vals, n)
    implicit none
    integer(kind=ki8), intent(inout) :: keys(:)
    real(kind=r8), intent(inout) :: vals(:)
    integer(kind=ki4), intent(in) :: n
    if (n <= 1) return
    call quicksort_key_values(keys, vals, 1_ki4, n)
  end subroutine sort_key_values

  recursive subroutine quicksort_key_values(keys, vals, lo, hi)
    implicit none
    integer(kind=ki8), intent(inout) :: keys(:)
    real(kind=r8), intent(inout) :: vals(:)
    integer(kind=ki4), intent(in) :: lo, hi
    integer(kind=ki4) :: i, j
    integer(kind=ki8) :: pivot, tkey
    real(kind=r8) :: tval

    i = lo
    j = hi
    pivot = keys((lo + hi) / 2)

    do
      do while (keys(i) < pivot)
        i = i + 1
      end do
      do while (keys(j) > pivot)
        j = j - 1
      end do

      if (i <= j) then
        tkey = keys(i)
        keys(i) = keys(j)
        keys(j) = tkey
        tval = vals(i)
        vals(i) = vals(j)
        vals(j) = tval
        i = i + 1
        j = j - 1
      end if

      if (i > j) exit
    end do

    if (lo < j) call quicksort_key_values(keys, vals, lo, j)
    if (i < hi) call quicksort_key_values(keys, vals, i, hi)
  end subroutine quicksort_key_values

  ! ------------------------------------------------------------------
  ! Exact tabular inbreeding + A-matrix builder (fallback)
  ! Computes full A (dense) by recursion and returns F.
  subroutine compute_inbreeding_tabular(n, sire, dam, F, writeA, filename)
    implicit none
    integer(kind=ki4), intent(in) :: n
    integer(kind=ki4), intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(out) :: F(n)
    logical, intent(in), optional :: writeA
    character(len=*), intent(in), optional :: filename

    real(kind=r8), allocatable :: A(:,:)
    integer :: i, j, si, di
    logical :: do_write
    character(len=256) :: fname
    integer :: unit, ios

    allocate(A(n,n))
    A = 0.0_r8
    F = 0.0_r8

    do i = 1, n
      si = sire(i)
      di = dam(i)
      if (si == 0 .and. di == 0) then
        ! founder
        A(i,i) = 1.0_r8
      else
        ! fill lower triangle A(i,1:i-1)
        do j = 1, i-1
          A(i,j) = 0.0_r8
          if (si /= 0) A(i,j) = A(i,j) + A(si,j)
          if (di /= 0) A(i,j) = A(i,j) + A(di,j)
          A(i,j) = 0.5_r8 * A(i,j)
          A(j,i) = A(i,j)
        end do
        if (si /= 0 .and. di /= 0) then
          F(i) = 0.5_r8 * A(si,di)
        else
          F(i) = 0.0_r8
        end if
        A(i,i) = 1.0_r8 + F(i)
      end if
    end do

    do_write = .false.
    if (present(writeA)) do_write = writeA
    if (do_write .and. present(filename)) then
      fname = filename
      open(newunit=unit, file=trim(fname), status='replace', action='write', iostat=ios)
      if (ios == 0) then
        do i = 1, n
          do j = 1, i
            if (abs(A(i,j)) > 0.0_r8) write(unit,'(I8,1x,I8,1x,ES16.8)') i, j, A(i,j)
          end do
        end do
        close(unit)
      else
        write(*,*) 'Cannot open triplet file:', trim(fname)
      end if
    end if

    deallocate(A)

  end subroutine compute_inbreeding_tabular

end module colleau_mod

