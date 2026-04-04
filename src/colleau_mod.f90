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
    integer(kind=ki4) :: stack_cap
    integer(kind=ki4) :: top, k, s, d, i
    integer(kind=ki4) :: dropped_push
    real(kind=r8) :: curr_w
    real(kind=r8), parameter :: MIN_PROP_WEIGHT = 1.0d-12

    stack_cap = max(1024_ki4, maxanc)
    allocate(stack(stack_cap))
    allocate(wstack(stack_cap))

    ! We'll traverse ancestors starting at `start`, but only record
    ! contributions for founders (animals with no parents). This builds
    ! k-vectors over founders which matches the tabular method.
    nanc = 0
    dropped_push = 0

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

      ! Prune negligible paths so ancestor propagation cannot explode in depth.
      if (abs(curr_w) <= MIN_PROP_WEIGHT) cycle

      ! If node has parents, push them for further expansion with half weight
      if (s > 0) then
        if (top < size(stack)) then
          top = top + 1
          stack(top) = s
          wstack(top) = 0.5_r8*curr_w
        else
          dropped_push = dropped_push + 1
        end if
      end if
      if (d > 0) then
        if (top < size(stack)) then
          top = top + 1
          stack(top) = d
          wstack(top) = 0.5_r8*curr_w
        else
          dropped_push = dropped_push + 1
        end if
      end if

    end do

    if (verbose .and. dropped_push > 0) then
      write(*,'(A,I0,2X,A,I0)') 'build_contrib_sparse_local: dropped pushes=', dropped_push, 'start=', start
    end if

  end subroutine build_contrib_sparse_local

  subroutine grow_local_stack(stack, wstack)
    implicit none
    integer(kind=ki4), allocatable, intent(inout) :: stack(:)
    real(kind=r8), allocatable, intent(inout) :: wstack(:)
    integer(kind=ki4), allocatable :: new_stack(:)
    real(kind=r8), allocatable :: new_wstack(:)
    integer(kind=ki4) :: old_n, new_n

    old_n = size(stack)
    new_n = max(2*old_n, old_n + 1024)

    allocate(new_stack(new_n))
    allocate(new_wstack(new_n))
    new_stack(1:old_n) = stack(1:old_n)
    new_wstack(1:old_n) = wstack(1:old_n)

    call move_alloc(new_stack, stack)
    call move_alloc(new_wstack, wstack)
  end subroutine grow_local_stack

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
    ! Inbreeding calculation with dependency-safe, level-wise parallelization.
    implicit none
    integer(kind=ki4), intent(in) :: n
    integer(kind=ki4), intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(out) :: F(n)

    integer(kind=ki4) :: i, n_s, n_d, maxanc
    integer(kind=ki4) :: k, jj
    integer(kind=ki4) :: anc_id, si, di
    integer(kind=ki4) :: lvl, max_level, n_curr, curr_ptr
    integer(kind=ki4) :: qh, qt, u
    integer(kind=ki4) :: edge_n, base
    integer(kind=ki4), allocatable :: indeg(:), level(:), order(:)
    integer(kind=ki4), allocatable :: child_count(:), child_start(:), child_fill(:), child_to(:)
    integer(kind=ki4), allocatable :: level_count(:), level_start(:), level_pos(:), level_nodes(:)
    real(kind=r8) :: dval
    integer(kind=ki4), allocatable :: anc_s(:), anc_d(:)
    real(kind=r8), allocatable :: val_s(:), val_d(:)
    real(kind=r8) :: a_sd
    real(kind=r8) :: t_start, t_end, t_elapsed, t_compute_total
    t_compute_total = 0.0_r8

    F = 0.0_r8

    maxanc = min(n, MAXANC_LIMIT)
    max_ns = 0
    max_nd = 0

    allocate(indeg(n), level(n), order(n), child_count(n), child_start(n), child_fill(n))
    indeg = 0
    level = 1
    child_count = 0

    do i = 1, n
      if (sire(i) > 0) then
        indeg(i) = indeg(i) + 1
        child_count(sire(i)) = child_count(sire(i)) + 1
      end if
      if (dam(i) > 0) then
        indeg(i) = indeg(i) + 1
        child_count(dam(i)) = child_count(dam(i)) + 1
      end if
    end do

    child_start(1) = 1
    do i = 2, n
      child_start(i) = child_start(i-1) + child_count(i-1)
    end do
    child_fill = 0
    edge_n = sum(child_count)
    allocate(child_to(max(1_ki4, edge_n)))

    do i = 1, n
      if (sire(i) > 0) then
        base = child_start(sire(i)) + child_fill(sire(i))
        child_to(base) = i
        child_fill(sire(i)) = child_fill(sire(i)) + 1
      end if
      if (dam(i) > 0) then
        base = child_start(dam(i)) + child_fill(dam(i))
        child_to(base) = i
        child_fill(dam(i)) = child_fill(dam(i)) + 1
      end if
    end do

    qh = 1
    qt = 0
    do i = 1, n
      if (indeg(i) == 0) then
        qt = qt + 1
        order(qt) = i
      end if
    end do

    n_curr = 0
    do while (qh <= qt)
      u = order(qh)
      qh = qh + 1
      n_curr = n_curr + 1

      if (child_count(u) > 0) then
        do k = child_start(u), child_start(u) + child_count(u) - 1
          i = child_to(k)
          level(i) = max(level(i), level(u) + 1)
          indeg(i) = indeg(i) - 1
          if (indeg(i) == 0) then
            qt = qt + 1
            order(qt) = i
          end if
        end do
      end if
    end do

    if (n_curr < n) then
      ! Fallback for malformed cyclic remnants: process remaining at highest level.
      max_level = maxval(level)
      do i = 1, n
        if (indeg(i) > 0) level(i) = max_level + 1
      end do
    end if
    max_level = maxval(level)

    allocate(level_count(max_level), level_start(max_level), level_pos(max_level), level_nodes(n))
    level_count = 0
    do i = 1, n
      level_count(level(i)) = level_count(level(i)) + 1
    end do

    level_start(1) = 1
    do lvl = 2, max_level
      level_start(lvl) = level_start(lvl-1) + level_count(lvl-1)
    end do
    level_pos = 0
    do i = 1, n
      lvl = level(i)
      k = level_start(lvl) + level_pos(lvl)
      level_nodes(k) = i
      level_pos(lvl) = level_pos(lvl) + 1
    end do

    do lvl = 1, max_level
      if (level_count(lvl) <= 0) cycle
!$omp parallel default(shared) &
!$omp private(i,n_s,n_d,anc_s,anc_d,val_s,val_d,a_sd,si,di,jj,k,anc_id,dval,t_start,t_end,t_elapsed,u,curr_ptr) &
!$omp reduction(+:t_compute_total) reduction(max:max_ns,max_nd)
      allocate(anc_s(maxanc), anc_d(maxanc))
      allocate(val_s(maxanc), val_d(maxanc))

!$omp do schedule(static)
      do curr_ptr = level_start(lvl), level_start(lvl) + level_count(lvl) - 1
        u = level_nodes(curr_ptr)
        if (sire(u) == 0 .or. dam(u) == 0) then
          F(u) = 0.0_r8
        else
          t_start = omp_get_wtime()
          n_s = 0
          n_d = 0
          call build_contrib_sparse_local(sire(u), sire, dam, anc_s, val_s, n_s, maxanc)
          call build_contrib_sparse_local(dam(u),  sire, dam, anc_d, val_d, n_d, maxanc)

          a_sd = 0.0_r8
          if (n_s > 0 .and. n_d > 0) then
            do jj = 1, n_s
              do k = 1, n_d
                if (anc_s(jj) == anc_d(k)) then
                  anc_id = anc_s(jj)
                  if (anc_id >= 1 .and. anc_id <= size(sire)) then
                    si = sire(anc_id)
                    di = dam(anc_id)
                    if (si >= 1 .and. si <= n .and. di >= 1 .and. di <= n) then
                      dval = 0.5_r8 - 0.25_r8*(F(si) + F(di))
                    else if (si >= 1 .and. si <= n) then
                      dval = 0.75_r8 - 0.25_r8*F(si)
                    else if (di >= 1 .and. di <= n) then
                      dval = 0.75_r8 - 0.25_r8*F(di)
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

          F(u) = 0.5_r8 * a_sd
          t_end = omp_get_wtime()
          t_elapsed = t_end - t_start
          t_compute_total = t_compute_total + t_elapsed

          if (n_s > max_ns) max_ns = n_s
          if (n_d > max_nd) max_nd = n_d

          if (verbose .and. mod(u, PROG_STEP) == 0) then
      !$omp critical(colleau_progress)
            if (verbose) then
              write(*,'(A,I8,2X,A,I8,2X,A,I8,2X,A,I8,A,F10.4)') &
                   'compute_inbreeding: processed', u, 'of', n, 'n_s=', n_s, 'n_d=', n_d, ' F=', F(u)
              if (n_s >= maxanc .or. n_d >= maxanc) then
                write(*,'(A)') 'WARNING: ancestor expansion reached maxanc limit; results may be truncated'
              end if
            end if
      !$omp end critical(colleau_progress)
          end if
        end if
      end do
!$omp end do

      deallocate(anc_s, anc_d, val_s, val_d)
!$omp end parallel
    end do

    deallocate(indeg, level, order, child_count, child_start, child_fill, child_to)
    deallocate(level_count, level_start, level_pos, level_nodes)
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
    integer(kind=ki4) :: j, i, ios, version
    integer(kind=ki4) :: row_out, col_out
    integer(kind=ki8) :: nnz, nnz_total, k
    real(kind=r8), allocatable :: x(:), y(:)
    integer :: unit, merge_u, nthreads, tid
    integer(kind=ki8), allocatable :: part_nnz(:)
    character(len=MAX_STR), allocatable :: part_files(:)
    real(kind=r8) :: val_in
    integer(kind=ki4) :: row_in, col_in
    character(len=4) :: magic
    character(len=MAX_STR) :: part_name
    real(kind=r8) :: eps
    character(len=32) :: eps_str
    integer :: ios_env

    eps = 0.0_r8
    eps_str = ''
    call get_environment_variable('RELPED_TRIPLET_EPS', eps_str, status=ios_env)
    if (ios_env == 0) then
      read(eps_str, *, iostat=ios) eps
      if (ios /= 0) eps = 0.0_r8
      if (eps < 0.0_r8) eps = 0.0_r8
    end if

    nthreads = max(1, omp_get_max_threads())
    allocate(part_files(0:nthreads-1), part_nnz(0:nthreads-1))
    part_nnz = 0_ki8

!$omp parallel default(shared) private(tid,unit,part_name,j,i,x,y,row_out,col_out,nnz)
    allocate(x(n), y(n))
    tid = omp_get_thread_num()
    write(part_name,'(A,".part.",I0)') trim(filename), tid
    part_files(tid) = trim(part_name)

    open(newunit=unit, file=trim(part_name), status='replace', action='write', &
      form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open temporary triplet file:', trim(part_name)
      stop 7
    end if

    nnz = 0_ki8
!$omp do schedule(static)
    do j = 1, n
      x = 0.0_r8
      x(j) = 1.0_r8
      call colleau_apply(n, sire, dam, F, x, y)
      do i = j, n
        if (abs(y(i)) > eps) then
          row_out = int(i, kind=ki4)
          col_out = int(j, kind=ki4)
          write(unit) row_out, col_out, y(i)
          nnz = nnz + 1_ki8
        end if
      end do
    end do
!$omp end do

    part_nnz(tid) = nnz
    close(unit)
    deallocate(x, y)
!$omp end parallel

    open(newunit=merge_u, file=filename, status='replace', action='readwrite', &
      form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open triplet file:', trim(filename)
      stop 7
    end if

    magic = 'ATRP'
    version = 1_ki4
    nnz_total = 0_ki8
    do tid = 0, nthreads - 1
      nnz_total = nnz_total + part_nnz(tid)
    end do
    write(merge_u) magic, version, int(n, kind=ki4), nnz_total

    do tid = 0, nthreads - 1
      open(newunit=unit, file=trim(part_files(tid)), status='old', action='read', &
        form='unformatted', access='stream', iostat=ios)
      if (ios /= 0) cycle
      do k = 1_ki8, part_nnz(tid)
        read(unit, iostat=ios) row_in, col_in, val_in
        if (ios /= 0) exit
        write(merge_u) row_in, col_in, val_in
      end do
      close(unit, status='delete')
    end do

    close(merge_u)
    deallocate(part_files, part_nnz)
  end subroutine write_triplets

  subroutine write_triplets_text_mod(n, sire, dam, F, filename)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(in) :: F(n)
    character(len=*), intent(in) :: filename
    integer(kind=ki4) :: j, i, ios
    integer(kind=ki4) :: row_out, col_out
    integer :: unit, merge_u, nthreads, tid
    integer(kind=ki8), allocatable :: part_nnz(:)
    character(len=MAX_STR), allocatable :: part_files(:)
    character(len=MAX_STR) :: part_name, line_buf
    real(kind=r8), allocatable :: x(:), y(:)

    nthreads = max(1, omp_get_max_threads())
    allocate(part_files(0:nthreads-1), part_nnz(0:nthreads-1))
    part_nnz = 0_ki8

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
        if (abs(y(i)) > 0.0_r8) then
          row_out = int(i, kind=ki4)
          col_out = int(j, kind=ki4)
          write(unit, '(I0,1X,I0,1X,ES24.16)') row_out, col_out, y(i)
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
    deallocate(part_files, part_nnz)
  end subroutine write_triplets_text_mod

  subroutine write_ainv_triplets(n, sire, dam, F, filename)
    implicit none
    integer(kind=ki4), intent(in) :: n
    integer(kind=ki4), intent(in) :: sire(n), dam(n)
    real(kind=r8), intent(in) :: F(n)
    character(len=*), intent(in) :: filename

    integer(kind=ki4) :: i, s, d, unit, ios, version
    integer(kind=ki4) :: nt, max_terms
    integer(kind=ki4), allocatable :: rows(:), cols(:)
    integer(kind=ki8), allocatable :: keys(:)
    real(kind=r8), allocatable :: vals(:)
    real(kind=r8) :: di, bi
    real(kind=r8) :: v_acc
    integer(kind=ki8) :: nnz
    integer(kind=ki8) :: key_acc
    integer(kind=ki4) :: row_acc, col_acc
    character(len=4) :: magic

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
      open(newunit=unit, file=filename, status='replace', action='readwrite', &
        form='unformatted', access='stream', iostat=ios)
      if (ios == 0) then
        magic = 'AINV'
        version = 1_ki4
        nnz = 0_ki8
        write(unit) magic, version, int(n, kind=ki4), nnz
        close(unit)
      end if
      deallocate(rows, cols, vals)
      return
    end if

    allocate(keys(nt))
    do i = 1, nt
      keys(i) = int(rows(i), kind=ki8) * int(n + 1, kind=ki8) + int(cols(i), kind=ki8)
    end do
    call sort_key_values(keys, vals, nt)

    open(newunit=unit, file=filename, status='replace', action='readwrite', &
      form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open inverse triplet file:', trim(filename)
      deallocate(rows, cols, vals, keys)
      return
    end if

    ! Binary stream header:
    ! [char4 magic="AINV"][int32 version][int32 n][int64 nnz]
    magic = 'AINV'
    version = 1_ki4
    nnz = 0_ki8
    write(unit) magic, version, int(n, kind=ki4), nnz

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
            write(unit) row_acc, col_acc, v_acc
            nnz = nnz + 1_ki8
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
        write(unit) row_acc, col_acc, v_acc
        nnz = nnz + 1_ki8
      end if
    end if

    ! Backfill final nnz in header.
    write(unit, pos=1) magic, version, int(n, kind=ki4), nnz

    close(unit)
    deallocate(rows, cols, vals, keys)
  end subroutine write_ainv_triplets

  subroutine write_ainv_triplets_text_mod(n, sire, dam, F, filename)
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
        if (abs(v_acc) > 1.0d-16) then
          row_acc = int(key_acc / int(n + 1, kind=ki8), kind=ki4)
          col_acc = int(mod(key_acc, int(n + 1, kind=ki8)), kind=ki4)
          if (row_acc >= 1 .and. col_acc >= 1) then
            write(unit, '(I0,1X,I0,1X,ES24.16)') row_acc, col_acc, v_acc
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
        write(unit, '(I0,1X,I0,1X,ES24.16)') row_acc, col_acc, v_acc
      end if
    end if

    close(unit)
    deallocate(rows, cols, vals, keys)
  end subroutine write_ainv_triplets_text_mod

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

