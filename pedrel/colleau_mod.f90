module colleau_mod
  use M_param
  implicit none
  integer(kind=ki4), save :: call_count = 0
  integer(kind=ki4), parameter :: PROG_STEP = 100
  integer(kind=ki4), save :: max_ns = 0
  integer(kind=ki4), save :: max_nd = 0
  logical, save :: verbose = .false.
  integer(kind=ki4), allocatable, save :: dbg_targets(:)
  integer(kind=ki4), save :: n_dbg_targets = 0
  
  ! Colleau 관련 서브루틴 모듈
  ! use real(kind=r8) for double precision
contains

  subroutine set_colleau_verbose(flag)
    logical, intent(in) :: flag
    verbose = flag
  end subroutine set_colleau_verbose

  subroutine set_colleau_debug_targets(n, targets)
    integer(kind=ki4), intent(in) :: n
    integer(kind=ki4), intent(in) :: targets(:)
    if (allocated(dbg_targets)) deallocate(dbg_targets)
    if (n > 0) then
      allocate(dbg_targets(n))
      dbg_targets = targets(1:n)
      n_dbg_targets = n
    else
      n_dbg_targets = 0
    end if
  end subroutine set_colleau_debug_targets


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

    nanc = 1
    anc(1) = start
    val(1) = 1.0_r8

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

      if (s > 0) then
        call accum_local(s, 0.5_r8*curr_w, anc, val, nanc)
        top = top + 1
        if (top > maxanc) stop 'stack overflow in build_contrib_sparse_local'
        stack(top) = s
        wstack(top) = 0.5_r8*curr_w
      end if

      if (d > 0) then
        call accum_local(d, 0.5_r8*curr_w, anc, val, nanc)
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
    integer :: outunit, ios
    character(len=128) :: fname
    integer(kind=ki4), allocatable :: anc_s(:), anc_d(:)
    real(kind=r8), allocatable :: val_s(:), val_d(:)
    real(kind=r8) :: a_sd

    F = 0.0_r8

    maxanc = min(n, MAXANC_LIMIT)

!$omp parallel do default(shared) private(i,n_s,n_d,anc_s,anc_d,val_s,val_d,a_sd)
    do i = 1, n
      if (sire(i) == 0 .or. dam(i) == 0) then
        F(i) = 0.0_r8
      else
        n_s = 0; n_d = 0
        allocate(anc_s(maxanc)); allocate(anc_d(maxanc))
        allocate(val_s(maxanc)); allocate(val_d(maxanc))
        anc_s = 0; anc_d = 0
        val_s = 0.0_r8; val_d = 0.0_r8
        ! normal execution (no per-target debug)
        call build_contrib_sparse_local(sire(i), sire, dam, anc_s, val_s, n_s, maxanc)
        call build_contrib_sparse_local(dam(i),  sire, dam, anc_d, val_d, n_d, maxanc)
        ! no pre-contrib debug dump
        a_sd = dot_sparse_local(anc_s, val_s, n_s, anc_d, val_d, n_d)
        F(i) = 0.5_r8 * a_sd
        ! If this index is among debug targets, print detailed info
      !$omp critical
        if (n_dbg_targets > 0) then
          ! declare locals outside inner loops
        end if
      !$omp end critical
      if (n_dbg_targets > 0) then
        !$omp critical
        do ttt = 1, n_dbg_targets
          if (i == dbg_targets(ttt)) then
            write(*,'(A,I8,2X,A,I8)') 'DBG: compute_inbreeding: index=', i, ' n_s=', n_s
            write(*,'(A,ES24.16,2X,A,ES24.16)') 'DBG: a_sd=', a_sd, ' F(i)=', F(i)
            write(*,'(A)') 'DBG: sample anc_s (idx,anc,val_s) up to 20:'
            do jj = 1, min(20,n_s)
              write(*,'(I6,1X,I8,1X,ES24.16)') jj, anc_s(jj), val_s(jj)
            end do
            write(*,'(A)') 'DBG: sample anc_d (idx,anc,val_d) up to 20:'
            do jj = 1, min(20,n_d)
              write(*,'(I6,1X,I8,1X,ES24.16)') jj, anc_d(jj), val_d(jj)
            end do
          end if
        end do
        !$omp end critical
      end if
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

end module colleau_mod
