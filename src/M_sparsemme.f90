module M_SparseMME
  !---------------------------------------------------------------------------
  ! MME(Mixed Model Equations) 구축·풀이 전용 모듈
  !
  !   add_ainv    – Henderson A-inverse 기여분 누적
  !   csr_matvec  – 대칭 CSR 행렬-벡터 곱 (y = A*x)
  !   solve_pcg   – 대각 전처리 PCG 풀이기
  !
  ! 의존: M_HashSort (hash_table_t, sparse_csr_t, hash_add 등)
  !---------------------------------------------------------------------------
  use M_Kinds,    only: r8
  use M_HashSort, only: hash_table_t, sparse_csr_t, hash_add
  implicit none
  private

  public :: add_ainv
  public :: csr_matvec
  public :: solve_pcg

contains

  !==========================================================================
  !  add_ainv — Henderson A-inverse 기여분을 해시 테이블에 누적
  !
  !  Quaas (1976):
  !    양친기지:  d_i = 1 / [0.5  - 0.25*(F_s + F_d)]
  !    부만기지:  d_i = 1 / [0.75 - 0.25*F_s]
  !    모만기지:  d_i = 1 / [0.75 - 0.25*F_d]
  !    양친불명:  d_i = 1
  !
  !  f_sire, f_dam은 optional (생략 시 F=0 가정)
  !==========================================================================
  subroutine add_ainv(ht, animal, sire, dam, f_sire, f_dam)
    type(hash_table_t), intent(inout) :: ht
    integer,  intent(in) :: animal, sire, dam
    real(r8), intent(in), optional :: f_sire, f_dam

    real(r8) :: fs, fd, di

    fs = 0.0_r8;  if (present(f_sire)) fs = f_sire
    fd = 0.0_r8;  if (present(f_dam))  fd = f_dam

    if (sire > 0 .and. dam > 0) then
      di = 1.0_r8 / (0.5_r8 - 0.25_r8 * (fs + fd))
    else if (sire > 0) then
      di = 1.0_r8 / (0.75_r8 - 0.25_r8 * fs)
    else if (dam > 0) then
      di = 1.0_r8 / (0.75_r8 - 0.25_r8 * fd)
    else
      di = 1.0_r8
    end if

    call hash_add(ht, animal, animal, di)

    if (sire > 0) then
      call hash_add(ht, animal, sire, -di * 0.5_r8)
      call hash_add(ht, sire,  sire,   di * 0.25_r8)
    end if

    if (dam > 0) then
      call hash_add(ht, animal, dam,  -di * 0.5_r8)
      call hash_add(ht, dam,   dam,    di * 0.25_r8)
    end if

    if (sire > 0 .and. dam > 0) then
      call hash_add(ht, sire, dam,     di * 0.25_r8)
    end if
  end subroutine add_ainv

  !==========================================================================
  !  csr_matvec — y = A*x  (대칭 상삼각 CSR)
  !==========================================================================
  subroutine csr_matvec(csr, x, y)
    type(sparse_csr_t), intent(in) :: csr
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    integer :: i, j, k
    real(r8) :: aij

    y(1:csr%n) = 0.0_r8
    do i = 1, csr%n
      do k = csr%ia(i), csr%ia(i + 1) - 1
        j = csr%ja(k);  aij = csr%a(k)
        y(i) = y(i) + aij * x(j)
        if (i /= j) y(j) = y(j) + aij * x(i)
      end do
    end do
  end subroutine csr_matvec

  !==========================================================================
  !  solve_pcg — Preconditioned Conjugate Gradient (대각 전처리)
  !==========================================================================
  subroutine solve_pcg(csr, rhs, sol, maxiter, tol, niter, resid)
    type(sparse_csr_t), intent(in) :: csr
    real(r8), intent(in)    :: rhs(:)
    real(r8), intent(inout) :: sol(:)
    integer,  intent(in)    :: maxiter
    real(r8), intent(in)    :: tol
    integer,  intent(out)   :: niter
    real(r8), intent(out)   :: resid

    real(r8), allocatable :: r(:), z(:), p(:), Ap(:), diag(:)
    real(r8) :: rz, rz_old, alpha, beta, bnorm
    integer  :: i, k, n

    n = csr%n
    allocate(r(n), z(n), p(n), Ap(n), diag(n))

    ! 대각 전처리기 추출
    do i = 1, n
      diag(i) = 0.0_r8
      do k = csr%ia(i), csr%ia(i + 1) - 1
        if (csr%ja(k) == i) then;  diag(i) = csr%a(k);  exit;  end if
      end do
      if (diag(i) == 0.0_r8) diag(i) = 1.0_r8
    end do

    ! r = b - A*x0
    call csr_matvec(csr, sol, Ap)
    r = rhs(1:n) - Ap
    bnorm = sqrt(dot_product(rhs(1:n), rhs(1:n)))
    if (bnorm == 0.0_r8) bnorm = 1.0_r8

    ! z = M^{-1} r,  p = z
    z = r / diag;  p = z;  rz = dot_product(r, z)

    do niter = 1, maxiter
      call csr_matvec(csr, p, Ap)
      alpha = rz / dot_product(p, Ap)
      sol(1:n) = sol(1:n) + alpha * p
      r = r - alpha * Ap

      resid = sqrt(dot_product(r, r)) / bnorm
      if (resid < tol) exit

      z = r / diag
      rz_old = rz;  rz = dot_product(r, z)
      beta = rz / rz_old
      p = z + beta * p
    end do

    deallocate(r, z, p, Ap, diag)
  end subroutine solve_pcg

end module M_SparseMME
