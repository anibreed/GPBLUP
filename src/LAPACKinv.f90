module LAPACKinv
implicit none
  ! External procedures defined in LAPACK
private
save
public LAPACK_inv
contains
  function LAPACK_inv(A) result(Ainv)
  double precision:: A(:,:)
  double precision:: Ainv(size(A,1),size(A,2))
  double precision,allocatable:: WORK(:)  ! work array for LAPACK
  integer,allocatable:: IPIV(:)       ! pivot indices
  integer :: N, M, INFO, LDA
  integer :: LWORK, i
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  N = size(A,1)
  if(N /= size(A,2)) then
        PRINT '(" Matrix is not square matrix ")'
        STOP
  endif
  LDA = N
  LWORK=N*N
  ALLOCATE (WORK(LWORK))
  ALLOCATE (IPIV(N))
  M=N
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  CALL DGETRF( M, N, Ainv, LDA, IPIV, INFO )

  select case(INFO)
         CASE ( : -1)
            PRINT '(" LU decomposition:  illegal value ")'
            STOP
       CASE (0)
            PRINT '(" LU decomposition successful ")'
       CASE (1 : )
            WRITE(*,35)INFO,INFO
            35       FORMAT( 'LU decomposition: U(',I4,',',I4,') = 0 ')
  END SELECT

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(N, Ainv, N, IPIV, WORK, LWORK, INFO)

  if (INFO /= 0) then
       stop 'Matrix inversion failed!'
  else
       PRINT '(" Inverse Successful ")'
  end if
  end function

end module
