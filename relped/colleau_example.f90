program colleau_example
  use colleau_mod
  implicit none
  integer :: n, i
  integer, allocatable :: sire(:), dam(:)
  real(selected_real_kind(15,307)), allocatable :: F(:), x(:), y(:)

  ! 간단한 예제: 작은 pedigree
  n = 8
  allocate(sire(n), dam(n), F(n), x(n), y(n))

  ! 1-based index, 0 = unknown
  ! 예: founding animals 1..2 are founders
  sire = (/ 0, 0, 1, 1, 3, 3, 5, 5 /)
  dam  = (/ 0, 0, 2, 2, 4, 4, 6, 6 /)

  print *, 'Pedigree (i:sire,dam):'
  do i = 1, n
    print '(I2, 2X, I2, 2X, I2)', i, sire(i), dam(i)
  end do

  call compute_inbreeding(n, sire, dam, F)

  print *, 'Inbreeding coefficients F:'
  do i = 1, n
    print '(I2, 2X, F10.6)', i, F(i)
  end do

  ! x 벡터 예시 (임의값)
  do i = 1, n
    x(i) = real(i, kind=selected_real_kind(15,307))
  end do

  call colleau_apply(n, sire, dam, F, x, y)

  print *, 'Result y = A * x (Colleau indirect):'
  do i = 1, n
    print '(I2, 2X, F12.6)', i, y(i)
  end do

  deallocate(sire, dam, F, x, y)
end program colleau_example
