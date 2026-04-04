module inbreed
  use M_Kinds, only: r8
  use link_heap
  implicit none

  type(Heap) :: Arg_AS
  type(Heap) :: Arg_AD

contains

  subroutine Inb_F(NPD, F)
    integer, intent(in) :: NPD(:,:)
    real(kind=r8), intent(inout) :: F(size(NPD,2))
    real(kind=r8) :: LS(size(NPD,2)), LD(size(NPD,2)), D(size(NPD,2))
    logical :: in_heap_s(size(NPD,2)), in_heap_d(size(NPD,2))
    integer :: i, j, k, si, di, sj, dj, sk, dk, N

    N = size(NPD,2)
    D = 0.0
    F = 0.0

    call heap_init(Arg_AS, max(128, N / 8))
    call heap_init(Arg_AD, max(128, N / 8))

    do i = 1, N
      LS = 0.0
      LD = 0.0
      in_heap_s = .false.
      in_heap_d = .false.
      Arg_AS%size = 0
      Arg_AD%size = 0

      si = NPD(1,i)
      if (si /= 0) then
        LS(si) = 1.0
        call insertL(Arg_AS, si)
        in_heap_s(si) = .true.
      end if

      di = NPD(2,i)
      if (di /= 0) then
        LD(di) = 1.0
        call insertL(Arg_AD, di)
        in_heap_d(di) = .true.
      end if

      if ((si > 0) .and. (di > 0)) then
        D(i) = 0.5 - 0.25 * (F(si) + F(di))
      else if (si > 0) then
        D(i) = 0.75 - 0.25 * F(si)
      else if (di > 0) then
        D(i) = 0.75 - 0.25 * F(di)
      else
        D(i) = 1.0
      end if

      do while (.not. heap_empty(Arg_AS) .or. .not. heap_empty(Arg_AD))
        j = getL(Arg_AS)
        k = getL(Arg_AD)

        if (j > k) then
          sj = NPD(1,j)
          if (sj > 0) then
            if (.not. in_heap_s(sj)) then
              call insertL(Arg_AS, sj)
              in_heap_s(sj) = .true.
            end if
            LS(sj) = LS(sj) + 0.5 * LS(j)
          end if
          dj = NPD(2,j)
          if (dj > 0) then
            if (.not. in_heap_s(dj)) then
              call insertL(Arg_AS, dj)
              in_heap_s(dj) = .true.
            end if
            LS(dj) = LS(dj) + 0.5 * LS(j)
          end if
          call deleteL(Arg_AS)
          in_heap_s(j) = .false.

        else if (k > j) then
          sk = NPD(1,k)
          if (sk > 0) then
            if (.not. in_heap_d(sk)) then
              call insertL(Arg_AD, sk)
              in_heap_d(sk) = .true.
            end if
            LD(sk) = LD(sk) + 0.5 * LD(k)
          end if
          dk = NPD(2,k)
          if (dk > 0) then
            if (.not. in_heap_d(dk)) then
              call insertL(Arg_AD, dk)
              in_heap_d(dk) = .true.
            end if
            LD(dk) = LD(dk) + 0.5 * LD(k)
          end if
          call deleteL(Arg_AD)
          in_heap_d(k) = .false.

        else
          sj = NPD(1,j)
          if (sj > 0) then
            if (.not. in_heap_s(sj)) then
              call insertL(Arg_AS, sj)
              in_heap_s(sj) = .true.
            end if
            if (.not. in_heap_d(sj)) then
              call insertL(Arg_AD, sj)
              in_heap_d(sj) = .true.
            end if
            LS(sj) = LS(sj) + 0.5 * LS(j)
            LD(sj) = LD(sj) + 0.5 * LD(j)
          end if

          dj = NPD(2,j)
          if (dj > 0) then
            if (.not. in_heap_s(dj)) then
              call insertL(Arg_AS, dj)
              in_heap_s(dj) = .true.
            end if
            if (.not. in_heap_d(dj)) then
              call insertL(Arg_AD, dj)
              in_heap_d(dj) = .true.
            end if
            LS(dj) = LS(dj) + 0.5 * LS(j)
            LD(dj) = LD(dj) + 0.5 * LD(j)
          end if

          F(i) = F(i) + LS(j) * LD(j) * 0.5 * D(j)
          call deleteL(Arg_AS)
          call deleteL(Arg_AD)
          in_heap_s(j) = .false.
          in_heap_d(j) = .false.
        end if
      end do

      if (mod(i, 1000) == 0) print *, i, ' th Animal calculated F ', F(i) * 100
    end do

    call heap_clear(Arg_AS)
    call heap_clear(Arg_AD)
  end subroutine Inb_F

end module inbreed
