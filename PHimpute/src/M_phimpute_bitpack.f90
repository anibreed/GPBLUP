module M_phimpute_bitpack
  use M_kinds
  implicit none

  integer(kind=ki1), parameter :: GENO_MISSING = 3_ki1

contains

  integer(kind=ki1) function get_geno_value(packed_row, snp_idx) result(val)
    integer(kind=ki4), intent(in) :: packed_row(:)
    integer, intent(in) :: snp_idx
    integer :: word_idx, bit_shift
    integer(kind=ki4) :: word_val

    if (snp_idx <= 0) then
      val = GENO_MISSING
      return
    end if
    word_idx = (snp_idx - 1) / 16 + 1
    if (word_idx > size(packed_row)) then
      val = GENO_MISSING
      return
    end if
    bit_shift = 2 * mod(snp_idx - 1, 16)
    word_val = packed_row(word_idx)
    val = int(iand(ishft(word_val, -bit_shift), 3_ki4), kind=ki1)
  end function get_geno_value

  subroutine set_geno_value(packed_row, snp_idx, val)
    integer(kind=ki4), intent(inout) :: packed_row(:)
    integer, intent(in) :: snp_idx
    integer(kind=ki1), intent(in) :: val
    integer :: word_idx, bit_shift
    integer(kind=ki4) :: mask

    if (snp_idx <= 0) return
    word_idx = (snp_idx - 1) / 16 + 1
    if (word_idx > size(packed_row)) return
    bit_shift = 2 * mod(snp_idx - 1, 16)
    mask = ishft(3_ki4, bit_shift)
    packed_row(word_idx) = iand(packed_row(word_idx), not(mask))
    packed_row(word_idx) = ior(packed_row(word_idx), ishft(int(val, ki4), bit_shift))
  end subroutine set_geno_value

  subroutine pack_geno_string(geno_str, n_snps, packed_row)
    character(len=*), intent(in) :: geno_str
    integer, intent(in) :: n_snps
    integer(kind=ki4), intent(inout) :: packed_row(:)
    integer :: s, word_idx, bit_shift
    integer(kind=ki4) :: val
    character(len=1) :: c

    packed_row = 0
    do s = 1, n_snps
      if (s > len_trim(geno_str)) then
        val = GENO_MISSING
      else
        c = geno_str(s:s)
        select case (c)
        case ('0')
          val = 0_ki4
        case ('1')
          val = 1_ki4
        case ('2')
          val = 2_ki4
        case default
          val = GENO_MISSING
        end select
      end if
      word_idx = (s - 1) / 16 + 1
      bit_shift = 2 * mod(s - 1, 16)
      packed_row(word_idx) = ior(packed_row(word_idx), ishft(val, bit_shift))
    end do
  end subroutine pack_geno_string

  subroutine pack_genotypes(geno, packed)
    integer(kind=ki1), intent(in) :: geno(:,:)
    integer(kind=ki4), allocatable, intent(out) :: packed(:,:)
    integer :: n_animals, n_snps, n_words

    n_animals = size(geno, 1)
    n_snps = size(geno, 2)
    n_words = (n_snps + 15) / 16

    allocate(packed(n_words, n_animals))
    packed = 0

    ! TODO: Implement 2-bit packing.
  end subroutine pack_genotypes

  subroutine unpack_genotypes(packed, n_snps, geno)
    integer(kind=ki4), intent(in) :: packed(:,:)
    integer, intent(in) :: n_snps
    integer(kind=ki1), allocatable, intent(out) :: geno(:,:)
    integer :: n_animals

    n_animals = size(packed, 2)
    allocate(geno(n_animals, n_snps))
    geno = GENO_MISSING

    ! TODO: Implement 2-bit unpacking.
  end subroutine unpack_genotypes

end module M_phimpute_bitpack
