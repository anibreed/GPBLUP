module M_phimpute_output
  use M_kinds
  use M_Variables, only: LEN_STR
  use M_phimpute_data, only: geno_matrix_t
  use M_phimpute_bitpack, only: get_geno_value, GENO_MISSING
  implicit none

  public :: write_geno_file

contains

  subroutine write_geno_file(file_name, geno_mat, animal_ids)
    character(len=*), intent(in) :: file_name
    type(geno_matrix_t), intent(in) :: geno_mat
    character(len=LEN_STR), intent(in) :: animal_ids(:)

    integer :: u_fd, i, s
    character(len=:), allocatable :: geno_str
    integer(kind=ki1) :: g

    if (geno_mat%n_animals <= 0 .or. geno_mat%n_snps <= 0) return
    if (size(animal_ids) < geno_mat%n_animals) return

    u_fd = 99
    open(unit=u_fd, file=trim(file_name), status='replace', action='write')

    allocate(character(len=geno_mat%n_snps) :: geno_str)

    do i = 1, geno_mat%n_animals
      do s = 1, geno_mat%n_snps
        g = get_geno_value(geno_mat%packed(:, i), s)
        select case (g)
        case (0)
          geno_str(s:s) = '0'
        case (1)
          geno_str(s:s) = '1'
        case (2)
          geno_str(s:s) = '2'
        case default
          geno_str(s:s) = '9'
        end select
      end do
      write(u_fd,'(A,1X,A)') trim(animal_ids(i)), geno_str
    end do

    deallocate(geno_str)
    close(u_fd)
  end subroutine write_geno_file

end module M_phimpute_output
