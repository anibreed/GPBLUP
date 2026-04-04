module M_QCGenoIO
  use M_IO_Geno, only: load_geno_file_mod, write_valid_geno_for_imputation_mod, &
      set_qcgenoio_debug_mod
  implicit none

  public :: load_geno_file_mod, write_valid_geno_for_imputation_mod, set_qcgenoio_debug_mod
end module M_QCGenoIO
