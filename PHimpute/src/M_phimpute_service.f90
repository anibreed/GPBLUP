module M_phimpute_service
  use M_phimpute_params, only: phimpute_params_t
  use M_readpar_phimpute, only: read_phimpute_params
  use M_phimpute_steps, only: run_phimpute_pipeline
  implicit none

  public :: run_phimpute_from_file, run_phimpute_params

contains

  subroutine run_phimpute_from_file(param_file)
    character(len=*), intent(in) :: param_file
    type(phimpute_params_t) :: params

    call read_phimpute_params(param_file, params)
    call run_phimpute_params(params)
  end subroutine run_phimpute_from_file

  subroutine run_phimpute_params(params)
    type(phimpute_params_t), intent(in) :: params

    call run_phimpute_pipeline(params)
  end subroutine run_phimpute_params

end module M_phimpute_service
