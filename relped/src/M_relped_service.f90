module M_relped_service
  use renPED, only: renumped
  use M_readpar_relped, only: relped_params_t, read_relped_params
  implicit none

  public :: run_relped_from_file, run_relped_params

contains

  subroutine run_relped_from_file(param_file, use_gpu, device_id)
    character(len=*), intent(in) :: param_file
    logical, intent(in) :: use_gpu
    integer, intent(in) :: device_id
    type(relped_params_t) :: params

    call read_relped_params(param_file, params)
    call run_relped_params(params, use_gpu, device_id)
  end subroutine run_relped_from_file

  subroutine run_relped_params(params, use_gpu, device_id)
    type(relped_params_t), intent(in) :: params
    logical, intent(in) :: use_gpu
    integer, intent(in) :: device_id


    if (params%has_reference) then
      if (params%refinfo(1) /= 0) then
        call renumped(PEDFile=params%pedfile, PEDinfo=params%pedinfo, REFFile=params%reffile, &
          REFinfo=params%refinfo, inb=params%inbreed, rel=params%relation, rel_inv=params%relation_inv, &
          rel_outfile=params%relation_file, relinv_outfile=params%relation_inv_file, use_gpu=use_gpu, &
          device_id=device_id, omp_threads_in=params%omp_threads)
      else
        call renumped(PEDFile=params%pedfile, PEDinfo=params%pedinfo, REFFile=params%reffile, &
          inb=params%inbreed, rel=params%relation, rel_inv=params%relation_inv, &
          rel_outfile=params%relation_file, relinv_outfile=params%relation_inv_file, use_gpu=use_gpu, &
          device_id=device_id, omp_threads_in=params%omp_threads)
      end if
    else
      call renumped(PEDFile=params%pedfile, PEDinfo=params%pedinfo, inb=params%inbreed, &
        rel=params%relation, rel_inv=params%relation_inv, rel_outfile=params%relation_file, &
        relinv_outfile=params%relation_inv_file, use_gpu=use_gpu, device_id=device_id, &
        omp_threads_in=params%omp_threads)
    end if
  end subroutine run_relped_params

end module M_relped_service
