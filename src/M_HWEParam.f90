module M_HWEParam
  use M_Kinds
  use M_ReadFile
  implicit none
  public :: get_hwe_qc_settings, extract_value
contains
  subroutine get_hwe_qc_settings(param_file, out_method, out_threshold)
    character(len=*), intent(in) :: param_file
    character(len=16), intent(out) :: out_method
    real(r8), intent(out) :: out_threshold
    character(len=256) :: line
    integer :: unit, found_method, found_threshold, eqpos
    out_method = "FDR"
    out_threshold = 1.0e-3_r8
    found_method = 0
    found_threshold = 0
    unit = fopen(trim(param_file))
    do
      read(unit, '(A)', end=100) line
      eqpos = index(line, "HWE_METHOD")
      if (eqpos > 0) then
        out_method = adjustl(extract_value(line))
        found_method = 1
      end if
      eqpos = index(line, "HWE_THRESHOLD")
      if (eqpos > 0) then
        read(line(index(line, "=") + 1:), *) out_threshold
        found_threshold = 1
      end if
    end do
100 continue
    close(unit)
    if (found_method == 0) out_method = "FDR"
    if (found_threshold == 0) out_threshold = 1.0e-3_r8
  end subroutine get_hwe_qc_settings

  function extract_value(line) result(val)
    character(len=*), intent(in) :: line
    character(len=16) :: val
    integer :: eqpos
    eqpos = index(line, "=")
    if (eqpos > 0) then
      val = adjustl(line(eqpos+1:))
    else
      val = ""
    end if
  end function extract_value
end module M_HWEParam
