program test_readpar
  use M_kinds
  use M_Variables
  use M_ReadFile
  use M_StrEdit
  use M_readpar
  implicit none
  
  character(len=512) :: param_file
  logical :: file_exists
  
  ! Parameter file path
  param_file = '/home/dhlee/DATA/FinalReport/param'
  
  ! Check if file exists
  inquire(file=trim(param_file), exist=file_exists)
  if (.not. file_exists) then
    write(*,'(A)') "ERROR: Parameter file not found: "//trim(param_file)
    stop 1
  end if
  
  write(*,'(A)') "=========================================="
  write(*,'(A)') "  M_readpar Module Test"
  write(*,'(A)') "=========================================="
  write(*,'(A)') ""
  write(*,'(A)') "Testing file: "//trim(param_file)
  write(*,'(A)') ""
  
  ! Call read_to_sentences to load into memory
  call read_to_sentences(param_file)

  ! Call legacy read_parameters (which still uses file reading for now)
  call read_parameters(param_file)
  
  ! Verify loaded parameters
  write(*,'(A)') ""
  write(*,'(A)') "=========================================="
  write(*,'(A)') "  Loaded Parameters Summary"
  write(*,'(A)') "=========================================="
  write(*,'(A)') ""
  
  call print_detail_info("PEDFile", PEDFile)
  call print_detail_info("DATAFile", DATAFile)
  call print_detail_info("MapFile", MapFile)
  call print_detail_info("FRFile", FRFile)
  call print_detail_info("GENOFile", GENOFile)
  
  write(*,'(A)') ""
  write(*,'(A)') "=========================================="
  write(*,'(A)') "  Quality Control Parameters"
  write(*,'(A)') "=========================================="
  write(*,'(A)') ""
  
  write(*,'(A)') "  FRQC (Frequency Quality Control):"
  write(*,'(A)') "    - Status: Initialized"
  write(*,'(A)') ""
  call print_frqc_info(FRQC)
  write(*,'(A)') "  POPQC (Population Quality Control):"
  write(*,'(A)') "    - Status: Initialized"
  write(*,'(A)') ""
  call print_popqc_info(POPQC)
  write(*,'(A)') "=========================================="
  write(*,'(A)') "  Test Complete - SUCCESS"
  write(*,'(A)') "=========================================="
  write(*,'(A)') ""
  write(*,'(A)') "=========================================="
  write(*,'(A)') "  output file Parameters"
  write(*,'(A)') "=========================================="
  write(*,'(A)') ""
  call print_ofile_info("FROUTFile", FROUTFile)
  call print_ofile_info("GENOUTFile", GENOUTFile)

  contains
    subroutine print_detail_info(name, info)
      character(len=*), intent(in) :: name
      type(FileInfo), intent(in) :: info
      integer :: i
      
      write(*,'(A)') "  "//trim(name)//":"
      if (len_trim(info%FileName) > 0) then
        write(*,'(A,A)') "    - File: ", trim(info%FileName)
        write(*,'(A,I0)') "    - Header lines: ", info%Header
        write(*,'(A,A)') "    - Delimiter: '", trim(info%Delim_char)//"'"
        write(*,'(A,I0)') "    - Variables: ", info%NVAR
        do i = 1, info%NVAR
          write(*,'(A,I2,A,A,A,I4)') "      [", i, "] ", trim(info%FieldName(i)), " at pos ", info%FieldLoc(i)
        end do
      else
        write(*,'(A)') "    - NOT CONFIGURED"
      end if
      write(*,'(A)') ""
    end subroutine print_detail_info

    subroutine print_frqc_info(qc_struct)
      type(FRQC_matrices), intent(in) :: qc_struct
      write(*,'(A,F10.4)') "    - GC Score:      ", qc_struct%gc_score
      write(*,'(A,F10.4,A,F10.4)') "    - R Intensity:   ", qc_struct%r_intensity(1), " to ", qc_struct%r_intensity(2)
      write(*,'(A,F10.4)') "    - GT Score:      ", qc_struct%gt_score
      write(*,'(A,F10.4)') "    - Cluster Sep:   ", qc_struct%cluster_sep
      write(*,'(A,F10.4)') "    - Animal CR:     ", qc_struct%animal_cr
      write(*,'(A,L2)')    "    - QC Passed:     ", qc_struct%qc_passed
      write(*,'(A)') ""
    end subroutine print_frqc_info

    subroutine print_popqc_info(qc_struct)
      type(POPQC_matrices), intent(in) :: qc_struct
      write(*,'(A,F10.4)') "    - SNP CR:        ", qc_struct%snp_cr
      write(*,'(A,F10.4)') "    - Animal CR:     ", qc_struct%animal_cr
      write(*,'(A,A)')     "    - HWE Method:    ", trim(qc_struct%HWE_method)
      write(*,'(A,E12.4)') "    - HWE Threshold: ", qc_struct%HWE_threshold
      write(*,'(A,F10.4)') "    - MAF:           ", qc_struct%maf
      write(*,'(A,L2)')    "    - Sex Check:     ", qc_struct%sex_check
      write(*,'(A,F10.4)') "    - Parent Check:  ", qc_struct%mendel_viol_threshold
      write(*,'(A,L2)')    "    - Parent Seek:   ", qc_struct%parent_seek
      write(*,'(A)') ""
    end subroutine print_popqc_info

    subroutine print_ofile_info(name, info)
      character(len=*), intent(in) :: name
      type(FileInfo), intent(in) :: info
      integer :: i
      write(*,'(A)') "  "//trim(name)//":"
      write(*,'(A,A)') "    - Prefix:    ", trim(info%FileName)
      write(*,'(A,A)') "    - Delimiter: ", trim(info%Delim_char)
      write(*,'(A,I0)') "    - Variables: ", info%NVAR
      do i = 1, info%NVAR
        write(*,'(A,I2,A,A,A,I4)') "      [", i, "] ", trim(info%FieldName(i)), " width ", info%FieldLoc(i)
      end do
      write(*,'(A)') ""
    end subroutine print_ofile_info
end program test_readpar
