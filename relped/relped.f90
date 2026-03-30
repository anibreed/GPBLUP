program relped
 use renPED
 use readpar_mod, only: relped_params_t, read_relped_params
 implicit none
 character(len=MAX_STR):: PARAM_File
 type(relped_params_t) :: params
 logical :: use_gpu
 integer :: device_id
 PARAM_File=missC
 use_gpu = .false.
 device_id = 0

 call getarg(1,PARAM_File)
 PARAM_File=trim(adjustl(PARAM_File))
 if((len_trim(PARAM_File)).lt.1) stop 'Input parameter file'
 print*,"Finish to read parameter file=  ",trim(PARAM_File)

 call read_relped_params(PARAM_File, params)

 if(params%has_reference) then
     if(params%refinfo(1) /= 0) then
         call renumped(PEDFile=params%pedfile, PEDinfo=params%pedinfo, REFFile=params%reffile, REFinfo=params%refinfo, inb=params%inbreed, rel=params%relation, rel_inv=params%relation_inv, rel_outfile=params%relation_file, relinv_outfile=params%relation_inv_file, use_gpu=use_gpu, device_id=device_id, omp_threads_in=params%omp_threads)
     else
         call renumped(PEDFile=params%pedfile, PEDinfo=params%pedinfo, REFFile=params%reffile, inb=params%inbreed, rel=params%relation, rel_inv=params%relation_inv, rel_outfile=params%relation_file, relinv_outfile=params%relation_inv_file, use_gpu=use_gpu, device_id=device_id, omp_threads_in=params%omp_threads)
     endif
 else
     call renumped(PEDFile=params%pedfile, PEDinfo=params%pedinfo, inb=params%inbreed, rel=params%relation, rel_inv=params%relation_inv, rel_outfile=params%relation_file, relinv_outfile=params%relation_inv_file, use_gpu=use_gpu, device_id=device_id, omp_threads_in=params%omp_threads)
 endif
 end program
