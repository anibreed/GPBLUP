module renPED
  use M_ped
  use M_ReadFile, only: fopen, N_recf, readline
  use M_StrEdit, only: NumToStr, to_upper
  use M_renped_config, only: ren_outfile_global
  use M_Sort, only: Qsort
  use colleau_mod
  use inbreed, only: Inb_F
  use omp_lib
  implicit none
  interface
    subroutine triplets_bin_to_text(bin_file, txt_file)
      character(len=*), intent(in) :: bin_file, txt_file
    end subroutine triplets_bin_to_text
  end interface
  integer, parameter :: FMT_BINARY = 1
  integer, parameter :: FMT_TEXT = 2
  integer, parameter :: FMT_BOTH = 3
  type(ped_hash_t) :: h
  integer(kind=ki4) :: stat_sex_conflict = 0
  integer(kind=ki4) :: stat_sex_corrected = 0
  integer :: rel_fmt, relinv_fmt
  logical :: write_rel_bin, write_rel_txt, write_relinv_bin, write_relinv_txt
  character(len=MAX_STR) :: rel_txt_file, relinv_txt_file, rel_bin_file, relinv_bin_file
  contains

    function strip_dir(path) result(name)
      character(len=*), intent(in) :: path
      character(len=MAX_STR) :: name
      integer :: p, lt

      name = trim(path)
      lt = len_trim(name)
      if (lt <= 0) return
      p = index(name, '/', back=.true.)
      if (p > 0 .and. p < lt) then
        name = name(p+1:lt)
      end if
    end function strip_dir

    subroutine renumped(PEDFile, PEDinfo, REFFile, REFinfo, inb, rel, rel_inv, &
      rel_outfile, relinv_outfile, use_gpu, device_id, omp_threads_in)
  implicit none
  character(LEN=*),intent(in):: PEDFile
  integer,intent(in):: PEDinfo(:)
  character(LEN=*),intent(in),optional:: REFFile
  character(LEN=*),intent(in),optional:: rel_outfile, relinv_outfile
  integer,intent(in),optional:: REFinfo(:)
  integer,intent(in),optional:: inb, rel, rel_inv
  logical, intent(in), optional :: use_gpu
  integer, intent(in), optional :: device_id
  integer, intent(in), optional :: omp_threads_in
  type(ped_hash_t) :: h_keep
  integer,allocatable:: order(:), ref_nid(:)
  character(LEN=LEN_STR),allocatable:: REFAnim(:)
  character(LEN=LEN_STR):: SV(MAX_VAR), string
  character(LEN=MAX_STR):: ofile, rel_file, relinv_file
  logical:: inb_cal=.false., rel_cal=.false., relinv_cal=.false.
  logical :: need_F
  logical :: ref_ok
  integer(kind=ki4) :: NREC, NA, i
  integer(kind=ki4) :: k
  integer(kind=ki4) :: ref_idx
  integer :: omp_nprocs, omp_threads
  real(kind=r8) :: t_trace_start, t_trace_end
  real(kind=r8) :: t_cache_start, t_cache_end
  real(kind=r8) :: t_lookup_parent_total
  real(kind=r8) :: t_lookup_anc_start, t_lookup_anc_end, t_lookup_anc_total
  integer(kind=ki4) :: n_deleted
  integer(kind=ki4) :: n_cycle_break
  integer(kind=ki4) :: n_cycle_self_break, n_cycle_deep_break, n_cycle_scan_limit
  integer(kind=ki4) :: n_lookup_parent_calls, n_lookup_anc_calls
  integer(kind=ki4) :: ref_valid, ref_found, ref_missing, ref_marked
  integer(kind=ki4) :: ref_unique, ref_dup, ref_seen_slot
  integer(kind=ki4) :: ref_seed_unique, ref_seed_dup
  integer(kind=ki4) :: ref_desc_added
  integer(kind=ki4) :: ref_dup_sex_conflict, ref_dup_sex_count
  integer(kind=ki4), allocatable :: ref_dup_sex_first(:), ref_dup_sex_dup(:)
  character(len=LEN_STR), allocatable :: ref_dup_sex_id(:)
  integer(kind=ki4) :: ped_dual_role_count
  character(len=MAX_STR) :: ped_dual_role_report
  logical :: ref_seen_found
  logical :: do_deep_cycle_scan
  type(str_hash_t) :: ref_seen
  integer :: log_unt, ios_log, log_env_status
  character(len=MAX_STR) :: log_path

  log_unt = -1
  log_path = ''
  call get_environment_variable('RELPED_LOG', log_path, status=log_env_status)
  if (log_env_status == 0 .and. len_trim(log_path) > 0) then
    open(newunit=log_unt, file=trim(log_path), status='unknown', position='append', action='write', iostat=ios_log)
    if (ios_log /= 0) log_unt = -1
  end if
  call log_stage(log_unt, '[renumped][01] Start')
  call log_kv(log_unt, '[renumped] PEDFile=', trim(PEDFile))
  if (present(REFFile)) then
    call log_kv(log_unt, '[renumped] REFFile(arg)=', trim(REFFile))
  else
    call log_stage(log_unt, '[renumped] REFFile(arg)=<not provided>')
  end if

  ref_ok = .false.
  if(present(REFFile)) then
   if(.not.present(REFinfo) .or. (present(REFinfo) .and. REFinfo(1) == missI)) then
     write(*,*) 'Warning: REFFile specified but REFINFO missing or invalid; ignoring REFFile'
   else
     ref_ok = .true.
   endif
  endif
  if(present(inb)) then
  if(inb == 1) inb_cal=.true.
  endif
  if(present(rel)) then
  if(rel == 1) rel_cal=.true.
  endif
  if(present(rel_inv)) then
  if(rel_inv == 1) relinv_cal=.true.
  endif
  rel_file = 'A.mat'
  relinv_file = 'Ainv.mat'
  if (present(rel_outfile)) then
    if (len_trim(rel_outfile) > 0) rel_file = trim(rel_outfile)
  end if
  if (present(relinv_outfile)) then
    if (len_trim(relinv_outfile) > 0) relinv_file = trim(relinv_outfile)
  end if
  ! Always write relationship outputs to the current working directory.
  rel_file = strip_dir(rel_file)
  relinv_file = strip_dir(relinv_file)
  rel_fmt = FMT_BINARY
  relinv_fmt = FMT_BINARY
  rel_bin_file = trim(rel_file)
  relinv_bin_file = trim(relinv_file)
  rel_txt_file = trim(rel_file) // '.txt'
  relinv_txt_file = trim(relinv_file) // '.txt'

  if (len_trim(rel_file) >= 4 .and. rel_file(len_trim(rel_file)-3:) == '.txt') then
    rel_fmt = FMT_TEXT
    rel_txt_file = trim(rel_file)
    rel_bin_file = trim(rel_file(1:len_trim(rel_file)-4))
  else if (len_trim(rel_file) >= 5 .and. rel_file(len_trim(rel_file)-4:) == '.both') then
    rel_fmt = FMT_BOTH
    rel_bin_file = trim(rel_file(1:len_trim(rel_file)-5))
    rel_txt_file = trim(rel_bin_file) // '.txt'
  end if

  if (len_trim(relinv_file) >= 4 .and. relinv_file(len_trim(relinv_file)-3:) == '.txt') then
    relinv_fmt = FMT_TEXT
    relinv_txt_file = trim(relinv_file)
    relinv_bin_file = trim(relinv_file(1:len_trim(relinv_file)-4))
  else if (len_trim(relinv_file) >= 5 .and. relinv_file(len_trim(relinv_file)-4:) == '.both') then
    relinv_fmt = FMT_BOTH
    relinv_bin_file = trim(relinv_file(1:len_trim(relinv_file)-5))
    relinv_txt_file = trim(relinv_bin_file) // '.txt'
  end if

  write_rel_bin = (rel_fmt == FMT_BINARY .or. rel_fmt == FMT_BOTH)
  write_rel_txt = (rel_fmt == FMT_TEXT .or. rel_fmt == FMT_BOTH)
  write_relinv_bin = (relinv_fmt == FMT_BINARY .or. relinv_fmt == FMT_BOTH)
  write_relinv_txt = (relinv_fmt == FMT_TEXT .or. relinv_fmt == FMT_BOTH)
  ! For BOTH, reuse binary output for text conversion to avoid recomputation.
  need_F = inb_cal .or. rel_cal .or. relinv_cal
  omp_nprocs = omp_get_num_procs()
  omp_threads = min(32, omp_nprocs)
  if (present(device_id)) then
    ! Reserved for optional GPU backend selection.
  end if
  if (present(omp_threads_in)) then
    omp_threads = omp_threads_in
  end if
  if (omp_threads < 1) omp_threads = 1
  if (omp_threads > omp_nprocs) omp_threads = omp_nprocs
  call omp_set_num_threads(omp_threads)
  write(*,'(A,L1,2X,A,L1,2X,A,L1)') '[renumped][02] Options inb=', inb_cal, 'rel=', rel_cal, 'rel_inv=', relinv_cal
  write(*,'(A,L1,2X,A,I0)') '[renumped] REF enabled=', ref_ok, 'OMP threads=', omp_threads
  if (log_unt > 0) then
    write(log_unt,'(A,L1,2X,A,L1,2X,A,L1)') '[renumped][02] Options inb=', inb_cal, 'rel=', rel_cal, 'rel_inv=', relinv_cal
    write(log_unt,'(A,L1,2X,A,I0)') '[renumped] REF enabled=', ref_ok, 'OMP threads=', omp_threads
  end if
  call flush(6)

  t_lookup_parent_total = 0.0_r8
  t_lookup_anc_total = 0.0_r8
  n_lookup_parent_calls = 0
  n_lookup_anc_calls = 0
  ref_valid = 0
  ref_found = 0
  ref_missing = 0
  ref_marked = 0
  ref_unique = 0
  ref_dup = 0
  ref_seed_unique = 0
  ref_seed_dup = 0
  ref_desc_added = 0
  ref_dup_sex_conflict = 0
  ref_dup_sex_count = 0
  ped_dual_role_count = 0
  ped_dual_role_report = missC
  if (allocated(ref_dup_sex_first)) deallocate(ref_dup_sex_first)
  if (allocated(ref_dup_sex_dup)) deallocate(ref_dup_sex_dup)
  if (allocated(ref_dup_sex_id)) deallocate(ref_dup_sex_id)
  call ref_seen%reset()
  stat_sex_conflict = 0
  stat_sex_corrected = 0

  call ingest_ped_records(PEDFile, PEDinfo, NA, ped_dual_role_count, ped_dual_role_report, log_unt)

  call log_stage(log_unt, '[renumped][04] PED ingest complete')

  NREC=NA
  print*,"end read PEDFile data", NREC
  print*,"Ped File: Hash  size = ", h%hash_mask,"    Hash filled(No. of animal on PED) = ", h%n_keys_stored
    write(*,'(A,I0,2X,A,I0)') 'PED merge summary: sex_conflict=', stat_sex_conflict, &
      'sex_corrected=', stat_sex_corrected
    write(*,'(A,I0,2X,A)') 'PED dual-role IDs=', ped_dual_role_count, 'report output disabled'
  if(ref_ok) then
   call ingest_ref_records(REFFile, REFinfo, NA, NREC, REFAnim, ref_dup_sex_first, ref_dup_sex_dup, ref_dup_sex_id, &
     ref_valid, ref_found, ref_missing, ref_unique, ref_dup, ref_dup_sex_conflict, ref_dup_sex_count, log_unt, ref_seen)
 else
  call log_stage(log_unt, '[renumped][05] REF phase skipped')
   do i = 0, h%n_buckets-1
     if (h%valid_index(i)) then
        h%vals(i)%GEN= 0
     endif
   enddo
 endif

 print*,"After reading Reference File: Hash size = ", h%hash_mask, &
   "    PED hash filled(total animals on PED) = ", h%n_keys_stored, &
   "    REF accepted(unique valid keys) = ", NREC
 if(ref_ok) then
   do i = 0, h%n_buckets-1
     if (.not. h%valid_index(i)) cycle
     if (h%vals(i)%OBS > 0) ref_marked = ref_marked + 1
   enddo
     write(*,'(A,I0,2X,A,I0,2X,A,I0,2X,A,I0)') 'REF lookup summary: valid=', ref_valid, &
       'found=', ref_found, 'missing=', ref_missing, 'OBS-marked=', ref_marked
     write(*,'(A,I0,2X,A,I0,2X,A,I0)') 'REF key summary: unique=', ref_unique, &
       'duplicate=', ref_dup, 'valid-unique=', ref_valid - ref_unique
     write(*,'(A,I0)') 'REF duplicate sex-conflict count=', ref_dup_sex_conflict
     if (ref_dup_sex_count > 0) write(*,'(A)') 'REF duplicate sex-conflict report output disabled'
     if (ref_unique /= ref_marked) then
    write(*,'(A)') 'WARNING: REF unique keys and OBS-marked differ. Check key normalization/ID truncation.'
     end if
 endif
 call cpu_time(t_cache_start)
 call ped_build_index_cache(h)
 call cpu_time(t_cache_end)
 call log_stage(log_unt, '[renumped][06] ped_build_index_cache complete')
 write(*,'(A,F10.4,A)') 'ped_build_index_cache(before trace) elapsed (CPU sec): ', t_cache_end - t_cache_start, ' s'

 ! Fast cycle guard first: clear direct self-parent links in O(N).
 n_cycle_self_break = 0
 do i = 0, h%n_buckets-1
   if (.not. h%valid_index(i)) cycle
   if (h%vals(i)%ASD(2) /= missA .and. trim(h%vals(i)%ASD(2)) == trim(h%vals(i)%ASD(1))) then
     h%vals(i)%ASD(2) = missA
     n_cycle_self_break = n_cycle_self_break + 1
   end if
   if (h%vals(i)%ASD(3) /= missA .and. trim(h%vals(i)%ASD(3)) == trim(h%vals(i)%ASD(1))) then
     h%vals(i)%ASD(3) = missA
     n_cycle_self_break = n_cycle_self_break + 1
   end if
 end do
 if (n_cycle_self_break > 0) then
   write(*,'(A,I0)') 'Cycle guard: direct self-parent links cleared=', n_cycle_self_break
 end if

 ! Deep recursive cycle scan is very expensive on large pedigrees.
 ! Keep it for small/medium data and skip for large runs to avoid startup stalls.
 n_cycle_scan_limit = 100000
 do_deep_cycle_scan = (h%n_keys_stored <= n_cycle_scan_limit)
 n_cycle_deep_break = 0
 if (do_deep_cycle_scan) then
   write(*,'(A,I0,A)') 'Cycle guard: deep scan enabled (records=', h%n_keys_stored, ')'
   do i = 0, h%n_buckets-1
     if (.not. h%valid_index(i)) cycle
     if (check_IDD(i)) n_cycle_deep_break = n_cycle_deep_break + 1
   end do
 else
   write(*,'(A,I0,A,I0,A)') 'Cycle guard: deep scan skipped (records=', h%n_keys_stored, &
     ' > limit=', n_cycle_scan_limit, ')'
 end if

 n_cycle_break = n_cycle_self_break + n_cycle_deep_break
 if (n_cycle_break > 0) then
   write(*,'(A,I0,2X,A,I0,2X,A,I0)') 'Cycle guard summary: total=', n_cycle_break, &
     'self=', n_cycle_self_break, 'deep=', n_cycle_deep_break
   call cpu_time(t_cache_start)
   call ped_build_index_cache(h)
   call cpu_time(t_cache_end)
   write(*,'(A,F10.4,A)') 'ped_build_index_cache(after cycle guard) elapsed (CPU sec): ', t_cache_end - t_cache_start, ' s'
 end if
      
 write(*,'(A)') 'trace_PED: enter'
 call flush(6)
 call log_stage(log_unt, '[renumped][07] trace_PED start')
 call cpu_time(t_trace_start)
 call trace_PED()
 call cpu_time(t_trace_end)
 call log_stage(log_unt, '[renumped][08] trace_PED end')
 write(*,'(A,F10.4,A)') 'trace_PED elapsed (CPU sec): ', t_trace_end - t_trace_start, ' s'
 call flush(6)

 ! Explicitly include REF animals and all their ancestors.
 if (ref_ok .and. allocated(REFAnim)) then
   call cpu_time(t_lookup_anc_start)
   ! trace_PED already propagates REF seeds to all ancestors.
   ! Keep dedup/accounting here only.
   call ref_seen%reset()

   do k = 1, size(REFAnim)
     if (len_trim(REFAnim(k)) == 0) cycle
      call ref_seen%get(trim(REFAnim(k)), ref_seen_slot, ref_seen_found)
      if (ref_seen_found) then
        ref_seed_dup = ref_seed_dup + 1
        cycle
      end if
      call ref_seen%put(trim(REFAnim(k)), 1)
      ref_seed_unique = ref_seed_unique + 1

      ref_idx = ped_slot_from_key(REFAnim(k))
      n_lookup_anc_calls = n_lookup_anc_calls + 1
      if (ref_idx < 0 .or. ref_idx >= h%n_buckets) cycle
      if (.not. h%valid_index(ref_idx)) cycle
      if (h%vals(ref_idx)%OBS <= 0) h%vals(ref_idx)%OBS = 1
   end do

   call cpu_time(t_lookup_anc_end)
   t_lookup_anc_total = t_lookup_anc_total + (t_lookup_anc_end - t_lookup_anc_start)
   write(*,'(A,I0,2X,A,I0)') 'REF seed summary: unique=', ref_seed_unique, 'duplicate_skipped=', ref_seed_dup
   write(*,'(A,F10.4,A,I0)') 'lookup(ref+ancestor) elapsed (CPU sec): ', &
     t_lookup_anc_total, ' calls=', n_lookup_anc_calls

   ! REF selection rule: ancestor tracing only (no descendant expansion).
   ref_desc_added = 0
   write(*,'(A,I0)') 'REF direct-descendant added=', ref_desc_added
 end if

  n_deleted = 0
  if(ref_ok) then
   ! Keep REF seeds and all traced ancestors (GEN >= 0).
   ! OBS>0 keeps only direct REF rows and drops ancestors needed for F/A/Ainv.
   call h_keep%reset()
   do i = 0, h%n_buckets-1
      if (.not. h%valid_index(i)) cycle
      if (h%vals(i)%GEN >= 0) then
        call h_keep%ustore_value(h%vals(i)%ASD(1), h%vals(i))
      else
        n_deleted = n_deleted + 1
      endif
   enddo
   call h%reset()
   h = h_keep
   call h_keep%reset()
   print*,"After tracing PED: Hash  size = ", h%hash_mask,"    Hash filled(No. of animal on PED) = ", h%n_keys_stored
 !    call sleep(10)
 endif
  if (n_deleted > 0) then
    call cpu_time(t_cache_start)
    call ped_build_index_cache(h)
    call cpu_time(t_cache_end)
    write(*,'(A,I0)') 'Deleted records after trace: ', n_deleted
    write(*,'(A,F10.4,A)') 'ped_build_index_cache(after delete) elapsed (CPU sec): ', t_cache_end - t_cache_start, ' s'
  end if

  ! Use all surviving animals after trace/deletion.
  ! The OBS-only filter can under-count reference-connected pedigrees.
  NREC = h%n_keys_stored
  write(string,'(I0)') NREC
  call log_kv(log_unt, '[renumped][09] Survivors after trace=', trim(string))
 
 allocate(ref_nid(NREC),order(NREC))

 ref_nid(:)=missI
 order(:)=missI
 NA=0
 do i = 0, h%n_buckets-1
    if(h%valid_index(i)) then
      NA=NA+1
      ref_nid(NA)=h%vals(i)%BYR*100+(100 - h%vals(i)%GEN)
      order(NA)=i
    endif
 enddo
 if(NREC/=NA) stop "Valid_index error"
 NREC=NA
 call Qsort(ref_nid,order)
 call stabilize_order_ties_by_id(order, NREC)
 do i=1, NREC
  h%vals(order(i))%NASD(1)=i
 enddo
 call build_parent_index_map(order, NREC, t_lookup_parent_total, n_lookup_parent_calls, omp_threads, omp_nprocs)
 write(*,'(A,F10.4,A,I0)') 'lookup(parent map) elapsed (CPU sec): ', &
      t_lookup_parent_total, ' calls=', n_lookup_parent_calls

 if(need_F) then
  call run_inbreeding_phase(order, NREC, rel_cal, relinv_cal, rel_file, relinv_file, log_unt)
 endif

 call log_stage(log_unt, '[renumped][12] Writing renumbered output')
  if (len_trim(ren_outfile_global) > 0) then
    ofile = trim(ren_outfile_global)
  else if (ref_ok) then
    ofile = trim(REFFile)//"_ren"
  else
    ofile = trim(PEDFile)//"_ren"
  end if
  ofile = strip_dir(ofile)
call write_ren_output(ofile, NREC, order)
 call log_ren_tail(ofile, 10, log_unt)
 call log_kv(log_unt, '[renumped][13] Done output=', trim(ofile))
call cleanup_renumped_state(order, ref_nid, REFAnim, ref_dup_sex_first, ref_dup_sex_dup, ref_dup_sex_id)
 if (log_unt > 0) then
   close(unit=log_unt)
 end if
 call flush(6)
 end subroutine

 subroutine run_inbreeding_phase(order, nrec, rel_cal, relinv_cal, rel_file, relinv_file, log_unt)
 integer, intent(in) :: order(:), nrec
 logical, intent(in) :: rel_cal, relinv_cal
 character(len=*), intent(in) :: rel_file, relinv_file
 integer, intent(in) :: log_unt
 real(8), allocatable :: F(:)
 integer(kind=ki4), allocatable :: sire_arr(:), dam_arr(:)
 character(len=LEN_STR), allocatable :: orig_id(:)
 integer :: i, ii
 integer :: cnt_gt1e6, cnt_gt1e8, cnt_gt1e12
 real(kind=r8) :: minpos, maxF, sumF
 real(kind=r8) :: t_inb_start, t_inb_end

 call log_stage(log_unt, '[renumped][10] Inbreeding/relationship phase start')

 allocate(F(nrec))
 allocate(sire_arr(nrec))
 allocate(dam_arr(nrec))
 do i = 1, nrec
   if (h%vals(order(i))%NASD(2) == missI) then
     sire_arr(i) = 0
   else
     sire_arr(i) = h%vals(order(i))%NASD(2)
   end if
   if (h%vals(order(i))%NASD(3) == missI) then
     dam_arr(i) = 0
   else
     dam_arr(i) = h%vals(order(i))%NASD(3)
   end if
 end do

 write(*,*) 'Prepared sire/dam arrays for', nrec, 'records'
 write(*,*) '>>> renumped: calling Inb_F (Luo-Meuwissen)'
 call timestamp()
 call set_colleau_verbose(.false.)

 ii = 0
 do i = 1, nrec
   if (sire_arr(i) == 0 .or. dam_arr(i) == 0) ii = ii + 1
 end do
 write(*,'(A,I8)') 'Number of records with at least one missing parent:', ii

 call cpu_time(t_inb_start)
 call inbreed_build_npd(nrec, sire_arr, dam_arr, F)
 call cpu_time(t_inb_end)
 write(*,'(A,F10.4,A)') 'Inb_F elapsed (CPU sec): ', t_inb_end - t_inb_start, ' s'

 ii = 0
 do i = 1, nrec
   if (abs(F(i)) > 1.0d-12) ii = ii + 1
 end do
 write(*,'(A,I8)') 'Number of non-zero F values:', ii

 cnt_gt1e6 = 0; cnt_gt1e8 = 0; cnt_gt1e12 = 0
 minpos = 1.0e30_r8; maxF = 0.0_r8; sumF = 0.0_r8
 do i = 1, nrec
   if (F(i) > 1.0d-6) cnt_gt1e6 = cnt_gt1e6 + 1
   if (F(i) > 1.0d-8) cnt_gt1e8 = cnt_gt1e8 + 1
   if (F(i) > 1.0d-12) cnt_gt1e12 = cnt_gt1e12 + 1
   if (F(i) > 0.0_r8) then
     if (F(i) < minpos) minpos = F(i)
     if (F(i) > maxF) maxF = F(i)
     sumF = sumF + F(i)
   end if
 end do
 if (minpos > 1.0e29_r8) minpos = 0.0_r8
 write(*,'(A,I8)') 'Count F>1e-6:', cnt_gt1e6
 write(*,'(A,I8)') 'Count F>1e-8:', cnt_gt1e8
 write(*,'(A,I8)') 'Count F>1e-12:', cnt_gt1e12
 write(*,'(A,1X,ES16.8)') 'min positive F=', minpos
 write(*,'(A,1X,ES16.8)') 'max F=', maxF
 write(*,'(A,1X,ES16.8)') 'sum F=', sumF

 if (rel_cal) then
   write(*,*) '>>> renumped: writing triplets ', trim(rel_file)
   if (write_rel_bin) then
     write(*,'(A)') 'write_triplets binary layout: [char4 magic="ATRP"][int32 version]' // &
                    '[int32 n][int64 nnz][nnz * (int32 row,int32 col,float64 value)]'
     if (log_unt > 0) then
       write(log_unt,'(A)') 'write_triplets binary layout: [char4 magic="ATRP"][int32 version]' // &
                            '[int32 n][int64 nnz][nnz * (int32 row,int32 col,float64 value)]'
     end if
   end if
   if (write_rel_txt) then
     write(*,'(A)') 'write_triplets text layout: row col value row_id col_id'
     if (log_unt > 0) write(log_unt,'(A)') 'write_triplets text layout: row col value row_id col_id'
   end if
   if ((write_rel_txt .and. .not. write_rel_bin) .or. (write_relinv_txt .and. .not. write_relinv_bin)) then
     allocate(orig_id(nrec))
     do i = 1, nrec
       orig_id(i) = h%vals(order(i))%ASD(1)
     end do
   end if

   call cpu_time(t_inb_start)
    if (write_rel_bin) call write_triplets(nrec, sire_arr, dam_arr, F, trim(rel_bin_file))
    if (write_rel_txt) then
      if (write_rel_bin) then
        call triplets_bin_to_text_with_ids(trim(rel_bin_file), trim(rel_txt_file), order, nrec)
      else
        call write_triplets_text_mod(nrec, sire_arr, dam_arr, F, orig_id, trim(rel_txt_file))
      end if
    end if
   call cpu_time(t_inb_end)
   write(*,'(A,F10.4,A)') 'write_triplets elapsed (CPU sec): ', t_inb_end - t_inb_start, ' s'
  if (write_rel_bin) call log_triplet_binary_block(trim(rel_bin_file), 'ATRP', 10, log_unt, 'A_triplets')
   write(*,*) '<<< renumped: triplets written'
 endif
 if (relinv_cal) then
   write(*,*) '>>> renumped: writing inverse triplets ', trim(relinv_file)
   if (write_relinv_bin) then
     write(*,'(A)') 'write_ainv_triplets binary layout: [char4 magic="AINV"][int32 version]' // &
                    '[int32 n][int64 nnz][nnz * (int32 row,int32 col,float64 value)]'
     if (log_unt > 0) then
       write(log_unt,'(A)') 'write_ainv_triplets binary layout: [char4 magic="AINV"][int32 version]' // &
                            '[int32 n][int64 nnz][nnz * (int32 row,int32 col,float64 value)]'
     end if
   end if
   if (write_relinv_txt) then
     write(*,'(A)') 'write_ainv_triplets text layout: row col value row_id col_id'
     if (log_unt > 0) write(log_unt,'(A)') 'write_ainv_triplets text layout: row col value row_id col_id'
   end if
   call cpu_time(t_inb_start)
    if (write_relinv_bin) call write_ainv_triplets(nrec, sire_arr, dam_arr, F, trim(relinv_bin_file))
    if (write_relinv_txt) then
      if (write_relinv_bin) then
        call triplets_bin_to_text_with_ids(trim(relinv_bin_file), trim(relinv_txt_file), order, nrec)
      else
        call write_ainv_triplets_text_mod(nrec, sire_arr, dam_arr, F, orig_id, trim(relinv_txt_file))
      end if
    end if
   call cpu_time(t_inb_end)
   write(*,'(A,F10.4,A)') 'write_ainv_triplets elapsed (CPU sec): ', t_inb_end - t_inb_start, ' s'
  if (write_relinv_bin) call log_triplet_binary_block(trim(relinv_bin_file), 'AINV', 10, log_unt, 'Ainv_triplets')
   write(*,*) '<<< renumped: inverse triplets written'
 endif

 write(*,*) '<<< renumped: returned from compute_inbreeding'
 call timestamp()

 if (allocated(orig_id)) deallocate(orig_id)
 do i = 1, nrec
   h%vals(order(i))%F = F(i)
 end do

 deallocate(F, sire_arr, dam_arr)
 call log_stage(log_unt, '[renumped][11] Inbreeding/relationship phase end')
 end subroutine

 subroutine build_parent_index_map(order, nrec, t_lookup_parent_total, n_lookup_parent_calls, omp_threads, omp_nprocs)
 integer, intent(in) :: order(:), nrec
 real(kind=r8), intent(inout) :: t_lookup_parent_total
 integer(kind=ki4), intent(inout) :: n_lookup_parent_calls
 integer, intent(in) :: omp_threads, omp_nprocs
 integer :: i, ix
 real(kind=r8) :: t_lookup_parent_start, t_lookup_parent_end

 write(*,'(A,I4,2X,A,I4)') 'OpenMP: using threads=', omp_threads, ' (procs=', omp_nprocs,')'
 t_lookup_parent_start = omp_get_wtime()
 !$omp parallel do default(shared) private(i,ix) reduction(+:n_lookup_parent_calls) schedule(static)
 do i=1,nrec
    ix = ped_slot_from_key(h%vals(order(i))%ASD(2))
    if(ix /= -1) then
      n_lookup_parent_calls = n_lookup_parent_calls + 1
      ! Always map sire index if found (avoid leaving parent index as missing)
      h%vals(order(i))%NASD(2) = h%vals(ix)%NASD(1)
    endif
    ix = ped_slot_from_key(h%vals(order(i))%ASD(3))
    if(ix /= -1) then
      n_lookup_parent_calls = n_lookup_parent_calls + 1
      ! Always map dam index if found
      h%vals(order(i))%NASD(3) = h%vals(ix)%NASD(1)
    endif
 enddo
 !$omp end parallel do
 t_lookup_parent_end = omp_get_wtime()
 t_lookup_parent_total = t_lookup_parent_total + (t_lookup_parent_end - t_lookup_parent_start)
 end subroutine

 subroutine ingest_ped_records(PEDFile, PEDinfo, NA, ped_dual_role_count, ped_dual_role_report, log_unt)
 character(LEN=*),intent(in):: PEDFile
 integer,intent(in):: PEDinfo(:)
 integer(kind=ki4), intent(out) :: NA
 integer(kind=ki4), intent(out) :: ped_dual_role_count
 character(len=*), intent(out) :: ped_dual_role_report
 integer, intent(in) :: log_unt
 type(CPED):: NX(3), OX(3)
 character(LEN=LEN_STR):: SV(MAX_VAR)
 integer:: IV(MAX_VAR)
 integer(kind=ki4) :: unt, n, i, stats(3)

 call log_stage(log_unt, '[renumped][03] Reading PED records')
 unt = fopen(PEDFile,MAX_RECL)
 NA=0

 do
   SV=missC
   IV=missI
   call initX(OX)
   call readline(unt,n,SV,IV)
   if(n <= 0) exit
   NA=NA+1
   call checkPED(SV(PEDinfo(1:3)),NX)
   if(PEDinfo(4)/=0) NX(1)%BYR=IV(PEDinfo(4))
   if(PEDinfo(5)/=0) NX(1)%SEX=IV(PEDinfo(5))
   do i=1,3
    if( .not. Valid(NX(i)%ASD(1))) cycle
    call h%get_value(NX(i)%ASD(1), OX(i), stats(i))
   enddo
   call update(OX, NX, stats)
   do i=1,3
    if( .not. Valid(NX(i)%ASD(1))) cycle
    call h%ustore_value(NX(i)%ASD(1), NX(i))
   enddo
 enddo
 close(unt)

 ! Normalize sex using pedigree parent roles: sire=>male, dam=>female.
 call enforce_parent_role_sex()
 ped_dual_role_count = 0
 ped_dual_role_report = 'disabled'
 end subroutine

 subroutine ingest_ref_records(REFFile, REFinfo, NA, NREC, REFAnim, ref_dup_sex_first, ref_dup_sex_dup, ref_dup_sex_id, &
   ref_valid, ref_found, ref_missing, ref_unique, ref_dup, ref_dup_sex_conflict, ref_dup_sex_count, log_unt, ref_seen)
 character(LEN=*), intent(in) :: REFFile
 integer, intent(in) :: REFinfo(:)
 integer(kind=ki4), intent(out) :: NA, NREC
 character(LEN=LEN_STR), allocatable, intent(inout) :: REFAnim(:)
 integer(kind=ki4), allocatable, intent(inout) :: ref_dup_sex_first(:), ref_dup_sex_dup(:)
 character(len=LEN_STR), allocatable, intent(inout) :: ref_dup_sex_id(:)
 integer(kind=ki4), intent(inout) :: ref_valid, ref_found, ref_missing, ref_unique, ref_dup
 integer(kind=ki4), intent(inout) :: ref_dup_sex_conflict, ref_dup_sex_count
 integer, intent(in) :: log_unt
 type(str_hash_t), intent(inout) :: ref_seen
 type(CPED) :: NX(3)
 character(LEN=LEN_STR):: SV(MAX_VAR)
 integer :: IV(MAX_VAR)
  integer(kind=ki4) :: unt, n, status, ref_curr_sex, ref_prev_sex, ref_seen_slot_local
  logical :: header_checked
  character(len=LEN_STR) :: ref_id_token, ref_sex_token
 logical :: ref_seen_found_local

 call log_stage(log_unt, '[renumped][05] Reading REF records')
 NREC= N_recf(REFFile)
 allocate(REFAnim(NREC))
 allocate(ref_dup_sex_first(NREC), ref_dup_sex_dup(NREC), ref_dup_sex_id(NREC))
 ref_dup_sex_first = missI
 ref_dup_sex_dup = missI
 ref_dup_sex_id = missC
 REFAnim=missC

 unt = fopen(REFFile,MAX_RECL)
 NA=0
 header_checked = .false.
 do
   SV=missC
   IV=missI
   call initX(NX)
   call readline(unt,n,SV,IV)
   if(n <= 0) exit
   if (.not. header_checked) then
     header_checked = .true.
     ref_id_token = to_upper(trim(SV(REFinfo(1))))
     ref_sex_token = ''
     if (size(REFinfo) >= 2 .and. REFinfo(2) /= 0) then
       ref_sex_token = to_upper(trim(SV(REFinfo(2))))
     end if
     if (ref_id_token == 'ANIMAL_ID' .or. ref_id_token == 'ID' .or. ref_id_token == 'REF_ID' .or. &
         ref_id_token == 'SAMPLE_ID' .or. ref_id_token == 'INDIVIDUAL_ID' .or. ref_id_token == 'INDIVIDUAL') then
       cycle
     end if
     if (len_trim(ref_sex_token) > 0) then
       if (ref_sex_token == 'SEX' .or. ref_sex_token == 'GENDER') cycle
     end if
   end if
   if(.not.Valid(SV(REFinfo(1)))) cycle
   ref_valid = ref_valid + 1
   ref_curr_sex = missI
   if (size(REFinfo) >= 2 .and. REFinfo(2) /= 0) ref_curr_sex = IV(REFinfo(2))

   call ref_seen%get(trim(SV(REFinfo(1))), ref_seen_slot_local, ref_seen_found_local)
   if (ref_seen_found_local) then
     ref_dup = ref_dup + 1
     ref_prev_sex = ref_seen_slot_local
     if (ref_prev_sex /= missI .and. ref_curr_sex /= missI .and. ref_prev_sex /= ref_curr_sex) then
       ref_dup_sex_conflict = ref_dup_sex_conflict + 1
       ref_dup_sex_count = ref_dup_sex_count + 1
       ref_dup_sex_id(ref_dup_sex_count) = trim(SV(REFinfo(1)))
       ref_dup_sex_first(ref_dup_sex_count) = ref_prev_sex
       ref_dup_sex_dup(ref_dup_sex_count) = ref_curr_sex
     end if
     cycle
   else
     ref_unique = ref_unique + 1
     call ref_seen%put(trim(SV(REFinfo(1))), ref_curr_sex)
   end if

   NA=NA+1

   call h%get_value(SV(REFinfo(1)), NX(1), status)
   if(status == -1) then
     ref_missing = ref_missing + 1
     print*,SV(REFinfo(1)),"not included on pedigree file"
     NX(1)%ASD(1)=SV(REFinfo(1))
     if(size(REFinfo) >= 2 .and. REFinfo(2) /= 0) NX(1)%SEX=IV(REFinfo(2))
   else
     ref_found = ref_found + 1
   endif
   NX(1)%GEN = 0
   NX(1)%OBS=NX(1)%OBS+1
   call h%ustore_value(NX(1)%ASD(1),NX(1))
   REFAnim(NA)=SV(REFinfo(1))
 enddo
 close(unt)
 NREC=NA
 print*,"end read REFFile data", NREC
 end subroutine

 subroutine write_ren_output(ofile, nrec, order)
 character(len=*), intent(in) :: ofile
 integer(kind=ki4), intent(in) :: nrec
 integer, intent(in) :: order(:)
 integer :: unt, i, strlen
 character(len=LEN_STR) :: string
 character(len=MAX_STR) :: ofmt

 unt=fopen(ofile,MAX_RECL)
 strlen = 0
 do i = 1, nrec
   if(len_trim(h%vals(order(i))%ASD(1)) > strlen) strlen=len_trim(h%vals(order(i))%ASD(1))
 end do

 string = NumToStr(strlen, '(i10)')
 ofmt="),i2,i9,i3,i5,f10.6,3(1x,a"//trim(string)//"))"
 strlen=floor(log10(real(nrec)) + 1)
 string = NumToStr(strlen, '(i10)')
 ofmt="(3(1x,i"//trim(string)//ofmt

 do i=1,nrec
   write(unt,fmt=ofmt) h%vals(order(i))%NASD(:), &
           h%vals(order(i))%SEX, &
           h%vals(order(i))%BYR, &
           h%vals(order(i))%GEN, &
           h%vals(order(i))%OBS, &
           h%vals(order(i))%F, &
           h%vals(order(i))%ASD(:)
 enddo
 close(unit=unt)
 end subroutine

 subroutine cleanup_renumped_state(order, ref_nid, REFAnim, ref_dup_sex_first, ref_dup_sex_dup, ref_dup_sex_id)
 integer, allocatable, intent(inout) :: order(:), ref_nid(:)
 character(len=LEN_STR), allocatable, intent(inout) :: REFAnim(:), ref_dup_sex_id(:)
 integer(kind=ki4), allocatable, intent(inout) :: ref_dup_sex_first(:), ref_dup_sex_dup(:)

 if (allocated(order)) deallocate(order)
 if (allocated(ref_nid)) deallocate(ref_nid)
 if (allocated(REFAnim)) deallocate(REFAnim)
 if (allocated(ref_dup_sex_first)) deallocate(ref_dup_sex_first)
 if (allocated(ref_dup_sex_dup)) deallocate(ref_dup_sex_dup)
 if (allocated(ref_dup_sex_id)) deallocate(ref_dup_sex_id)

 call h%reset()
 call ped_clear_index_cache()
 end subroutine

 subroutine log_stage(log_unt, msg)
 integer, intent(in) :: log_unt
 character(len=*), intent(in) :: msg
 write(*,'(A)') trim(msg)
 if (log_unt > 0) write(log_unt,'(A)') trim(msg)
 call flush(6)
 end subroutine

 subroutine log_kv(log_unt, label, value)
 integer, intent(in) :: log_unt
 character(len=*), intent(in) :: label, value
 write(*,'(A,1X,A)') trim(label), trim(value)
 if (log_unt > 0) write(log_unt,'(A,1X,A)') trim(label), trim(value)
 call flush(6)
 end subroutine

 subroutine log_ren_tail(ren_file, tail_n, log_unt)
 character(len=*), intent(in) :: ren_file
 integer, intent(in) :: tail_n, log_unt
 character(len=MAX_RECL), allocatable :: ring(:)
 character(len=MAX_RECL) :: line
 integer :: u, ios, cnt, idx, start_i, n_show, i

 if (tail_n <= 0) return
 allocate(ring(tail_n))
 ring = ''
 cnt = 0

 open(newunit=u, file=trim(ren_file), status='old', action='read', iostat=ios)
 if (ios /= 0) then
   write(*,'(A,1X,A)') 'ren tail log skipped (open failed):', trim(ren_file)
   if (log_unt > 0) write(log_unt,'(A,1X,A)') 'ren tail log skipped (open failed):', trim(ren_file)
   deallocate(ring)
   return
 end if

 do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   cnt = cnt + 1
   idx = mod(cnt - 1, tail_n) + 1
   ring(idx) = trim(line)
 end do
 close(u)

 n_show = min(cnt, tail_n)
 write(*,'(A,I0,A)') 'ren tail(', n_show, ') from file: '//trim(ren_file)
 if (log_unt > 0) write(log_unt,'(A,I0,A)') 'ren tail(', n_show, ') from file: '//trim(ren_file)
 if (n_show <= 0) then
   deallocate(ring)
   return
 end if

 start_i = mod(cnt - n_show, tail_n) + 1
 do i = 1, n_show
   idx = mod(start_i + i - 2, tail_n) + 1
   write(*,'(A)') trim(ring(idx))
   if (log_unt > 0) write(log_unt,'(A)') trim(ring(idx))
 end do
 deallocate(ring)
 end subroutine

 subroutine log_triplet_binary_block(file_path, expected_magic, block_n, log_unt, title)
 character(len=*), intent(in) :: file_path, expected_magic, title
 integer, intent(in) :: block_n, log_unt
 integer :: u, ios, i, j
 integer(kind=ki4) :: ver, nfile, row, col, bs, start_idx
 integer(kind=ki8) :: nnz, k
 real(kind=r8) :: val
 character(len=4) :: magic
 real(kind=r8), allocatable :: blk(:,:)

 open(newunit=u, file=trim(file_path), status='old', action='read', form='unformatted', access='stream', iostat=ios)
 if (ios /= 0) then
   write(*,'(A,1X,A)') trim(title)//' block log skipped (open failed):', trim(file_path)
   if (log_unt > 0) write(log_unt,'(A,1X,A)') trim(title)//' block log skipped (open failed):', trim(file_path)
   return
 end if

 read(u, iostat=ios) magic, ver, nfile, nnz
 if (ios /= 0) then
   write(*,'(A,1X,A)') trim(title)//' block log skipped (invalid header):', trim(file_path)
   if (log_unt > 0) write(log_unt,'(A,1X,A)') trim(title)//' block log skipped (invalid header):', trim(file_path)
   close(u)
   return
 end if
 if (trim(magic) /= trim(expected_magic)) then
     write(*,'(A,1X,A,1X,A,1X,A)') trim(title)//' block log skipped (magic mismatch):', 'expected', &
       trim(expected_magic), 'got '//trim(magic)
     if (log_unt > 0) write(log_unt,'(A,1X,A,1X,A,1X,A)') trim(title)//' block log skipped (magic mismatch):', &
       'expected', trim(expected_magic), 'got '//trim(magic)
   close(u)
   return
 end if

 bs = min(block_n, nfile)
 if (bs <= 0) then
   close(u)
   return
 end if
 start_idx = nfile - bs + 1
 allocate(blk(bs, bs))
 blk = 0.0_r8

 do k = 1, nnz
   read(u, iostat=ios) row, col, val
   if (ios /= 0) exit
   if (row >= start_idx .and. row <= nfile .and. col >= start_idx .and. col <= nfile) then
     blk(row-start_idx+1, col-start_idx+1) = val
     blk(col-start_idx+1, row-start_idx+1) = val
   end if
 end do
 close(u)

 write(*,'(A,1X,A,1X,A,I0,A,I0,A)') trim(title)//' lower block from', trim(file_path), &
   '[', start_idx, ':', nfile, ']'
 if (log_unt > 0) write(log_unt,'(A,1X,A,1X,A,I0,A,I0,A)') trim(title)//' lower block from', trim(file_path), &
   '[', start_idx, ':', nfile, ']'
 do i = 1, bs
   write(*,'(10(1X,F10.6))') (blk(i,j), j=1,bs)
   if (log_unt > 0) write(log_unt,'(10(1X,F10.6))') (blk(i,j), j=1,bs)
 end do

 deallocate(blk)
 end subroutine

 subroutine checkPED(SV,X)
 character(len=*),intent(in):: SV(:)
 type(CPED),intent(out):: X(:)
 character(len=len(SV(1))):: TV(size(SV))
 integer:: i
 call initX(X)
  if(.not.Valid(SV(1))) return
  TV=SV
  ! check same ID between Sire and dam
  if(SV(2)==SV(3)) then
     TV(2)=missA
     TV(3)=missA
  endif
  ! check same ID between animal and parents
  do i = 2, 3
       if(SV(1) == SV(i))  TV(i)=missA
  enddo
  do i =1, 3
      X(1)%ASD(i)=TV(i)
  enddo
  X(2)%ASD(1)=TV(2)
  X(2)%SEX=Male
  X(3)%ASD(1)=TV(3)
  X(3)%SEX=Female
 end subroutine

    subroutine enforce_parent_role_sex()
    integer(kind=ki4) :: i, idx
    integer(kind=ki4), allocatable :: role_mask(:)

    allocate(role_mask(0:h%n_buckets-1))
    role_mask = 0

    do i = 0, h%n_buckets-1
      if (.not. h%valid_index(i)) cycle

      if (h%vals(i)%ASD(2) /= missA) then
        idx = ped_slot_from_key(h%vals(i)%ASD(2))
        if (idx /= -1) role_mask(idx) = ior(role_mask(idx), 1)
      end if

      if (h%vals(i)%ASD(3) /= missA) then
        idx = ped_slot_from_key(h%vals(i)%ASD(3))
        if (idx /= -1) role_mask(idx) = ior(role_mask(idx), 2)
      end if
    end do

    do i = 0, h%n_buckets-1
      if (.not. h%valid_index(i)) cycle

      if (role_mask(i) == 1) then
        if (h%vals(i)%SEX /= Male) then
          h%vals(i)%SEX = Male
          stat_sex_corrected = stat_sex_corrected + 1
        end if
      else if (role_mask(i) == 2) then
        if (h%vals(i)%SEX /= Female) then
          h%vals(i)%SEX = Female
          stat_sex_corrected = stat_sex_corrected + 1
        end if
      end if
    end do

    deallocate(role_mask)
    end subroutine enforce_parent_role_sex

    subroutine write_ref_duplicate_sex_report(reffile, n_conflict, ids, first_sex, dup_sex, report_file)
    character(len=*), intent(in) :: reffile
    integer(kind=ki4), intent(in) :: n_conflict
    character(len=*), intent(in) :: ids(:)
    integer(kind=ki4), intent(in) :: first_sex(:), dup_sex(:)
    character(len=*), intent(out) :: report_file
    integer :: unt, ios
    integer(kind=ki4) :: i
    character(len=MAX_STR) :: base

    base = strip_dir(reffile)
    report_file = trim(base)//'.duplicate_sex_conflicts.txt'
    open(newunit=unt, file=trim(report_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      report_file = 'report_open_failed'
      return
    end if

    write(unt,'(A)') '# Duplicate REF IDs with conflicting sex codes'
    write(unt,'(A)') '# id first_seen_sex duplicate_seen_sex'
    do i = 1, n_conflict
      write(unt,'(A,1X,I0,1X,I0)') trim(ids(i)), first_sex(i), dup_sex(i)
    end do

    close(unt)
    end subroutine write_ref_duplicate_sex_report

  subroutine update(A,B, stat)  ! modify A to B
  type(CPED),intent(inout):: A(:), B(:)
  integer,intent(in):: stat(:)
  integer:: ii, jj
  do  ii = 1, 3
     if(stat(ii) == -1) cycle 
    do jj=2, 3
       if(A(ii)%ASD(jj) /= missA .and. B(ii)%ASD(jj) == missA) then
         B(ii)%ASD(jj) = A(ii)%ASD(jj)
       elseif(A(ii)%ASD(jj) /= missA .and. B(ii)%ASD(jj) /= missA) then
         if(A(ii)%ASD(jj) /= B(ii)%ASD(jj)) B(ii)%ASD(jj)= missA
       endif
    enddo
    if(A(ii)%BYR /= missI .and. B(ii)%BYR == missI) B(ii)%BYR=A(ii)%BYR
    if(A(ii)%SEX /= missI) then
      if(B(ii)%SEX == missI) then
        B(ii)%SEX=A(ii)%SEX
      elseif(B(ii)%SEX /= missI) then
        if(A(ii)%SEX /= B(ii)%SEX ) then
          stat_sex_conflict = stat_sex_conflict + 1
          ! Keep record and parent links; correct sex using PED role.
          if (ii == 2) then
            B(ii)%SEX = Male
          else if (ii == 3) then
            B(ii)%SEX = Female
          end if
       endif
      endif
    endif
  enddo
  return
  end subroutine

    subroutine trace_PED()
    integer:: i, k, nvalid, genA, updates, parent_slot, parent_k
    integer:: q_head, q_tail, q_count, q_cap, qk
    integer, allocatable :: valid_slots(:), sire_slot(:), dam_slot(:), slot_to_k(:), queue_k(:)
    logical, allocatable :: in_queue(:)
    integer :: trace_processed, trace_log_step, trace_next_log
    integer :: max_gen_cap, new_gen
    real(kind=r8) :: t_lookup_start, t_lookup_end, t_lookup_total
    real(kind=r8) :: t_trace_loop_start, t_trace_loop_now
    integer :: n_lookup_calls

    t_lookup_total = 0.0_r8
    n_lookup_calls = 0

    nvalid = 0
    do i = 0, h%n_buckets-1
      if (h%valid_index(i)) nvalid = nvalid + 1
    enddo
    if (nvalid <= 0) return

    trace_log_step = max(50000, nvalid / 20)
    trace_processed = 0
    trace_next_log = trace_log_step
    max_gen_cap = nvalid
    call cpu_time(t_trace_loop_start)
    write(*,'(A,I0,A,I0,A,I0)') 'trace_PED: start nvalid=', nvalid, ' log_step=', trace_log_step, ' gen_cap=', max_gen_cap
    call flush(6)

    allocate(valid_slots(nvalid), sire_slot(nvalid), dam_slot(nvalid))
    allocate(slot_to_k(0:h%n_buckets-1))
    allocate(queue_k(nvalid), in_queue(nvalid))

    slot_to_k = 0
    in_queue = .false.

    k = 0
    do i = 0, h%n_buckets-1
      if (.not. h%valid_index(i)) cycle
      k = k + 1
      valid_slots(k) = i
      slot_to_k(i) = k
    enddo

    do k = 1, nvalid
      i = valid_slots(k)
      sire_slot(k) = -1
      dam_slot(k) = -1
      call cpu_time(t_lookup_start)
      if (h%vals(i)%ASD(2) /= missA) then
        sire_slot(k) = ped_slot_from_key(h%vals(i)%ASD(2))
        n_lookup_calls = n_lookup_calls + 1
      end if
      if (h%vals(i)%ASD(3) /= missA) then
        dam_slot(k) = ped_slot_from_key(h%vals(i)%ASD(3))
        n_lookup_calls = n_lookup_calls + 1
      end if
      call cpu_time(t_lookup_end)
      t_lookup_total = t_lookup_total + (t_lookup_end - t_lookup_start)
    enddo

    ! Skip recursive cycle-scan (very expensive on large pedigrees).
    ! Bounded generation propagation below guarantees termination.
    write(*,'(A)') 'trace_PED: cycle-scan skipped, using bounded propagation'
    call flush(6)

    q_cap = nvalid
    q_head = 1
    q_tail = 1
    q_count = 0
    do k = 1, nvalid
      i = valid_slots(k)
      if (h%vals(i)%GEN >= 0) then
        queue_k(q_tail) = k
        q_tail = q_tail + 1
        if (q_tail > q_cap) q_tail = 1
        q_count = q_count + 1
        in_queue(k) = .true.
      endif
    enddo
    write(*,'(A,I0)') 'trace_PED: queue initialized, count=', q_count
    call flush(6)

    updates = 0
    do while (q_count > 0)
      qk = queue_k(q_head)
      q_head = q_head + 1
      if (q_head > q_cap) q_head = 1
      q_count = q_count - 1
      in_queue(qk) = .false.
      trace_processed = trace_processed + 1

      if (trace_processed >= trace_next_log) then
        call cpu_time(t_trace_loop_now)
        write(*,'(A,I0,A,I0,A,I0,A,F10.4,A)') 'trace_PED: processed=', trace_processed, &
          '/', nvalid, ' queue=', q_count, ' elapsed=', t_trace_loop_now - t_trace_loop_start, ' s'
        call flush(6)
        trace_next_log = trace_next_log + trace_log_step
      end if

      i = valid_slots(qk)
      genA = h%vals(i)%GEN

      parent_slot = sire_slot(qk)
      if (parent_slot /= -1) then
        new_gen = min(max_gen_cap, genA + 1)
        if (h%vals(parent_slot)%GEN < new_gen) then
          h%vals(parent_slot)%GEN = new_gen
          updates = updates + 1
          parent_k = slot_to_k(parent_slot)
          if (parent_k > 0 .and. .not. in_queue(parent_k)) then
            queue_k(q_tail) = parent_k
            q_tail = q_tail + 1
            if (q_tail > q_cap) q_tail = 1
            q_count = q_count + 1
            in_queue(parent_k) = .true.
          endif
        endif
      endif

      parent_slot = dam_slot(qk)
      if (parent_slot /= -1) then
        new_gen = min(max_gen_cap, genA + 1)
        if (h%vals(parent_slot)%GEN < new_gen) then
          h%vals(parent_slot)%GEN = new_gen
          updates = updates + 1
          parent_k = slot_to_k(parent_slot)
          if (parent_k > 0 .and. .not. in_queue(parent_k)) then
            queue_k(q_tail) = parent_k
            q_tail = q_tail + 1
            if (q_tail > q_cap) q_tail = 1
            q_count = q_count + 1
            in_queue(parent_k) = .true.
          endif
        endif
      endif
    enddo

    call cpu_time(t_trace_loop_now)
    write(*,'(A,I0,A,F10.4,A)') 'trace_PED: done processed=', trace_processed, &
      ' elapsed=', t_trace_loop_now - t_trace_loop_start, ' s'
    call flush(6)
    print*, 'trace_PED: generation updates =', updates
    write(*,'(A,F10.4,A,I0)') 'trace_PED lookup elapsed (CPU sec): ', &
      t_lookup_total, ' calls=', n_lookup_calls

    deallocate(valid_slots, sire_slot, dam_slot, slot_to_k, queue_k, in_queue)
  end subroutine

   subroutine stabilize_order_ties_by_id(order, nrec)
   integer, intent(inout) :: order(:)
   integer, intent(in) :: nrec
   integer :: i, j, k, key_nid, curr_nid, tmp

   if (nrec <= 1) return

   i = 1
   do while (i <= nrec)
     key_nid = h%vals(order(i))%BYR*100 + (100 - h%vals(order(i))%GEN)
     j = i + 1
     do while (j <= nrec)
       curr_nid = h%vals(order(j))%BYR*100 + (100 - h%vals(order(j))%GEN)
       if (curr_nid /= key_nid) exit
       j = j + 1
     end do

     if (j - i > 1) then
       do k = i + 1, j - 1
         tmp = order(k)
         curr_nid = k - 1
         do while (curr_nid >= i)
           if (trim(h%vals(order(curr_nid))%ASD(1)) <= trim(h%vals(tmp)%ASD(1))) exit
           order(curr_nid + 1) = order(curr_nid)
           curr_nid = curr_nid - 1
         end do
         order(curr_nid + 1) = tmp
       end do
     end if

     i = j
   end do
   end subroutine stabilize_order_ties_by_id

 logical function check_IDD(step) result(check)
 integer,intent(in):: step

 check = .false.
 if (.not. h%valid_index(step)) return

 if (h%vals(step)%ASD(2) /= missA) then
   if (check_ID(h%vals(step)%ASD(2), h%vals(step)%ASD(1), 1)) then
     h%vals(step)%ASD(2) = missA
     check = .true.
   end if
 end if

 if (h%vals(step)%ASD(3) /= missA) then
   if (check_ID(h%vals(step)%ASD(3), h%vals(step)%ASD(1), 1)) then
     h%vals(step)%ASD(3) = missA
     check = .true.
   end if
 end if
 end function check_IDD

 recursive logical function check_ID(curr_key, root_key, depth) result(tangled)
 character(len=*),intent(in) :: curr_key, root_key
 integer,intent(in) :: depth
 integer :: slot

 tangled = .false.
 if (curr_key == missA) return
 if (trim(curr_key) == trim(root_key)) then
   tangled = .true.
   return
 end if

 if (depth > h%n_keys_stored) then
   ! Defensive cutoff against malformed deep loops.
   tangled = .true.
   return
 end if

 slot = h%get_index(curr_key)
 if (slot == -1) return

 if (h%vals(slot)%ASD(2) /= missA) then
   if (check_ID(h%vals(slot)%ASD(2), root_key, depth + 1)) then
     tangled = .true.
     return
   end if
 end if

 if (h%vals(slot)%ASD(3) /= missA) then
   if (check_ID(h%vals(slot)%ASD(3), root_key, depth + 1)) then
     tangled = .true.
     return
   end if
 end if
 end function check_ID


elemental subroutine initX(X)
 TYPE(CPED), intent(inout) :: X
 integer:: ii
 do ii=1,3
    X%ASD(ii)=missA
    X%NASD(ii)=missI
 enddo
 X%SEX=missI
 X%GEN= -1
 X%BYR=missI
 X%OBS=missI
 X%F=missR
end subroutine

subroutine triplets_bin_to_text_with_ids(bin_file, txt_file, order, nrec)
  use M_param
  implicit none
  character(len=*), intent(in) :: bin_file, txt_file
  integer, intent(in) :: order(:), nrec
  real(kind=r8) :: eps
  character(len=32) :: eps_str
  integer :: ios_env
  integer :: u_bin, u_txt, ios
  integer(kind=ki4) :: ver, nfile, row, col
  integer(kind=ki8) :: nnz, k
  real(kind=r8) :: val
  character(len=4) :: magic
  character(len=LEN_STR), allocatable :: orig_id(:)

  eps = 1.0d-6
  eps_str = ''
  call get_environment_variable('RELPED_TRIPLET_EPS', eps_str, status=ios_env)
  if (ios_env == 0) then
    read(eps_str, *, iostat=ios) eps
    if (ios /= 0) eps = 1.0d-6
    if (eps < 0.0_r8) eps = 0.0_r8
  end if

  allocate(orig_id(nrec))
  do k = 1_ki8, nrec
    orig_id(k) = h%vals(order(k))%ASD(1)
  end do

  open(newunit=u_bin, file=trim(bin_file), status='old', action='read', &
    form='unformatted', access='stream', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open triplet binary file:', trim(bin_file)
    deallocate(orig_id)
    return
  end if

  read(u_bin, iostat=ios) magic, ver, nfile, nnz
  if (ios /= 0) then
    write(*,*) 'Invalid triplet binary header:', trim(bin_file)
    close(u_bin)
    deallocate(orig_id)
    return
  end if

  open(newunit=u_txt, file=trim(txt_file), status='replace', action='write', &
    form='formatted', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Cannot open triplet text file:', trim(txt_file)
    close(u_bin)
    deallocate(orig_id)
    return
  end if

  do k = 1, nnz
    read(u_bin, iostat=ios) row, col, val
    if (ios /= 0) exit
    if (abs(val) > eps) then
      if (row >= 1 .and. row <= nrec .and. col >= 1 .and. col <= nrec) then
        write(u_txt, '(I0,1X,I0,1X,F12.6,1X,A,1X,A)') row, col, val, &
          trim(orig_id(row)), trim(orig_id(col))
      else
        write(u_txt, '(I0,1X,I0,1X,F12.6,1X,A,1X,A)') row, col, val, 'NA', 'NA'
      end if
    end if
  end do

  close(u_bin)
  close(u_txt)
  deallocate(orig_id)
end subroutine triplets_bin_to_text_with_ids

 subroutine inbreed_build_npd(nrec, sire_arr, dam_arr, F)
 integer, intent(in) :: nrec
 integer(kind=ki4), intent(in) :: sire_arr(:), dam_arr(:)
 real(kind=r8), intent(inout) :: F(:)
 integer, allocatable :: npd(:,:)

 allocate(npd(2, nrec))
 npd(1, 1:nrec) = sire_arr(1:nrec)
 npd(2, 1:nrec) = dam_arr(1:nrec)

 call Inb_F(npd, F)

 deallocate(npd)
 end subroutine inbreed_build_npd

   logical function Valid(XC)
   character(len=*),intent(in):: XC
   character(len=len(XC)):: XX
   integer:: N
   Valid = .true.
   XX=trim(adjustl(XC))
   if(iachar(XX(1:1)) < 33 .or. iachar(XX(1:1)) > 126 ) then
    Valid = .false.
   else
    N = len_trim(XX)
    select case (N) 
     case( :0 )
       Valid = .false.
     case( 1 )
       if(XX(1:1) == '0') Valid = .false.
     case( 2: )
       Valid = .true.
    end select
   endif
   end function

end module
