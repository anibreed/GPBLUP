module M_QCMendel
  use M_Kinds
  use M_Variables, only: LEN_STR
  use M_ped, only: CPED
  use M_Hashtable
  use H_PedRead
  implicit none

  private
  public :: calculate_mendelian_inconsistency_mod
  public :: recheck_mendel_after_parent_seek_mod
  public :: build_geno_row_map
  public :: find_geno_row_map
  public :: is_duo_mendel_error
  public :: is_trio_mendel_error

contains

  subroutine calculate_mendelian_inconsistency_mod(in_n_animals, in_n_snps, in_snp_error_threshold, in_error_threshold, &
      in_animal_cr, in_geno_matrix, in_geno_animal_ids, in_pht, in_sex_qc_status, out_action)
    integer, intent(in)           :: in_n_animals
    integer, intent(in)           :: in_n_snps
    real,    intent(in)           :: in_snp_error_threshold
    real,    intent(in)           :: in_error_threshold
    real(r8), intent(in)          :: in_animal_cr
    integer(kind=ki1), intent(in) :: in_geno_matrix(in_n_animals, in_n_snps)
    character(len=*), intent(in)  :: in_geno_animal_ids(in_n_animals)
    type(PEDHashTable), intent(in) :: in_pht
    integer, intent(in), optional :: in_sex_qc_status(:)
    integer, intent(out)          :: out_action(in_n_animals)

    real(r8), parameter :: DELETE_FACTOR = 5.0_r8

    integer :: i, j
    integer :: sire_row, dam_row, parent_row, child_mode
    integer :: g_child, g_sire, g_dam
    integer :: n_non_missing, n_snp_removed
    integer :: n_trio, n_duo, n_no_parent
    integer :: n_child_not_in_ped, n_parent_id_missing, n_parent_not_genotyped
    integer :: n_eval, n_cleared, n_deleted, n_sex_switch
    real(r8) :: err_rate, cr_rate
    logical :: sex_also_fail, hard_delete
    character(len=LEN_STR) :: mode_str, action_str
    character(len=24) :: animal_disp
    character(len=4) :: mode_disp
    character(len=20) :: reason_str
    type(CPED) :: ped_child
    logical :: found_child
    character(len=LEN_STR), allocatable :: map_keys(:)
    integer,                allocatable :: map_rows(:)
    logical,                allocatable :: map_used(:)
    integer,                allocatable :: valid_cnt(:), error_cnt(:), mode_cnt(:)
    integer,                allocatable :: snp_valid_cnt(:), snp_error_cnt(:)
    logical,                allocatable :: snp_mask(:)

    allocate(valid_cnt(in_n_animals), error_cnt(in_n_animals), mode_cnt(in_n_animals))
    allocate(snp_valid_cnt(in_n_snps), snp_error_cnt(in_n_snps), snp_mask(in_n_snps))
    valid_cnt     = 0;  error_cnt  = 0;  mode_cnt  = 0
    snp_valid_cnt = 0;  snp_error_cnt = 0;  snp_mask = .true.
    out_action    = 0

    call build_geno_row_map(in_n_animals, in_geno_animal_ids, map_keys, map_rows, map_used)

    print *, "  Loop 1: SNP scan - counting Mendel errors per SNP..."

    do i = 1, in_n_animals
      if (len_trim(in_geno_animal_ids(i)) == 0) cycle
      found_child = pht_search(in_pht, trim(in_geno_animal_ids(i)), ped_child)
      if (.not. found_child) cycle

      sire_row = 0;  dam_row = 0
      if (len_trim(ped_child%SIRE) > 0 .and. trim(ped_child%SIRE) /= '0') &
          sire_row = find_geno_row_map(trim(ped_child%SIRE), map_keys, map_rows, map_used)
      if (len_trim(ped_child%DAM)  > 0 .and. trim(ped_child%DAM)  /= '0') &
          dam_row  = find_geno_row_map(trim(ped_child%DAM),  map_keys, map_rows, map_used)

      if      (sire_row > 0 .and. dam_row > 0) then;  child_mode = 3
      else if (sire_row > 0 .or.  dam_row > 0) then;  child_mode = 2
      else;  cycle
      end if
      mode_cnt(i) = child_mode
      parent_row  = sire_row;  if (parent_row == 0) parent_row = dam_row

      do j = 1, in_n_snps
        g_child = int(in_geno_matrix(i, j));  if (g_child == 9) cycle
        if (child_mode == 3) then
          g_sire = int(in_geno_matrix(sire_row, j));  g_dam = int(in_geno_matrix(dam_row, j))
          if (g_sire == 9 .or. g_dam == 9) cycle
          snp_valid_cnt(j) = snp_valid_cnt(j) + 1
          if (is_trio_mendel_error(g_child, g_sire, g_dam)) snp_error_cnt(j) = snp_error_cnt(j) + 1
        else
          g_sire = int(in_geno_matrix(parent_row, j));  if (g_sire == 9) cycle
          snp_valid_cnt(j) = snp_valid_cnt(j) + 1
          if (is_duo_mendel_error(g_child, g_sire)) snp_error_cnt(j) = snp_error_cnt(j) + 1
        end if
      end do
    end do

    n_snp_removed = 0
    do j = 1, in_n_snps
      if (snp_valid_cnt(j) > 0) then
        if (real(snp_error_cnt(j), r8) / real(snp_valid_cnt(j), r8) > real(in_snp_error_threshold, r8)) then
          snp_mask(j) = .false.
          n_snp_removed = n_snp_removed + 1
        end if
      end if
    end do

    print '(a, i8)', "    Total SNPs            : ", in_n_snps
    print '(a, i8)', "    SNPs evaluable (>=1 pair): ", count(snp_valid_cnt > 0)
    print '(a, i8)', "    SNPs not evaluable    : ", count(snp_valid_cnt == 0)
    print '(a, a)',  "      (no genotyped parent-offspring pair at this SNP - kept by default)"
    print '(a, i8)', "    SNPs removed (err>thr): ", n_snp_removed
    print '(a, i8)', "    SNPs retained         : ", count(snp_mask)
    print '(a, i8)', "      = evaluable passed  : ", count(snp_valid_cnt > 0) - n_snp_removed
    print '(a, i8)', "      + not evaluable     : ", count(snp_valid_cnt == 0)

    print *, "  Loop 2: Sample scan - computing per-animal Mendel error rates..."

    n_trio = 0;  n_duo = 0;  n_no_parent = 0
    n_child_not_in_ped = 0;  n_parent_id_missing = 0;  n_parent_not_genotyped = 0
    n_eval = 0;  n_cleared = 0;  n_deleted = 0;  n_sex_switch = 0

    do i = 1, in_n_animals
      if (len_trim(in_geno_animal_ids(i)) == 0) then
        n_no_parent = n_no_parent + 1
        cycle
      end if

      found_child = pht_search(in_pht, trim(in_geno_animal_ids(i)), ped_child)
      if (.not. found_child) then
        n_child_not_in_ped = n_child_not_in_ped + 1
        n_no_parent        = n_no_parent + 1
        cycle
      end if

      sire_row = 0; dam_row = 0
      if (len_trim(ped_child%SIRE) > 0 .and. trim(ped_child%SIRE) /= '0') &
          sire_row = find_geno_row_map(trim(ped_child%SIRE), map_keys, map_rows, map_used)
      if (len_trim(ped_child%DAM) > 0 .and. trim(ped_child%DAM) /= '0') &
          dam_row = find_geno_row_map(trim(ped_child%DAM), map_keys, map_rows, map_used)

      if ((len_trim(ped_child%SIRE) == 0 .or. trim(ped_child%SIRE) == '0') .and. &
          (len_trim(ped_child%DAM)  == 0 .or. trim(ped_child%DAM)  == '0')) then
        n_parent_id_missing = n_parent_id_missing + 1
      else if (sire_row == 0 .and. dam_row == 0) then
        n_parent_not_genotyped = n_parent_not_genotyped + 1
      end if

      if      (sire_row > 0 .and. dam_row > 0) then;  child_mode = 3;  n_trio = n_trio + 1
      else if (sire_row > 0 .or.  dam_row > 0) then;  child_mode = 2;  n_duo  = n_duo + 1
      else;  n_no_parent = n_no_parent + 1; cycle
      end if
      parent_row = sire_row; if (parent_row == 0) parent_row = dam_row

      do j = 1, in_n_snps
        if (.not. snp_mask(j)) cycle
        g_child = int(in_geno_matrix(i, j)); if (g_child == 9) cycle
        if (child_mode == 3) then
          g_sire = int(in_geno_matrix(sire_row, j)); g_dam = int(in_geno_matrix(dam_row, j))
          if (g_sire == 9 .or. g_dam == 9) cycle
          valid_cnt(i) = valid_cnt(i) + 1
          if (is_trio_mendel_error(g_child, g_sire, g_dam)) error_cnt(i) = error_cnt(i) + 1
        else
          g_sire = int(in_geno_matrix(parent_row, j)); if (g_sire == 9) cycle
          valid_cnt(i) = valid_cnt(i) + 1
          if (is_duo_mendel_error(g_child, g_sire)) error_cnt(i) = error_cnt(i) + 1
        end if
      end do

      if (valid_cnt(i) <= 0) cycle
      n_eval = n_eval + 1
      err_rate = real(error_cnt(i), r8) / real(valid_cnt(i), r8)
      if (err_rate <= real(in_error_threshold, r8)) cycle

      n_non_missing = 0
      do j = 1, in_n_snps
        if (int(in_geno_matrix(i, j)) /= 9) n_non_missing = n_non_missing + 1
      end do
      cr_rate = real(n_non_missing, r8) / real(in_n_snps, r8)

      sex_also_fail = .false.
      if (present(in_sex_qc_status)) then
        if (i <= size(in_sex_qc_status)) sex_also_fail = (in_sex_qc_status(i) == 1)
      end if

      hard_delete = (in_animal_cr > 0.0_r8 .and. cr_rate < in_animal_cr)

      if (hard_delete) then
        out_action(i) = 2
        n_deleted = n_deleted + 1
      else if (sex_also_fail) then
        out_action(i) = 3
        n_sex_switch = n_sex_switch + 1
      else
        out_action(i) = 1
        n_cleared = n_cleared + 1
      end if
    end do

    print *, ""
    print *, "  Mendelian inconsistency summary (two-pass):"
    print '(a, f8.4)', "    SNP error threshold:           ", real(in_snp_error_threshold, r8)
    print '(a, f8.4)', "    Animal error threshold:        ", real(in_error_threshold, r8)
    print '(a, f8.4)', "    Deferred delete threshold(x5): ", DELETE_FACTOR * real(in_error_threshold, r8)
    print '(a, i8)',   "    SNPs removed   (Loop 1):       ", n_snp_removed
    print '(a, i8)',   "    Trio animals:                  ", n_trio
    print '(a, i8)',   "    Duo animals:                   ", n_duo
    print '(a, i8)',   "    No evaluable parent:           ", n_no_parent
    print '(a, i8)',   "      - Child not in PED:          ", n_child_not_in_ped
    print '(a, i8)',   "      - Parent ID missing in PED:  ", n_parent_id_missing
    print '(a, i8)',   "      - Parent not genotyped:      ", n_parent_not_genotyped
    print '(a, i8)',   "    Evaluated animals:             ", n_eval
    print '(a, i8)',   "    PASS:                          ", n_eval - n_cleared - n_sex_switch - n_deleted
    print '(a, i8)',   "    CLEAR_PARENT (action=1):       ", n_cleared
    print '(a, i8)',   "    SEX_SWITCH   (action=3):       ", n_sex_switch
    print '(a, i8)',   "    DELETE        (action=2):       ", n_deleted
    if (n_sex_switch > 0) print '(a)', "    (SEX_SWITCH animals: sex will be flipped in main then re-enter parent-seek)"
    print '(a)', "    (High Mendel-error deletes are deferred to post-parent-seek recheck)"

    if (n_cleared + n_sex_switch + n_deleted > 0) then
      print *, ""
      print *, "  Flagged animal details:"
      print *, "  Note: ValidCmpSNPs = Mendel-comparable SNP count after Loop-1 SNP mask"
      print *, "        (NOT the parent-seek window/bin sampled SNP count)."
      print *, "  Animal                   Mode  Errors  ValidCmpSNPs  ErrRate  CallRate SexFail  Action        Reason"
      print *, "  -------------------------------------------------------------------------------------------------------"
      do i = 1, in_n_animals
        if (out_action(i) == 0) cycle
        select case (out_action(i))
        case (3)
          action_str = 'SEX_SWITCH  '
          reason_str = 'SEX_MISMATCH_ONLY'
        case (2)
          action_str = 'DELETE      '
          n_non_missing = 0
          do j = 1, in_n_snps
            if (int(in_geno_matrix(i, j)) /= 9) n_non_missing = n_non_missing + 1
          end do
          cr_rate = real(n_non_missing, r8) / real(in_n_snps, r8)
          if (in_animal_cr > 0.0_r8 .and. cr_rate < in_animal_cr) then
            reason_str = 'LOW_CALL_RATE'
          else
            reason_str = 'POSTSEEK_DELETE'
          end if
        case default
          action_str = 'CLEAR_PARENT'
          reason_str = 'ERR_GT_THRESHOLD'
        end select

        mode_str = merge('TRIO', 'DUO ', mode_cnt(i) == 3)
        sex_also_fail = .false.
        if (present(in_sex_qc_status)) then
          if (i <= size(in_sex_qc_status)) sex_also_fail = (in_sex_qc_status(i) == 1)
        end if
        n_non_missing = 0
        do j = 1, in_n_snps
          if (int(in_geno_matrix(i, j)) /= 9) n_non_missing = n_non_missing + 1
        end do
        cr_rate = real(n_non_missing, r8) / real(in_n_snps, r8)
        err_rate = real(error_cnt(i), r8) / real(max(1, valid_cnt(i)), r8)
        animal_disp = adjustl(trim(in_geno_animal_ids(i)))
        mode_disp = adjustl(trim(mode_str))
        print '(2x, a24, 1x, a4, 2x, i6, 1x, i10, 2x, f7.4, 2x, f7.4, 3x, l1, 2x, a12, 2x, a20)', &
            animal_disp, mode_disp, error_cnt(i), valid_cnt(i), err_rate, cr_rate, sex_also_fail, trim(action_str), trim(reason_str)
      end do
    end if

    if (allocated(map_keys)) deallocate(map_keys)
    if (allocated(map_rows)) deallocate(map_rows)
    if (allocated(map_used)) deallocate(map_used)
    if (allocated(snp_valid_cnt)) deallocate(snp_valid_cnt)
    if (allocated(snp_error_cnt)) deallocate(snp_error_cnt)
    if (allocated(snp_mask)) deallocate(snp_mask)
    deallocate(valid_cnt, error_cnt, mode_cnt)
  end subroutine calculate_mendelian_inconsistency_mod

  subroutine build_geno_row_map(in_n, in_ids, out_keys, out_rows, out_used)
    integer, intent(in) :: in_n
    character(len=*), intent(in) :: in_ids(:)
    character(len=LEN_STR), allocatable, intent(out) :: out_keys(:)
    integer, allocatable, intent(out) :: out_rows(:)
    logical, allocatable, intent(out) :: out_used(:)
    integer :: i, map_size

    map_size = htab_get_prime_size(max(101, in_n * 2 + 1))
    allocate(out_keys(map_size), out_rows(map_size), out_used(map_size))
    out_keys = ''
    out_rows = 0
    out_used = .false.

    do i = 1, in_n
      if (len_trim(in_ids(i)) == 0) cycle
      call insert_geno_row_map(trim(in_ids(i)), i, out_keys, out_rows, out_used)
    end do
  end subroutine build_geno_row_map

  subroutine insert_geno_row_map(in_key, in_row, io_keys, io_rows, io_used)
    character(len=*), intent(in) :: in_key
    integer, intent(in) :: in_row
    character(len=LEN_STR), intent(inout) :: io_keys(:)
    integer, intent(inout) :: io_rows(:)
    logical, intent(inout) :: io_used(:)
    integer :: idx, start_idx, nslot

    nslot = size(io_keys)
    start_idx = htab_hash_string(trim(in_key), nslot) + 1
    idx = start_idx

    do
      if (.not. io_used(idx)) then
        io_used(idx) = .true.
        io_keys(idx) = trim(in_key)
        io_rows(idx) = in_row
        return
      end if
      if (trim(io_keys(idx)) == trim(in_key)) then
        io_rows(idx) = in_row
        return
      end if
      idx = idx + 1
      if (idx > nslot) idx = 1
      if (idx == start_idx) exit
    end do
  end subroutine insert_geno_row_map

  integer function find_geno_row_map(in_key, in_keys, in_rows, in_used) result(out_row)
    character(len=*), intent(in) :: in_key
    character(len=LEN_STR), intent(in) :: in_keys(:)
    integer, intent(in) :: in_rows(:)
    logical, intent(in) :: in_used(:)
    integer :: idx, start_idx, nslot

    out_row = 0
    if (len_trim(in_key) == 0) return

    nslot = size(in_keys)
    start_idx = htab_hash_string(trim(in_key), nslot) + 1
    idx = start_idx

    do
      if (.not. in_used(idx)) return
      if (trim(in_keys(idx)) == trim(in_key)) then
        out_row = in_rows(idx)
        return
      end if
      idx = idx + 1
      if (idx > nslot) idx = 1
      if (idx == start_idx) return
    end do
  end function find_geno_row_map

  logical function is_duo_mendel_error(in_child, in_parent) result(out_err)
    integer, intent(in) :: in_child, in_parent
    out_err = .false.
    if (in_parent == 0 .and. in_child == 2) out_err = .true.
    if (in_parent == 2 .and. in_child == 0) out_err = .true.
  end function is_duo_mendel_error

  logical function is_trio_mendel_error(in_child, in_sire, in_dam) result(out_err)
    integer, intent(in) :: in_child, in_sire, in_dam
    integer :: s_min, s_max, d_min, d_max
    integer :: c_min, c_max

    if (in_sire == 0) then
      s_min = 0; s_max = 0
    else if (in_sire == 1) then
      s_min = 0; s_max = 1
    else
      s_min = 1; s_max = 1
    end if

    if (in_dam == 0) then
      d_min = 0; d_max = 0
    else if (in_dam == 1) then
      d_min = 0; d_max = 1
    else
      d_min = 1; d_max = 1
    end if

    c_min = s_min + d_min
    c_max = s_max + d_max
    out_err = (in_child < c_min .or. in_child > c_max)
  end function is_trio_mendel_error

  subroutine recheck_mendel_after_parent_seek_mod(in_n_animals, in_n_snps, in_error_threshold, &
                                                  in_animal_cr, in_geno_matrix, in_geno_animal_ids, &
                                                  in_pht, io_action)
    integer, intent(in) :: in_n_animals, in_n_snps
    real, intent(in) :: in_error_threshold
    real(r8), intent(in) :: in_animal_cr
    integer(kind=ki1), intent(in) :: in_geno_matrix(in_n_animals, in_n_snps)
    character(len=*), intent(in) :: in_geno_animal_ids(in_n_animals)
    type(PEDHashTable), intent(in) :: in_pht
    integer, intent(inout) :: io_action(:)

    real(r8), parameter :: DELETE_FACTOR = 5.0_r8
    integer, parameter :: MIN_VALID_SNPS = 100

    integer :: i, j
    integer :: sire_row, dam_row, parent_row, g_child, g_sire, g_dam
    integer :: n_valid, n_err, n_non_missing
    integer :: n_eval, n_del_now, n_keep
    real(r8) :: err_rate, cr_rate, del_thr
    type(CPED) :: ped_child
    logical :: found_child
    character(len=LEN_STR), allocatable :: map_keys(:)
    integer, allocatable :: map_rows(:)
    logical, allocatable :: map_used(:)

    if (in_error_threshold <= 0.0) then
      print *, "  INFO: Mendel threshold <= 0. Recheck skipped."
      return
    end if

    call build_geno_row_map(in_n_animals, in_geno_animal_ids, map_keys, map_rows, map_used)
    del_thr = DELETE_FACTOR * real(in_error_threshold, r8)
    n_eval = 0
    n_del_now = 0
    n_keep = 0

    do i = 1, in_n_animals
      if (i > size(io_action)) exit
      if (io_action(i) /= 1) cycle
      if (io_action(i) == 2) cycle
      if (len_trim(in_geno_animal_ids(i)) == 0) cycle

      found_child = pht_search(in_pht, trim(in_geno_animal_ids(i)), ped_child)
      if (.not. found_child) cycle

      sire_row = 0
      dam_row = 0
      if (len_trim(ped_child%SIRE) > 0 .and. trim(ped_child%SIRE) /= '0') then
        sire_row = find_geno_row_map(trim(ped_child%SIRE), map_keys, map_rows, map_used)
      end if
      if (len_trim(ped_child%DAM) > 0 .and. trim(ped_child%DAM) /= '0') then
        dam_row = find_geno_row_map(trim(ped_child%DAM), map_keys, map_rows, map_used)
      end if

      if (sire_row == 0 .and. dam_row == 0) cycle
      parent_row = sire_row
      if (parent_row == 0) parent_row = dam_row

      n_valid = 0
      n_err = 0
      do j = 1, in_n_snps
        g_child = int(in_geno_matrix(i, j))
        if (g_child == 9) cycle
        if (sire_row > 0 .and. dam_row > 0) then
          g_sire = int(in_geno_matrix(sire_row, j))
          g_dam  = int(in_geno_matrix(dam_row, j))
          if (g_sire == 9 .or. g_dam == 9) cycle
          n_valid = n_valid + 1
          if (is_trio_mendel_error(g_child, g_sire, g_dam)) n_err = n_err + 1
        else
          g_sire = int(in_geno_matrix(parent_row, j))
          if (g_sire == 9) cycle
          n_valid = n_valid + 1
          if (is_duo_mendel_error(g_child, g_sire)) n_err = n_err + 1
        end if
      end do

      if (n_valid < MIN_VALID_SNPS) cycle
      n_eval = n_eval + 1
      err_rate = real(n_err, r8) / real(n_valid, r8)

      n_non_missing = 0
      do j = 1, in_n_snps
        if (int(in_geno_matrix(i, j)) /= 9) n_non_missing = n_non_missing + 1
      end do
      cr_rate = real(n_non_missing, r8) / real(in_n_snps, r8)

      if (err_rate > del_thr .or. &
          (in_animal_cr > 0.0_r8 .and. cr_rate < in_animal_cr)) then
        io_action(i) = 2
        n_del_now = n_del_now + 1
      else
        n_keep = n_keep + 1
      end if
    end do

    print '(a, f8.4)', "  Deferred delete threshold (x5):      ", del_thr
    print '(a, i8)',   "  Action=1 animals re-evaluated:        ", n_eval
    print '(a, i8)',   "  Newly deleted after recheck:          ", n_del_now
    print '(a, i8)',   "  Retained after recheck:               ", n_keep

    if (allocated(map_keys)) deallocate(map_keys)
    if (allocated(map_rows)) deallocate(map_rows)
    if (allocated(map_used)) deallocate(map_used)
  end subroutine recheck_mendel_after_parent_seek_mod

end module M_QCMendel
