program popQC_main
  ! Population Quality Control Program
  ! Analyzes genomic SNP data from ReadFR for quality metrics
  !
  ! Input: GENO file (from ReadFR), PED file, MAP file
  ! Output: QC statistics, filtered data, quality reports
  !  
  use M_Kinds
  use M_Stamp
  use M_Variables
  use M_ReadPar
  use M_ReadFile
  use M_Hashtable
  use H_PedRead
  use H_MapRead
  use M_StrEdit
  use M_QCStepRunner
  use M_QCSNP
  use M_QCSexCheck
  use M_QCIDMerge, only: set_id_merge_context, merge_highly_identical_animals, finalize_deferred_identical_pairs
  use M_QCGenoIO, only: load_geno_file_mod, write_valid_geno_for_imputation_mod, set_qcgenoio_debug_mod
  use M_QCSetup, only: setup_parameters_mod, resolve_animal_key_fields_mod, load_pedigree_info_mod, set_qcsetup_debug_mod
  use M_QCMendel, only: calculate_mendelian_inconsistency_mod, recheck_mendel_after_parent_seek_mod
  use M_QCParentSeek, only: set_parent_seek_context, identify_correct_parentage_mod
  implicit none

  character(len=256) :: Par_File, version_str
  character(len=64) :: mode_arg2, mode_arg3
  logical :: parent_seek_only
  logical :: setup_debug_logs
  logical :: geno_debug_logs
  logical :: file_exists
  
  ! Variables for QC analysis
  character(len=256) :: ped_file
  integer :: n_animals, n_map_snps, n_geno_snps, n_x_snps
  integer, allocatable :: ped_sex(:)  ! 1=Male, 2=Female
  integer, allocatable :: x_snp_indices(:)  ! X chromosome SNP indices
  real(r8), allocatable :: sex_het_expected(:), sex_het_observed(:)
  integer, allocatable :: sex_qc_status(:)  ! 0=Pass, 1=Mismatch, 2=Uncertain
  integer, allocatable :: mendel_action(:)   ! 0=pass, 1=clear_parent, 2=delete
  integer, allocatable :: id_merge_seek_flag(:) ! 1=deferred identical-ID pair -> force parent-seek
  
  ! Variables for SNP-level QC
  ! NOTE: n_map_snps = MAP file total records (e.g., 76756, includes duplicates)
  !       n_geno_snps = GENO unique positions (e.g., 76614, indexed by Array_All)
  !       Output arrays sized by n_geno_snps and indexed by GENO position (1 to n_geno_snps)
  real(r8), allocatable :: snp_call_rate(:)  ! Per-SNP call rate (0.0 to 1.0)
  real(r8), allocatable :: snp_maf(:)  ! Minor Allele Frequency
  real(r8), allocatable :: snp_allele_freq(:,:)  ! Allele frequencies (2 columns)
  integer, allocatable :: snp_monomorphic(:)  ! 0=Polymorphic, 1=Monomorphic
  real(r8), allocatable :: snp_hwe_pvalue(:)  ! Hardy-Weinberg p-value
  integer, allocatable :: snp_qc_status(:)  ! 0=Pass, 1=Filter (MAF), 2=Monomorphic, 3=HWE_fail, 5=Non-autosome
  logical, allocatable :: snp_is_autosome(:)  ! .true. for Chr 1-18; .false. for sex/unknown chr
  real(r8), allocatable :: animal_call_rate(:)  ! Per-animal call rate
  integer :: n_snps_pass_qc, n_snps_fail_maf, n_snps_fail_low_maf, n_snps_fail_low_callrate, n_snps_monomorphic, n_snps_hwe
  integer :: k_map, j_col  ! loop vars for building snp_is_autosome mask

  ! HWE QC threshold
  real(r8) :: hwe_threshold
  

  ! Global PED/MAP records for shared access
  type(PEDHashTable), save :: PHT
  type(SNPInfo), allocatable :: map_records(:)
  type(MAPHashTable), save :: MHT
  integer(kind=ki1), allocatable :: geno_matrix(:,:)
  character(len=LEN_STR), save :: ped_animal_key_field
  character(len=LEN_STR), save :: geno_animal_key_field
  character(len=LEN_STR), allocatable :: geno_animal_ids(:)

  ! Initialize variables

  ! Print header
  version_str = '0.1 - Population Quality Control (popQC) Pipeline'
  call version(trim(version_str))
  print '(a)'
  call timestamp()
  print '(a)'

  ! Get parameter file from command line
  call getarg(1, Par_File)
  mode_arg2 = ''
  mode_arg3 = ''
  call getarg(2, mode_arg2)
  call getarg(3, mode_arg3)
  mode_arg2 = trim(to_upper(mode_arg2))
  mode_arg3 = trim(to_upper(mode_arg3))

  parent_seek_only = (mode_arg2 == 'PARENT_SEEK_ONLY' .or. mode_arg3 == 'PARENT_SEEK_ONLY')
  setup_debug_logs = (mode_arg2 == 'DEBUG_SETUP' .or. mode_arg3 == 'DEBUG_SETUP' .or. &
                      mode_arg2 == 'DEBUG' .or. mode_arg3 == 'DEBUG')
  geno_debug_logs = (mode_arg2 == 'DEBUG_GENO' .or. mode_arg3 == 'DEBUG_GENO' .or. &
                     mode_arg2 == 'DEBUG' .or. mode_arg3 == 'DEBUG')
  call set_qcsetup_debug_mod(setup_debug_logs)
  call set_qcgenoio_debug_mod(geno_debug_logs)
  
  if (len_trim(Par_File) < 1) then
    print *, "ERROR: Parameter file required"
    print *, ""
    print *, "Usage: popQC <parameter_file> [PARENT_SEEK_ONLY] [DEBUG_SETUP|DEBUG_GENO|DEBUG]"
    print *, ""
    print *, "Example: popQC parameter.txt"
    print *, "Example: popQC parameter.txt DEBUG_SETUP"
    print *, "Example: popQC parameter.txt PARENT_SEEK_ONLY DEBUG"
    stop
  end if

  ! Check if parameter file exists
  inquire(file=trim(Par_File), exist=file_exists)
  if (.not. file_exists) then
    print *, "ERROR: Parameter file not found: ", trim(Par_File)
    print *, "Please check the file path and try again."
    stop
  end if

  print *, ""
  print *, "========================================="
  print *, "popQC - Population Quality Control"
  print *, "========================================="
  print *, ""
  print *, "Reading parameter file: ", trim(Par_File)
  if (parent_seek_only) then
    print *, "Run mode: PARENT_SEEK_ONLY (execute through Step 8 only)"
  end if
  if (setup_debug_logs) then
    print *, "Setup logging: DEBUG_SETUP enabled"
  end if
  if (geno_debug_logs) then
    print *, "GENO logging: DEBUG_GENO enabled"
  end if

  ! ===========================================
  ! Main QC Workflow
  ! ===========================================
  ! 1. Read parameter file using M_readpar module
  !    Sets: PEDFile, MAPFile, GENOFile with FieldName and FieldLoc
  ! 2. Load pedigree information using PEDFile metadata
  ! 3. Load marker/map information
  ! 4. SNP-level QC analysis
  !    - Call rate (per animal, per SNP)
  !    - Allele frequencies and MAF
  !    - Hardy-Weinberg equilibrium test
  !    - Monomorphism detection
  ! 5. Calculate sex concordance QC
  ! 5.5 Merge highly identical IDs (>=95% concordance on QC-pass SNPs)
  ! 6. Mendelian inconsistency testing
  ! 7. Generate QC reports
  ! 8. Output filtered genotype data

  print *, ""
  print *, "========================================="
  print *, "[SETUP PHASE]"
  print *, "========================================="
  
  print *, ""
  print *, "Step 1: Reading parameters information from parameter file..."
 ! Load parameters from file - this populates global FileInfo structures
  ! (PEDFile, MAPFile, GENOFile with FieldName and FieldLoc)
  call setup_parameters_mod(trim(Par_File))
  call resolve_animal_key_fields_mod(ped_animal_key_field, geno_animal_key_field)
  
  ! Step 2: Load pedigree information into HashTable
  print *, ""
  print *, "Step 2: Loading pedigree information into HashTable..."
  call load_pedigree_info_mod(PHT, ped_file, n_animals, ped_animal_key_field)
  
  ! Print all records from HashTable for verification
  print *, ""
  print *, "--- Verification: All PED Records from HashTable ---"
  
  ! Since we use HashTable now, n_animals is mostly for GENO Matrix allocation
  ! Actual total animal count from GENO file (instead of PED)
  n_animals = N_recf(GENOFile%FileName) - GENOFile%Header
  print *, "  Total Animals in GENO file for QC: ", n_animals
  
  ! Step 3: Load marker/map information
  print *, ""
  print *, "Step 3: Loading marker information..."
  n_x_snps = 0
  call mht_load_map_file(MHT, map_records, n_map_snps)  ! n_map_snps = Total MAP records
  call detect_x_chromosome_snps(map_records, n_map_snps, x_snp_indices, n_x_snps)
  print *, "  X chromosome (CHR=20) SNPs: ", n_x_snps
  
  ! n_geno_snps = GENO field length
  call load_geno_file_mod(n_animals, map_records, geno_animal_key_field, n_geno_snps, geno_matrix, geno_animal_ids)
  print *, "  Unique GENO positions:      ", n_geno_snps, " (indexed by Array_All)"
  print *, "  Duplicate MAP SNPs:         ", n_map_snps - n_geno_snps
  
  print *, ""
  print *, "========================================="
  print *, "[QC ANALYSIS PHASE]"
  print *, "========================================="
 
  ! Step 4: SNP-level QC analysis
  print *, ""
  print *, "Step 4: Performing SNP-level quality control..."
  if (n_geno_snps > 0) then
    allocate(snp_call_rate(n_geno_snps))
    allocate(snp_maf(n_geno_snps))
    allocate(snp_allele_freq(n_geno_snps, 2))
    allocate(snp_monomorphic(n_geno_snps))
    allocate(snp_hwe_pvalue(n_geno_snps))
    allocate(snp_qc_status(n_geno_snps))
    allocate(animal_call_rate(n_animals))

    ! Build autosome mask: Chr 1-18 = valid autosome; Chr>18 = sex/unknown
    allocate(snp_is_autosome(n_geno_snps))
    snp_is_autosome = .false.
    do k_map = 1, n_map_snps
      j_col = map_records(k_map)%Array_All
      if (j_col >= 1 .and. j_col <= n_geno_snps) then
        if (map_records(k_map)%Chr >= 1 .and. map_records(k_map)%Chr <= 18) &
          snp_is_autosome(j_col) = .true.
      end if
    end do
    print *, "  Autosome SNPs  (Chr  1-18): ", count(snp_is_autosome)
    print *, "  Sex/Unk SNPs   (Chr >18):  ", n_geno_snps - count(snp_is_autosome)

    print *, "  DEBUG: About to call calculate_snp_statistics with:"
    print *, "    n_geno_snps=", n_geno_snps
    print *, "    size(snp_call_rate)=", size(snp_call_rate)
    print *, "    size(snp_maf)=", size(snp_maf)
    print *, "    size(snp_allele_freq,1)=", size(snp_allele_freq, 1)
    
    call calculate_snp_statistics(trim(Par_File), n_geno_snps, n_animals, geno_matrix, &
            snp_call_rate, snp_allele_freq, snp_maf, animal_call_rate, snp_hwe_pvalue, hwe_threshold)
    
    call identify_monomorphic_snps(n_geno_snps, snp_maf, snp_monomorphic)
    
    if (hwe_threshold > 0.0_r8) then
      print *, "  HWE p-values calculated in calculate_snp_statistics"
      print *, "  HWE threshold: ", hwe_threshold
    else
      print *, "  HWE test disabled (HWE_THRESHOLD=0.0)"
    end if
    

    call determine_snp_qc_status(n_geno_snps, snp_call_rate, snp_maf, snp_monomorphic, &
             snp_hwe_pvalue, hwe_threshold, snp_qc_status, n_snps_pass_qc, &
             n_snps_fail_low_maf, n_snps_fail_low_callrate, n_snps_fail_maf, n_snps_monomorphic, n_snps_hwe, &
             snp_is_autosome)

    call report_snp_qc(n_geno_snps, n_animals, snp_call_rate, snp_maf, snp_monomorphic, &
           snp_hwe_pvalue, snp_qc_status, animal_call_rate, &
           n_snps_pass_qc, n_snps_fail_low_maf, n_snps_fail_low_callrate, n_snps_fail_maf, n_snps_monomorphic, n_snps_hwe, &
           snp_is_autosome)
    
    ! Keep SNP QC arrays alive until Step 8.
    ! They are reused by ID-merge (Step 6) and parent-seek G filtering.
  else
    print *, "WARNING: No SNPs loaded. SNP QC skipped."
  end if
  
  ! Step 5: Calculate sex concordance QC metrics
  print *, ""
  print *, "Step 5: Calculating sex concordance metrics (X-chromosome based)..."
  if (POPQC%sex_concordence_threshold <= 0.0) then
    print *, "  INFO: SEX_CONCORDENCE <= 0. Sex concordance check skipped."
  else if (n_x_snps == 0) then
    print *, "  WARNING: No X chromosome (CHR=20) SNPs found in MAP. Sex concordance check skipped."
  else
    allocate(ped_sex(n_animals))
    allocate(sex_het_expected(n_animals))
    allocate(sex_het_observed(n_animals))
    allocate(sex_qc_status(n_animals))
    call fill_ped_sex_from_hashtable(PHT, n_animals, geno_animal_ids, ped_sex)
    call calculate_sex_concordance(n_animals, n_x_snps, x_snp_indices, geno_matrix, ped_sex, &
                          real(POPQC%sex_concordence_threshold, r8), sex_het_observed, sex_qc_status)
    call report_sex_concordance(n_animals, geno_animal_ids, ped_sex, sex_het_observed, &
                        sex_qc_status, real(POPQC%sex_concordence_threshold, r8))
    ! NOTE: sex arrays kept alive so Step 6 can consider sex concordance in Mendel decision
  end if

  ! Step 6: Merge highly identical IDs before Mendelian/parent seek
  print *, ""
  print *, "Step 6: Merging highly identical IDs (>=95% SNP concordance)..."
  if (.not. allocated(id_merge_seek_flag)) allocate(id_merge_seek_flag(n_animals))
  id_merge_seek_flag = 0
  if (allocated(snp_qc_status)) then
    call set_id_merge_context(PHT, map_records, geno_matrix, geno_animal_ids, ped_animal_key_field, &
        snp_call_rate, snp_maf, snp_hwe_pvalue, hwe_threshold, id_merge_seek_flag)
    call merge_highly_identical_animals(n_animals, n_geno_snps, snp_qc_status, 0.95_r8)
  else
    print *, "  WARNING: SNP QC status unavailable. ID merge skipped."
  end if

  ! Step 7: Perform Mendelian inconsistency testing (Parental check)
  !   Two-pass algorithm:
  !     Loop 1 – SNP scan: count parent(AA)↔child(BB) errors per SNP;
  !              SNPs with error rate > mendel_snp_threshold are removed.
  !     Loop 2 – Sample scan: recompute per-animal error rate using only
  !              clean SNPs; flag animals above mendel_viol_threshold.
  !   Decision per flagged animal:
  !     sex_fail OR err > 5×threshold OR call_rate < animal_cr → DELETE
  !     otherwise → CLEAR_PARENT (sire/dam set to 0; Step 8 will seek parents)
  print *, ""
  print *, "Step 7: Performing Mendelian inconsistency testing (two-pass)..."
  if (POPQC%mendel_viol_threshold > 0.0) then
    allocate(mendel_action(n_animals))
    mendel_action = 0
    if (allocated(sex_qc_status)) then
      call calculate_mendelian_inconsistency_mod(n_animals, n_geno_snps, POPQC%mendel_snp_threshold, POPQC%mendel_viol_threshold, &
          real(POPQC%animal_cr, r8), geno_matrix, geno_animal_ids, PHT, in_sex_qc_status=sex_qc_status, out_action=mendel_action)
    else
      call calculate_mendelian_inconsistency_mod(n_animals, n_geno_snps, POPQC%mendel_snp_threshold, POPQC%mendel_viol_threshold, &
          real(POPQC%animal_cr, r8), geno_matrix, geno_animal_ids, PHT, out_action=mendel_action)
    end if
  else
    print *, "  INFO: MENDIAN_INCONSISTENCY animal threshold <= 0. Step 7 skipped."
  end if

  ! Post-process action=3 (SEX_SWITCH): flip registered sex, reset QC status,
  ! then reclassify as CLEAR_PARENT (action=1) so parent-seek runs with corrected sex.
  if (allocated(mendel_action) .and. allocated(ped_sex) .and. allocated(sex_qc_status)) then
    block
      integer :: i_sw
      logical :: any_sw
      any_sw = .false.
      do i_sw = 1, n_animals
        if (mendel_action(i_sw) /= 3) cycle
        if (.not. any_sw) then
          print *, ""
          print *, "  SEX_SWITCH post-processing (sex corrected before parent-seek):"
          print *, "  Animal                    Old_Sex  -> New_Sex"
          print *, "  -----------------------------------------------"
        end if
        any_sw = .true.
        if (ped_sex(i_sw) == MALE) then
          ped_sex(i_sw) = FEMALE
          print '(4x, a24, 2x, a, a)', trim(geno_animal_ids(i_sw)), 'MALE   ', '-> FEMALE'
        else if (ped_sex(i_sw) == FEMALE) then
          ped_sex(i_sw) = MALE
          print '(4x, a24, 2x, a, a)', trim(geno_animal_ids(i_sw)), 'FEMALE ', '-> MALE'
        else
          print '(4x, a24, 2x, a)',    trim(geno_animal_ids(i_sw)), 'UNKNOWN (no flip, parent-seek with both slots)'
        end if
        sex_qc_status(i_sw) = 0  ! clear mismatch flag → valid parent candidate
        mendel_action(i_sw) = 1  ! CLEAR_PARENT → enters Step 8
      end do
    end block
  end if

  ! Step 8: Find correct parentage for animals whose parents were cleared (action=1)
  print *, ""
  print *, "Step 8: Identifying correct parentage for flagged animals..."
  if (POPQC%parent_seek_threshold > 0.0) then
    ! Ensure mendel_action is allocated (all-zero if Mendelian check was skipped)
    if (.not. allocated(mendel_action)) then
      allocate(mendel_action(n_animals))
      mendel_action = 0
    end if

    ! Deferred identical-ID pairs: force into parent-seek unless already DELETE.
    if (allocated(id_merge_seek_flag)) then
      block
        integer :: i_force, n_forced
        n_forced = 0
        do i_force = 1, n_animals
          if (id_merge_seek_flag(i_force) /= 1) cycle
          if (len_trim(geno_animal_ids(i_force)) == 0) cycle
          if (mendel_action(i_force) == 2) cycle
          if (mendel_action(i_force) /= 1) n_forced = n_forced + 1
          mendel_action(i_force) = 1
        end do
        if (n_forced > 0) then
          print '(a, i8)', "  Deferred identical-ID rows forced to parent-seek: ", n_forced
        end if
      end block
    end if

    call set_parent_seek_context(PHT, map_records, geno_matrix, geno_animal_ids, ped_animal_key_field, &
        snp_call_rate, snp_maf, snp_hwe_pvalue, hwe_threshold)
    if (allocated(sex_qc_status)) then
      call identify_correct_parentage_mod(n_animals, n_geno_snps, mendel_action, real(POPQC%parent_seek_threshold, r8), POPQC, &
          in_sex_qc_status=sex_qc_status, in_corrected_sex=ped_sex)
    else
      call identify_correct_parentage_mod(n_animals, n_geno_snps, mendel_action, real(POPQC%parent_seek_threshold, r8), POPQC)
    end if
  else
    print *, "  INFO: PARENT_SEEK <= 0. Step 8 skipped."
  end if

  if (.not. parent_seek_only) then
    print *, ""
    print *, "Step 9: Finalizing deferred identical-ID pairs..."
    call finalize_deferred_identical_pairs(n_animals, n_geno_snps, 0.95_r8)

    print *, ""
    print *, "Step 10: Post-parent-seek Mendelian recheck (deferred delete)..."
    if (allocated(mendel_action)) then
      call recheck_mendel_after_parent_seek_mod(n_animals, n_geno_snps, &
          POPQC%mendel_viol_threshold, real(POPQC%animal_cr, r8), geno_matrix, geno_animal_ids, PHT, mendel_action)
    else
      print *, "  INFO: Mendelian action state unavailable. Recheck skipped."
    end if

    ! Step 11: Write valid genotype records for downstream imputation
    print *, ""
    print *, "Step 11: Writing valid genotype output for imputation..."
    if (allocated(mendel_action)) then
      call write_valid_geno_for_imputation_mod(n_animals, n_geno_snps, geno_matrix, geno_animal_ids, PHT, &
          n_snps_pass_qc, snp_qc_status, snp_is_autosome, mendel_action)
    else
      call write_valid_geno_for_imputation_mod(n_animals, n_geno_snps, geno_matrix, geno_animal_ids, PHT, &
          n_snps_pass_qc, snp_qc_status, snp_is_autosome)
    end if
  else
    print *, ""
    print *, "Step 9-11 skipped by PARENT_SEEK_ONLY mode."
  end if

  if (allocated(mendel_action)) deallocate(mendel_action)

  ! Release sex QC arrays (no longer needed after Step 8)
  if (allocated(ped_sex))          deallocate(ped_sex)
  if (allocated(sex_het_expected)) deallocate(sex_het_expected)
  if (allocated(sex_het_observed)) deallocate(sex_het_observed)
  if (allocated(sex_qc_status))    deallocate(sex_qc_status)

  ! Release SNP QC arrays (no longer needed after Step 8)
  if (allocated(snp_call_rate))   deallocate(snp_call_rate)
  if (allocated(snp_maf))         deallocate(snp_maf)
  if (allocated(snp_allele_freq)) deallocate(snp_allele_freq)
  if (allocated(snp_monomorphic)) deallocate(snp_monomorphic)
  if (allocated(snp_hwe_pvalue))  deallocate(snp_hwe_pvalue)
  if (allocated(snp_qc_status))   deallocate(snp_qc_status)
  if (allocated(snp_is_autosome)) deallocate(snp_is_autosome)
  if (allocated(animal_call_rate)) deallocate(animal_call_rate)
  if (allocated(id_merge_seek_flag)) deallocate(id_merge_seek_flag)
  
  ! Step 12: Finalize and cleanup
  print *, ""
  print *, "Quality control analysis completed successfully."
  print *, ""
  call timestamp()
  end program popQC_main
