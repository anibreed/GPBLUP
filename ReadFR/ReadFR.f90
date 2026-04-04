program ReadFR
  use M_Kinds
  use M_Stamp
  use M_Variables
  use M_StrEdit
  use M_ReadFile
  use M_readpar, only: read_parameters
   use M_ped, only: CPED, init_cped_record
   use M_IO_Ped
   use M_IO_Map
  implicit none

  character(len=MAX_RECL) :: Par_File
  type(PEDHashTable) :: PED_BY_ID, PED_BY_ARN
  type(MAPHashTable) :: MHT_MAP  
  type(SNPInfo), allocatable :: MapInfo(:)
  
  integer :: unitF, nSNP, NREC, i, j, ios, n_vars
  integer :: Col_AniID, Col_SNPID, Col_Allele1, Col_Allele2, MaxArraySize
  integer :: Col_R, Col_GC, Col_GT, Col_Clus
  integer :: out_seq
  integer :: next_progress_pct, current_progress_pct
  real :: Val_R, Val_GC, Val_GT, Val_Clus
  integer :: SNP_Counter, Actual_Map_Idx
  character(len=LEN_STR) :: fr_animal_key_field
  character(len=256) :: version_str
  character(len=MAX_STR) :: out_file
  character(len=MAX_RECL) :: line_buffer
  character(len=MAX_STR) :: base_out_file
  character(len=MAX_STR) :: out_value
  character(len=MAX_STR) :: XC(MAX_VAR)
  character(len=1) :: out_delim
  logical :: out_exists
  
  ! Data storage for genotypes
  integer(kind=ki1), allocatable :: Genotypes(:,:)
  character(len=LEN_STR), allocatable :: FR_AnimalIDs(:)
   type(CPED), allocatable :: FR_AnimalPeds(:)
  integer :: nFR_Animals, current_animal_idx, current_snp_idx
  
  ! Statistics gathering
  integer, allocatable :: Count_Valid_Ani(:)
  integer, allocatable :: Count_Valid_Total(:)
  integer, allocatable :: Count_Valid_Chr(:,:)
  integer, allocatable :: Chr_Max_Idx(:)
  integer, allocatable :: Count_Fail_R(:), Count_Fail_GC(:), Count_Fail_GT(:), Count_Fail_Clus(:)
  integer :: Max_Chr_Num, chr_num
  logical :: qc_fail_r, qc_fail_gc, qc_fail_gt, qc_fail_clus, qc_failed
   logical :: allele1_is_topf, allele2_is_topf, use_atgc_to_ab
  
   character(len=LEN_STR) :: animal_id, snp_id, genotype_str, prev_id, raw_animal_key
   character(len=2) :: genotype_code
  character(len=MAX_RECL) :: raw_file_name
   character(len=MAX_STR) :: date_prefix, fr_basename
   integer :: slash_pos
   type(CPED) :: ani_ped

  version_str = '1.5 - Numericalization for FinalReport file to Compatible GENO Output'
  call version(trim(version_str))
  print '(a)'
  call timestamp()
  print '(a)'

  ! 1. 파라미터 파일 로드
  call getarg(1, Par_File)
  if (len_trim(Par_File) == 0) Par_File = "param"
  
  print *, "Reading parameter file: ", trim(Par_File)
  call read_parameters(Par_File)

   ! QC threshold summary (log once at parameter load time)
   print *
   print *, "=================== FRQC Thresholds (from param) ==================="
   write(*, '(A, F10.4, A, F10.4)') "R_INTENSITY   : ", FRQC%r_intensity(1), " ~ ", FRQC%r_intensity(2)
   write(*, '(A, F10.4)')          "GC_SCORE      : ", FRQC%gc_score
   write(*, '(A, F10.4)')          "GT_SCORE      : ", FRQC%gt_score
   write(*, '(A, F10.4)')          "CLUSTER_SEP   : ", FRQC%cluster_sep
   write(*, '(A, F10.4)')          "ANIMAL_CR     : ", FRQC%animal_cr
   print *, "===================================================================="
  
  ! Identify column positions from FRFILE Section in Par_File
  Col_AniID = 0; Col_SNPID = 0; Col_Allele1 = 0; Col_Allele2 = 0
  Col_R = 0; Col_GC = 0; Col_GT = 0; Col_Clus = 0
  fr_animal_key_field = ''
   allele1_is_topf = .false.
   allele2_is_topf = .false.
  do i = 1, FRFile%NVAR
     select case(trim(FRFile%FieldName(i)))
     case("ANIMAL_ID") ; Col_AniID = FRFile%FieldLoc(i)
     case("ANIMAL_ARN")
        Col_AniID = FRFile%FieldLoc(i)
        fr_animal_key_field = "ANIMAL_ARN"
     case("ARN")
        if (Col_AniID == 0) Col_AniID = FRFile%FieldLoc(i)
        if (len_trim(fr_animal_key_field) == 0) fr_animal_key_field = "ANIMAL_ARN"
     case("ID")
        if (Col_AniID == 0) Col_AniID = FRFile%FieldLoc(i)
        if (len_trim(fr_animal_key_field) == 0) fr_animal_key_field = "ANIMAL_ID"
     case("SNP_ID")    ; Col_SNPID = FRFile%FieldLoc(i)
     case("ALLELE1_AB")
        Col_Allele1 = FRFile%FieldLoc(i)
        allele1_is_topf = .false.
     case("ALLELE2_AB")
        Col_Allele2 = FRFile%FieldLoc(i)
        allele2_is_topf = .false.
     case("ALLELE1_TOP", "ALLELE1_FORWARD")
        Col_Allele1 = FRFile%FieldLoc(i)
        allele1_is_topf = .true.
     case("ALLELE2_TOP", "ALLELE2_FORWARD")
        Col_Allele2 = FRFile%FieldLoc(i)
        allele2_is_topf = .true.
     case("R_INTENSITY"); Col_R = FRFile%FieldLoc(i)
     case("GC_SCORE")   ; Col_GC = FRFile%FieldLoc(i)
     case("GT_SCORE")   ; Col_GT = FRFile%FieldLoc(i)
     case("CLUSTER_SEP"); Col_Clus = FRFile%FieldLoc(i)
     ! ANIMAL_CR is a threshold-only parameter; no column to read from file
     end select

     if (index(trim(to_upper(FRFile%FieldName(i))), "ANIMAL_") == 1 .and. len_trim(fr_animal_key_field) == 0) then
        Col_AniID = FRFile%FieldLoc(i)
        fr_animal_key_field = trim(FRFile%FieldName(i))
     end if
  end do

  ! Default values to prevent 0-index error if not found in param
  if (Col_AniID == 0)   Col_AniID = 1
  if (Col_SNPID == 0)   Col_SNPID = 5
  if (Col_Allele1 == 0) Col_Allele1 = 11
  if (Col_Allele2 == 0) Col_Allele2 = 12
  ! QC columns often at these positions if not specified
  if (Col_R == 0) Col_R = 25
  if (Col_GC == 0) Col_GC = 27
  if (Col_GT == 0) Col_GT = 30
  if (Col_Clus == 0) Col_Clus = 31

  if (len_trim(fr_animal_key_field) == 0) fr_animal_key_field = "ANIMAL_ARN"
  use_atgc_to_ab = allele1_is_topf .and. allele2_is_topf
  print *, "FR key field: ", trim(fr_animal_key_field)
  if (use_atgc_to_ab) then
     print *, "Genotype mode: TOP/FORWARD (A,T->A ; G,C->B)"
  else
     print *, "Genotype mode: AB (direct AA/AB/BB parsing)"
  end if
  
  if (len_trim(FRFile%FileName) == 0) then
     print *, "ERROR: FRFILE not specified."
     stop
  end if

  ! Pointer 연결 (FROUTFile 기반)
  call init_pointers(FROUTFile)

  ! 2. 출력 파일명 생성: GENO_[YYMMDD].txt
  raw_file_name = trim(adjustl(FRFile%FileName))
  slash_pos = scan(raw_file_name, '/', back=.true.)
  if (slash_pos > 0) then
     fr_basename = raw_file_name(slash_pos+1:)
  else
     fr_basename = raw_file_name
  end if

  if (len_trim(fr_basename) >= 6) then
     date_prefix = fr_basename(1:6)
  else
     date_prefix = "FRDATA"
  end if

  base_out_file = "GENO_" // trim(date_prefix)
  out_file = trim(base_out_file) // ".txt"
  inquire(file=trim(out_file), exist=out_exists)
  if (out_exists) then
     out_seq = 1
     do
        write(out_file, '(A, "_", I0, ".txt")') trim(base_out_file), out_seq
        inquire(file=trim(out_file), exist=out_exists)
        if (.not. out_exists) exit
        out_seq = out_seq + 1
     end do
  end if
  print *, "Output file: ", trim(out_file)

  inquire(file=trim(MAPFile%FileName), exist=out_exists)
  if (.not. out_exists) then
     print *, "ERROR: MAP file not found: ", trim(MAPFile%FileName)
     stop 2
  end if

  inquire(file=trim(PEDFile%FileName), exist=out_exists)
  if (.not. out_exists) then
     print *, "ERROR: PED file not found: ", trim(PEDFile%FileName)
     stop 2
  end if

  inquire(file=trim(raw_file_name), exist=out_exists)
  if (.not. out_exists) then
     print *, "ERROR: FR file not found: ", trim(raw_file_name)
     stop 2
  end if

  ! 3. MAP & PED 로드
  call mht_load_map_file(MHT_MAP, MapInfo, nSNP, quiet=.true.)
  print *, "MAP loaded: ", trim(MAPFile%FileName), " (", nSNP, " SNPs)"
  
  ! Determine the maximum array address for memory allocation
  MaxArraySize = 0
  Max_Chr_Num = 0
  do i = 1, nSNP
     if (MapInfo(i)%Array_All > MaxArraySize) MaxArraySize = MapInfo(i)%Array_All
     if (MapInfo(i)%Chr > Max_Chr_Num)        Max_Chr_Num  = MapInfo(i)%Chr
  end do
  if (MaxArraySize == 0) MaxArraySize = nSNP ! Fallback
  print *, "Max Array Address: ", MaxArraySize
  print *, "Max Chromosome:    ", Max_Chr_Num
  
  ! Statistics memory allocation
  allocate(Chr_Max_Idx(Max_Chr_Num))
  Chr_Max_Idx = 0
  do i = 1, nSNP
     chr_num = MapInfo(i)%Chr
     if (chr_num > 0 .and. chr_num <= Max_Chr_Num) then
        if (MapInfo(i)%Array_Chr > Chr_Max_Idx(chr_num)) then
           Chr_Max_Idx(chr_num) = MapInfo(i)%Array_Chr
        end if
     end if
  end do
  
  call pht_load_ped_file(PED_BY_ID, NREC, unitF, "ID", quiet=.true.)
  call pht_load_ped_file(PED_BY_ARN, NREC, unitF, "ARN", quiet=.true.)
  print *, "PED loaded for ID+ARN: ", NREC, " records"
  
  ! 4. FinalReport 1차 스캔: 개체 리스트 확보
  unitF = fopen(raw_file_name)
  if (unitF < 0) then
     print *, "ERROR: Failed to open FR file: ", trim(raw_file_name)
     stop 2
  end if
  do i = 1, FRFile%Header; read(unitF, '(A)'); end do
  
  nFR_Animals = 0
  allocate(FR_AnimalIDs(10000)) 
  prev_id = ""
  do
     read(unitF, '(A)', iostat=ios) line_buffer
     if (ios /= 0) exit
     if (len_trim(line_buffer) == 0) cycle
     call split_string(line_buffer, XC, n_vars, CHAR(9))
     raw_animal_key = adjustl(trim(XC(Col_AniID)))
     if (pht_search(PED_BY_ID, raw_animal_key, ani_ped)) then
        animal_id = trim(ani_ped%ID)
     else if (pht_search(PED_BY_ARN, raw_animal_key, ani_ped)) then
        animal_id = trim(ani_ped%ID)
     else
        animal_id = raw_animal_key
     end if
     
     if (animal_id /= prev_id) then
        nFR_Animals = nFR_Animals + 1
        if (nFR_Animals > size(FR_AnimalIDs)) then
           block
             character(len=LEN_STR), allocatable :: tmp(:)
             allocate(tmp(size(FR_AnimalIDs)*2))
             tmp(1:size(FR_AnimalIDs)) = FR_AnimalIDs
             call move_alloc(tmp, FR_AnimalIDs)
           end block
        end if
        FR_AnimalIDs(nFR_Animals) = animal_id
        prev_id = animal_id
     end if
  end do
  close(unitF)
  print *, "Animals in FR: ", nFR_Animals

  ! 5. 유전자형 데이터 채우기 (AlphaImpute2 요청: Missing=9)
  allocate(Genotypes(MaxArraySize, nFR_Animals))
  allocate(FR_AnimalPeds(nFR_Animals))
  allocate(Count_Valid_Ani(nFR_Animals))
  allocate(Count_Valid_Total(MaxArraySize))
  allocate(Count_Valid_Chr(Max_Chr_Num, nFR_Animals))
  allocate(Count_Fail_R(nFR_Animals))
  allocate(Count_Fail_GC(nFR_Animals))
  allocate(Count_Fail_GT(nFR_Animals))
  allocate(Count_Fail_Clus(nFR_Animals))
  
  Genotypes = 9 
  Count_Valid_Ani = 0
  Count_Valid_Total = 0
  Count_Valid_Chr = 0
  Count_Fail_R    = 0
  Count_Fail_GC   = 0
  Count_Fail_GT   = 0
  Count_Fail_Clus = 0
  
  unitF = fopen(raw_file_name)
  do i = 1, FRFile%Header; read(unitF, '(A)'); end do
  
  current_animal_idx = 0
  prev_id = ""
  SNP_Counter = 0
  do
     read(unitF, '(A)', iostat=ios) line_buffer
     if (ios /= 0) exit
     if (len_trim(line_buffer) == 0) cycle
     call split_string(line_buffer, XC, n_vars, CHAR(9))
     
     raw_animal_key = adjustl(trim(XC(Col_AniID)))
     if (pht_search(PED_BY_ID, raw_animal_key, ani_ped)) then
        animal_id = trim(ani_ped%ID)
     else if (pht_search(PED_BY_ARN, raw_animal_key, ani_ped)) then
        animal_id = trim(ani_ped%ID)
     else
        animal_id = raw_animal_key
     end if
     snp_id    = trim(XC(Col_SNPID))
   call build_genotype_code(XC(Col_Allele1), XC(Col_Allele2), use_atgc_to_ab, genotype_code)
   genotype_str = genotype_code
     
     if (animal_id /= prev_id) then
        current_animal_idx = current_animal_idx + 1
        prev_id = animal_id
        SNP_Counter = 0
        
        ! 패스 간 개체 수 불일치 방지 (경계 밖 쓰기 → 세그폴트 예방)
        if (current_animal_idx > nFR_Animals) then
           print *, "WARNING: More animals in 2nd pass than 1st; extra=", trim(animal_id)
           current_animal_idx = nFR_Animals
           cycle
        end if

        ! 현재 개체의 가계 정보 미리 로드 (XC(2) 참조를 위해)
        if (.not. pht_search(PED_BY_ID, animal_id, ani_ped)) then
           call init_cped_record(ani_ped)
           ani_ped%ID = animal_id
           ani_ped%ARN = raw_animal_key
        end if
        FR_AnimalPeds(current_animal_idx) = ani_ped
     end if
     
     SNP_Counter = SNP_Counter + 1
     
     ! Capture QC metrics; on parse failure keep default (1.0 = passes QC)
     Val_R = 1.0; Val_GC = 1.0; Val_GT = 1.0; Val_Clus = 1.0
     if (Col_R > 0 .and. Col_R <= n_vars)      call safe_read_real(XC(Col_R),    1.0, Val_R)
     if (Col_GC > 0 .and. Col_GC <= n_vars)     call safe_read_real(XC(Col_GC),   1.0, Val_GC)
     if (Col_GT > 0 .and. Col_GT <= n_vars)     call safe_read_real(XC(Col_GT),   1.0, Val_GT)
     if (Col_Clus > 0 .and. Col_Clus <= n_vars) call safe_read_real(XC(Col_Clus), 1.0, Val_Clus)

     ! MapFile의 순서와 FinalReport의 개체내 Record 순서가 동일함을 이용
     ! nSNP(76756) 범위를 벗어나지 않도록 안전장치 추가
     if (SNP_Counter > 0 .and. SNP_Counter <= nSNP) then
        Actual_Map_Idx = SNP_Counter
        i = MapInfo(Actual_Map_Idx)%Array_All
        
        if (i > 0 .and. i <= MaxArraySize .and. current_animal_idx > 0 .and. current_animal_idx <= nFR_Animals) then
           ! Apply per-SNP QC Thresholds (ANIMAL_CR is applied at animal level after all SNPs)
           call check_snp_qc(Val_R, Val_GC, Val_GT, Val_Clus, &
                              qc_fail_r, qc_fail_gc, qc_fail_gt, qc_fail_clus, qc_failed)
           if (qc_fail_r)    Count_Fail_R(current_animal_idx)    = Count_Fail_R(current_animal_idx)    + 1
           if (qc_fail_gc)   Count_Fail_GC(current_animal_idx)   = Count_Fail_GC(current_animal_idx)   + 1
           if (qc_fail_gt)   Count_Fail_GT(current_animal_idx)   = Count_Fail_GT(current_animal_idx)   + 1
           if (qc_fail_clus) Count_Fail_Clus(current_animal_idx) = Count_Fail_Clus(current_animal_idx) + 1

           if (qc_failed) then
              Genotypes(i, current_animal_idx) = 9 ! Failed per-SNP QC
           else
              select case(trim(genotype_str))
                 case("AA") ; Genotypes(i, current_animal_idx) = 0
                 case("AB") ; Genotypes(i, current_animal_idx) = 1
                 case("BB") ; Genotypes(i, current_animal_idx) = 2
                 case default ; Genotypes(i, current_animal_idx) = 9
              end select
              
              ! Statistics: Count valid genotypes (0, 1, 2)
              if (Genotypes(i, current_animal_idx) /= 9) then
                 Count_Valid_Ani(current_animal_idx) = Count_Valid_Ani(current_animal_idx) + 1
                 Count_Valid_Total(i) = Count_Valid_Total(i) + 1
                 chr_num = MapInfo(Actual_Map_Idx)%Chr
                 if (chr_num > 0 .and. chr_num <= Max_Chr_Num) then
                    Count_Valid_Chr(chr_num, current_animal_idx) = Count_Valid_Chr(chr_num, current_animal_idx) + 1
                 end if
              end if
           end if
        end if
     end if
  end do
  close(unitF)

  ! ANIMAL_CR 기준 미달 개체는 출력 시 cycle로 제외 (아래 출력 루프 참조)
  ! Genotypes 배열은 그대로 유지 (통계 log는 per-SNP QC 기준으로 출력)

  ! 6. 동적 필드 출력 (OUTPUT NO_VARIABLES 기준)
  ! PREFIX (FROUTFile%FileName)와 Fieldwidth (FROUTFile%FieldLoc) 사용
  open(unit=40, file=trim(out_file), status='replace')
  out_delim = FindDelim(trim(FROUTFile%Delim_char))

  ! OUTPUT HEADER 값에 따라 헤더 출력 여부 결정
  if (FROUTFile%Header > 0) then
     do j = 1, FROUTFile%NVAR
        select case(trim(FROUTFile%FieldName(j)))
        case("ANIMAL_ID", "ANIMAL_ARN", "ID", "ARN")
           out_value = 'ANIMAL_ID'
           call write_fixed(40, out_value, Fieldwidth(j))
        case("GENO")
           write(40, '(A)', advance='no') "GENO"
        case default
           out_value = trim(FROUTFile%FieldName(j))
           call write_fixed(40, out_value, Fieldwidth(j))
        end select

        if (j < FROUTFile%NVAR) then
           write(40, '(A1)', advance='no') out_delim
        end if
     end do
     write(40, '()')
  end if

   next_progress_pct = 10
  print *, "Processing animals (simple mode):", nFR_Animals, " total"

  do i = 1, nFR_Animals
     ! ANIMAL_CR 기준 미달 개체는 출력에서 완전히 제외
     if (FRQC%animal_cr > 0.0) then
        if (real(Count_Valid_Ani(i)) / real(max(1, nSNP)) < FRQC%animal_cr) cycle
     end if

     ! 해당 개체의 Pedigree 정보 검색 (First Pass에서 저장됨)
     ani_ped = FR_AnimalPeds(i)
     animal_id = FR_AnimalIDs(i)

     current_progress_pct = int(real(i) * 100.0 / real(max(1, nFR_Animals)))
     if (current_progress_pct >= next_progress_pct .or. i == nFR_Animals) then
        write(*, '(A, I0, A, I0, A, I0, A)') "  Progress: ", i, " / ", nFR_Animals, " (", &
           min(100, current_progress_pct), "%)"
        do while (next_progress_pct <= current_progress_pct)
           next_progress_pct = next_progress_pct + 10
        end do
     end if
     
     ! OUTPUT 섹션의 NO_VARIABLES 만큼 필드 출력
     do j = 1, FROUTFile%NVAR
        ! 각 변수명에 따른 데이터 매칭
        select case(trim(FROUTFile%FieldName(j)))
        case("ANIMAL_ID", "ID")
           out_value = trim(ani_ped%ID)
           call write_fixed(40, out_value, Fieldwidth(j))
        case("ANIMAL_ARN", "ARN")
           out_value = trim(ani_ped%ARN)
           call write_fixed(40, out_value, Fieldwidth(j))
        case("BREED")
           out_value = trim(ani_ped%BREED)
           call write_fixed(40, out_value, Fieldwidth(j))
        case("SIRE")
           out_value = trim(ani_ped%SIRE)
           call write_fixed(40, out_value, Fieldwidth(j))
        case("DAM")
           out_value = trim(ani_ped%DAM)
           call write_fixed(40, out_value, Fieldwidth(j))
        case("SEX")
           write(out_value, '(I0)') ani_ped%SEX
           out_value = trim(adjustl(out_value))
           call write_fixed(40, out_value, Fieldwidth(j))
        case("BIRTH", "BDATE", "BIRTH_DATE", "BIRTHDATE")
           write(out_value, '(I0)') ani_ped%BDate
           out_value = trim(adjustl(out_value))
           call write_fixed(40, out_value, Fieldwidth(j))
        case("LOC")
           out_value = trim(ani_ped%LOC)
           call write_fixed(40, out_value, Fieldwidth(j))
        case("GENO")
           ! GENO는 공백 없이 이어서 출력
           do current_snp_idx = 1, MaxArraySize
              write(40, '(I1)', advance='no') Genotypes(current_snp_idx, i)
           end do
        end select

        if (j < FROUTFile%NVAR) then
           write(40, '(A1)', advance='no') out_delim
        end if
     end do
     write(40, '()') ! 행 종료
  end do
  close(40)

  ! 7. QC 통계 보고서: 개체별 조건별 Invalid SNP 수 및 ANIMAL_CR 출력
  ! (중복 포함 - 동일 SNP이 여러 조건에 동시 실패 가능)
  print *
  print *, "================================================================================================================"
  print *, " Per-Animal SNP QC Yield Log  [per-SNP conditions; ANIMAL_CR applied at animal level; overlaps included]  "
  print *, "================================================================================================================"
  write(*, '(A15, 4(1X,A8), 1X,A9, 1X,A10, 1X,A6)') &
     "Animal_ID", "Fail_R", "Fail_GC", "Fail_GT", "Fail_Sep", &
     "Valid_N", "CallRate%", "CR_Pass"
  print *, "-------------------------------------------------------------------"
  do i = 1, nFR_Animals
     block
        real :: call_rate
        character(len=6) :: cr_status
        call_rate = real(Count_Valid_Ani(i)) / real(max(1, nSNP)) * 100.0
        if (real(Count_Valid_Ani(i)) / real(max(1, nSNP)) >= FRQC%animal_cr) then
           cr_status = "PASS"
        else
           cr_status = "FAIL"
        end if
        write(*, '(A15, 4(1X,I8), 1X,I9, 1X,F10.2, 1X,A6)') &
           trim(FR_AnimalIDs(i)), &
           Count_Fail_R(i), Count_Fail_GC(i), Count_Fail_GT(i), Count_Fail_Clus(i), &
           Count_Valid_Ani(i), call_rate, trim(cr_status)
     end block
  end do
  print *, "================================================================================================================"
  print *
  ! Per-animal × Per-chromosome valid SNP table suppressed (verbose log)

  call mht_free(MHT_MAP)
   call pht_free(PED_BY_ID)
   call pht_free(PED_BY_ARN)
  print *, "ReadFR Process Finished Successfully."

contains

   subroutine build_genotype_code(in_allele1, in_allele2, in_use_atgc_to_ab, out_code)
      character(len=*), intent(in) :: in_allele1, in_allele2
      logical, intent(in) :: in_use_atgc_to_ab
      character(len=2), intent(out) :: out_code

      character(len=1) :: a1, a2, m1, m2

      a1 = first_allele_char(in_allele1)
      a2 = first_allele_char(in_allele2)
      m1 = map_allele_to_ab(a1, in_use_atgc_to_ab)
      m2 = map_allele_to_ab(a2, in_use_atgc_to_ab)

      out_code = '??'
      if (m1 == '?' .or. m2 == '?') return

      if (m1 == 'A' .and. m2 == 'A') then
         out_code = 'AA'
      else if (m1 == 'B' .and. m2 == 'B') then
         out_code = 'BB'
      else if ((m1 == 'A' .and. m2 == 'B') .or. (m1 == 'B' .and. m2 == 'A')) then
         out_code = 'AB'
      end if
   end subroutine build_genotype_code

   character(len=1) function first_allele_char(in_text) result(out_char)
      character(len=*), intent(in) :: in_text
      character(len=LEN_STR) :: tmp

      tmp = trim(adjustl(in_text))
      out_char = '?'
      if (len_trim(tmp) <= 0) return

      out_char = to_upper(tmp(1:1))
   end function first_allele_char

   character(len=1) function map_allele_to_ab(in_char, in_use_atgc_to_ab) result(out_char)
      character(len=1), intent(in) :: in_char
      logical, intent(in) :: in_use_atgc_to_ab

      out_char = '?'

      if (in_use_atgc_to_ab) then
         select case(in_char)
         case('A', 'T')
            out_char = 'A'
         case('G', 'C')
            out_char = 'B'
         end select
      else
         select case(in_char)
         case('A')
            out_char = 'A'
         case('B')
            out_char = 'B'
         end select
      end if
   end function map_allele_to_ab

   subroutine write_fixed(unit_no, value, width)
      integer, intent(in) :: unit_no, width
      character(len=*), intent(in) :: value
      integer :: w, n
      character(len=MAX_STR) :: tmp

      tmp = adjustl(value)
      w = max(0, width)
      n = len_trim(tmp)

      if (w == 0) then
         write(unit_no, '(A)', advance='no') trim(tmp)
      else if (n >= w) then
         write(unit_no, '(A)', advance='no') tmp(1:w)
      else
         write(unit_no, '(A)', advance='no') trim(tmp) // repeat(' ', w - n)
      end if
   end subroutine write_fixed

   ! Per-SNP QC condition checker (ANIMAL_CR excluded; handled at animal level)
   ! 각 조건의 실패 여부를 개별 반환 (중복 포함: 동일 SNP이 여러 조건 동시 실패 가능)
   subroutine check_snp_qc(v_r, v_gc, v_gt, v_clus, &
                            fail_r, fail_gc, fail_gt, fail_clus, fail_any)
      real, intent(in)  :: v_r, v_gc, v_gt, v_clus
      logical, intent(out) :: fail_r, fail_gc, fail_gt, fail_clus, fail_any

      fail_r    = v_r < FRQC%r_intensity(1) .or. v_r > FRQC%r_intensity(2)
      fail_gc   = v_gc < FRQC%gc_score
      fail_gt   = FRQC%gt_score    > 0.0 .and. v_gt   < FRQC%gt_score
      fail_clus = FRQC%cluster_sep > 0.0 .and. v_clus < FRQC%cluster_sep
      fail_any  = fail_r .or. fail_gc .or. fail_gt .or. fail_clus
   end subroutine check_snp_qc

   subroutine safe_read_real(in_text, default_val, out_val)
      character(len=*), intent(in) :: in_text
      real, intent(in) :: default_val
      real, intent(out) :: out_val
      integer :: ios_local
      real :: tmp

      out_val = default_val
      if (len_trim(in_text) <= 0) return

      tmp = default_val
      read(in_text, *, iostat=ios_local) tmp
      if (ios_local == 0) out_val = tmp
   end subroutine safe_read_real

end program ReadFR
