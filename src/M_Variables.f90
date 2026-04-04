module M_Variables
use M_Kinds
implicit none   
! define VARIABLE for PROGRAMS
INTEGER,PARAMETER,PRIVATE:: lock_y=2028, lock_m=12, lock_d=31
INTEGER, PARAMETER, PUBLIC:: missI = 0
INTEGER, PARAMETER, PUBLIC:: ZERO=missI
INTEGER, PARAMETER, PUBLIC:: LEN_KEY = 32
INTEGER, PARAMETER, PUBLIC:: LEN_STR = LEN_KEY
INTEGER, PARAMETER, PUBLIC:: MAX_STR = 128
INTEGER, PARAMETER, PUBLIC:: MAX_VAR = MAX_STR
INTEGER, PARAMETER, PUBLIC:: MAX_RECL = 512
REAL, PARAMETER, PUBLIC:: missR = 0.d0
REAL,PUBLIC:: STARTIME, FINISHTIME
integer, parameter, PUBLIC:: UNKNOWN=0,MALE=2,FEMALE=1
CHARACTER(LEN=MAX_RECL), PUBLIC :: RECORD
CHARACTER(len=1), PARAMETER, PUBLIC:: missC = achar(32)  ! " "
CHARACTER(len=1), PARAMETER, PUBLIC:: nullC = missC      ! " "
CHARACTER(len=1), PARAMETER, PUBLIC:: SPACE = missC      ! " "
CHARACTER(len=1), PARAMETER, PUBLIC:: missA = achar(48)  ! "0"
CHARACTER(len=1), PARAMETER, PUBLIC:: TAB = achar(9)
LOGICAL, PARAMETER, PUBLIC :: TRUE = .true., FALSE = .false.
LOGICAL, PARAMETER, PUBLIC ::  Desc=.true., Ordered=.true.
  
type, PUBLIC :: FileInfo
    character(len=MAX_RECL) :: FileName
    character(len=LEN_STR) :: FieldName(MAX_VAR)
    integer :: FieldLoc(MAX_VAR)
    character(len=LEN_STR) :: Delim_char
    integer :: Header
    integer :: NVAR
end type FileInfo

type, PUBLIC :: SNPInfo
    character(len=LEN_STR) :: SNP_ID
    integer :: Chr
    integer :: Pos
    integer :: Array_All
    integer :: Array_Chr
end type SNPInfo

! =========================================================================
! QC 데이터 구조체 (SNP 품질 관리용)
! =========================================================================

type, PUBLIC :: FRQC_matrices
   character(len=LEN_STR) :: SNP_ID ! SNP identifier 
   integer :: NVAR
   real :: gc_score                 ! GC Score value
   real :: r_intensity(2)           ! R - Intensity value
   real :: gt_score                 ! GT Score value
   real :: cluster_sep              ! Cluster Separation value
   real :: animal_cr                ! Animal call rate
   logical :: qc_passed            ! QC pass/fail flag
end type FRQC_matrices

type, PUBLIC :: POPQC_matrices
   integer :: NVAR
   real :: snp_cr                       ! SNP call rate
   real :: animal_cr                    ! Animal call rate
   real :: maf                          ! Minor Allele Frequency
   character(len=LEN_STR) :: HWE_method ! HWE method (e.g., "FDR", "Bonferroni")
   real :: HWE_threshold                ! FDR: 0.05, Bonferroni: 0.05/NVAR
    integer :: HWE_min_n                 ! Minimum non-missing N required for HWE test
    integer :: HWE_min_mac               ! Minimum minor-allele count required for HWE test
    logical :: sex_check                 ! legacy sex check flag (compat)
    real :: sex_concordence_threshold    ! female F threshold for sex concordance
    real :: mendel_snp_threshold         ! Mendelian inconsistency threshold (per SNP)
    real :: mendel_viol_threshold        ! Mendelian inconsistency threshold (per animal)
    logical :: parent_seek               ! legacy parent seek flag (compat)
    real :: parent_seek_threshold        ! candidate parent selection threshold
    integer :: parent_seek_nsnps         ! target SNP count for parent-seek panel
    real :: parent_seek_g_tol            ! G relationship tolerance for parent-seek
    real :: id_merge_sample_fraction     ! sampling fraction for ID_MERGE contexts
    real :: id_merge_confirm_match_threshold ! full-autosome confirm threshold for ID_MERGE
    real :: id_swap_compat_min           ! Minimum Mendel compatibility for ID swap detection
    integer :: id_swap_max_depth          ! Maximum recursion depth for ID swap chain
    integer :: multi_resolved_action      ! 0=update PED, 1=delete record
end type POPQC_matrices

! =========================================================================
! QC Threshold 설정 (Parameter 파일에서 읽어옴)
! =========================================================================

! Module variables declaration with PUBLIC attribute
type(FileInfo), public :: PEDFile, DATAFile, MAPFile, FRFile, GENOFile, SNPFile
type(FileInfo), target, save :: FROUTFile, GENOUTFile
character(len=MAX_RECL), pointer, save :: PREFIX => null()
integer, pointer, save :: Fieldwidth(:) => null()  ! Pointer alias for FieldLoc in output context
type(FRQC_matrices), public :: FRQC
type(POPQC_matrices), public :: POPQC

! Public subroutines
public :: init_File, init_FRQC, init_POPQC

contains

subroutine init_pointers(OUTFile)
    type(FileInfo), intent(in), target :: OUTFile 
    
    nullify(PREFIX)
    nullify(Fieldwidth)
    
    PREFIX => OUTFile%FileName
    Fieldwidth => OUTFile%FieldLoc 
end subroutine init_pointers

subroutine init_File(File)
    type(FileInfo), intent(out) :: File
    File%FileName = ''
    File%FieldName = ''
    File%FieldLoc = 0
    File%Delim_char = ''
    File%Header = 0
    File%NVAR = 0
end subroutine init_File

subroutine init_FRQC(FR_QC)
   type(FRQC_matrices), intent(out) :: FR_QC
   FR_QC%SNP_ID = ''
   FR_QC%NVAR = 0
   FR_QC%gc_score = 0.15
   FR_QC%r_intensity(1) = 0.20
   FR_QC%r_intensity(2) = 2.00
   FR_QC%gt_score = 0.15
   FR_QC%cluster_sep = 0.30
   FR_QC%animal_cr = 0.70
   FR_QC%qc_passed = .false.
end subroutine init_FRQC

subroutine init_POPQC(POP_QC)
   type(POPQC_matrices), intent(out) :: POP_QC
   POP_QC%NVAR = 0
   POP_QC%maf = 0.01
   POP_QC%snp_cr = 0.70
   POP_QC%animal_cr = 0.70
   POP_QC%HWE_method = 'BONFERRONI'
   POP_QC%HWE_threshold = 0.05
    POP_QC%HWE_min_n = 20
    POP_QC%HWE_min_mac = 5
    POP_QC%sex_check = .true.
    POP_QC%sex_concordence_threshold = 0.02
    POP_QC%mendel_snp_threshold = 0.01
    POP_QC%mendel_viol_threshold = 0.01
    POP_QC%parent_seek = .true.
    POP_QC%parent_seek_threshold = 0.95
    POP_QC%parent_seek_nsnps = 500
    POP_QC%parent_seek_g_tol = 0.15
    POP_QC%id_merge_sample_fraction = 0.20
    POP_QC%id_merge_confirm_match_threshold = 0.995
    POP_QC%id_swap_compat_min = 0.95
    POP_QC%id_swap_max_depth = 3
    POP_QC%multi_resolved_action = 0
end subroutine init_POPQC

end module M_Variables

