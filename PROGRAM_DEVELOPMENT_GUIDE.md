# 새로운 프로그램 개발 & 통합 가이드

## 개요

GPBLUP 프로젝트에 새로운 단위 프로그램을 개발하고 통합하는 완전한 워크플로우를 설명합니다.

---

## 📋 개발 체크리스트

프로그램 개발 시 다음 순서를 따릅니다:

- [ ] 1. 프로그램 템플릿 생성
- [ ] 2. 프로그램 로직 구현
- [ ] 3. CMakeLists.txt에 등록
- [ ] 4. 빌드 및 테스트
- [ ] 5. 문서 작성

---

## 🚀 Step 1: 프로그램 템플릿 생성

### 자동 생성 (권장)

```bash
cd /home/dhlee/GPBLUP
./create_program.sh MyProgram "My program description"
```

결과:
```
MyProgram/
├── MyProgram.f90                    ← 메인 소스코드
├── README.md                        ← 프로그램 문서
└── test/
    └── parameter_example.txt        ← 파라미터 템플릿
```

### 수동 생성

프로그램 디렉토리 구조:
```
MyProgram/
├── MyProgram.f90              ← 메인 프로그램
├── README.md                  ← 프로젝트 설명
└── test/                      ← 테스트 데이터
    ├── parameter_example.txt  ← 파라미터 파일
    └── test_data/             ← 테스트 입력 파일 (선택)
```

---

## 🔧 Step 2: 프로그램 로직 구현

### 기본 구조

생성된 템플릿의 기본 구조:

```fortran
program MyProgram
  use M_Kinds
  use M_Variables
  use M_readpar
  implicit none
  
  character(len=256) :: param_file
  
  ! 1. 사용자 입력 (파라미터 파일) 받기
  call get_command_argument(1, param_file)
  
  ! 2. 파라미터 파일 파싱 (M_readpar 사용)
  call setup_parameters(param_file)
  
  ! 3. 메인 처리
  call main_process()
  
  ! 4. 결과 출력
  call generate_output()
  
contains
  ! subroutine 정의
end program MyProgram
```

### 주요 전역 변수 (M_Variables에서)

```fortran
! 파일 정보
type(FileInfo), public :: PEDFile      ! 혈통 파일
type(FileInfo), public :: MAPFile      ! 마커 맵 파일
type(FileInfo), public :: SNPFile      ! SNP 파일
type(FileInfo), public :: GENOFile     ! 유전형 파일
type(FileInfo), public :: DATAFile     ! 데이터 파일

! 분석 파라미터
type(QC_Thresholds), public :: QCThresholds
character(len=MAX_STR), public :: OutputPrefix

! FileInfo 구조체 필드
type :: FileInfo
    character(len=MAX_STR) :: FileName
    character(len=MAX_STR) :: FieldName(MAX_VAR)    ! 필드 이름
    integer :: FieldLoc(MAX_VAR)                    ! 필드 위치
    character(len=LEN_STR) :: Delim_char           ! 구분자
    integer :: Header                              ! 헤더라인 수
    integer :: NVAR                                ! 변수 개수
end type FileInfo
```

### 예: PED 파일 읽기

```fortran
subroutine main_process()
  integer :: iounit, iostat, n_animals
  character(len=MAX_RECL) :: line
  character(len=MAX_STR) :: fields(MAX_VAR)
  integer :: n_fields
  
  ! PED 파일이 정의되어 있는지 확인
  if (len_trim(PEDFile%FileName) == 0) then
    print *, "ERROR: PED file not specified"
    stop
  end if
  
  print *, "Reading PED file: ", trim(PEDFile%FileName)
  
  ! 파일 열기
  open(unit=10, file=trim(PEDFile%FileName), status='old', action='read')
  
  ! 헤더 라인 스킵
  do i = 1, PEDFile%Header
    read(10, *)
  end do
  
  ! 데이터 읽기
  n_animals = 0
  do
    read(10, '(A)', iostat=iostat) line
    if (iostat /= 0) exit
    
    ! 라인을 필드로 분리 (구분자 사용)
    call split_string_by_delimiter(line, PEDFile%Delim_char, &
                                   fields, n_fields)
    
    ! 필드 처리 예)
    ! ANIMAL_ID_idx = find_field_index(PEDFile, "ANIMAL_ID")
    ! print *, "Animal: ", trim(fields(ANIMAL_ID_idx))
    
    n_animals = n_animals + 1
  end do
  
  close(unit=10)
  
  print *, "Read ", n_animals, " animals"
  
end subroutine main_process

subroutine split_string_by_delimiter(str, delim, fields, n_fields)
  character(len=*), intent(in) :: str, delim
  character(len=*), intent(out) :: fields(:)
  integer, intent(out) :: n_fields
  
  ! TODO: Implement string splitting
  ! (M_StrEdit 모듈 사용 가능)
end subroutine split_string_by_delimiter
```

---

## 📝 Step 3: CMakeLists.txt에 등록

### CMakeLists.txt 수정 (마지막에 추가)

```cmake
# ==========================================================================
# MyProgram
# ==========================================================================
message(STATUS "Configuring MyProgram...")

set(MYPROGRAM_SOURCES
    MyProgram/MyProgram.f90
)

add_executable(MyProgram ${MYPROGRAM_SOURCES})
target_link_libraries(MyProgram PRIVATE gpblup_modules OpenMP::OpenMP_Fortran)
set_target_properties(MyProgram PROPERTIES
    OUTPUT_NAME MyProgram
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# MyProgram depends on all common modules
add_dependencies(MyProgram gpblup_modules)

# Install to bin directory
install(TARGETS MyProgram
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)
```

### CMakeLists.txt 최종 구조

```cmake
# .... 기존 코드 ....

# ReadFR 설정
add_executable(ReadFR ReadFR/ReadFR.f90)
# ...

# popQC 설정
add_executable(popQC popQC/popQC.f90)
# ...

# 👇 새 프로그램들은 여기에 추가
# MyProgram 설정
add_executable(MyProgram MyProgram/MyProgram.f90)
# ...

# Program2 설정
add_executable(Program2 Program2/Program2.f90)
# ...
```

---

## 🔨 Step 4: 빌드 및 테스트

### 빌드

```bash
cd /home/dhlee/GPBLUP

# 전체 재빌드 (새 프로그램 포함)
make rebuild

# 또는
cmake -B build -DCMAKE_BUILD_TYPE=Release
cd build
make -j$(nproc)
```

### 빌드 확인

```bash
# 바이너리 생성 확인
ls -lh bin/MyProgram

# 또는 빌드 상태 확인
make status
```

### 테스트

```bash
# 프로그램 실행
./bin/MyProgram MyProgram/test/parameter_example.txt

# 결과 확인
ls -la output/
```

---

## 📚 Step 5: 문서 작성

### MyProgram/README.md 작성

```markdown
# MyProgram

## 설명
프로그램의 목적과 기능 설명

## 입력 파일 형식
- PED 파일: 형식 설명
- MAP 파일: 형식 설명
- ...

## 파라미터
파라미터 파일의 의미 설명

## 출력
생성되는 출력 파일 설명

## 사용 예
./bin/MyProgram parameter.txt
```

---

## 🔄 워크플로우 요약

### 빠른 개발 사이클

```bash
# 1️⃣ 프로젝트 초기화
./ create_program.sh NewProgram "My description"

# 2️⃣ 소스코드 구현
vim NewProgram/NewProgram.f90

# 3️⃣ CMakeLists.txt 수정 (NewProgram 추가)
vim CMakeLists.txt

# 4️⃣ 빌드
make rebuild

# 5️⃣ 테스트
./bin/NewProgram NewProgram/test/parameter_example.txt

# 6️⃣ 반복 (파일 수정 후)
make              # 자동으로 NewProgram만 재컴파일
```

### 공통 모듈 수정 시

```bash
# M_Variables.f90 수정
vim source/M_Variables.f90

# 빌드 (모든 프로그램 자동 재컴파일)
make

# 결과: ReadFR, popQC, MyProgram, ... 모두 자동 업데이트
```

---

## 💡 모범 사례

### 1. 파라미터 파일 설계

```bash
# 좋은 예
PEDFILE: data/pedigree.txt
  HEADER: 1
  DELIM: SPACE
  NO_VARIABLES: 8
  1 ID, 2 BREED, 3 SIRE, 4 DAM, ...

# 프로그램별 파라미터 섹션
MYPROGRAM:
  THRESHOLD: 0.95
  NUM_ITERATIONS: 1000
```

### 2. 에러 처리

```fortran
subroutine main_process()
  integer :: iostat
  logical :: file_exists
  
  ! 파일 존재 확인
  inquire(file=trim(PEDFile%FileName), exist=file_exists)
  if (.not. file_exists) then
    print *, "ERROR: File not found"
    stop 1
  end if
  
  ! I/O 에러 확인
  open(unit=10, file=trim(PEDFile%FileName), &
       status='old', action='read', iostat=iostat)
  if (iostat /= 0) then
    print *, "ERROR: Cannot open file, code:", iostat
    stop 1
  end if
end subroutine main_process
```

### 3. 로깅/출력

```fortran
print *, ""
print *, "========== MyProgram =========="
print *, ""
print *, "Input files:"
print *, "  PED: ", trim(PEDFile%FileName)
print *, "  MAP: ", trim(MAPFile%FileName)
print *, ""
print *, "Analysis parameters:"
print *, "  Threshold: ", my_threshold
print *, ""
print *, "Processing..."
print *, ""
print *, "✓ Completed successfully"
print *, ""
```

---

## 🎯 프로젝트 통합 체크리스트

프로그램을 완성하기 전에 확인할 사항:

- [ ] 프로그램이 M_readpar로부터 파라미터 읽음
- [ ] 공통 모듈들 (M_Kinds, M_Variables 등) 사용
- [ ] 파라미터 파일 템플릿 포함
- [ ] README.md 작성
- [ ] CMakeLists.txt에 등록
- [ ] 에러 처리 구현
- [ ] 테스트 데이터 준비
- [ ] `make` 명령으로 성공적으로 빌드됨
- [ ] bin/ 에 실행파일 생성됨

---

## 📦 배포 준비

프로그램 완성 후:

```bash
# 전체 프로젝트 확인
make status

# 모든 프로그램이 표시되는지 확인
ls -lh bin/

# git에 추가 (Optional)
cd /home/dhlee/GPBLUP
git add NewProgram/
git add CMakeLists.txt
git commit -m "Add NewProgram"
```

---

## 🤝 새 프로그램 추가 예시

### 예: 유전형 필터링 프로그램 ("GenoFilter")

```bash
# 1. 템플릿 생성
./create_program.sh GenoFilter "Genotype filtering program"

# 2. 소스코드 구현
vim GenoFilter/GenoFilter.f90

# 3. CMakeLists.txt 수정
# (GenoFilter 추가)
vim CMakeLists.txt

# 4. 빌드
make rebuild

# 5. 테스트
./bin/GenoFilter GenoFilter/test/parameter_example.txt

# 결과
# ✓ ReadFR (89K)
# ✓ popQC (124K)
# ✓ GenoFilter (NEW)
# ...
```

---

## ✅ 완료

이제 다음과 같은 장점을 갖게 됩니다:

- ✅ 일관된 프로그램 구조
- ✅ 자동 dependency 추적
- ✅ 공통 모듈 공유
- ✅ 간편한 빌드 시스템
- ✅ 확장 가능한 아키텍처
- ✅ 통합된 파라미터 시스템

새로운 프로그램을 계속 추가해도 기존 프로그램에 영향 없고, 공통 모듈 수정 시 모든 프로그램이 자동으로 동기화됩니다! 🎉
