# 빠른 시작: 새 프로그램 추가하기 (Quick Start)

## 30초 가이드

### 1️⃣ 프로그램 템플릿 생성
```bash
cd /home/dhlee/GPBLUP
./create_program.sh MyProgram "My program description"
```

**생성되는 파일:**
```
MyProgram/
├── MyProgram.f90                    ← 메인 소스코드 (구현하기)
├── README.md                        ← 문서
└── test/
    └── parameter_example.txt        ← 파라미터 템플릿
```

---

### 2️⃣ 프로그램 로직 구현
```bash
vim MyProgram/MyProgram.f90
```

`main_process()` 서브루틴에 당신의 로직 작성:

```fortran
subroutine main_process()
  print *, "Starting main process..."
  
  ! TODO: 당신의 코드를 여기에 작성하세요
  ! M_Variables의 전역 변수 사용 가능:
  !   - PEDFile, MAPFile, SNPFile, GENOFile
  !   - OutputPrefix, QCThresholds
  
  print *, "✓ Main process completed"
end subroutine main_process
```

---

### 3️⃣ CMakeLists.txt에 등록

**파일:** `/home/dhlee/GPBLUP/CMakeLists.txt`

마지막 `add_dependencies` 항목 다음에 다음 코드 추가:

```cmake
# ============================================================================
# MyProgram
# ============================================================================
message(STATUS "Configuring MyProgram program...")

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

# MyProgram depends on all modules
add_dependencies(MyProgram gpblup_modules)

# Install to bin directory
install(TARGETS MyProgram
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)
```

> ⚠️ **중요:** 변수명은 대문자로, 프로그램명은 원래 케이스로
> ```
> set(MYPROGRAM_SOURCES ...)      ← 대문자
> add_executable(MyProgram ...)   ← 원래 케이스
> ```

---

### 4️⃣ 빌드
```bash
make rebuild
# 또는
make             # 기존 프로그램 유지하면서 새로운 것만 빌드
```

**결과:**
```
✅ MyProgram (새로 생성됨)
✅ ReadFR (89K)
✅ popQC (124K)
```

---

### 5️⃣ 확인 및 테스트
```bash
# 바이너리 확인
ls -lh bin/MyProgram

# 실행
./bin/MyProgram MyProgram/test/parameter_example.txt
```

---

## 🔄 사이클: 프로그램 수정 후 재빌드

프로그램 파일 수정 → 자동 재빌드:

```bash
# 1. 파일 수정
vim MyProgram/MyProgram.f90

# 2. 재빌드 (MyProgram만 재컴파일)
make

# 3. 완료 ✓
```

---

## 🔗 공통 모듈 수정 후

src/ 파일 수정 → **모든 프로그램 자동 재컴파일:**

```bash
# 1. 공통 모듈 수정
vim src/M_Variables.f90

# 2. 재빌드 하면 모든 프로그램이 자동 재컴파일
make

# 결과: ReadFR, popQC, MyProgram, ... 모두 업데이트
```

---

## 📝 파라미터 파일 형식

`MyProgram/test/parameter_example.txt`:

```bash
# 필수 입력 파일
PEDFILE: pedigree.txt
  HEADER: 1
  DELIM: SPACE
  NO_VARIABLES: 8
  1 ANIMAL_ID, 2 BREED, 3 SIRE, 4 DAM, 5 SEX, 6 BDATE, 7 LOC, 8 ARN

MAPFILE: map.txt
  HEADER: 0
  DELIM: SPACE
  NO_VARIABLES: 4
  1 SNP_ID, 2 CHR, 3 POS, 4 ALLELE

# 출력 디렉토리
OUTPUTPREFIX: output/MyProgram

# 프로그램 특정 파라미터
MY_PARAM: value
```

---

## 💡 팁

### 파일 필드 접근 예

```fortran
! M_readpar로부터 파일 정보 읽음
! PEDFile%FieldName(:)  - 필드 이름 배열
! PEDFile%FieldLoc(:)   - 필드 위치 배열
! PEDFile%Delim_char    - 구분자
! PEDFile%Header        - 헤더 라인 수
! PEDFile%NVAR          - 변수 개수

! 예: ANIMAL_ID 필드 위치 찾기
integer :: i, id_idx
do i = 1, PEDFile%NVAR
  if (trim(adjustl(PEDFile%FieldName(i))) == "ANIMAL_ID") then
    id_idx = i
    exit
  end if
end do
```

### 파일I/O 패턴

```fortran
integer :: iounit, iostat

! 파일 열기
open(unit=10, file=trim(filename), status='old', action='read', iostat=iostat)
if (iostat /= 0) then
  print *, "ERROR: Cannot open file"
  stop 1
end if

! 헤더 스킵
do i = 1, PEDFile%Header
  read(10, *)
end do

! 데이터 읽기
do
  read(10, *, iostat=iostat) data_variables
  if (iostat /= 0) exit
  ! 처리
end do

close(10)
```

---

## ✅ 체크리스트

새 프로그램 완성 전 확인:

```
□ MyProgram/MyProgram.f90 작성
□ README.md 업데이트
□ CMakeLists.txt 등록
□ make rebuild 성공
□ bin/MyProgram 생성됨
□ bin/MyProgram test/parameter_example.txt 실행 성공
□ 파라미터 파일 테스트 완료
```

---

## 🚨 문제 해결

### "add_executable" 에러
```
CMake Error: No CMAKE_Fortran_COMPILER could be found.
```
→ gfortran 설치: `sudo apt-get install gfortran`

### 컴파일 에러
```bash
rm -rf build/      # CMake 캐시 초기화
make rebuild       # 재빌드
```

### 파일을 찾을 수 없음
```
ERROR: File not found
```
→ 파라미터 파일의 경로 확인 (상대 경로 vs 절대 경로)

---

## 📚 더 자세한 정보

- [PROGRAM_DEVELOPMENT_GUIDE.md](PROGRAM_DEVELOPMENT_GUIDE.md) - 전체 개발 가이드
- [BUILD_GUIDE.md](BUILD_GUIDE.md) - 빌드 시스템 가이드
- ReadFR/README.md - ReadFR 참고 예제
- popQC/README.md - popQC 참고 예제

---

## 예: Hello World 프로그램

### 1. 생성
```bash
./create_program.sh HelloWorld "Simple hello world program"
```

### 2. 구현 (HelloWorld/HelloWorld.f90)
```fortran
subroutine main_process()
  print *, "Hello from HelloWorld!"
  print *, "PED file: ", trim(PEDFile%FileName)
end subroutine main_process
```

### 3. CMakeLists.txt 추가
```cmake
add_executable(HelloWorld HelloWorld/HelloWorld.f90)
target_link_libraries(HelloWorld PRIVATE gpblup_modules OpenMP::OpenMP_Fortran)
set_target_properties(HelloWorld PROPERTIES OUTPUT_NAME HelloWorld RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
add_dependencies(HelloWorld gpblup_modules)
```

### 4. 빌드
```bash
make rebuild
```

### 5. 실행
```bash
./bin/HelloWorld HelloWorld/test/parameter_example.txt
```

---

**완료! 이제 통합 GPBLUP 프로젝트의 일부가 되었습니다.**

프로그램을 계속 추가하면서 모든 프로그램이 최신 공통 모듈을 사용하도록 자동으로 동기화됩니다! 🎉
