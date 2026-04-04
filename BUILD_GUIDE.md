# GPBLUP 빌드 시스템 가이드

## 개요

GPBLUP 프로젝트는 CMake 기반의 자동 빌드 시스템을 사용합니다. 이 시스템은 공통 모듈의 변경 시 자동으로 dependency가 있는 실행파일들을 재컴파일합니다.

## 프로젝트 구조

```
GPBLUP/
├── src/                  # 공통 모듈 (모든 프로그램에서 공유)
│   ├── M_Kinds.f90
│   ├── M_Variables.f90
│   ├── M_readpar.f90        # 파라미터 파일 파서
│   ├── M_ReadFile.f90
│   ├── M_StrEdit.f90
│   ├── M_Stamp.f90
│   ├── M_HashTable.f90
│   ├── M_PEDHashTable.f90
│   ├── M_param.f90
│   └── Qsort4.f90
│
├── popQC/                   # popQC 프로그램
│   └── popQC.f90
│
├── ReadFR/                  # ReadFR 프로그램
│   └── ReadFR.f90
│
├── bin/                     # 실행파일 저장 위치 (자동 생성)
│   ├── ReadFR              # ReadFR 실행파일
│   └── popQC               # popQC 실행파일
│
├── build/                   # CMake 빌드 디렉토리 (자동 생성)
│   ├── modules/            # .mod 파일 저장 위치
│   └── libdkblupf90.a       # 정적 라이브러리
│
├── CMakeLists.txt          # CMake 설정 파일
├── Makefile                # 루트 Makefile (간단한 인터페이스)
└── build.sh                # 빌드 스크립트 (선택적)
```

## 빌드 방법

### 1. 간단한 방법 - Makefile 사용

#### One-batch 설치 (권장)
```bash
cd /home/dhlee/GPBLUP
./install.sh
# 또는
./install.sh /opt/gpblup
```
ReadFR, popQC, relped, relgeno를 한번에 빌드하여 PREFIX/bin에 설치합니다.

#### Release 빌드 (최적화, 권장)
```bash
cd /home/dhlee/GPBLUP
make              # 또는 make release
```

#### Debug 빌드 (디버깅 심볼 포함)
```bash
make debug
```

#### 전체 재빌드
```bash
make rebuild      # clean + rebuild
```

#### 빌드 결과 정리
```bash
make clean        # 모든 빌드 결과 삭제
```

#### 현재 빌드 상태 확인
```bash
make status
```

### 2. 스크립트 방법 - build.sh 사용

```bash
cd /home/dhlee/GPBLUP
./build.sh                  # Release 빌드 (기본값)
./build.sh release          # Release 빌드
./build.sh debug            # Debug 빌드
./build.sh rebuild          # Clean + 빌드
./build.sh clean            # 모든 빌드 결과 정리
```

### 3. CMake 직접 사용

```bash
cd /home/dhlee/GPBLUP
mkdir -p build
cd build

# Release 빌드
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# 또는 Debug 빌드
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
```

## 자동 Dependency 처리

### 핵심 메커니즘

1. **모든 공통 모듈은 정적 라이브러리로 컴파일됨** (`libgpblup_modules.a`)
   - src/ 디렉토리의 모든 .f90 파일이 컴파일됨
   - 모듈들간의 dependency 자동 추적

2. **각 프로그램은 정적 라이브러리에 링크됨**
   - `ReadFR` ← libgpblup_modules.a
   - `popQC` ← libgpblup_modules.a

3. **CMake가 dependency 자동 감지**
   - src/ 파일 수정 → libgpblup_modules.a 재컴파일
   - libgpblup_modules.a 재컴파일 → ReadFR, popQC 자동 재링크

### 예시: 모듈 수정 후 빌드

```bash
# M_Variables.f90 수정
vim /home/dhlee/GPBLUP/src/M_Variables.f90

# 다시 빌드하면 자동으로:
make
# 1. M_Variables.f90 재컴파일
# 2. libgpblup_modules.a 재컴파일
# 3. ReadFR 재링크
# 4. popQC 재링크
# 5. bin/ReadFR, bin/popQC 자동 업데이트
```

## 새 프로그램 추가 방법

혹시 새로운 프로그램을 추가하려면 CMakeLists.txt를 수정하면 됩니다:

### 예: `NewProgram` 추가

**CMakeLists.txt에 추가:**
```cmake
# ============================================================================
# NewProgram Program
# ============================================================================
message(STATUS "Configuring NewProgram program...")

set(NEWPROGRAM_SOURCES
    NewProgram/NewProgram.f90
)

add_executable(NewProgram ${NEWPROGRAM_SOURCES})
target_link_libraries(NewProgram PRIVATE gpblup_modules OpenMP::OpenMP_Fortran)
set_target_properties(NewProgram PROPERTIES
    OUTPUT_NAME NewProgram
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# NewProgram depends on all modules
add_dependencies(NewProgram gpblup_modules)
```

그 다음 재빌드하면 자동으로 `bin/NewProgram`이 생성됩니다.

## 빌드 성능 최적화

### 병렬 컴파일 활용
```bash
# 시스템의 모든 CPU 코어 사용
make -j$(nproc)

# 또는 특정 개수 지정
make -j4           # 4개 코어 사용
```

### 컴파일 플래그

**Release 빌드** (기본값):
- `-O2` : 최적화 레벨 2
- `-march=native` : CPU별 최적화

**Debug 빌드**:
- `-g` : 디버깅 심볼 포함
- `-O0` : 최적화 비활성화
- `-fcheck=all` : 런타임 체크 활성화
- `-fbacktrace` : 백트레이스 지원

## 주요 특징

✅ **자동 Dependency 추적**
- CMake가 파일 변경 감지
- 필요한 부분만 재컴파일

✅ **병렬 컴파일**
- `make -j$(nproc)`로 모든 CPU 코어 활용
- 빌드 시간 단축

✅ **모든 프로그램이 동일한 모듈 사용**
- readpar 모듈 수정 시 ReadFR과 popQC 모두 자동 재컴파일
- 호환성 보장

✅ **깔끔한 빌드 디렉토리 분리**
- src/ : 소스 코드
- build/ : 빌드 아티팩트 (.o, .mod)
- bin/ : 최종 실행파일

## 문제 해결

### 빌드가 실패하는 경우

1. **완전히 정리한 후 재빌드**
   ```bash
   make clean
   make rebuild
   ```

2. **CMake 캐시 초기화**
   ```bash
   rm -rf build/
   make
   ```

3. **Fortran 컴파일러 버전 확인**
   ```bash
   gfortran --version
   ```

### 특정 모듈만 재컴파일하려면

```bash
# build 디렉토리에서 직접 make 실행
cd /home/dhlee/GPBLUP/build
make gpblup_modules -j$(nproc)  # 모듈만 재컴파일
make ReadFR -j$(nproc)          # ReadFR만 재컴파일
make popQC -j$(nproc)           # popQC만 재컴파일
```

## 빌드 결과 확인

```bash
# 바이너리 크기 확인
ls -lh /home/dhlee/GPBLUP/bin/

# 최근 수정 시간 확인
ls -lh /home/dhlee/GPBLUP/bin/ | grep Feb
```

## 제거 및 정리

모든 빌드 결과를 제거하고 소스만 남기려면:

```bash
make clean          # bin/과 build/ 정리
# 또는
rm -rf build/
rm -f bin/ReadFR bin/popQC
```

---

**최종 구성:**
- ✅ CMakeLists.txt : CMake 빌드 시스템 정의
- ✅ Makefile : 간단한 빌드 인터페이스
- ✅ build.sh : 선택적 빌드 스크립트
