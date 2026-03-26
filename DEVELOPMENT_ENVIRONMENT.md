# 📦 GPBLUP 프로그램 개발 환경 구성 완료

**작성일:** Feb 19, 2026  
**프로젝트:** GPBLUP (Genomic Prediction with BLUP)

---

## 🎯 목표 달성 현황

✅ **통합 프로젝트 구조 완성**
- 공통 모듈 (source/) 중앙화
- 단위 프로그램 (ReadFR, popQC) 독립 개발
- 자동 dependency 추적 및 재컴파일

✅ **자동 빌드 시스템 구축**
- CMake 기반 통합 빌드
- Makefile 인터페이스 제공
- 병렬 컴파일 지원

✅ **새 프로그램 개발 자동화**
- 프로그램 템플릿 생성 스크립트
- 완전한 개발 가이드 제공
- CMakeLists.txt 자동 등록

---

## 📚 개발 가이드 & 도구

| 파일 | 설명 | 대상 |
|------|------|------|
| **QUICK_START.md** | 🚀 30초 시작 가이드 | 신입 개발자 |
| **PROGRAM_DEVELOPMENT_GUIDE.md** | 📖 상세 개발 가이드 | 모든 개발자 |
| **BUILD_GUIDE.md** | 🔨 빌드 시스템 가이드 | 빌드 관리자 |
| **create_program.sh** | 🤖 템플릿 생성 도구 | 자동화 |
| **add_to_cmake.sh** | 🔧 CMake 등록 헬퍼 | 자동화 |

---

## 🚀 빠른 시작 (3단계)

### Step 1: 프로그램 생성
```bash
cd /home/dhlee/GPBLUP
./create_program.sh MyProgram "My program description"
```

### Step 2: 로직 구현
```bash
vim MyProgram/MyProgram.f90
# main_process() 서브루틴에 코드 작성
```

### Step 3: 빌드
```bash
# CMakeLists.txt에 프로그램 추가 (템플릿 제공)
vim CMakeLists.txt

# 빌드
make rebuild
```

---

## 📁 프로젝트 최종 구조

```
GPBLUP/
├── 📚 개발 가이드
│   ├── QUICK_START.md                    ← 30초 시작 (필독!)
│   ├── PROGRAM_DEVELOPMENT_GUIDE.md      ← 상세 가이드
│   ├── BUILD_GUIDE.md                    ← 빌드 시스템
│   └── README.md
│
├── 🤖 개발 도구
│   ├── create_program.sh                 ← 템플릿 생성
│   └── add_to_cmake.sh                   ← CMake 등록
│
├── 🔨 빌드 설정
│   ├── CMakeLists.txt                    ← 모든 프로그램 정의
│   ├── Makefile                          ← 간편한 빌드
│   └── build.sh                          ← 빌드 스크립트
│
├── 📦 공통 모듈 라이브러리
│   └── source/                           ← 10개 공통 모듈
│       ├── M_Kinds.f90
│       ├── M_Variables.f90               ← 전역변수 정의
│       ├── M_readpar.f90                 ← 파라미터 파서
│       ├── M_ReadFile.f90
│       ├── M_StrEdit.f90
│       ├── M_Stamp.f90
│       ├── M_HashTable.f90
│       ├── M_PEDHashTable.f90
│       ├── M_param.f90
│       └── Qsort4.f90
│
├── 📦 프로그램 A: ReadFR
│   ├── ReadFR.f90                        ← 메인 소스
│   ├── 문서 (5개)
│   └── test/                             ← 테스트 데이터
│
├── 📦 프로그램 B: popQC
│   ├── popQC.f90                         ← 메인 소스
│   ├── 문서 (4개)
│   └── test/                             ← 테스트 파라미터
│
├── 📦 프로그램 C: [앞으로 추가할 프로그램]
│   ├── [Program].f90
│   ├── README.md
│   └── test/
│
├── 📁 bin/                               ← 최종 실행파일
│   ├── ReadFR (89K)
│   ├── popQC (124K)
│   └── [앞으로 추가될 프로그램들]
│
├── 📁 build/                             ← CMake 빌드 디렉토리
│   ├── modules/                          ← .mod 파일들
│   └── ...
│
└── 📁 lib/                               ← 라이브러리
```

---

## 🔄 개발 워크플로우

### 일반적인 사이클

```
1. 프로그램 생성
   ./create_program.sh NewProgram

2. 로직 구현
   vim NewProgram/NewProgram.f90
   ↓
   (반복: 수정 → make → 테스트)

3. CMakeLists.txt 등록
   vim CMakeLists.txt
   (템플릿 제공)

4. 최종 빌드
   make rebuild

5. 배포
   bin/NewProgram 자동 생성됨
```

### 공통 모듈 수정 시

```
1. 모듈 수정
   vim source/M_Variables.f90

2. 재빌드 (자동으로 모든 프로그램 재컴파일)
   make

결과: ReadFR, popQC, [모든 프로그램] 자동 업데이트 ✓
```

---

## 💡 핵심 개념

### 1. 계층적 모듈 설계

```
┌─────────────────────────────────────┐
│  개별 프로그램                       │
│  (ReadFR, popQC, MyProgram, ...)   │
└──────────────┬──────────────────────┘
               │ 링크
┌──────────────▼──────────────────────┐
│  공통 모듈 라이브러리               │
│  (gpblup_modules.a)                │
│  - M_readpar (파라미터 파서)       │
│  - M_Variables (전역 변수)         │
│  - 기타 유틸리티 모듈들             │
└──────────────┬──────────────────────┘
               │ 사용 (use)
┌──────────────▼──────────────────────┐
│  Fortran 표준 라이브러리             │
└─────────────────────────────────────┘
```

### 2. 자동 Dependency 추적

```cmake
# CMakeLists.txt의 메커니즘

add_library(gpblup_modules ${MODULE_SOURCES})
# source 디렉토리의 모든 .f90 → 정적 라이브러리로 컴파일

add_executable(MyProgram MyProgram/MyProgram.f90)
target_link_libraries(MyProgram PRIVATE gpblup_modules)
add_dependencies(MyProgram gpblup_modules)
# MyProgram은 gpblup_modules에 완전히 의존

# 결과:
# source/ 파일 수정 → gpblup_modules 재컴파일
# → MyProgram, ReadFR, popQC 모두 자동 재링크 ✓
```

### 3. 공유 파라미터 시스템

```
파라미터 파일
    ↓
M_readpar.read_parameters()
    ↓
전역 변수 채우기:
  - PEDFile%FileName, FieldName, FieldLoc, ...
  - MAPFile%...
  - SNPFile%...
  - GENOFile%...
  - OutputPrefix
  - QCThresholds
    ↓
모든 프로그램이 공통으로 사용
```

---

## 📊 현재 상태

| 항목 | 상태 | 세부사항 |
|------|------|---------|
| **ReadFR** | ✅ 완성 | 89K, Feb 19 17:47 |
| **popQC** | ✅ 완성 | 124K, Feb 19 17:47 |
| **공통 모듈** | ✅ 10개 | M_Kinds ~ Qsort4 |
| **빌드 시스템** | ✅ CMake | Makefile 인터페이스 |
| **개발 도구** | ✅ 자동화 | create_program.sh 제공 |
| **문서** | ✅ 완전 | 5개 개발 가이드 |

---

## 🎓 학습 경로

### 신입 개발자

1. **QUICK_START.md 읽기** (5분)
2. **예제 프로그램 생성 & 실행** (10분)
3. **PROGRAM_DEVELOPMENT_GUIDE.md 읽기** (15분)
4. **자신의 프로그램 개발 시작** (30분+)

### 경험 개발자

- **create_program.sh** → 템플릿 자동 생성
- **CMakeLists.txt 템플릿** 복사 + 수정
- **make rebuild** → 완료

---

## 🔐 모범 사례

### ✅ DO (해야 할 것)

- ✓ 공통 기능은 source/ 모듈에 작성
- ✓ 프로그램별 로직은 각 디렉토리에 작성
- ✓ M_readpar로 파라미터 읽기
- ✓ 전역 변수 (M_Variables) 활용
- ✓ 에러 처리 구현

### ❌ DON'T (하지 말 것)

- ✗ 각 프로그램에서 파라미터 파싱 직접 구현
- ✗ 겹치는 기능을 여러 프로그램에서 재작성
- ✗ source/ 모듈을 직접 수정 없이 개발
- ✗ 하드코딩된 파일명 사용
- ✗ 에러 처리 생략

---

## 🚀 향후 계획

### Phase 1: ✅ 완료
- [x] 통합 프로젝트 구조 구축
- [x] CMake 빌드 시스템
- [x] ReadFR, popQC 통합
- [x] 자동 dependency 추적

### Phase 2: 🔄 진행 중
- [ ] 새 프로그램들 지속적 개발
- [ ] 프로그램별 add_dependencies 추가
- [ ] bin/ 에 자동 업데이트

### Phase 3: 🎯 계획 중
- [ ] 통합 테스트 스위트
- [ ] CI/CD 파이프라인
- [ ] 성능 프로파일링
- [ ] 병렬 처리 최적화

---

## 🆘 도움말

### 문서 찾기
```bash
cd /home/dhlee/GPBLUP
ls -lh *.md                  # 모든 가이드 문서
cat QUICK_START.md           # 빠른 시작 (필독!)
cat PROGRAM_DEVELOPMENT_GUIDE.md  # 상세 개발 가이드
```

### 스크립트 도움말
```bash
./create_program.sh          # 사용법 표시
./add_to_cmake.sh            # 사용법 표시
make help                    # 빌드 명령 설명
```

### 빌드 확인
```bash
make status                  # 현재 빌드 상태
ls -lh bin/                  # 생성된 바이너리
make clean && make           # 전체 정리 후 재빌드
```

---

## 📈 통계

| 항목 | 수량 |
|------|------|
| **소스 모듈** | 10개 |
| **프로그램** | 2개 (ReadFR, popQC) |
| **개발 가이드** | 5개 |
| **자동화 스크립트** | 2개 |
| **총 소스 라인** | ~3000 라인 |
| **프로젝트 크기** | 2.1GB |

---

## ✨ 특징 요약

```
🎯 통합 프로젝트
  • 공통 모듈 중앙화
  • 자동 dependency 추적
  • 호환성 보장

🤖 자동화
  • 템플릿 생성 (create_program.sh)
  • 자동 재컴파일
  • CMake 통합

📚 완전한 문서
  • QUICK_START (필독!)
  • 개발 가이드
  • 빌드 시스템 설명

🚀 효율적 개발
  • 프로그램 추가 3단계
  • 병렬 컴파일
  • 빠른 피드백 루프
```

---

## 🎉 결론

GPBLUP 프로젝트는 이제 **지속적인 프로그램 개발에 최적화된 통합 플랫폼**입니다.

### 다음 단계
1. [QUICK_START.md](QUICK_START.md) 읽기
2. `./create_program.sh` 로 첫 프로그램 생성
3. 계속 프로그램 추가하기

### 연락처 및 지원
- 개발 가이드: [PROGRAM_DEVELOPMENT_GUIDE.md](PROGRAM_DEVELOPMENT_GUIDE.md)
- 빌드 문제: [BUILD_GUIDE.md](BUILD_GUIDE.md)
- 빠른 시작: [QUICK_START.md](QUICK_START.md)

---

**행운을 빕니다! Happy Coding! 🎯**

---

*마지막 업데이트: Feb 19, 2026*
