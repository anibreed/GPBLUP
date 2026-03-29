# PED 해시 테이블 기반 검색 프로그램 - 사용 가이드

## 개요
PEDInfo 타입의 **ARN(Animal Registration Number)** 변수를 해시 키로 사용하여 O(1) 시간 복잡도의 고속 검색을 지원합니다.
가변 크기의 테이블을 지원하며, 자동으로 소수 크기로 최적화됩니다.

## 특징
- **O(1)** 시간 복잡도의 검색 (평균)
- **가변 테이블 크기**: 자동 소수 조정으로 충돌 최소화
- **64비트 ARN 지원**: 큰 숫자의 동물 ID 처리 가능
- **체이닝 충돌 처리**: 메모리 효율적인 해시 충돌 관리
- **동적 메모리 관리**: 필요한 만큼만 할당
- **전체 PED 데이터 저장**: ID, BREED, SIRE, DAM, SEX, BDate, LOC 등 모든 정보 저장

## 파일 구조
```
src/M_Variables.f90         - PEDInfo 타입 정의
src/M_PEDHashTable.f90      - PED 해시 테이블 구현
test_ped_hashtable.f90         - 테스트/사용 예제
```

## PEDInfo 타입 정의
```fortran
type, PUBLIC :: PEDInfo
    character(len=LEN_STR) :: BREED    ! 품종
    character(len=LEN_STR) :: ID       ! 동물 ID
    integer(kind=ki8) :: ARN           ! 동물 등록 번호 (해시 키)
    character(len=LEN_STR) :: SIRE     ! 부 ID
    character(len=LEN_STR) :: DAM      ! 모 ID
    integer :: SEX                     ! 성별 (1=수, 2=암)
    integer :: BDate                   ! 출생 날짜 (YYYYMMDD)
    character(len=LEN_STR) :: LOC      ! 위치
end type PEDInfo
```

## 핵심 함수

### 1. pht_create - 해시 테이블 생성
```fortran
call pht_create(ped_ht, table_size)
```
- **table_size**: 요청한 테이블 크기 (자동으로 소수로 조정)
- 예: 50 요청 → 101(소수)로 자동 조정

### 2. pht_insert - PED 정보 삽입
```fortran
type(PEDInfo) :: ped
ped%ARN = 2190500974_ki8
ped%ID = 'P001'
ped%BREED = 'Pig'
ped%SIRE = 'S001'
ped%DAM = 'D001'
ped%SEX = 2
ped%BDate = 20220115
ped%LOC = 'Farm1'

call pht_insert(ped_ht, ped)
```
- **기존 ARN 재삽입**: 자동으로 데이터 업데이트

### 3. pht_search - PED 정보 검색 (O(1) 속도)
```fortran
type(PEDInfo) :: found_ped
logical :: found

found = pht_search(ped_ht, 2190500974_ki8, found_ped)
if (found) then
    print *, "찾음: ", found_ped%ID
else
    print *, "찾지 못함"
end if
```
- **반환값**: .true. (찾음), .false. (못 찾음)

### 4. pht_delete - PED 정보 삭제
```fortran
logical :: deleted

deleted = pht_delete(ped_ht, 2190500974_ki8)
if (deleted) then
    print *, "삭제 완료"
end if
```

### 5. pht_free - 메모리 해제
```fortran
call pht_free(ped_ht)
```
- 프로그램 종료 전 반드시 호출 (메모리 누수 방지)

### 6. pht_print_stats - 통계 출력
```fortran
call pht_print_stats(ped_ht)
```

### 7. pht_print_ped_info - PED 정보 출력
```fortran
call pht_print_ped_info(found_ped)
```

## 사용 예제

### 기본 사용법
```fortran
program UsePEDHashTable
  use M_Kinds
  use M_Variables
  use M_PEDHashTable
  implicit none
  
  type(PEDHashTable) :: ped_ht
  type(PEDInfo) :: my_ped, found_ped
  logical :: found
  integer :: table_size
  
  ! 1. 테이블 생성 (크기 1000 요청)
  table_size = 1000
  call pht_create(ped_ht, table_size)
  
  ! 2. 데이터 준비 및 삽입
  my_ped%BREED = 'Duroc'
  my_ped%ID = 'DUROC001'
  my_ped%ARN = 2190500974_ki8
  my_ped%SIRE = 'S001'
  my_ped%DAM = 'D001'
  my_ped%SEX = 2
  my_ped%BDate = 20220115
  my_ped%LOC = 'Farm1'
  
  call pht_insert(ped_ht, my_ped)
  
  ! 3. 검색
  found = pht_search(ped_ht, 2190500974_ki8, found_ped)
  if (found) then
    call pht_print_ped_info(found_ped)
  end if
  
  ! 4. 통계
  call pht_print_stats(ped_ht)
  
  ! 5. 정리
  call pht_free(ped_ht)
end program UsePEDHashTable
```

### 파일에서 데이터 읽고 해시 테이블에 저장
```fortran
! 파일에서 읽은 데이터를 PED 해시 테이블에 저장
do i = 1, NREC
  ped%BREED = trim(XC(1))
  ped%ID = trim(XC(2))
  ped%ARN = int(XI(3), ki8)   ! ARN을 64비트로 변환
  ped%SIRE = trim(XC(4))
  ped%DAM = trim(XC(5))
  ped%SEX = XI(6)
  ped%BDate = XI(7)
  ped%LOC = trim(XC(8))
  
  call pht_insert(ped_ht, ped)  ! 해시 테이블에 저장
end do
```

## 성능 분석

| 작업 | 평균 | 최악 | 비고 |
|------|------|------|------|
| 삽입 | O(1) | O(n) | 로드 팩터 < 0.75일 때 |
| 검색 | O(1) | O(n) | 체이닝 충돌 가능 |
| 삭제 | O(1) | O(n) | 체이닝 충돌 가능 |

## 권장 설정

### 테이블 크기 선택
```fortran
! 데이터 100개 예상 시
table_size = 150  ! 자동으로 151(소수)로 조정

! 데이터 1000개 예상 시
table_size = 1500  ! 자동으로 1511(소수)로 조정
```

### 로드 팩터 모니터링
- **로드 팩터 = 저장된 항목 / 테이블 크기**
- 최적: 0.5 ~ 0.75
- 초과 시: 테이블 크기 확대 권장

## 테스트 실행

```bash
# 프로젝트 빌드
make clean && make

# PED 해시 테이블 테스트 컴파일
gfortran -O2 -g -Ibuild -Iinclude test_ped_hashtable.f90 \
  -Llib -ldkblupf90 -o test_ped_hashtable -Wl,-rpath,$(pwd)/lib

# 테스트 실행
./test_ped_hashtable
```

## 예상 출력
```
===== PED 해시 테이블 테스트 =====

[1] PED 해시 테이블 생성 (가변 크기)
PED 해시 테이블 생성: 요청 크기=        50 실제 크기(소수)=       101

[2] PED 데이터 삽입
  삽입: ARN=2190500974 ID=      P001
  삽입: ARN=2190500952 ID=      P002
  ... (생략)

[3] PED 해시 테이블 통계
===== PED 해시 테이블 통계 =====
테이블 크기:        101
저장된 개체:          8
로드 팩터:     0.0792
평균 체인 길이:     0.0792
최대 체인 길이:          1

[4] PED 검색 테스트 (성공)
--- ARN=2190500974 검색 ---
 검색 성공!
ARN:  2190500974 | ID:      P001 | BREED:         Pig | ...

[5] PED 검색 테스트 (실패)
--- ARN=9999999999 검색 ---
 검색 실패 (예상된 결과 - 존재하지 않음)

[6] PED 데이터 수정
... (생략)

[7] PED 데이터 삭제
--- ARN=2190500927 삭제 ---
 삭제 성공

[9] 메모리 해제
   PED 해시 테이블 메모리 해제 완료
```

## 주의사항

1. **ARN은 64비트(ki8)**: ARN 값 할당 시 `_ki8` 접미사 사용
   ```fortran
   ped%ARN = 2190500974_ki8  ! O
   ped%ARN = 2190500974      ! X (컴파일 오류)
   ```

2. **메모리 해제 필수**: 프로그램 종료 전 `pht_free()` 호출
   ```fortran
   call pht_free(ped_ht)  ! 필수!
   ```

3. **소수 자동 조정**: 사용자가 요청한 크기와 실제 크기가 다를 수 있음
   ```fortran
   call pht_create(ped_ht, 100)  ! 요청 100 → 실제 101(소수)
   ```

4. **포인터 관리**: 테이블은 자동으로 메모리 관리 (수동 deallocate 금지)

## 비교: 선형 탐색 vs 해시 테이블

| 데이터 수 | 선형 탐색 | 해시 테이블 | 개선율 |
|----------|---------|----------|-------|
| 100      | 50회    | 1회      | 50배  |
| 1,000    | 500회   | 1회      | 500배 |
| 10,000   | 5,000회 | 1회      | 5,000배 |
| 100,000  | 50,000회| 1회      | 50,000배|

## API 레퍼런스

### 데이터 타입
```fortran
type(PEDHashTable)  ! 해시 테이블
type(PEDInfo)       ! PED 레코드
```

### 공개 함수
```fortran
subroutine pht_create(pht, table_size)
subroutine pht_insert(pht, ped)
logical function pht_search(pht, arn_key, ped)
logical function pht_delete(pht, arn_key)
subroutine pht_free(pht)
subroutine pht_print_stats(pht)
subroutine pht_print_ped_info(ped)
integer function pht_get_prime_size(min_size)
```

## 라이선스
내부 프로젝트용

---

**문의**: 프로젝트 관리자
**최종 수정**: 2026-02-09
