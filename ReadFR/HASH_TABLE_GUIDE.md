# 해시 테이블 기반 숫자 문자열 검색 알고리즘

## 개요
숫자로만 이루어진 문자형 변수에 대한 해시 기반 검색 알고리즘입니다.
체이닝(Chaining)을 사용하여 해시 충돌을 처리합니다.

## 특징
- **O(1)** 시간 복잡도의 검색 (평균)
- 동적 메모리 할당으로 유연한 크기 조정
- 충돌 처리: 체이닝 방식
- 로드 팩터 모니터링 가능
- 삽입, 검색, 삭제 모두 지원

## 파일 구조
```
src/M_HashTable.f90    - 해시 테이블 구현
test_hashtable.f90        - 테스트/사용 예제
```

## 사용 방법

### 1. 모듈 사용
```fortran
use M_HashTable
```

### 2. 데이터 타입 선언
```fortran
type(HashTable) :: my_table
```

### 3. 테이블 생성
```fortran
integer :: table_size
table_size = 101  ! 소수 권장 (충돌 최소화)
call ht_create(my_table, table_size)
```

### 4. 데이터 삽입
```fortran
call ht_insert(my_table, "21905009744", 1)
call ht_insert(my_table, "21905009529", 2)
```

### 5. 데이터 검색
```fortran
integer :: value
logical :: found

found = ht_search(my_table, "21905009744", value)
if (found) then
    print *, "찾음: ", value
else
    print *, "찾지 못함"
end if
```

### 6. 데이터 삭제
```fortran
logical :: deleted
deleted = ht_delete(my_table, "21905009744")
if (deleted) then
    print *, "삭제 완료"
end if
```

### 7. 통계 확인
```fortran
call ht_print_stats(my_table)
```

### 8. 메모리 해제
```fortran
call ht_free(my_table)
```

## 성능 비교

| 작업 | 평균 시간 | 최악 시간 |
|------|---------|---------|
| 삽입 | O(1) | O(n) |
| 검색 | O(1) | O(n) |
| 삭제 | O(1) | O(n) |

**n** = 테이블의 항목 수

## 해시 함수
숫자 문자열에 최적화된 다항 해시(Polynomial Rolling Hash):
```
hash = (문자1 × 31^(k-1) + 문자2 × 31^(k-2) + ... + 문자k) mod table_size
```

## 권장 사항
1. **테이블 크기**: 소수 사용 (101, 503, 1009, 10007, 100003)
2. **로드 팩터**: 0.75 미만일 때 최적 성능
3. **충돌 관리**: 로드 팩터가 높으면 테이블 크기 증가 고려

## 테스트 실행
```bash
./test_hashtable
```

## 출력 예제
```
===== 해시 테이블 통계 =====
테이블 크기:        101
저장된 항목:         10
로드 팩터:   0.0990
평균 체인 길이:   0.0990
최대 체인 길이:          2
```

## API 함수

### ht_create
해시 테이블을 생성합니다.
```fortran
subroutine ht_create(ht, table_size)
    type(HashTable), intent(out) :: ht
    integer, intent(in) :: table_size
end subroutine
```

### ht_insert
키-값 쌍을 삽입합니다.
```fortran
subroutine ht_insert(ht, key, value)
    type(HashTable), intent(inout) :: ht
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
end subroutine
```

### ht_search
키를 검색합니다.
```fortran
logical function ht_search(ht, key, value) result(found)
    type(HashTable), intent(in) :: ht
    character(len=*), intent(in) :: key
    integer, intent(out) :: value
end function
```

### ht_delete
키를 삭제합니다.
```fortran
logical function ht_delete(ht, key) result(deleted)
    type(HashTable), intent(inout) :: ht
    character(len=*), intent(in) :: key
end function
```

### ht_free
메모리를 해제합니다.
```fortran
subroutine ht_free(ht)
    type(HashTable), intent(inout) :: ht
end subroutine
```

### ht_print_stats
해시 테이블 통계를 출력합니다.
```fortran
subroutine ht_print_stats(ht)
    type(HashTable), intent(in) :: ht
end subroutine
```

## 예제 코드 (ReadFR에 적용)

아래는 Animal ID를 해시 테이블로 관리하는 예제입니다:

```fortran
program ReadFRWithHashTable
  use M_HashTable
  use M_Variables
  use M_ReadFile
  implicit none
  
  type(HashTable) :: animal_table
  integer :: value
  logical :: found
  
  ! 테이블 생성 (동물 수만큼 소수 크기 선택)
  call ht_create(animal_table, 1009)
  
  ! 데이터 읽으면서 삽입
  do i = 1, NREC
    call ht_insert(animal_table, trim(Animal_ID), i)
  end do
  
  ! 검색 (O(1) 속도)
  if (ht_search(animal_table, "21905009744", value)) then
    print *, "동물 찾음: 인덱스 =", value
  end if
  
  ! 정리
  call ht_free(animal_table)
end program
```

## 노트
- 숫자 문자열만 입력 (예: "12345", "21905009744")
- 해시 함수는 자동으로 숫자만 처리
- 같은 키로 재삽입하면 값이 업데이트됨
- 메모리 누수 방지를 위해 프로그램 종료 전 `ht_free` 반드시 호출
