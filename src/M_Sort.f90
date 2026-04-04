module M_Sort
implicit none

interface UniqX
   module procedure uniq_I, uniq_S
end interface

interface Search
   module procedure Search_S, Search_I, Search_R
end interface

interface Qsort
   module procedure Qsort_S, Qsort_I, Qsort_R, Qsort_S2
end interface

save
private
public UniqX, Search, Qsort, SearchC
contains
!-------------------------------------------------------------------------------------------------------
! Unique subroutine to rearrange by deleting duplicate value on array elements
!-------------------------------------------------------------------------------------------------------
 subroutine uniq_I(A, N)
  integer,intent(inout):: A(:)
  integer,intent(out):: N
  integer,allocatable:: swapA(:)
  integer::i, kk, missI = 0
  if(size(A).le.1) then
    N=1
    return
  endif
  allocate(swapA(size(A)))
  swapA=missI
  kk = 1
  do while(A(kk) == missI)
     kk = kk + 1
  enddo
  swapA(1) = A(kk)
  N=1
  do i=kk+1,size(A)
        if(swapA(N) == A(i)) cycle
        N = N + 1
        swapA(N) = A(i)
  enddo
  A(1:N)=swapA(1:N)
  A(N+1:)=missI
  return
end subroutine

subroutine uniq_S(A, N)
  character(LEN=*),intent(inout)::A(:)
  integer,intent(out):: N
  character(LEN=len(A(1))),allocatable:: swapA(:)
  character(len=1):: missC=achar(32), missA=achar(48)
  integer:: i, kk
  if(size(A) <= 1) then
    N=1
    return
  endif
  allocate(swapA(size(A)))
  swapA=missC
  kk = 1
  do while(A(kk) == missC .or. A(kk) == missA )
     kk = kk + 1
  enddo
  swapA(1) = A(kk)
  N=1
  do i=kk+1,size(A)
        if(swapA(N) == A(i)) cycle
        N = N + 1
        swapA(N) = A(i)
  enddo
  A(1:N)=swapA(1:N)
  A(N+1:)=missC
  return
end subroutine

recursive function Search_S (a, value) result (bsresult)
    character(len=*), intent(in) :: a(:), value
    integer          :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
        bsresult = 0        ! not found
    else if ( LGT(a(mid), value) ) then
        bsresult= Search_S(a(:mid-1), value)
    else if ( LLT(a(mid), value) ) then
        bsresult = Search_S(a(mid+1:), value)
        if (bsresult /= 0) then
            bsresult = mid + bsresult
        end if
    else
        bsresult = mid      ! SUCCESS!!
    end if
end function

recursive function Search_I (a, value) result (bsresult)
    INTEGER, intent(in) :: a(:), value
    integer          :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
        bsresult = 0        ! not found
    else if (a(mid) > value ) then
        bsresult= Search_I(a(:mid-1), value)
    else if ( a(mid) < value ) then
        bsresult = Search_I(a(mid+1:), value)
        if (bsresult /= 0) then
            bsresult = mid + bsresult
        end if
    else
        bsresult = mid      ! SUCCESS!!
    end if
end function

recursive function Search_R (a, value) result (bsresult)
    REAL, intent(in) :: a(:), value
    integer          :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
        bsresult = 0        ! not found
    else if (a(mid) > value ) then
        bsresult= Search_R(a(:mid-1), value)
    else if ( a(mid) < value ) then
        bsresult = Search_R(a(mid+1:), value)
        if (bsresult /= 0) then
            bsresult = mid + bsresult
        end if
    else
        bsresult = mid      ! SUCCESS!!
    end if
end function

RECURSIVE SUBROUTINE Qsort_S(list, order, desc)
character(len=*), INTENT(INOUT):: list(:)
INTEGER, INTENT(INOUT),optional:: order(:)
logical, INTENT(IN),optional:: desc
CALL quick_sort_S(1, SIZE(list))
CONTAINS
 RECURSIVE SUBROUTINE quick_sort_S(left_end, right_end)
 INTEGER, INTENT(IN):: left_end, right_end
 INTEGER:: i, j, itemp
 character(len=len(list)):: reference, temp
 INTEGER, PARAMETER  :: max_simple_sort_size = 6
 IF (right_end < left_end + max_simple_sort_size) THEN
   CALL swap_S(left_end, right_end)
 ELSE
   reference = list((left_end + right_end)/2)
   i = left_end - 1; j = right_end + 1
   DO
     DO
        i = i + 1
        if(present(desc)) then
            IF (LLE(list(i),reference)) EXIT
        else
            IF (LGE(list(i),reference)) EXIT
        endif
     END DO
     DO
        j = j - 1
        if(present(desc)) then
            IF (LGE(list(j),reference)) EXIT
        else
            IF (LLE(list(j),reference)) EXIT
        endif
     END DO
     IF (i < j) THEN
        temp = list(i); list(i) = list(j); list(j) = temp
        if(present(order)) then
            itemp = order(i); order(i) = order(j); order(j) = itemp
        endif
     ELSE IF (i == j) THEN
       i = i + 1
       EXIT
     ELSE
       EXIT
     END IF
   END DO
   IF (left_end < j) CALL quick_sort_S(left_end, j)
   IF (i < right_end) CALL quick_sort_S(i, right_end)
 END IF
 END SUBROUTINE quick_sort_S
 SUBROUTINE swap_S(left_end, right_end)
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 character(len=len(list)):: temp
 DO i = left_end, right_end - 1
    DO j = i+1, right_end
        IF (LGT(list(i),list(j))) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             if(present(order)) then
                 itemp = order(i); order(i) = order(j); order(j) = itemp
             endif
        END IF
    END DO
 END DO
 END SUBROUTINE swap_S
END SUBROUTINE Qsort_S

RECURSIVE SUBROUTINE Qsort_I(list, order,desc)
Integer,INTENT(INOUT):: list(:)
INTEGER,INTENT(INOUT), optional:: order(:)
logical, INTENT(IN),optional:: desc
CALL quick_sort_I(1, SIZE(list))
CONTAINS
 RECURSIVE SUBROUTINE quick_sort_I(left_end, right_end)
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 INTEGER             :: reference, temp
 INTEGER, PARAMETER  :: max_simple_sort_size = 6
 IF (right_end < left_end + max_simple_sort_size) THEN
   CALL swap_I(left_end, right_end)
 ELSE
   reference = list((left_end + right_end)/2)
   i = left_end - 1; j = right_end + 1
   DO
     DO
       i = i + 1
       if(present(desc)) then
           IF (list(i) <= reference) EXIT
       else
           IF (list(i) >= reference) EXIT
       endif
     END DO
     DO
       j = j - 1
       if(present(desc)) then
          IF (list(j) >= reference) EXIT
       else
          IF (list(j) <= reference) EXIT
       endif
     END DO
     IF (i < j) THEN
         temp = list(i); list(i) = list(j); list(j) = temp
        if(present(order)) then
            itemp = order(i); order(i) = order(j); order(j) = itemp
        endif
     ELSE IF (i == j) THEN
       i = i + 1
       EXIT
     ELSE
       EXIT
     END IF
   END DO
   IF (left_end < j) CALL quick_sort_I(left_end, j)
   IF (i < right_end) CALL quick_sort_I(i, right_end)
 END IF
 END SUBROUTINE quick_sort_I
 SUBROUTINE swap_I(left_end, right_end)
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 INTEGER             :: temp
 DO i = left_end, right_end - 1
   DO j = i+1, right_end
     IF (list(i) > list(j)) THEN
        temp = list(i); list(i) = list(j); list(j) = temp
        if(present(order)) then
            itemp = order(i); order(i) = order(j); order(j) = itemp
        endif
     END IF
   END DO
 END DO
 END SUBROUTINE swap_I
END SUBROUTINE Qsort_I

RECURSIVE SUBROUTINE Qsort_R(list, order,desc)
REAL,INTENT(INOUT):: list(:)
INTEGER,INTENT(INOUT),optional:: order(:)
logical, INTENT(IN),optional:: desc
CALL quick_sort_R(1, SIZE(list))
CONTAINS
 RECURSIVE SUBROUTINE quick_sort_R(left_end, right_end)
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 REAL                    :: reference, temp
 INTEGER, PARAMETER  :: max_simple_sort_size = 6
 IF (right_end < left_end + max_simple_sort_size) THEN
   CALL swap_R(left_end, right_end)
 ELSE
   reference = list((left_end + right_end)/2)
   i = left_end - 1; j = right_end + 1
   DO
     DO
       i = i + 1
       if(present(desc)) then
           IF (list(i) <= reference) EXIT
       else
           IF (list(i) >= reference) EXIT
       endif
     END DO
     DO
       j = j - 1
       if(present(desc)) then
           IF (list(j) >= reference) EXIT
       else
           IF (list(j) <= reference) EXIT
       endif
     END DO
     IF (i < j) THEN
        temp = list(i); list(i) = list(j); list(j) = temp
        if(present(order)) then
            itemp = order(i); order(i) = order(j); order(j) = itemp
        endif
     ELSE IF (i == j) THEN
       i = i + 1
       EXIT
     ELSE
       EXIT
     END IF
   END DO
   IF (left_end < j) CALL quick_sort_R(left_end, j)
   IF (i < right_end) CALL quick_sort_R(i, right_end)
 END IF
 END SUBROUTINE quick_sort_R
 SUBROUTINE swap_R(left_end, right_end)
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 REAL                     :: temp
 DO i = left_end, right_end - 1
   DO j = i+1, right_end
     IF (list(i) > list(j)) THEN
         temp = list(i); list(i) = list(j); list(j) = temp
        if(present(order)) then
            itemp = order(i); order(i) = order(j); order(j) = itemp
        endif
     END IF
   END DO
 END DO
 END SUBROUTINE swap_R
END SUBROUTINE Qsort_R

RECURSIVE SUBROUTINE Qsort_S2(list, key,order)
implicit none
character(len=*), INTENT(INOUT):: list(:,:)
INTEGER, INTENT(IN):: key
INTEGER, INTENT(INOUT):: order(:)
CALL quick_sort_S2(1, SIZE(list,1))
CONTAINS

 RECURSIVE SUBROUTINE quick_sort_S2(left_end, right_end)
 implicit none
 INTEGER, INTENT(IN):: left_end, right_end
 INTEGER:: i, j, itemp
 character(len=len(list)):: reference, temp(size(list,2))
 INTEGER, PARAMETER  :: max_simple_sort_size = 6
 IF (right_end < left_end + max_simple_sort_size) THEN
   CALL swap_S2(left_end, right_end)
 ELSE
   reference = list((left_end + right_end)/2, key)
   i = left_end - 1
   j = right_end + 1
   DO
     DO
       i = i + 1
       IF ( LGE(list(i,key), reference) ) EXIT
     END DO
     DO
       j = j - 1
       IF ( LLE(list(j,key), reference) ) EXIT
     END DO
     IF (i < j) THEN
        temp = list(i,:); list(i,:) = list(j,:); list(j,:) = temp
        itemp = order(i); order(i) = order(j); order(j) = itemp
     ELSE IF (i == j) THEN
         i = i + 1
         EXIT
     ELSE
         EXIT
     END IF
   END DO
   IF (left_end < j) CALL quick_sort_S2(left_end, j)
   IF (i < right_end) CALL quick_sort_S2(i, right_end)
 END IF
 END SUBROUTINE quick_sort_S2

 SUBROUTINE swap_S2(left_end, right_end)
 implicit none
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 character(len=len(list)):: temp(size(list,2))
 DO i = left_end, right_end - 1
   DO j = i+1, right_end
     IF ( LGT(list(i,key), list(j,key)) ) THEN
        temp = list(i,:); list(i,:) = list(j,:); list(j,:) = temp
        itemp = order(i); order(i) = order(j); order(j) = itemp
     END IF
   END DO
 END DO
 END SUBROUTINE swap_S2
END SUBROUTINE Qsort_S2

recursive function SearchC(a,value) result(idx)
   character(len=*),intent(in):: a(:), value
   integer:: idx
   idx=binary_search_left(a,value,1,size(a))
   if(a(idx)/=value) idx= -1
end function
recursive function binary_search_left(a,value,low,high) result(idx)
   character(len=*),intent(in):: a(:), value
   integer,intent(in):: low, high
   integer:: idx
   if(high<=low) then
     idx=low
   else
     idx=low+(high-low)/2
     if(a(idx) < value) then
        idx=binary_search_left(a,value,idx+1,high)
     else
        idx=binary_search_left(a,value,low,idx)
     endif
   endif
end function

end module M_Sort
