module Qsort4
implicit none
interface Qsort
   module procedure Qsort_S, Qsort_S2, Qsort_R
end interface

! Quick sort routine from: Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990)
! "Programmer's Guide to Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.

save
private
public Qsort, SearchC
contains
 
RECURSIVE SUBROUTINE Qsort_S(list, order)
implicit none
character(len=*), INTENT(INOUT):: list(:)
INTEGER, INTENT(INOUT):: order(:)
CALL quick_sort_S(1, SIZE(list))
CONTAINS

 RECURSIVE SUBROUTINE quick_sort_S(left_end, right_end)
 implicit none
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
       IF (LGE(list(i),reference)) EXIT
     END DO
     DO
       j = j - 1
       IF (LLE(list(j),reference)) EXIT
     END DO
     IF (i < j) THEN
       temp = list(i); list(i) = list(j); list(j) = temp
       itemp = order(i); order(i) = order(j); order(j) = itemp
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
 implicit none
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 character(len=len(list)):: temp
 DO i = left_end, right_end - 1
    DO j = i+1, right_end
        IF (LGT(list(i),list(j))) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
        END IF
    END DO
 END DO
 END SUBROUTINE swap_S
END SUBROUTINE Qsort_S

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

RECURSIVE SUBROUTINE Qsort_R(list, order)
implicit none
REAL,INTENT(INOUT):: list(:)
INTEGER,INTENT(INOUT):: order(:)
CALL quick_sort_R(1, SIZE(list))
CONTAINS

 RECURSIVE SUBROUTINE quick_sort_R(left_end, right_end)
 implicit none
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
       IF (list(i) >= reference) EXIT
     END DO
     DO
       j = j - 1
       IF (list(j) <= reference) EXIT
     END DO
     IF (i < j) THEN
       temp = list(i); list(i) = list(j); list(j) = temp
       itemp = order(i); order(i) = order(j); order(j) = itemp
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
 implicit none
 INTEGER, INTENT(IN) :: left_end, right_end
 INTEGER             :: i, j, itemp
 REAL                     :: temp
 DO i = left_end, right_end - 1
   DO j = i+1, right_end
     IF (list(i) > list(j)) THEN
       temp = list(i); list(i) = list(j); list(j) = temp
       itemp = order(i); order(i) = order(j); order(j) = itemp
     END IF
   END DO
 END DO
 END SUBROUTINE swap_R
END SUBROUTINE Qsort_R

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

end module Qsort4

