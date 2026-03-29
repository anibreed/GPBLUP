module M_StrEdit
! ============================================================================
! M_StrEdit - String Editing and Manipulation Module
! ============================================================================
! Purpose:
!   Comprehensive string manipulation functions for Fortran programs
!   Optimized for performance with frequently-used functions
!
! Core Functions (Most Used):
!   - to_upper, to_lower: Case conversion
!   - split_string: Parse strings by delimiter
!   - FindDelim: Convert delimiter names to characters
!   - is_numeric, is_letter, is_digit: Type checking
!   - str_to_*, NumToStr: Type conversions
!   - resolve_file_path, is_absolute_path: Path utilities
!   - to_logical: String to boolean conversion
!
! Author: Dr. DEUKMIN LEE
! Updated: March 2026 - Performance optimizations
! ============================================================================
use M_Kinds
use M_Variables
implicit none
private

! Public Interface (alphabetically organized)
public :: concat, to_upper, to_lower, to_logical
public :: remove_special_chars, remove_all_special_chars, clean_string
public :: parse, compact, removesp
public :: shiftstr, insertstr, delsubstr, delall, match, trimzero
public :: NumToStr
public :: str_to_i4, str_to_i8, str_to_r4, str_to_r8
public :: strvec_to_r4, strvec_to_i4, strvec_to_num
public :: remove_duplicate_spaces, split_string
public :: is_str, is_letter, is_digit, is_numeric
public :: FindDelim
public :: get_directory_path, resolve_file_path, is_absolute_path
public :: yyyymmdd_to_date, bdate_to_jdn
public :: get_field_pos
public :: operator(+)

INTERFACE OPERATOR (+)
 MODULE PROCEDURE concat
END INTERFACE

interface NumToStr
   module procedure i4_to_str
   module procedure i8_to_str
   module procedure r4_to_str
   module procedure r8_to_str
end interface

interface strvec_to_num
   module procedure strvec_to_r4
   module procedure strvec_to_i4
end interface

contains
!**********************************************************************
 FUNCTION concat(cha,chb)
 IMPLICIT NONE
 CHARACTER (LEN=*), INTENT(IN) :: cha, chb 
 CHARACTER (LEN=(LEN_TRIM(cha) + LEN_TRIM(chb))) :: concat
 concat = TRIM(cha)//TRIM(chb)
 END FUNCTION concat
!**********************************************************************
subroutine remove_special_chars(str, keep_chars)
! 문자열에서 특수문자를 제거하는 subroutine
! 알파벳(a-z, A-Z), 숫자(0-9), 공백만 남기고 나머지 제거
! 선택적으로 keep_chars에 지정된 문자는 유지 가능
! Arguments:
!   str        - (inout) 처리할 문자열
!   keep_chars - (optional, in) 유지할 추가 문자들 (예: ".-_")

character(len=*), intent(inout) :: str
character(len=*), intent(in), optional :: keep_chars
character(len=len(str)) :: temp_str
integer :: i, j
logical :: keep_char

temp_str = ' '
j = 0

do i = 1, len_trim(str)
   keep_char = .false.
   
   if (is_letter(str(i:i)) .or. is_digit(str(i:i))) keep_char = .true.
   if (str(i:i) == ' ') keep_char = .true.
   
   ! 추가로 유지할 문자 확인
   if (present(keep_chars)) then
      if (index(keep_chars, str(i:i)) > 0) keep_char = .true.
   end if
   
   ! 문자 유지
   if (keep_char) then
      j = j + 1
      temp_str(j:j) = str(i:i)
   end if
end do

str = temp_str

end subroutine remove_special_chars

!**********************************************************************

subroutine remove_all_special_chars(str)
! 문자열에서 모든 특수문자를 공백으로 대체
! 알파벳(a-z, A-Z)과 숫자(0-9)만 유지, 나머지는 공백으로 변환
!
! Arguments:
!   str - (inout) 처리할 문자열

character(len=*), intent(inout) :: str
integer :: i

do i = 1, len(str)
   ! 숫자가 아니고, 대문자가 아니고, 소문자가 아니면 공백으로 변경
   if (.not. (is_letter(str(i:i)) .or. is_digit(str(i:i)))) then
      str(i:i) = ' '
   end if
end do

end subroutine remove_all_special_chars

!**********************************************************************

subroutine clean_string(str, mode)
! 다양한 모드로 문자열 정리
!
! Arguments:
!   str  - (inout) 처리할 문자열
!   mode - (optional, in) 정리 모드
!          'alphanum'   : 알파벳과 숫자만 유지 (기본값)
!          'alpha'      : 알파벳만 유지
!          'numeric'    : 숫자만 유지
!          'printable'  : 출력 가능한 ASCII만 유지 (32-126)
!          'identifier' : 식별자 형식 (알파벳, 숫자, 언더스코어)

character(len=*), intent(inout) :: str
character(len=*), intent(in), optional :: mode
character(len=len(str)) :: temp_str
character(len=20) :: clean_mode
integer :: i, j, ichar_val
logical :: keep_char

! 기본 모드 설정
clean_mode = 'alphanum'
if (present(mode)) then
   clean_mode = trim(adjustl(mode))
   clean_mode = to_lower(clean_mode)
end if

temp_str = ' '
j = 0

do i = 1, len_trim(str)
   ichar_val = iachar(str(i:i))
   keep_char = .false.
   
   select case (trim(clean_mode))
   case ('alphanum')
      ! 알파벳과 숫자, 공백 유지
      if ((ichar_val >= 48 .and. ichar_val <= 57) .or. &   ! 0-9
          (ichar_val >= 65 .and. ichar_val <= 90) .or. &   ! A-Z
          (ichar_val >= 97 .and. ichar_val <= 122) .or. &  ! a-z
          (ichar_val == 32)) then                          ! space
         keep_char = .true.
      end if
      
   case ('alpha')
      ! 알파벳과 공백만 유지
      if ((ichar_val >= 65 .and. ichar_val <= 90) .or. &   ! A-Z
          (ichar_val >= 97 .and. ichar_val <= 122) .or. &  ! a-z
          (ichar_val == 32)) then                          ! space
         keep_char = .true.
      end if
      
   case ('numeric')
      ! 숫자와 공백만 유지
      if ((ichar_val >= 48 .and. ichar_val <= 57) .or. &   ! 0-9
          (ichar_val == 32)) then                          ! space
         keep_char = .true.
      end if
      
   case ('printable')
      ! 출력 가능한 ASCII 문자만 유지 (32-126)
      if (ichar_val >= 32 .and. ichar_val <= 126) then
         keep_char = .true.
      end if
      
   case ('identifier')
      ! 식별자(형식: 알파벳, 숫자, 언더스코어) 유지
      if ((ichar_val >= 48 .and. ichar_val <= 57) .or. &   ! 0-9
          (ichar_val >= 65 .and. ichar_val <= 90) .or. &   ! A-Z
          (ichar_val >= 97 .and. ichar_val <= 122) .or. &  ! a-z
          (ichar_val == 95)) then                          ! underscore
         keep_char = .true.
      end if
      
   case default
      ! 기본(alphanum)만 유지
      if ((ichar_val >= 48 .and. ichar_val <= 57) .or. &   ! 0-9
          (ichar_val >= 65 .and. ichar_val <= 90) .or. &   ! A-Z
          (ichar_val >= 97 .and. ichar_val <= 122) .or. &  ! a-z
          (ichar_val == 32)) then                          ! space
         keep_char = .true.
      end if
   end select
   
   if (keep_char) then
      j = j + 1
      temp_str(j:j) = str(i:i)
   end if
end do

str = temp_str

end subroutine clean_string


!**********************************************************************

function to_upper(str) result(outstr)
! Convert string to upper case (function)

character(len=*), intent(in) :: str
character(len=len(str)) :: outstr
integer :: i, iav, ioffset

outstr = str
ioffset = iachar('A') - iachar('a')
do i = 1, len(str)
   iav = iachar(str(i:i))
   if (iav >= iachar('a') .and. iav <= iachar('z')) then
      outstr(i:i) = achar(iav + ioffset)
   end if
end do

end function to_upper

!**********************************************************************

function to_lower(str) result(outstr)
! Convert string to lower case (function)

character(len=*), intent(in) :: str
character(len=len(str)) :: outstr
integer :: i, iav, ioffset

outstr = str
ioffset = iachar('A') - iachar('a')
do i = 1, len(str)
   iav = iachar(str(i:i))
   if (iav >= iachar('A') .and. iav <= iachar('Z')) then
      outstr(i:i) = achar(iav - ioffset)
   end if
end do

end function to_lower

function i4_to_str(num, fmt) result(str)
    integer(kind=ki4), intent(in) :: num
    character(len=*), intent(in), optional :: fmt
    character(len=64) :: str
    integer :: ios
    if (present(fmt)) then
       write(str, fmt=fmt, iostat=ios) num
    else
       write(str, fmt='(i0)', iostat=ios) num
    end if
    if (ios /= 0) str = ' '
    str = trim(adjustl(str))
end function

function i8_to_str(num, fmt) result(str)
    integer(kind=ki8), intent(in) :: num
    character(len=*), intent(in), optional :: fmt
    character(len=64) :: str
    integer :: ios
    if (present(fmt)) then
       write(str, fmt=fmt, iostat=ios) num
    else
       write(str, fmt='(i0)', iostat=ios) num
    end if
    if (ios /= 0) str = ' '
    str = trim(adjustl(str))
end function

function r4_to_str(num, fmt) result(str)
    real(kind=r4), intent(in) :: num
    character(len=*), intent(in), optional :: fmt
    character(len=64) :: str
    integer :: ios
    if (present(fmt)) then
       write(str, fmt=fmt, iostat=ios) num
    else
       write(str, fmt='(g0)', iostat=ios) num
    end if
    if (ios /= 0) str = ' '
    str = trim(adjustl(str))
end function

function r8_to_str(num, fmt) result(str)
    real(kind=r8), intent(in) :: num
    character(len=*), intent(in), optional :: fmt
    character(len=64) :: str
    integer :: ios
    if (present(fmt)) then
       write(str, fmt=fmt, iostat=ios) num
    else
       write(str, fmt='(g0)', iostat=ios) num
    end if
    if (ios /= 0) str = ' '
    str = trim(adjustl(str))
end function

!**********************************************************************

function str_to_i4(str) result(num)
! Convert string to single precision integer
character(len=*), intent(in) :: str
integer(kind=ki4) :: num
real(kind=r8) :: rnum

if(.not.is_numeric(str)) then
   num=0_ki4
   return
end if   
rnum = str_to_r8(str)
if(abs(rnum) > real(huge(num), r8)) then
  num = 0
  return
end if
num = nint(rnum, ki4)
end function str_to_i4
!**********************************************************************

function str_to_i8(str) result(num)
! Convert string to double precision integer
character(len=*), intent(in) :: str
integer(kind=ki8) :: num
real(r8) :: rnum
if(.not.is_numeric(str)) then
   num=0_ki8
   return
end if   
rnum = str_to_r8(str)
if(abs(rnum) > real(huge(num), r8)) then
  num=0
  return
end if
num=nint(rnum,ki8)
end function str_to_i8

!**********************************************************************

function str_to_r4(str) result(num)
! Convert string to single precision real
character(len=*), intent(in) :: str
real(kind=r4) :: num
real(kind=r8) :: rnum 

! str_to_r8_internal을 호출하여 재귀 호출 방지
if(.not.is_numeric(str)) then
   num=0.0_r4
   return
end if   
rnum = str_to_r8(str)
if( abs(rnum) > huge(num) ) then
  num=0.0_r4
  return
end if
if( abs(rnum) > 0.0_r8 .and. abs(rnum) < tiny(num) ) then
  num=0.0_r4
  return
end if
num=real(rnum,r4)
end function str_to_r4

!**********************************************************************

function str_to_r8(str) result(num)
! Convert string to double precision real (internal function)
character(len=*), intent(in) :: str
real(kind=r8) :: num
integer :: ilen, ipos, ios

ilen=len_trim(str)
! 빈 문자열 처리
if(ilen == 0) then
   num=0.0_r8
   return
end if
ipos=scan(str,'EeDd')
! 과학적 표기법이 있는데 마지막 문자가 숫자가 아닌 경우 (형식 오류)
if(ipos /= 0 .and. .not.is_numeric(str)) then
   num=0.0_r8
   return
end if
read(str,*,iostat=ios) num
if(ios/=0) num=0.0_r8
end function str_to_r8

!**********************************************************************

subroutine strvec_to_r4(strs, nums)
! Convert string vector to r4 vector
! 숫자이면 실수로 변환하고, 문자이면 0.0_r4로 변환
character(len=*), intent(in) :: strs(:)
real(kind=r4), intent(out) :: nums(:)
integer :: i, n

nums = 0.0_r4
n = min(size(strs), size(nums))
do i = 1, n
   nums(i) = str_to_r4(strs(i))
end do

end subroutine strvec_to_r4

!**********************************************************************

subroutine strvec_to_i4(strs, nums)
! Convert string vector to i4 vector
! 숫자이면 정수로 변환하고, 문자이면 0으로 변환
character(len=*), intent(in) :: strs(:)
integer(kind=ki4), intent(out) :: nums(:)
integer :: i, n

nums = 0_ki4
n = min(size(strs), size(nums))
do i = 1, n
   nums(i) = str_to_i4(strs(i))
end do

end subroutine strvec_to_i4

!**********************************************************************

subroutine parse(str,delims,args,nargs)
! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

character(len=*) :: str,delims
character(len=len_trim(str)) :: strsav
character(len=*),dimension(:) :: args
integer :: nargs,i,na,k,lenstr

strsav=str
if(delims /= achar(9)) call compact(str)
na=size(args)
do i=1,na
  args(i)=' '
end do  
nargs=0
lenstr=len_trim(str)
if(lenstr==0) return
k=0

do
   if(len_trim(str) == 0) exit
   nargs=nargs+1
   call split(str,delims,args(nargs))
   call removebksl(args(nargs))
end do   
str=strsav

end subroutine parse

!**********************************************************************

subroutine compact(str)
! 1. Converts multiple spaces and tabs to single spaces
! 2. Deletes control characters
! 3. Removes initial spaces.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str)):: outstr
integer :: i,lenstr,ich,k,isp

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  
  select case(ich)
  
    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1
      
    case(33:)      ! not a space, quote, or control character(over space)
      k=k+1
      outstr(k:k)=ch
      isp=0
      
  end select
  
end do

str=adjustl(outstr)

end subroutine compact

!**********************************************************************

subroutine removesp(str)

! Removes spaces, tabs, and control characters in string str

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr
integer :: i,lenstr,ich, k

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)    
    case(0:32)  ! space, tab, or control character
         cycle       
    case(33:)  
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!**********************************************************************

subroutine shiftstr(str,n)
 
! Shifts characters in in the string 'str' n positions
! Positive values denote a right shift
! Negative values denote a left shift
! Characters that are shifted off the end are lost
! Positions opened up by the shift are replaced by spaces.

character(len=*):: str
integer :: n
integer :: lenstr,nabs  

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
return

end subroutine shiftstr

!**********************************************************************

subroutine insertstr(str,strins,loc)

! Inserts the string 'strins' into the string 'str' at position 'loc'. 
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. 
! Trailing spaces of 'strins' are removed prior to insertion

character(len=*):: str,strins
character(len=len(str))::tempstr
integer:: loc,lenstrins

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!**********************************************************************

subroutine delsubstr(str,substr)

! Deletes first occurrence of substring 'substr' from string 'str' 
! Shifts characters left to fill hole
! Trailing spaces or blanks are not considered part of 'substr'.

character(len=*):: str,substr
integer :: ipos,lensubstr

lensubstr=len_trim(substr)
ipos=index(str,substr)
if(ipos==0) return
if(ipos == 1) then
   str=str(lensubstr+1:)
else
   str=str(:ipos-1)//str(ipos+lensubstr:)
end if   
return

end subroutine delsubstr

!**********************************************************************

subroutine delall(str,substr)

! Deletes all occurrences of substring 'substr' from string 'str' 
! Shifts characters left to fill holes

character(len=*):: str,substr
integer :: ipos,lensubstr

lensubstr=len_trim(substr)
do
   ipos=index(str,substr)
   if(ipos == 0) exit
   if(ipos == 1) then
      str=str(lensubstr+1:)
   else
      str=str(:ipos-1)//str(ipos+lensubstr:)
   end if
end do   
return

end subroutine delall

!**********************************************************************


subroutine match(str,ipos,imatch)

! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.

character(len=*) :: str
character :: delim1,delim2,ch
integer :: ipos,imatch,lenstr,istart,iend,inc,idelim2,isum
integer :: i

lenstr=len_trim(str)
delim1=str(ipos:ipos)
select case(delim1)
   case('(')
      idelim2=iachar(delim1)+1
      istart=ipos+1
      iend=lenstr
      inc=1
   case(')')
      idelim2=iachar(delim1)-1
      istart=ipos-1
      iend=1
      inc=-1
   case('[','{','<')
      idelim2=iachar(delim1)+2
      istart=ipos+1
      iend=lenstr
      inc=1
   case(']','}','>')
      idelim2=iachar(delim1)-2
      istart=ipos-1
      iend=1
      inc=-1
   case default
      write(*,*) delim1,' is not a valid delimiter'
      return
end select
if(istart < 1 .or. istart > lenstr) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if
delim2=achar(idelim2) ! matching delimiter

isum=1
do i=istart,iend,inc
   ch=str(i:i)
   if(ch /= delim1 .and. ch /= delim2) cycle
   if(ch == delim1) isum=isum+1
   if(ch == delim2) isum=isum-1
   if(isum == 0) exit
end do
if(isum /= 0) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if   
imatch=i

return

end subroutine match

!***********************************************************************

subroutine trimzero(str)

! Deletes nonsignificant trailing zeroes from number string str. 
! If number string ends in a decimal point, one trailing zero is added.

character(len=*) :: str
character :: ch
character(len=10) :: exp
integer:: i,lstr,ipos

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle          
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!**********************************************************************

function is_str(str) result(res)
! Returns .true. if str contains at least one alphabet or Hangul character
character(len=*), intent(in) :: str
logical :: res
integer :: idx
res = .false.
if(len_trim(str) == 0) then
   return
else
   do idx=1,len_trim(str)
      if(is_letter(str(idx:idx)) .or. is_hangul(str(idx:idx))) then
         res = .true.
         return
      end if
   end do
end if
end function is_str

!**********************************************************************

elemental function is_letter(ch) result(res)
! Returns .true. if ch is an English alphabet letter (A-Z or a-z)
character, intent(in) :: ch
logical :: res
integer :: iav
iav = iachar(ch)
res = ((iav >= iachar('A') .and. iav <= iachar('Z')) .or. &
       (iav >= iachar('a') .and. iav <= iachar('z')))
end function is_letter

elemental function is_digit(ch) result(res)
! Returns .true. if ch is a digit (0-9)
character, intent(in) :: ch
logical :: res
integer :: iav
iav = iachar(ch)
res = (iav >= iachar('0') .and. iav <= iachar('9'))
end function is_digit

!**********************************************************************

elemental function is_hangul(ch) result(res)
! Returns .true. if ch is likely a Hangul character (UTF-8 byte > 127)
character, intent(in) :: ch
logical :: res

res = (iachar(ch) > 127)
return

end function is_hangul

!**********************************************************************

function is_numeric(str) result(res)
! Returns .true. if str can be read as a numeric value (integer or real)
! OPTIMIZED: Single validation pass - removed redundant character checking
!
! Performance: Uses direct READ test which is most reliable and efficient
character(len=*), intent(in) :: str
logical :: res
character(len=len(str)) :: s
integer :: lenstr, ios
real(kind=r8) :: tmp

s = adjustl(str)
lenstr = len_trim(s)

if (lenstr == 0) then
   res = .false.
   return
end if

! Direct READ test - most reliable and efficient method
! No need for pre-validation since READ handles all valid numeric formats
read(s(1:lenstr), *, iostat=ios) tmp
res = (ios == 0)

end function is_numeric

!**********************************************************************

subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

character(len=*) :: str,delims,before
character,optional :: sep
logical :: pres
character :: ch,cha
integer :: lenstr,ipos,iposa,k,ibsl,i

pres=present(sep)
str=adjustl(str)
if(delims /= achar(9)) call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\') then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if
   ipos=index(delims,ch)         
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)
   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do
if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return

end subroutine split

!**********************************************************************

subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr
integer :: i,k,ibsl,lenstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive

do i=1,lenstr
  ch=str(i:i)
  if(ibsl == 1) then          ! backslash active
   k=k+1
   outstr(k:k)=ch
   ibsl=0
   cycle
  end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl


!**********************************************************************
subroutine remove_duplicate_spaces(input_str, output_str)
! 문자열에서 중복된 space를 하나로 만드는 서브루틴
! 연속된 여러 개의 space를 하나의 space로 축약
! Arguments:
!   input_str  - (in)  입력 문자열
!   output_str - (out) 처리된 출력 문자열

character(len=*), intent(in) :: input_str
character(len=*), intent(out) :: output_str
integer :: i, j, len_input
logical :: prev_space

len_input = len_trim(input_str)
output_str = ''
j = 1
prev_space = .false.

do i = 1, len_input
    if (input_str(i:i) == ' ') then
        if (.not. prev_space) then
            output_str(j:j) = ' '
            j = j + 1
            prev_space = .true.
        end if
    else
        output_str(j:j) = input_str(i:i)
        j = j + 1
        prev_space = .false.
    end if
end do

end subroutine remove_duplicate_spaces

!**********************************************************************
subroutine split_string(str, parts, num_parts, delimiter)
! 문자열을 delimiter를 중심으로 분리하여 배열에 저장
! OPTIMIZED: Uses allocatable array for better memory efficiency
!
! Arguments:
!   str       - (in)  입력 문자열
!   delimiter - (in)  delimiter 문자 (CHAR(9)=tab, CHAR(32)=space 등)
!   parts     - (out) 분리된 문자열들의 배열
!   num_parts - (out) 분리된 부분의 개수
!
! Performance: O(n) time complexity, memory-efficient with allocatable arrays

character(len=*), intent(in) :: str
character(len=1), intent(in) :: delimiter
character(len=*), intent(out) :: parts(:)
integer, intent(out) :: num_parts

integer :: i, str_len, delim_count
integer, allocatable :: delim_pos(:)  ! OPTIMIZED: allocatable for exact memory usage

parts(:) = ''
str_len = len_trim(str)

if (str_len == 0) then
    num_parts = 1
    parts(1) = ''
    return
end if

! Count delimiters (first pass)
delim_count = 0
do i = 1, str_len
    if (str(i:i) == delimiter) delim_count = delim_count + 1
end do

num_parts = delim_count + 1

! Handle no delimiter case
if (delim_count == 0) then
    parts(1) = trim(adjustl(str))
    return
end if

! Allocate exact size needed
allocate(delim_pos(delim_count))

! Store delimiter positions (second pass)
delim_count = 0
do i = 1, str_len
    if (str(i:i) == delimiter) then
        delim_count = delim_count + 1
        delim_pos(delim_count) = i
    end if
end do

! Extract parts
parts(1) = trim(adjustl(str(1:delim_pos(1)-1)))

do i = 1, size(delim_pos) - 1
    parts(i+1) = trim(adjustl(str(delim_pos(i)+1:delim_pos(i+1)-1)))
end do

parts(num_parts) = trim(adjustl(str(delim_pos(size(delim_pos))+1:str_len)))

deallocate(delim_pos)
end subroutine split_string

!**********************************************************************
 character(len=1) function FindDelim(str) result(dlm)
  character(len=*):: str
  character(len=len_trim(str)) :: str1_temp
  integer :: Len
  ! 문자열을 대문자로 변환
  str1_temp = to_upper(str)
  str1_temp = adjustl(str1_temp)
  Len = len_trim(str1_temp)
  select case(str1_temp(1:Len))
     case("SPACE","BLANK",CHAR(32))
        dlm=CHAR(32)
     case("TAB",CHAR(9))
        dlm=CHAR(9)
     case("COMMA",",")
        dlm=char(44)
     case("SEMICOLON",";")
        dlm=char(59)
     case("COLON",":")
        dlm=char(58)
     case("PIPE","|")
        dlm=char(124)
     case("PERIOD",".")
        dlm=char(46)
     case("-")
        dlm=char(45)    
     case("_")  
        dlm=char(95)
     case("/")
        dlm=char(47)
     case("\")
        dlm=char(92)
     case default
        ! 만약 첫 문자가 이미 delimiter라면 그것을 반환
        if (len_trim(str1_temp) > 0) then
           dlm = str(1:1)
        else
           dlm = CHAR(32)
        end if
  end select  
  end function FindDelim

! ========================================
! 파일 경로 처리 함수들
! ========================================

function is_absolute_path(path) result(is_abs)
  implicit none
  character(len=*), intent(in) :: path
  logical :: is_abs
  
  ! Unix/Linux 절대 경로: /로 시작
  if (path(1:1) == '/') then
    is_abs = .true.
  else if (len_trim(path) > 1) then
    ! Windows 절대경로: C:\, D:\ 등 (드물지만 지원)
    if (path(2:2) == ':') then
      is_abs = .true.
    else
      is_abs = .false.
    end if
  else
    is_abs = .false.
  end if
end function is_absolute_path

function get_directory_path(filepath) result(dirpath)
  implicit none
  character(len=*), intent(in) :: filepath
  character(len=MAX_STR) :: dirpath
  integer :: i, last_slash
  
  last_slash = 0
  do i = len_trim(filepath), 1, -1
    if (filepath(i:i) == '/') then
      last_slash = i
      exit
    end if
  end do
  
  if (last_slash > 0) then
    dirpath = filepath(1:last_slash-1)
  else
    ! 디렉토리 부분이 없으면 현재 디렉토리 반환
    dirpath = '.'
  end if
end function get_directory_path

function resolve_file_path(base_dir, file_path) result(resolved_path)
  implicit none
  character(len=*), intent(in) :: base_dir, file_path
  character(len=MAX_STR) :: resolved_path
  
  ! 절대 경로이면 그대로 반환
  if (is_absolute_path(file_path)) then
    resolved_path = trim(adjustl(file_path))
  else
    ! 상대 경로이면 기본 디렉토리를 기준으로 경로 구성
    ! base_dir이 '.'이면 현재 디렉토리
    if (trim(base_dir) == '.') then
      resolved_path = trim(adjustl(file_path))
    else
      resolved_path = trim(adjustl(base_dir))//'/'//trim(adjustl(file_path))
    end if
  end if
end function resolve_file_path

function to_logical(str) result(res)
    implicit none
    character(len=*), intent(in) :: str
    logical :: res
    character(len=len_trim(adjustl(str))) :: tmp_str, temp_str

    ! 앞뒤 공백 제거 및 소문자 변환 준비
    tmp_str = trim(adjustl(str))
    
    ! 대소문자 구분 없애기 (소문자로 통일)
    temp_str = to_lower(tmp_str)
    
    ! 논리값 판별
    select case (temp_str)
    case ('true', 't', 'yes', 'y', '1')
        res = .true.
    case ('false', 'f', 'no', 'n', '0')
        res = .false.
    case default
        res = .false. ! 기본값 (또는 에러 처리)
        print *, "Warning: Invalid logical string '", trim(str), "', returning .false."
    end select
end function to_logical

!**********************************************************************
! yyyymmdd 정수를 (year, month, day)로 분해하는 서브루틴
! bdate_int : yyyymmdd 형식의 정수 (예: 20210315)
! year, month, day : 분해된 연, 월, 일
subroutine yyyymmdd_to_date(bdate_int, year, month, day)
  integer, intent(in)  :: bdate_int
  integer, intent(out) :: year, month, day
  year  = bdate_int / 10000
  month = mod(bdate_int, 10000) / 100
  day   = mod(bdate_int, 100)
end subroutine yyyymmdd_to_date

!**********************************************************************
! yyyymmdd 정수를 율리우스 일수(Julian Day Number)로 변환하는 함수
! 그레고리안력 기준의 JDN 반환 (날짜 간격 계산에 사용)
! bdate_int <= 0이면 0 반환 (날짜 불명)
function bdate_to_jdn(bdate_int) result(jdn)
  integer, intent(in) :: bdate_int
  integer :: jdn
  integer :: yy, mm, dd, aa
  if (bdate_int <= 0) then
    jdn = 0
    return
  end if
  call yyyymmdd_to_date(bdate_int, yy, mm, dd)
  aa  = (14 - mm) / 12
  yy  = yy + 4800 - aa
  mm  = mm + 12 * aa - 3
  jdn = dd + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045
end function bdate_to_jdn

  ! ============================================================
  ! get_field_pos  — locate start/end of a delimited field
  !   Returns out_start/out_end byte positions (1-based) within
  !   in_line for the in_col_pos-th field.  No memory allocation.
  ! ============================================================
  subroutine get_field_pos(in_line, in_line_len, in_delim, in_col_pos, &
                            out_start, out_end, out_found)
    character(len=*), intent(in)  :: in_line
    integer,          intent(in)  :: in_line_len
    character(len=1), intent(in)  :: in_delim
    integer,          intent(in)  :: in_col_pos
    integer,          intent(out) :: out_start, out_end
    logical,          intent(out) :: out_found

    integer :: i, col_idx, start_pos
    logical :: collapse_ws

    out_found = .false.; out_start = 0; out_end = -1
    if (in_col_pos <= 0 .or. in_line_len <= 0) return

    collapse_ws = (in_delim == ' ' .or. in_delim == achar(9))

    if (collapse_ws) then
      i = 1; col_idx = 0
      do while (i <= in_line_len)
        do while (i <= in_line_len .and. &
                  (in_line(i:i) == ' ' .or. in_line(i:i) == achar(9)))
          i = i + 1
        end do
        if (i > in_line_len) exit
        col_idx   = col_idx + 1
        start_pos = i
        do while (i <= in_line_len .and. &
                  in_line(i:i) /= ' ' .and. in_line(i:i) /= achar(9))
          i = i + 1
        end do
        if (col_idx == in_col_pos) then
          out_start = start_pos; out_end = i - 1; out_found = .true.; return
        end if
      end do
    else
      col_idx = 1; start_pos = 1
      do i = 1, in_line_len + 1
        if (i > in_line_len .or. in_line(i:i) == in_delim) then
          if (col_idx == in_col_pos) then
            out_start = start_pos; out_end = i - 1; out_found = .true.; return
          end if
          col_idx   = col_idx + 1
          start_pos = i + 1
        end if
      end do
    end if
  end subroutine get_field_pos

end module M_StrEdit


