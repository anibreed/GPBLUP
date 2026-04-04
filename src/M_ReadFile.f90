module M_ReadFile
use M_Kinds
use M_Variables
use M_StrEdit
implicit none
public :: fopen, N_recf, readline

interface readline
  module procedure readline_int, readline_real
end interface

contains

 SUBROUTINE readline_int(unt,n,XXC,XX,dlm_str)
  INTEGER:: unt,n
  CHARACTER(LEN=*) :: XXC(MAX_VAR)
  integer(kind=ki4) :: XX(MAX_VAR)
  CHARACTER(LEN=*),intent(in),optional:: dlm_str
  CHARACTER(LEN=1):: dlm1
  INTEGER :: iio, ii, ic, ichar_val
  CHARACTER(len=MAX_RECL) :: a
  CHARACTER(len=MAX_RECL) :: a_compact
  INTEGER :: len_a
   dlm1 = CHAR(32)
   if(present(dlm_str)) dlm1=FindDelim(dlm_str)
   n=0
   DO WHILE(n == 0)
     READ(unt,'(a)',IOSTAT=iio) a
     IF (iio /= 0) THEN
         n=-1
         RETURN
     ENDIF
     a = trim(ADJUSTL(a))
     
     ! 제어 문자 제거 (CR, LF 등만) - Record 끝에서만 연속으로 제거
     ! 현재 delimiter를 제외한 제어문자만 제거
     len_a = len_trim(a)
     do ic = len_a, 1, -1
        if (len_a > 0) then
           ichar_val = ichar(a(ic:ic))
           ! CR(13), LF(10), NUL(0) 등 제거, 하지만 delimiter는 유지
           if ((ichar_val == 13 .or. ichar_val == 10 .or. ichar_val == 0) &
               .and. char(ichar_val) /= dlm1) then
              len_a = ic - 1
           else
              exit  ! 제거할 제어 문자가 아니면 중단
           end if
        end if
     end do
     if(len_a > 0) then
        a = a(1:len_a)
     else
        a = ''
     end if
     
     ! # 주석 처리
     ii = SCAN(a,'#')
     IF (ii == 1 .or. len_a == 0) then
       n = 0
     else  
       if(ii > 0) then
         len_a = ii - 1
         a = a(1:len_a)
       end if
       n = len_a
     endif  
   ENDDO
   
   ! delimiter 처리
   if(dlm1 == CHAR(32)) then
      call remove_duplicate_spaces(a, a_compact)
      a = a_compact
   end if
   
   ! 디버그 출력 (필요시 주석 해제)
   ! print '(A,I0,A,I0)', 'readline_int: n=', n, ', dlm=', ichar(dlm1)
   
   ! 문자열 분리 및 변환
   call split_string(a, XXC, n, dlm1)
   if(n > 0) then
      call strvec_to_I4(XXC(1:n), XX(1:n))
   else
      XX = 0_ki4
   end if
 END SUBROUTINE readline_int

 SUBROUTINE readline_real(unt,n,XXC,XX,dlm_str)
   INTEGER:: unt,n
   CHARACTER(LEN=*) :: XXC(MAX_VAR)
   real(kind=r4) :: XX(MAX_VAR)
   CHARACTER(LEN=*),intent(in),optional:: dlm_str
   CHARACTER(LEN=1):: dlm1
   INTEGER :: iio, ii, ic, ichar_val
   CHARACTER(len=MAX_RECL) :: a
   CHARACTER(len=MAX_RECL) :: a_compact
   INTEGER :: len_a
   dlm1 = CHAR(32)
   
   if(present(dlm_str)) dlm1=FindDelim(dlm_str)
   n=0
   DO WHILE(n == 0)
     READ(unt,'(a)',IOSTAT=iio) a
     IF (iio /= 0) THEN
         n=-1
         RETURN
     ENDIF
     a = trim(ADJUSTL(a))
     
     ! 제어 문자 제거 (CR, LF 등만) - Record 끝에서만 연속으로 제거
     ! 현재 delimiter를 제외한 제어문자만 제거
     len_a = len_trim(a)
     do ic = len_a, 1, -1
        if (len_a > 0) then
           ichar_val = ichar(a(ic:ic))
           ! CR(13), LF(10), NUL(0) 등 제거, 하지만 delimiter는 유지
           if ((ichar_val == 13 .or. ichar_val == 10 .or. ichar_val == 0) &
               .and. char(ichar_val) /= dlm1) then
              len_a = ic - 1
           else
              exit  ! 제거할 제어 문자가 아니면 중단
           end if
        end if
     end do
     if(len_a > 0) then
        a = a(1:len_a)
     else
        a = ''
     end if
     
     ! # 주석 처리
     ii = SCAN(a,'#')
     IF (ii == 1 .or. len_a == 0) then
       n = 0
     else  
       if(ii > 0) then
         len_a = ii - 1
         a = a(1:len_a)
       end if
       n = len_a
     endif  
   ENDDO
   
   ! delimiter 처리
   if(dlm1 == CHAR(32)) then
     call remove_duplicate_spaces(a, a_compact)
     a = a_compact
   end if
   
   ! 디버그 출력 (필요시 주석 해제)
   ! print '(A,I0,A,I0)', 'readline_real: n=', n, ', dlm=', ichar(dlm1)
   
   ! 문자열 분리 및 변환
   call split_string(a, XXC, n, dlm1)
   if(n > 0) then
      call strvec_to_r4(XXC(1:n), XX(1:n))
   else
      XX = 0.0_r4
   end if
 END SUBROUTINE readline_real

integer function N_recf(filename, rec_size) result(rec_N)
 character(LEN=*),intent(in)::filename
 integer,intent(in),optional:: rec_size
 character(LEN=512):: REC
 integer::ios,unitF

 if(present(rec_size)) then
    unitF=fopen(filename, rec_size)
 else
    unitF=fopen(filename)
 endif
 rec_N=0
 if (unitF < 0) return
      ! Loop through input file to count records
      do
         read(unitF, '(A)', iostat=ios) REC
         if (ios /= 0) exit

         ! Skip empty lines and comment lines (starting with #)
         if (len_trim(REC) > 0 .and. REC(1:1) /= '#') then
            rec_N = rec_N + 1
         end if
      end do
         close(unit=unitF)
 end function N_recf

integer function fopen(filename, recl)
      character(len=*), intent(in) :: filename
      integer, intent(in), optional :: recl
      logical :: file_exists, file_opened
      integer :: rec_len

      rec_len = MAX_RECL
      if (present(recl)) rec_len = recl

      inquire(file=trim(filename), exist=file_exists, opened=file_opened, number=fopen)
      if (file_opened) close(unit=fopen)

      fopen = 10
      do
         fopen = fopen + 1
         inquire(unit=fopen, opened=file_opened)
         if (.not. file_opened) exit
      end do

      if (file_exists) then
         open(unit=fopen, file=trim(filename), recl=rec_len, status='unknown')
      else
         open(unit=fopen, file=trim(filename), recl=rec_len, status='new')
      endif
end function fopen

end module M_ReadFile
