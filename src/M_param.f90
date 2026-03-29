module M_param
INTEGER, PARAMETER :: ki1      = SELECTED_INT_KIND(1)  ! 1 Byte  integer : range( -127 ~ 127)
INTEGER, PARAMETER :: ki2      = SELECTED_INT_KIND(4)  ! 2 Byte integer :  range( -32767 ~ 32767)
INTEGER, PARAMETER :: ki4      = SELECTED_INT_KIND(9)  ! single precision integer : range(-2147483647 ~ 2147483647)
INTEGER, PARAMETER :: ki8      = SELECTED_INT_KIND(18) ! double precision integer
INTEGER, PARAMETER :: r4      = SELECTED_REAL_KIND( 6, 37 )
INTEGER, PARAMETER :: r8      = SELECTED_REAL_KIND( 15, 307 )
INTEGER, PARAMETER :: r16    = SELECTED_REAL_KIND( 18, 4931 )

! define parameter for raw data and pedigree
INTEGER,PARAMETER,PRIVATE:: lock_y=2028, lock_m=12, lock_d=31
INTEGER, PARAMETER, PUBLIC:: missI = 0
INTEGER, PARAMETER, PUBLIC:: ZERO=missI
INTEGER, PARAMETER, PUBLIC:: LEN_KEY = 20
INTEGER, PARAMETER, PUBLIC:: LEN_STR = LEN_KEY
INTEGER, PARAMETER, PUBLIC:: MAX_SNP = 100000
INTEGER, PARAMETER, PUBLIC:: MAX_RECL = 512
INTEGER, PARAMETER, PUBLIC:: MAX_NT = 20
INTEGER, PARAMETER, PUBLIC:: MAX_VAR = 128
INTEGER, PARAMETER, PUBLIC:: MAX_INTVAR = MAX_VAR
INTEGER, PARAMETER, PUBLIC:: MAX_STR=MAX_VAR
REAL, PARAMETER, PUBLIC:: missR = 0.d0
REAL,PUBLIC:: STARTIME, FINISHTIME
CHARACTER(LEN=MAX_RECL), PUBLIC :: RECORD
CHARACTER(len=1), PARAMETER, PUBLIC:: missC = achar(32)  ! " "
CHARACTER(len=1), PARAMETER, PUBLIC:: nullC = missC      ! " "
CHARACTER(len=1), PARAMETER, PUBLIC:: SPACE = missC      ! " "
CHARACTER(len=1), PARAMETER, PUBLIC:: missA = achar(48)  ! "0"
CHARACTER(len=1), PARAMETER, PUBLIC:: TAB = achar(9)     ! "TAB"

 contains
  subroutine timestamp (io_u)
  integer,optional::io_u
  integer::io_t
  character ( len = 40 ) string
  io_t = 6
  if(present(io_u)) io_t = io_u
  call timestring ( string )
  write (io_t, '(a)' ) trim ( string )
  end subroutine

  subroutine timestring ( string )
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if
  write ( string, '(1x,a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, trim ( ampm )
  end subroutine

  logical function timelock()
  character ( len = 8 ) date
  character ( len = 10 ) time
  character ( len = 5 ) zone
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer y,m,d, leftday
  integer values(8)
  real(kind=8):: expired, currentd

  call timestamp()
  timelock=.false.

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)

  expired=JULDAY(lock_y,lock_m,lock_d)
  currentd=JULDAY(y,m,d)

  leftday=int(expired-currentd)

  if(leftday.LE.60) then
    print*
    print '(a33,a,a1,i3,a1,i5)','This program will be expired at ', &
           trim(month(lock_m)),',',lock_d,',',lock_y    
    if(leftday.GE.0) then
       print '(i4,a20)',leftday,' days will be left'
    else
       print '(i4,a20)',abs(leftday),' days is passed over'
    endif
    print*
    print '(a)', ' Please update new version - contact with Dr. DEUKHWAN LEE'
    print '(a)', '                    (Hankyong National Univ.)'
    print '(a)', '                    (TEL: 031-670-5091, Mobile: 010-9320-5091)'
    print '(a)', '                    (E-mail: dhlee@hknu.ac.kr)'
    print*

    if(leftday.le.0) timelock=.true.
  endif
 end function

 FUNCTION JULDAY(year,month,day)
 integer:: year, month, day
 integer:: iflg
 real(KIND=r8):: JULDAY,jdn
 JULDAY=0.d0
 call TP_JDN(year,month,day,jdn,iflg)
 if(iflg.eq.-1) then
    print*,'Month was wrong'
 elseif(iflg.eq.-2) then
    print*,'Day was wrong'
 else
    JULDAY=jdn
 endif
 end function

 SUBROUTINE TP_JDN(year,month,day,jdn,iflg)
 !     INPUT:
 !         year: 4-digit number                  [I4]
 !        month: month of the year (1-12)        [I4]
 !          day: day of month (1-31)             [I4]
 !
 !     OUTPUT:
 !          jdn: julian day number               [R8]
 !         iflg: status flag                     [I4]
 !               0: no problems encountered
 !              -1: invalid month
 !              -2: invalid day
 !TLU = AINT ((153*M-457)/5) where M = the month (1-12)

 INTEGER :: TLU(12)=(/306,337,0,31,61,92,122,153,184,214,245,275/)
 INTEGER,INTENT(IN):: year,month,day
 REAL(kind=r8),INTENT(OUT):: jdn
 INTEGER,INTENT(OUT):: IFLG
 INTEGER:: yyyy
 IFLG=0
 IF ((month.LT.1).OR.(month.GT.12)) THEN
   IFLG=-1
   RETURN
 ELSEIF ((day.LT.1).OR.(day.GT.31)) THEN
   IFLG=-2
   RETURN
 ENDIF
 yyyy=year
 IF (month.LT.3) THEN
   yyyy=yyyy-1
 ENDIF
 jdn = Day + TLU(month)  + 365.0D+00*yyyy            &
                         + FLOOR(0.25D+00*yyyy)      &
                         - FLOOR(0.01D+00*yyyy)      &
                         + FLOOR(0.0025D+00*yyyy)    &
                         + 1721118.5D+00
 RETURN
 END SUBROUTINE

end module

