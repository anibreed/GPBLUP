module M_param
  use M_Kinds, only: ki1, ki2, ki4, ki8, r4, r8, r16
  use M_Variables, only: missI, ZERO, LEN_KEY, LEN_STR, MAX_RECL, MAX_VAR, MAX_STR, &
                         missR, STARTIME, FINISHTIME, RECORD, missC, nullC, missA, TAB
  use M_Stamp, only: timestamp, timelock
  implicit none

  integer, parameter, public :: MAX_INTVAR = MAX_VAR
  integer, parameter, public :: MAX_NT = 20
  integer, parameter, public :: MAXANC_LIMIT = 20000

  public :: JULDAY, TP_JDN

contains

  function JULDAY(year,month,day) result(jul)
    integer, intent(in) :: year, month, day
    integer :: iflg
    real(kind=r8) :: jul, jdn

    jul = 0.0_r8
    call TP_JDN(year, month, day, jdn, iflg)
    if (iflg == 0) jul = jdn
  end function JULDAY

  subroutine TP_JDN(year,month,day,jdn,iflg)
    integer :: TLU(12)=(/306,337,0,31,61,92,122,153,184,214,245,275/)
    integer, intent(in) :: year,month,day
    real(kind=r8), intent(out) :: jdn
    integer, intent(out) :: iflg
    integer :: yyyy

    iflg = 0
    if ((month < 1) .or. (month > 12)) then
      iflg = -1
      return
    elseif ((day < 1) .or. (day > 31)) then
      iflg = -2
      return
    endif

    yyyy = year
    if (month < 3) yyyy = yyyy - 1

    jdn = day + TLU(month) + 365.0_r8*yyyy &
        + floor(0.25_r8*yyyy) &
        - floor(0.01_r8*yyyy) &
        + floor(0.0025_r8*yyyy) + 1721118.5_r8
  end subroutine TP_JDN

end module M_param

