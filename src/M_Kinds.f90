module M_Kinds
INTEGER, PARAMETER :: ki1      = SELECTED_INT_KIND(1)  ! 1 Byte  integer : range( -127 ~ 127)
INTEGER, PARAMETER :: ki2      = SELECTED_INT_KIND(4)  ! 2 Byte integer :  range( -32767 ~ 32767)
INTEGER, PARAMETER :: ki4      = SELECTED_INT_KIND(9)  ! single precision integer : range(-2147483647 ~ 2147483647)
INTEGER, PARAMETER :: ki8      = SELECTED_INT_KIND(18) ! double precision integer
INTEGER, PARAMETER :: r4      = SELECTED_REAL_KIND( 6, 37 )
INTEGER, PARAMETER :: r8      = SELECTED_REAL_KIND( 15, 307 )
INTEGER, PARAMETER :: r16    = SELECTED_REAL_KIND( 18, 4931 )
END MODULE M_Kinds