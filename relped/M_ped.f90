module M_ped
use M_param
implicit none
 integer, public, parameter:: Male=2, Female=1

type, public :: CPED
     character(len=LEN_STR):: ASD(3)
     integer::  NASD(3)
     integer:: SEX, BYR, BDate, GEN, CG, OBS, NMat
     real(kind=r8):: F
     real(kind=r8):: EBV
end type CPED

type, public :: VPED
     character(len=LEN_STR):: ID
     integer::  NID
end type VPED

#define FFH_KEY_TYPE character(len=LEN_STR)
#define FFH_KEY_IS_STRING
#define FFH_VAL_TYPE type(CPED)
#include "./fhash/hash_inc.f90"
end module M_ped

