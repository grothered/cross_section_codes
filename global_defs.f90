MODULE global_defs
! This module contains various constants which we would like all other modules
! to see --- i.e. things which are not going to change. 
IMPLICIT NONE

! Define the type of real used by other routines
INTEGER, PARAMETER, PUBLIC  :: dp = SELECTED_REAL_KIND(12, 60)

! Define pi
REAL(dp), PARAMETER, PUBLIC ::  pi=atan(1._dp)*4._dp

! Define a standard length of character strings
INTEGER, PARAMETER, PUBLIC :: char_len=20

END MODULE global_defs

