! Convert a real value to a string with given number of digits
SUBROUTINE real2string(value,ndigits_1,ndigits_2,string,error)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: value      ! real value
INTEGER,          INTENT(IN)  :: ndigits_1  ! number of digits BEFORE decimal point
INTEGER,          INTENT(IN)  :: ndigits_2  ! number of digits AFTER  decimal point
CHARACTER(LEN=*), INTENT(OUT) :: string     ! character string
LOGICAL,          INTENT(OUT) :: error      ! conversion error

INTEGER :: required_length,lf,nd2
CHARACTER(LEN=20) :: format
!EXTERNAL :: rmsp

error=.true.
string=' '
nd2=MAX(ndigits_2,0)

required_length = ndigits_1 + nd2 + 1
IF(required_length > LEN(string)) RETURN

WRITE(format,100,ERR=10) required_length,nd2
100 FORMAT('(F',I3,'.',I3,')')
CALL rmsp(format,lf)
WRITE(string,FMT=format,ERR=10) value
IF(INDEX(string,'*') == 0) error=.false.

10 CONTINUE

END SUBROUTINE real2string
