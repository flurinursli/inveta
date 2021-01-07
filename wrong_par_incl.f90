flag = .false.

IF (PRESENT(strict)) flag = strict

IF (PRESENT(v1) .and. PRESENT(v2)) THEN
  IF (flag) THEN
    wrong_par = ANY(par .le. v1) .or. ANY(par .ge. v2)
  ELSE
    wrong_par = ANY(par .lt. v1) .or. ANY(par .gt. v2)
  ENDIF

ELSEIF (PRESENT(v1) .and. .not.PRESENT(v2)) THEN
  IF (flag) THEN
    wrong_par = ANY(par .le. v1)
  ELSE
    wrong_par = ANY(par .lt. v1)
  ENDIF

ELSE
  IF (flag) THEN
    wrong_par = par(1) .ge. par(2)
  ELSE
    wrong_par = par(1) .gt. par(2)
  ENDIF
ENDIF
