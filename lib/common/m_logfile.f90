MODULE m_logfile

  ! Purpose:
  !
  !
  !   Routines are not threadsafe.
  !
  !   To compile:
  !       gfortran -c -O2 -march=native precisions.f90 llsq.f90 -cpp (-DDOUBLE_PREC, -DMKL) (-fopenmp)
  !       ifort -c -O3 -xHost precisions.f90 llsq.f90 -cpp (-DOUBLE_PREC, -DMKL) (-qopenmp)
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   02/09/20                  original version
  !

  USE, INTRINSIC     :: iso_fortran_env, stdout => output_unit, stderr => error_unit
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: set_log_module, report_error, update_log, char_cr

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CHARACTER(LEN=1),                           PARAMETER :: char_cr = ACHAR(13)

  INTEGER(i32),     ALLOCATABLE, DIMENSION(:)           :: log_units

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --


  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE set_log_module(ok, disk, screen)

      ! Purpose:
      !   to setup the log-file environment. Log messages can be reported to disk and/or to screen.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                       INTENT(OUT) :: ok
      LOGICAL,                  OPTIONAL, INTENT(IN)  :: disk, screen
      CHARACTER(8)                                    :: date
      CHARACTER(10)                                   :: time
      CHARACTER(:), ALLOCATABLE                       :: fo, timestamp

      !-----------------------------------------------------------------------------------------------------------------------------

      ! open a formatted log-file and add it to "log_units"
      IF (PRESENT(disk)) THEN

         IF (disk .eqv. .true.) THEN

            CALL DATE_AND_TIME(date, time)

            timestamp = date(7:8) // date(5:6) // date(3:4) // '_' // time(1:2) // time(3:4) // time(5:6)

            fo = 'logfile_' // timestamp // '.txt'

            OPEN(stderr, file = fo, status = 'unknown', form = 'formatted', iostat = ok)

            IF (ok .ne. 0) THEN
              WRITE(stdout, *) 'Error: it was not possible to open log file ' // fo
              RETURN
            ENDIF

            log_units = [stderr]

         ENDIF

      ENDIF

      ! add screen to "log_units"
      IF (PRESENT(screen)) THEN
         IF (screen .eqv. .true.) THEN
            IF (ALLOCATED(log_units)) THEN
              log_units = [stdout, log_units]
            ELSE
              log_units = [stdout]
            ENDIF
         ENDIF
      ENDIF

    END SUBROUTINE set_log_module

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE report_error(msg)

      ! Purpose:
      !   to convert an error code into an error message based on the provided external funtion and to report content to selected
      !   units.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*), INTENT(IN) :: msg
      INTEGER(i32)             :: i, n

      !-----------------------------------------------------------------------------------------------------------------------------

      n = 0

      IF (ALLOCATED(log_units)) n = SIZE(log_units)

      DO i = 1, n
        WRITE(log_units(i), '(1X, A)') msg
      ENDDO

    END SUBROUTINE report_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE update_log(text, sameline, blankline)

      ! Purpose:
      !   to add a text to the log output units. By default every call includes an introductory blank line (use "blankline=.false."
      !   otherwise). Input text is split upon detection of a carriage return character and each single message is then reported on
      !   a new line (use "sameline=.true." otherwise).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),                       INTENT(IN) :: text
      LOGICAL,                  OPTIONAL, INTENT(IN) :: sameline, blankline
      CHARACTER(:), ALLOCATABLE                      :: adv, buffer
      INTEGER(i32)                                   :: i, j, n

      !-----------------------------------------------------------------------------------------------------------------------------

      adv = 'yes'

      IF (PRESENT(sameline)) THEN
        IF (sameline .eqv. .true.) adv = 'no'
      ENDIF

      n = 0

      IF (ALLOCATED(log_units)) n = SIZE(log_units)

      ! update all log units
      DO i = 1, n

        IF (PRESENT(blankline)) THEN
          IF (blankline .eqv. .true.) WRITE(log_units(i), '(A)')
        ELSE
          WRITE(log_units(i), '(A)')
        ENDIF

        j = 1

        ! loop over input text, split based on carriage return character
        DO

          buffer = split(text, char_cr, j)

          ! exit if no more chunks are available
          IF (buffer .eq. "") EXIT

          WRITE(log_units(i), '(A)', ADVANCE = adv) buffer

          j = j + 1

        ENDDO

      ENDDO

      DEALLOCATE(adv)

    END SUBROUTINE update_log

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_logfile
