MODULE m_inveta

  ! Purpose:
  !   To provide a series of routines used by program INVETA for the inversion of observed seismogram envelopes, including function
  !   "misfit" called by the inversion algorithm NA and representing the objective function.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_rtt
  USE, NON_INTRINSIC :: m_llsq_r32
  USE, NON_INTRINSIC :: m_logfile
  USE                :: mpi
  USE                :: omp_lib

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE wrong_par
    MODULE PROCEDURE wrong_par_i32, wrong_par_r32
  END INTERFACE wrong_par

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CHARACTER(:), ALLOCATABLE                 :: acf                                        !< autocorrelation function (scalar RTT)
  INTEGER(i32)                              :: world_size, world_rank
  INTEGER(i32)                              :: comm1, comm2, comm3
  INTEGER(i32)                              :: mode                                       !< inversion mode
  INTEGER(i32)                              :: threshold                                  !< minimum number of receivers for event
  INTEGER(i32)                              :: errors
  INTEGER(i32)                              :: nsi, ns, nr, itermax, seed                 !< NA parameters
  INTEGER(i32), ALLOCATABLE, DIMENSION(:)   :: iobs, nobs                                 !< observations per recording
  INTEGER(i32), ALLOCATABLE, DIMENSION(:)   :: pprank1, pprank2, pprank3
  INTEGER(i32), ALLOCATABLE, DIMENSION(:,:) :: gs, ge                                     !< global indices for cartesian topology
  LOGICAL                                   :: elastic                                    !< scalar or elastic RTT
  LOGICAL                                   :: noweight                                   !< disable weighted linear regression
  REAL(r32)                                 :: drespl
  REAL(r32)                                 :: beta
  REAL(r32)                                 :: fwin
  REAL(r32)                                 :: tlim
  REAL(r32)                                 :: pdwindow, sdwindow, pcwindow, scwindow     !< windows for direct and coda P-/S-waves
  REAL(r32),                 DIMENSION(7)   :: error_params
  REAL(r32),    ALLOCATABLE, DIMENSION(:)   :: etass, etass2pp, etaps2pp                  !< scattering parameters
  REAL(r32),    ALLOCATABLE, DIMENSION(:)   :: nu, hurst
  REAL(r32),    ALLOCATABLE, DIMENSION(:)   :: tobs, envobs, tpobs, tsobs                 !< observables used during inversion
  REAL(r32),    ALLOCATABLE, DIMENSION(:,:) :: fbands                                     !< frequency bands for inversion
  REAL(r64)                                 :: coefficient
  REAL(r64),                 DIMENSION(2)   :: cf

  TYPE :: info
    CHARACTER(200)                           :: folder
    CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: code, network, channel
    CHARACTER(30), ALLOCATABLE, DIMENSION(:) :: event
    REAL(r32),     ALLOCATABLE, DIMENSION(:) :: tp, ts
  END TYPE info

  TYPE(info), ALLOCATABLE, DIMENSION(:) :: recvr

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_search(f, code, inverted, mpar, fit)

      ! Purpose:
      !   to write to disk the "mpar" parameters space explored by NA for the "f"-th frequency band, including misfit values "fit"
      !   averaged over the number of observations.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      INTEGER(i32),                           INTENT(IN) :: f
      CHARACTER(10),            DIMENSION(:), INTENT(IN) :: code
      CHARACTER(30),            DIMENSION(:), INTENT(IN) :: inverted
      REAL(r32),                DIMENSION(:), INTENT(IN) :: mpar, fit
      CHARACTER(:), ALLOCATABLE                          :: fo
      INTEGER(i32)                                       :: nd, n, i, lu, ok, lf, hf
      REAL(r32)                                          :: avg

      !-----------------------------------------------------------------------------------------------------------------------------

      lf = NINT(fbands(1, f))
      hf = NINT(fbands(2, f))

      SELECT CASE (mode)
        CASE(0)
          fo = 'nasearch_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(1)) + '_' + TRIM(inverted(1)) + '.txt'
        CASE(1)
          fo = 'nasearch_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(1)) + '.txt'
        CASE(2)
          fo = 'nasearch_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(inverted(1)) + '.txt'
      END SELECT

      OPEN(newunit = lu, file = fo, status = 'unknown', form = 'formatted', access = 'sequential', action = 'write', iostat = ok)

      nd = SIZE(mpar) / SIZE(fit)
      n  = 0

      avg = SUM(nobs)

      DO i = 1, SIZE(mpar), nd
        n = n + 1
        WRITE(lu, *) mpar(i:i + nd - 1), fit(n) / avg            !< add average misfit value
      ENDDO

      CLOSE(lu)

    END SUBROUTINE write_search

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE bestfit(f, code, inverted, bestmodel)

      ! Purpose:
      !   to write to disk the best-fitting model as obtained from the inversion. Attenuation parameter "b" and amplitude factors for
      !   each set of observables are computed by solving the associated linear system.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      INTEGER(i32),                                                 INTENT(IN) :: f
      CHARACTER(10),              DIMENSION(:),                     INTENT(IN) :: code
      CHARACTER(30),              DIMENSION(:),                     INTENT(IN) :: inverted
      REAL(r32),                  DIMENSION(:),                     INTENT(IN) :: bestmodel
      CHARACTER(:),  ALLOCATABLE                                               :: fo
      INTEGER(i32)                                                             :: lu, lf, hf, i, j, k, l, n, p, is, ie, j0, j1, ok
      INTEGER(i32)                                                             :: rank, ierr, req
      INTEGER(i32),               DIMENSION(0:SIZE(pprank2)-1)                 :: displs, pprank
      REAL(r32)                                                                :: gpp, gps, gsp, gss, gi, aksq, kappa, t, const
      REAL(r32),                                                    PARAMETER  :: tau = 0.25_r32, wp = 1._r32, ws = 23.4_r32
      REAL(r32),     ALLOCATABLE, DIMENSION(:)                                 :: time, envelope
      REAL(r32),                  DIMENSION(SUM(nobs))                         :: delta, weight, b, tobs
      REAL(r32),                  DIMENSION(SUM(nobs),SIZE(nobs)+1)            :: a
      REAL(r64),                  DIMENSION(2)                                 :: tictoc

      !-----------------------------------------------------------------------------------------------------------------------------

      lf = NINT(fbands(1, f))
      hf = NINT(fbands(2, f))

      IF (elastic) THEN
        gss = bestmodel(1)
        gpp = gss / bestmodel(2)
        gps = gpp * bestmodel(3)
        gsp = gps / 6._r32
      ELSE
        gss   = bestmodel(1)
        aksq  = bestmodel(2)          !< nu
        kappa = bestmodel(3)          !< hurst
      ENDIF

      gi = 0._r32

#ifdef PERF
      CALL watch_start(tictoc(1), comm2)
#endif

#include "linsys_incl.f90"

      n = SIZE(nobs)

      ! divide "b" by velocity to obtain intrinsic attenuation
      gi = b(n + 1) / beta

      IF (world_rank .eq. 0) THEN

        ! report best fitting parameters to a file
        IF (mode .le. 1) THEN

          ! fo = 'bestpar_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(1))  + '.txt'
          fo = 'bestpar_' + TRIM(code(1)) + '.txt'

          OPEN(newunit = lu, file = fo, status = 'unknown', form = 'formatted', access = 'sequential', position = 'append',    &
               action = 'write', iostat = ok)

          IF (elastic) THEN
            WRITE(lu, *) mean(fbands(:,f)), gpp, gps, gsp, gss, gi * beta
          ELSE
            IF ( (nu(1) .eq. nu(2)) .and. (hurst(1) .eq. hurst(2)) ) THEN
              WRITE(lu, *) mean(fbands(:,f)), gss, gi * beta, gm(gss, aksq, kappa)
            ELSEIF (hurst(1) .eq. hurst(2)) THEN
              WRITE(lu, *) mean(fbands(:,f)), gss, aksq, gi * beta, gm(gss, aksq, kappa)
            ELSE
              WRITE(lu, *) mean(fbands(:,f)), gss, aksq, kappa, gi * beta, gm(gss, aksq, kappa)
            ENDIF
          ENDIF

          CLOSE(lu)

        ELSE

          DO j = 1, SIZE(nobs)

            ! fo = 'bestpar_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(j)) + '.txt'
            fo = 'bestpar_' + TRIM(code(j)) + '.txt'

            OPEN(newunit = lu, file = fo, status = 'unknown', form = 'formatted', access = 'sequential', position = 'append',    &
                 action = 'write', iostat = ok)

            IF (elastic) THEN
              WRITE(lu, *) mean(fbands(:,f)), gpp, gps, gsp, gss, gi * beta
            ELSE
              IF ( (nu(1) .eq. nu(2)) .and. (hurst(1) .eq. hurst(2)) ) THEN
                WRITE(lu, *) mean(fbands(:,f)), gss, gi * beta, gm(gss, aksq, kappa)
              ELSEIF (hurst(1) .eq. hurst(2)) THEN
                WRITE(lu, *) mean(fbands(:,f)), gss, aksq, gi * beta, gm(gss, aksq, kappa)
              ELSE
                WRITE(lu, *) mean(fbands(:,f)), gss, aksq, kappa, gi * beta, gm(gss, aksq, kappa)
              ENDIF
            ENDIF

            CLOSE(lu)

          ENDDO

        ENDIF

      ENDIF

#ifdef PERF
      CALL watch_stop(tictoc(1), comm2)
      CALL watch_start(tictoc(2), comm2)
#endif

      CALL mpi_comm_rank(comm2, rank, ierr)

      displs(0) = 0

      DO i = 1, SIZE(pprank2) - 1
        displs(i) = displs(i - 1) + pprank2(i - 1)
      ENDDO

      j0 = displs(rank) + 1
      j1 = j0 + pprank2(rank) - 1

      CALL mpi_comm_rank(comm1, rank, ierr)

      ! write also envelopes associated to best fitting parameters
      DO j = j0, j1

        n = iobs(j)

        ALLOCATE(time(n), envelope(n))

        DO i = 1, n
          time(i) = (i - 1) * drespl
        ENDDO

        IF (elastic) THEN
          CALL rtt(comm3, pprank3, time, tpobs(j) + tau, tsobs(j) + tau, gpp, gps, gsp, gss, gi, beta, wp, ws, tau, envelope, ok)
        ELSE
          CALL rtt(comm3, pprank3, time, tsobs(j) + tau, gss, gi, beta, acf, kappa, aksq, tau, envelope, ok)
        ENDIF

        IF (mode .le. 1) THEN
          fo = 'bestfit_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(1)) + '_' + TRIM(inverted(j)) + '.txt'
        ELSE
          fo = 'bestfit_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(j)) + '_' + TRIM(inverted(1)) + '.txt'
        ENDIF

        IF (rank .eq. 0) THEN

          OPEN(newunit = lu, file = fo, status = 'unknown', form = 'formatted', access = 'sequential', action = 'write',  &
               iostat = ok)

          const = EXP(b(j))

          p = SUM(iobs(1:j)) - iobs(j)

          DO i = 1, n
            WRITE(lu, *) time(i), envobs(i + p), const * envelope(i)
          ENDDO

          CLOSE(lu)

        ENDIF

        DEALLOCATE(time, envelope)

      ENDDO

#ifdef PERF
      CALL watch_stop(tictoc(2), comm2)
      CALL mpi_allreduce(mpi_in_place, tictoc, 2, mpi_double, mpi_max, comm2, ierr)

      IF (world_rank .eq. 0) THEN
        CALL update_log(num2char('Max Exe&IO Time', width=29, fill='.') +  &
                        num2char('[' + num2char(tictoc(1), notation='s', width=10, precision=3) + ',' +   &
                        num2char(tictoc(2), notation='s', width=11, precision=3) + ']', width=36, justify='r'),blankline = .false.)
      ENDIF
#endif

    END SUBROUTINE bestfit

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION misfit(nd, mpar, iter, imod)

      ! Purpose:
      !   to return misfit between observed and simulated coda, where the computation of the latter is based on the RTT. The misfit
      !   is defined as "sum((observed - simulated)**2)"
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      INTEGER(i32),                                                INTENT(IN) :: nd
      REAL(r32),                 DIMENSION(nd),                    INTENT(IN) :: mpar
      INTEGER(i32),                                                INTENT(IN) :: iter, imod
      INTEGER(i32)                                                            :: i, j, j0, j1, k, l, n, p, is, ie, ok, rank, ierr
      INTEGER(i32)                                                            :: req
      INTEGER(i32),              DIMENSION(2)                                 :: dbgrank
      INTEGER(i32),              DIMENSION(0:SIZE(pprank2)-1)                 :: displs, pprank
      REAL(r32)                                                               :: gss, gsp, gpp, gps, aksq, kappa, t
      REAL(r32),                                                   PARAMETER  :: gi = 0._r32, tau = 0.25_r32
      REAL(r32),                                                   PARAMETER  :: wp = 1._r32, ws = 23.4_r32
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                                 :: time, envelope
      REAL(r32),                 DIMENSION(SUM(nobs))                         :: delta, weight, b, tobs
      REAL(r32),                 DIMENSION(SUM(nobs),SIZE(nobs)+1)            :: a
      REAL(r64),                 DIMENSION(2)                                 :: tictoc

      !-----------------------------------------------------------------------------------------------------------------------------

      misfit = 0._r32

#ifdef ERROR_TRAP
      ! return immediately if an important error was raised
      IF (errors .gt. 1) RETURN
#endif

      IF (elastic) THEN
        gss = mpar(1)
        gpp = gss / mpar(2)
        gps = gpp * mpar(3)
        gsp = gps / 6._r32
      ELSE
        gss   = mpar(1)
        aksq  = mpar(2)          !< nu
        kappa = mpar(3)          !< hurst
      ENDIF

#ifdef PERF
      CALL watch_start(tictoc(1), comm2)
#endif

#include "linsys_incl.f90"

#ifdef PERF
      CALL watch_stop(tictoc(1), comm2)
      CALL mpi_allreduce(mpi_in_place, tictoc, 2, mpi_double, mpi_max, comm2, ierr)

      CALL mpi_comm_rank(comm2, dbgrank(1), ierr)
      CALL mpi_comm_rank(comm3, dbgrank(2), ierr)

      IF (ALL(dbgrank .eq. 0)) THEN
        CALL update_log(num2char('Iter, Model, Exe, Wait', width=29, fill='.')         +     &
                        num2char('[' + num2char(iter, width=4)                   + ',' +     &
                        num2char(imod, width=4)                                  + ',' +     &
                        num2char(tictoc(1), notation='s', width=10, precision=3) + ',' +     &
                        num2char(tictoc(2), notation='s', width=11, precision=3) + ']', width=36, justify='r'), blankline=.false.)
      ENDIF
#endif

      n = SIZE(nobs)

      ! misfit is defined as "LN(obs) - (LN(syn) + LN(a) -b*t)"
      ! remember that "b = Qsi*beta = Qpi*alpha"
      DO j = 1, n
        ie = SUM(nobs(1:j))
        is = ie - nobs(j) + 1
        DO i = is, ie
          misfit = misfit + (delta(i) - b(j) + b(n + 1) * tobs(i))**2
        ENDDO
      ENDDO

#ifdef ERROR_TRAP
      CALL share_error(comm2)
#endif

    END FUNCTION misfit

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION gm(g0, aksq, kappa)

      ! Purpose:
      !   To compute the tranport (momentum transfer) scattering coefficient as defined in Eq. 4.28 of Sato et al., 2012.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      REAL(r32),                     INTENT(IN) :: g0, aksq, kappa
      INTEGER(i32)                              :: neval, ier
      INTEGER(i32)                              :: last
      INTEGER(i32),                  PARAMETER  :: limit = 500
      INTEGER(i32),                  PARAMETER  :: lenw = limit * 4
      INTEGER(i32), DIMENSION(limit)            :: iwork
      REAL(r64)                                 :: abserr, x
      REAL(r64),                     PARAMETER  :: pi = 3.14159265358979323846_r64
      REAL(r64),                     PARAMETER  :: epsabs = 1.0E-09_r64
      REAL(r64),                     PARAMETER  :: epsrel = 1.0E-07_r64
      REAL(r64),    DIMENSION(lenw)             :: work

      !-----------------------------------------------------------------------------------------------------------------------------

      ! value afet epsrel can be between 1 and 6 (increase accuracy)

      ! integrate scattering pattern function between 0 and pi
      IF (acf .eq. 'vk') THEN

        cf(1) = 4._r64 * aksq
        cf(2) = kappa + 1.5_r64

        IF (aksq .eq. 0._r64) THEN
          coefficient = 1._r64
        ELSE
          coefficient = cf(1) * (kappa + 0.5_r64) / (1._r64 - (cf(1) + 1)**(-0.5_r64 - kappa))
        ENDIF

        CALL dqag(vkfun, 0._r64, pi, epsabs, epsrel, 2, x, abserr, neval, ier, limit, lenw, last, iwork, work)

      ELSE

        cf(1) = aksq

        IF (aksq .eq. 0._r64) THEN
          coefficient = 1._r64
        ELSE
          coefficient = cf(1) / (1._r64 - EXP(-cf(1)))
        ENDIF

        CALL dqag(gsfun, 0._r64, pi, epsabs, epsrel, 2, x, abserr, neval, ier, limit, lenw, last, iwork, work)

      ENDIF

      gm = 0.5_r32 * g0 * x

    END FUNCTION

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r64) FUNCTION vkfun(psi)

      ! Purpose:
      !   To return the weighted scattering pattern function at angle "psi" for Von Karman media.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      REAL(r64), INTENT(IN) :: psi

      !-----------------------------------------------------------------------------------------------------------------------------

      vkfun = coefficient / (1._r64 + cf(1) * SIN(psi / 2._r64)**2)**cf(2)

      vkfun = vkfun * SIN(psi) * (1._r64 - COS(psi))

    END FUNCTION vkfun

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r64) FUNCTION gsfun(psi)

      ! Purpose:
      !   To return the weighted scattering pattern function at angle "psi" for Gaussian media.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      REAL(r64), INTENT(IN) :: psi

      !-----------------------------------------------------------------------------------------------------------------------------

      gsfun = coefficient * EXP(-cf(1) * SIN(psi / 2._r64)**2)

      gsfun = gsfun * SIN(psi) * (1._r64 - COS(psi))

    END FUNCTION gsfun

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) PURE FUNCTION mean(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      REAL(r32),   DIMENSION(:), INTENT(IN) :: r
      INTEGER(i32)                          :: i
      INTEGER(i32)                          :: n
      REAL(r64)                             :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      v = 0._r64
      c = 1._r64

      DO i = 1, n
        v = v + (r(i) - v) / c
        c = c + 1._r64
      ENDDO

      mean = v

    END FUNCTION mean

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION show(str) RESULT(msg)

      ! Purpose:
      !   To concatenate a list of strings using comma as separator.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      CHARACTER(*), DIMENSION(:), INTENT(IN) :: str
      CHARACTER(:), ALLOCATABLE              :: msg
      INTEGER(i32)                           :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      msg = ''

      DO i = 1, SIZE(str)
        msg = msg + TRIM(str(i)) + ', '
      ENDDO

    END FUNCTION show

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION wrong_par_i32(par, v1, v2, strict) RESULT(wrong_par)

      ! Purpose:
      !   To determine if variable "par" belongs to range [v1, v2]. Square brackets are replaced by round brackets if "strict=.true."
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      INTEGER(i32), DIMENSION(:),           INTENT(IN) :: par
      INTEGER(i32),               OPTIONAL, INTENT(IN) :: v1, v2
      LOGICAL,                    OPTIONAL, INTENT(IN) :: strict
      LOGICAL                                          :: flag, wrong_par

      !-----------------------------------------------------------------------------------------------------------------------------

#include "wrong_par_incl.f90"

    END FUNCTION wrong_par_i32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION wrong_par_r32(par, v1, v2, strict) RESULT(wrong_par)

      ! Purpose:
      !   To determine if variable "par" belongs to range [v1, v2]. Square brackets are replaced by round brackets if "strict=.true."
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      REAL(r32), DIMENSION(:),           INTENT(IN) :: par
      REAL(r32),               OPTIONAL, INTENT(IN) :: v1, v2
      LOGICAL,                 OPTIONAL, INTENT(IN) :: strict
      LOGICAL                                       :: flag, wrong_par

      !-----------------------------------------------------------------------------------------------------------------------------

#include "wrong_par_incl.f90"

    END FUNCTION wrong_par_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE build_comms()

      INTEGER(i32)                              :: cartopo, ierr, nthreads, rank, ntasks
      INTEGER(i32),                   PARAMETER :: ndims = 3
      INTEGER(i32), DIMENSION(ndims)            :: dims, npts, coords
      LOGICAL,                        PARAMETER :: reorder = .true.
      LOGICAL,      DIMENSION(ndims), PARAMETER :: isperiodic = [.false., .false., .false.]

      !-----------------------------------------------------------------------------------------------------------------------------

      nthreads = 1
      !$ nthreads = omp_get_max_threads()

      ! logic grid of points to be shared amongst cpus
      npts = [ns, SIZE(nobs), query_fft_size() / nthreads]

      ! return optimal cpu configuration. In any case, all cpus are used.
      CALL best_cpu_grid(world_size, npts, dims)

      IF (world_rank .eq. 0) THEN
        CALL update_log(num2char('Processors grid', width=29, fill='.') +    &
                        num2char('[' + num2char(dims(1)) + ', ' + num2char(dims(2)) + ', ' + num2char(dims(3)) + ']', width=18,  &
                        justify='r') + num2char('(' + num2char(world_size - PRODUCT(dims)) + ' left out)', width=18, justify='r'))
      ENDIF

      ! create topology
      CALL mpi_cart_create(mpi_comm_world, ndims, dims, isperiodic, reorder, cartopo, ierr)

      ! set new communicators to null because some cpus may be left out
      comm1 = mpi_comm_null
      comm2 = mpi_comm_null
      comm3 = mpi_comm_null

      IF (cartopo .ne. mpi_comm_null) THEN

        CALL mpi_comm_rank(cartopo, rank, ierr)              !< this call is necessary if "reorder = .true."
        CALL mpi_comm_size(cartopo, ntasks, ierr)

        ALLOCATE(gs(ndims, 0:ntasks-1), ge(ndims, 0:ntasks-1))

        ! return process coordinates in current topology
        CALL mpi_cart_coords(cartopo, rank, ndims, coords, ierr)

        ! remove any multithreading to assign indices
        npts(3) = query_fft_size()

        ! return first/last index along each direction for the calling process. Note: first point has first index equal to 1.
        CALL coords2index(npts, dims, coords, gs(:, rank), ge(:, rank))

        ! make all processes aware of global indices
        CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, ndims, mpi_integer, cartopo, ierr)
        CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, ndims, mpi_integer, cartopo, ierr)

        ! build communicators
        CALL build_pencil(cartopo, 0, comm1, pprank1)
        CALL build_pencil(cartopo, 1, comm2, pprank2)
        CALL build_pencil(cartopo, 2, comm3, pprank3)

        ! release cartesian topology
        CALL mpi_comm_free(cartopo, ierr)

        DEALLOCATE(gs, ge)

      ENDIF

    END SUBROUTINE build_comms

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE best_cpu_grid(cpus, npts, dims)

      ! Purpose:
      ! To return an optimal grid of process according to some specific cost function.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      INTEGER(i32),                            INTENT(IN)  :: cpus
      INTEGER(i32),              DIMENSION(3), INTENT(IN)  :: npts
      INTEGER(i32),              DIMENSION(3), INTENT(OUT) :: dims
      INTEGER(i32)                                         :: i, j, k, l, c
      INTEGER(i32)                                         :: n2, n3
      INTEGER(i32)                                         :: a, b
      INTEGER(i32),              DIMENSION(3)              :: v1, v2
      INTEGER(i32), ALLOCATABLE, DIMENSION(:,:)            :: fact3, fact2
      INTEGER(i32), ALLOCATABLE, DIMENSION(:,:)            :: list
      REAL(r32)                                            :: val, minimum

      !-----------------------------------------------------------------------------------------------------------------------------

      minimum = HUGE(1._r32)

      c = 0

      ! factorise the number of available processes (return pairs)
      CALL factorization(world_size, fact3)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! factorise each pair member (thus resulting in a triplet of numbers), evaluate cost function and eventually store triplet
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! number of pairs found
      n3 = SIZE(fact3, 2)

      ! allocate a list where to store all the factorized triplets. Factor of 10 seems to be large enough for 200k processes
      ALLOCATE(list(3, n3 * 10))

      list(:,:) = 0

      ! loop over factorized processes
      DO l = 1, n3

        ! loop over pair members
        DO k = 1, 2

          IF (k .eq. 1) THEN
            a = fact3(1, l)
            b = fact3(2, l)
          ELSE
            b = fact3(1, l)
            a = fact3(2, l)
          ENDIF

          ! factorization of a number returns a new pair
          CALL factorization(b, fact2)

          n2 = SIZE(fact2, 2)

          ! loop over new pairs
          DO j = 1, n2

            ! build candidate triplet
            v1 = [a, fact2(:, j)]

            ! skip to next pair if current triplet already analysed ("v1" is already in "list")
            IF (match() .eqv. .true.) CYCLE

            c = c + 1

            ! triplet is new: add it to the list
            list(:, c) = v1

            ! evaluate cost function for all three possible arrangements (permutations) of the triplet
            DO i = 0, 2

              v1 = CSHIFT(v1, 1)

              ! evaluate cost function
              val = cost_fun(npts, v1)

              IF (val .lt. minimum) THEN
                dims    = v1
                minimum = val
              ENDIF

              v2 = [v1(1), v1(3), v1(2)]

              val = cost_fun(npts, v2)

              IF (val .lt. minimum) THEN
                dims    = v2
                minimum = val
              ENDIF

            ENDDO
            ! end permutations

          ENDDO
          ! end loop over factor pairs for "a/b"

        ENDDO
        ! end loop over "a/b"

      ENDDO
      ! end loop over factor pairs for "world_size"

      DEALLOCATE(fact3, fact2, list)

      CONTAINS

        LOGICAL FUNCTION match()

          INTEGER(i32) :: i

          !-------------------------------------------------------------------------------------------------------------------------

          match = .false.

          DO i = 1, c
            match = ANY(v1(1) .eq. list(:, i)) .and. ANY(v1(2) .eq. list(:, i)) .and. ANY(v1(3) .eq. list(:, i))
            IF (match .eqv. .true.) EXIT
          ENDDO

        END FUNCTION match

    END SUBROUTINE best_cpu_grid

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE factorization(n, factors)

      ! Purpose:
      ! To compute the factorization of an integer number based on the trial division method. For example, for "n=12" the output looks
      ! like "[1 12; 2 6; 3 4]"
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      INTEGER(i32),                              INTENT(IN)  :: n
      INTEGER(i32), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: factors
      INTEGER(i32)                                           :: i, c, s
      INTEGER(i32)                                           :: x
      INTEGER(i32), ALLOCATABLE, DIMENSION(:,:)              :: buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      ! max possible number of factors
      s = FLOOR(SQRT(REAL(n, r32)))

      ALLOCATE(buffer(2, s))

      buffer(:,:) = 0

      ! test factors
      DO i = 1, s

        x = n / i

        IF (MOD(n, i) .eq. 0) THEN
          buffer(1, i) = i                          !< add factor ...
          buffer(2, i) = x                          !< ... and its companion
        ENDIF

      ENDDO

      ! actual factors found
      i = COUNT(buffer(1, :) .ne. 0)

      ALLOCATE(factors(2, i))

      ! copy factors to output array
      c = 0
      DO i = 1, s
        IF (buffer(1, i) .ne. 0) THEN
          c = c + 1
          factors(:, c) = buffer(:, i)
        ENDIF
      ENDDO

      DEALLOCATE(buffer)

    END SUBROUTINE factorization

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION cost_fun(npts, dims)

      ! Purpose:
      ! To return a cost estimate based on the size of a grid block assigned to each process, where the cost is represented by differences
      ! along each side (the goal is to having blocks in shape of cubes).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      INTEGER(i32), DIMENSION(3), INTENT(IN)  :: npts, dims           !< number of points and processes in each direction
      REAL(r32),    DIMENSION(3)              :: side                 !< grid nodes in a block along each direction

      !-----------------------------------------------------------------------------------------------------------------------------

      side = REAL(npts, r32) / REAL(dims, r32)

      !cost_fun = ABS(side(1) - side(2)) + ABS(side(1) - side(3)) + ABS(side(2) - side(3))

      IF (ANY(side .lt. 1._r32)) THEN
        cost_fun = HUGE(0._r32)                                      !< we cannot have more cpus than points
      ELSE
        IF (ANY(nu .gt. 0._r32)) THEN
          cost_fun = side(1) + side(2) + side(3) * 1.0E+03_r32       !< favor more cpus on RTT for fwd scattering
        ELSE
          cost_fun = side(1) + side(2) + side(3) * 1.0E-03_r32       !< hinder cpus on RTT for iso scattering
        ENDIF
      ENDIF

    END FUNCTION cost_fun

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE destroy_comms()

      INTEGER(i32) :: ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      DEALLOCATE(pprank1, pprank2, pprank3)

      CALL mpi_comm_free(comm1, ierr)
      CALL mpi_comm_free(comm2, ierr)
      CALL mpi_comm_free(comm3, ierr)

    END SUBROUTINE destroy_comms

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE build_pencil(comm, dir, newcomm, n)

      ! Purpose:
      ! To group processes falling inside the same pencil oriented along a specific direction "dir" into new communicator "newcomm",
      ! whose size is "ntasks". The rank of the calling process in "newcomm" is "rank" and the number of points for each of them is
      ! assigned to "n".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      INTEGER(i32),                                        INTENT(IN)  :: comm, dir
      INTEGER(i32),                                        INTENT(OUT) :: newcomm
      INTEGER(i32), ALLOCATABLE, DIMENSION(:),             INTENT(OUT) :: n
      INTEGER(i32)                                                     :: i, ierr, rank, ntasks, orank, ontasks
      INTEGER(i32),              DIMENSION(0:world_size-1)             :: color
      LOGICAL,                   DIMENSION(2)                          :: bool

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL mpi_comm_rank(comm, rank, ierr)
      CALL mpi_comm_size(comm, ntasks, ierr)

      ! group processes into pencils
      DO i = 0, ntasks - 1

        color(i) = 0

        ! pencil oriented along x-axis
        IF (dir .eq. 0) THEN
          bool(1) = (gs(2, i) .eq. gs(2, rank)) .and. (ge(2, i) .eq. ge(2, rank))
          bool(2) = (gs(3, i) .eq. gs(3, rank)) .and. (ge(3, i) .eq. ge(3, rank))
        ! pencil oriented along y-axis
        ELSEIF (dir .eq. 1) THEN
          bool(1) = (gs(1, i) .eq. gs(1, rank)) .and. (ge(1, i) .eq. ge(1, rank))
          bool(2) = (gs(3, i) .eq. gs(3, rank)) .and. (ge(3, i) .eq. ge(3, rank))
        ! pencil oriented along z-axis
        ELSEIF (dir .eq. 2) THEN
          bool(1) = (gs(1, i) .eq. gs(1, rank)) .and. (ge(1, i) .eq. ge(1, rank))
          bool(2) = (gs(2, i) .eq. gs(2, rank)) .and. (ge(2, i) .eq. ge(2, rank))
        ENDIF

        IF (ALL(bool .eqv. .true.)) color(i) = i + 1

      ENDDO

      ! process belonging to the same pencil have same color
      color(rank) = MAXVAL(color, dim = 1)

      ! create communicator subgroup
      CALL mpi_comm_split(comm, color(rank), rank, newcomm, ierr)

      ! process id and communicator size
      CALL mpi_comm_rank(newcomm, orank, ierr)
      CALL mpi_comm_size(newcomm, ontasks, ierr)

      ALLOCATE(n(0:ontasks - 1))

      ! number of points along pencil direction for calling process
      n(orank) = ge(dir + 1, rank) - gs(dir + 1, rank) + 1

      ! make whole communicator aware of points for each process
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, n, 1, mpi_integer, newcomm, ierr)

    END SUBROUTINE build_pencil

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE explain_error()

      ! Purpose:
      ! To describe specific error codes.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   03/02/21                  original version
      !

      !---------------------------------------------------------------------------------------------------------------------------------

      IF (errors .eq. 1) THEN
        CALL update_log('Possible loss of accuracy detected')
      ELSEIF (errors .eq. 2) THEN
        CALL update_log('Could not find harmonic coefficients')
      ELSEIF (errors .eq. 3) THEN
        CALL update_log('Could not expand scattering pattern function accurately')
      ELSEIF (errors .eq. 4) THEN
        CALL update_log('Detected negative values in synthetic envelope')
        IF (elastic) THEN
          CALL update_log(num2char('Triggered for', width=29, fill='.') + num2char('Max Time', width=10, justify='r') + '|' +  &
                          num2char('DT', width=10, justify='r') + '|' +  num2char('Tp', width=10, justify='r') + '|' +  &
                          num2char('Ts', width=10, justify='r') + '|' +  num2char('EtaSS', width=10, justify='r') + '|' +  &
                          num2char('EtaSS/PP', width=10, justify='r') + '|' + num2char('EtaPS/PP', width=10, justify='r') + '|' )

          CALL update_log(num2char('', width=29, fill='.') +     &
                          num2char(error_params(1), notation='f', width=10, precision=2, justify='r') + '|' +  &
                          num2char(error_params(2), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(3), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(4), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(5), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(6), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(7), notation='f', width=10, precision=4, justify='r') + '|', blankline = .false.)
        ELSE
          CALL update_log(num2char('Triggered for', width=29, fill='.') + num2char('Max Time', width=10, justify='r') + '|' +  &
                          num2char('DT', width=10, justify='r') + '|' +  num2char('Ts', width=10, justify='r') + '|' +  &
                          num2char('EtaSS', width=10, justify='r') + '|' + num2char('Nu', width=10, justify='r') + '|' + &
                          num2char('Hurst', width=10, justify='r') + '|')

          CALL update_log(num2char('', width=29, fill='.') +   &
                          num2char(error_params(1), notation='f', width=10, precision=2, justify='r') + '|' +  &
                          num2char(error_params(2), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(3), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(4), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(5), notation='f', width=10, precision=4, justify='r') + '|' +  &
                          num2char(error_params(6), notation='f', width=10, precision=4, justify='r') + '|' , blankline = .false.)
        ENDIF
      ELSEIF (errors .gt. 10) THEN
        CALL update_log('Linear regression failed with code ' + num2char(errors - 10))
      ENDIF

    END SUBROUTINE explain_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
    !===================================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

    SUBROUTINE share_error(comm)

      ! Purpose:
      !   To share error code and parameters that triggered the error among members of "comm".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   11/01/21                  original version
      !

      INTEGER(i32),                           INTENT(IN) :: comm
      INTEGER(i32)                                       :: rank, ntasks, ierr
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)            :: list

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL mpi_comm_size(comm, ntasks, ierr)
      CALL mpi_comm_rank(comm, rank, ierr)

      ALLOCATE(list(0:ntasks - 1))

      list(rank) = errors

      ! exchange error codes inside "comm"
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, list, 1, mpi_int, comm, ierr)

      ! all members of "comm" receive parameters that triggered error
      IF (ANY(list .eq. 4)) THEN
        rank   = MINLOC(ABS(list - 4), DIM=1) - 1
        errors = 4
        CALL mpi_bcast(error_params, 7, mpi_real, rank, comm, ierr)
      ELSE
        errors = MAXVAL(list)
      ENDIF

      DEALLOCATE(list)

    END SUBROUTINE share_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_inveta

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

PROGRAM main

  ! Purpose:
  !   to invert coda waves and return best fitting model parameters based on RTT. Sampling of parameters space relies on the NA by
  !   Sambridge. The following inversion setups are possible:
  !
  !     a) inversion of each single coda recorded by a receiver
  !     b) single inversion of all coda recorded by a receiver
  !     c) single inversion of coda generated by the same event and recorded by a group of receivers
  !
  !   The code read recordings in miniseed ASCII format, filter and downsample them, extract envelope and then run the inversion. It
  !   is assumed that recordings first time sample corresponds to earthquake origin time, the latter being set to zero (P/S-wave
  !   arrival times must be based on this convention).
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !

  USE, INTRINSIC     :: iso_fortran_env, only: compiler_version
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_fft_cmplx
  USE, NON_INTRINSIC :: m_filter_r32
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_na
  USE, NON_INTRINSIC :: m_inveta
  USE                :: mpi

  IMPLICIT none

  CHARACTER(30)                              :: current
  CHARACTER(:),  ALLOCATABLE                 :: fo
  CHARACTER(10), ALLOCATABLE, DIMENSION(:)   :: code
  CHARACTER(30), ALLOCATABLE, DIMENSION(:)   :: inverted, blacklisted
  COMPLEX(r32),  ALLOCATABLE, DIMENSION(:)   :: analytic, spectrum
  INTEGER(i32)                               :: ierr, ok, lu, n, i, j, k, p, l, is, ie, pts, rank
  INTEGER(i32),               DIMENSION(2)   :: v
  INTEGER(i32),  ALLOCATABLE, DIMENSION(:)   :: lrank
  REAL(r32)                                  :: dt, t, gss, gpp, gsp, gps, aksq, kappa
  REAL(r64),                  DIMENSION(2)   :: tictoc
  REAL(r32),     ALLOCATABLE, DIMENSION(:)   :: time, trespl, h, envlp, respl, bestmodel, lrange, hrange, sampled, fit
  REAL(r32),     ALLOCATABLE, DIMENSION(:,:) :: timeseries

  !---------------------------------------------------------------------------------------------------------------------------------

  CALL mpi_init(ierr)

  CALL watch_start(tictoc(1), mpi_comm_world)

  CALL mpi_comm_size(mpi_comm_world, world_size, ierr)
  CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

  ! IF (world_rank .eq. 0) CALL set_log_module(ok, screen = .true.)
  CALL set_log_module(ok, screen = .true.)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! --------------------------------- read main input file and stop execution if error is raised -----------------------------------

  IF (world_rank .eq. 0) CALL read_input_file(ok)

  CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

  IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------- echo input parameters  -------------------------------------------------------

  IF (world_rank .eq. 0) THEN

    CALL update_log('**********************************************************************************************')
    CALL update_log('INVETA v1.0')
    CALL update_log('Compiled with ' + COMPILER_VERSION())
    CALL update_log('----------------------------------------------------------------------------------------------')
    CALL update_log('Summary of input parameters', blankline = .false.)
    CALL update_log(num2char('Inversion mode', width=29, fill='.') + num2char(mode, width=17, justify='r') + '|')
    IF (mode .eq. 2) &
      CALL update_log(num2char('Receivers threshold', width=29, fill='.') + num2char(threshold), blankline = .false.)

    CALL update_log(num2char('NA parameters', width=29, fill='.') + num2char('InitialModels', width=17, justify='r') + '|' +    &
                    num2char('Models', width=13, justify='r') + '|' + num2char('Resampled', width=13, justify='r') + '|' +  &
                    num2char('Iterations', width=13, justify='r') + '|' + num2char('Seed', width=13, justify='r') + '|',    &
                    blankline = .false.)

    CALL update_log(num2char('', width=29) + num2char(nsi, width=17, justify='r') + '|' +                           &
                    num2char(ns, width=13, justify='r') + '|' + num2char(nr, width=13, justify='r') + '|' +         &
                    num2char(itermax, width=13, justify='r') + '|' + num2char(seed, width=13, justify='r') + '|',   &
                    blankline = .false.)

    CALL update_log(num2char('Shear-wave speed', width=29, fill='.') + num2char(beta, notation='f', precision = 3, width=17, &
                    justify='r') + '|', blankline = .false.)
    CALL update_log(num2char('Frequency band', width=29, fill='.') + num2char('low (Hz)', width=17, justify = 'r') + '|' +   &
                    num2char('high (Hz)', width=13, justify = 'r') + '|', blankline = .false.)
    DO i = 1, SIZE(fbands, 2)
      CALL update_log(num2char('', width=29) + num2char(fbands(1,i), notation='f', width=17, precision=2) + '|' +   &
                      num2char(fbands(2,i), notation='f', width=13, precision=2) + '|', blankline = .false.)
    ENDDO
    CALL update_log(num2char('Number of receivers found', width = 29, fill='.') + num2char(SIZE(recvr), width=17, justify='r') +  &
                    '|', blankline = .false.)
    CALL update_log(num2char('Weighting of direct waves', width=29, fill='.')+ num2char(.not.noweight, width=17, justify='r') +  &
                    '|', blankline = .false.)
    CALL update_log(num2char('Elastic RTT', width=29, fill='.') + num2char(elastic, width=17, justify='r') + '|',  &
                    blankline = .false.)

    IF (elastic) THEN
      CALL update_log(num2char('Window width direct waves', width=29, fill='.') + num2char('P', width=17, justify='r') + '|' +  &
                      num2char('S', width=13, justify='r') + '|' + num2char('Factor (%)', width=13, justify='r') + '|',         &
                      blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(pdwindow, notation = 'f', width=17, precision=3, justify='r') + '|' +  &
                      num2char(sdwindow, notation='f', width=13, precision=3, justify='r') + '|' + num2char(fwin*100._r32,     &
                      notation='f', width=13, precision=1, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('Window width coda waves', width=29, fill='.') + num2char('P', width=17, justify='r') + '|' +  &
                      num2char('S', width=13, justify='r') + '|' + num2char('Tlim', width=13, justify='r') + '|',    &
                      blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(pcwindow, notation = 'f', width=17, precision=3, justify='r') + '|' +  &
                      num2char(scwindow, notation='f', width=13, precision=3, justify='r') + '|' + num2char(tlim, notation='f',  &
                      width=13, precision=1, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('Parameters search range', width=29, fill='.') + num2char('EtaSS', width=17, justify='r') + '|' + &
                      num2char('EtaSS/PP', width=13, justify='r') + '|' + num2char('EtaPS/PP', width=13, justify='r') + '|',  &
                      blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(etass(1),    notation='s', width=8, precision=1, justify='r') + ',' + &
                                               num2char(etass(2),    notation='s', width=8, precision=1) + '|' + &
                                               num2char(etass2pp(1), notation='f', width=6, precision=1, justify='r') + ',' + &
                                               num2char(etass2pp(2), notation='f', width=6, precision=1) + '|' + &
                                               num2char(etaps2pp(1), notation='f', width=6, precision=1, justify='r') + ',' + &
                                               num2char(etaps2pp(2), notation='f', width=6, precision=1) + '|', blankline = .false.)

    ELSE
      CALL update_log(num2char('Window width direct S-wave', width=29, fill='.') + num2char('S', width=17, justify='r') + '|' +  &
                      num2char('Factor (%)', width=13, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(sdwindow, notation='f', width=17, precision=3, justify='r') + '|' +  &
                      num2char(fwin*100._r32, notation='f', width=13, precision=1, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('Window width S-coda waves', width=29, fill='.') + num2char('S', width=17, justify='r') + '|' +  &
                      num2char('Tlim', width=13, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(scwindow, notation = 'f', width=17, precision=3, justify='r') + '|' +   &
                      num2char(tlim, notation='f', width=13, precision=1) + '|', blankline = .false.)
      CALL update_log(num2char('Parameters search range', width=29, fill='.') + num2char('EtaSS', width=17, justify='r') + '|' +  &
                      num2char('Nu', width=13, justify='r') + '|' + num2char('Hurst', width=13, justify='r') + '|',               &
                      blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(etass(1), notation='s', width=8, precision=1, justify='r') + ',' +  &
                      num2char(etass(2),    notation='s', width=8, precision=1) + '|' +                      &
                      num2char(nu(1), notation='f', width=6, precision=2, justify='r') + ',' +         &
                      num2char(nu(2), notation='f', width=6, precision=2) + '|' +                      &
                      num2char(hurst(1), notation='f', width=6, precision=2, justify='r') + ',' +      &
                      num2char(hurst(2), notation='f', width=6, precision=2) + '|', blankline = .false.)
      IF (ANY(nu .ne. 0._r32)) THEN
        CALL update_log(num2char('acf of choice', width=29, fill='.') + num2char(acf, width=17, justify='r') + '|',    &
                        blankline = .false.)
      ENDIF
    ENDIF

    CALL update_log('----------------------------------------------------------------------------------------------')

  ENDIF

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------- broadcast relevant input parameters  -----------------------------------------------

  CALL mpi_bcast(mode, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(threshold, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nsi, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(ns, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nr, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(itermax, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(seed, 1, mpi_int, 0, mpi_comm_world, ierr)

  CALL mpi_bcast(beta, 1, mpi_real, 0, mpi_comm_world, ierr)

  n = SIZE(fbands, 2)

  CALL mpi_bcast(n, 1, mpi_int, 0, mpi_comm_world, ierr)

  IF (world_rank .ne. 0) ALLOCATE(fbands(2, n))
  CALL mpi_bcast(fbands, n*2, mpi_real, 0, mpi_comm_world, ierr)

  n = SIZE(recvr)

  CALL mpi_bcast(n, 1, mpi_int, 0, mpi_comm_world, ierr)

  IF (world_rank .ne. 0) ALLOCATE(recvr(n))

  DO i = 1, SIZE(recvr)

    n = SIZE(recvr(i)%ts)
    CALL mpi_bcast(n, 1, mpi_int, 0, mpi_comm_world, ierr)

    IF (world_rank .ne. 0) ALLOCATE(recvr(i)%event(n), recvr(i)%code(n), recvr(i)%tp(n), recvr(i)%ts(n))

    CALL mpi_bcast(recvr(i)%event, LEN(recvr(i)%event(1))*n, mpi_character, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(recvr(i)%code, LEN(recvr(i)%code)*n, mpi_character, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(recvr(i)%tp, n, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(recvr(i)%ts, n, mpi_real, 0, mpi_comm_world, ierr)

  ENDDO

  CALL mpi_bcast(noweight, 1, mpi_logical, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(elastic, 1, mpi_logical, 0, mpi_comm_world, ierr)

  CALL mpi_bcast(fwin, 1, mpi_real, 0, mpi_comm_world, ierr)

  IF (world_rank .ne. 0) ALLOCATE(etass(2))
  CALL mpi_bcast(etass, 2, mpi_real, 0, mpi_comm_world, ierr)

  CALL mpi_bcast(sdwindow, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(scwindow, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(tlim, 1, mpi_real, 0, mpi_comm_world, ierr)

  IF (elastic) THEN
    IF (world_rank .ne. 0) ALLOCATE(etass2pp(2), etaps2pp(2))
    CALL mpi_bcast(etass2pp, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(etaps2pp, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(pdwindow, 1, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(pcwindow, 1, mpi_real, 0, mpi_comm_world, ierr)
  ELSE
    IF (world_rank .ne. 0) ALLOCATE(nu(2), hurst(2))
    IF (world_rank .ne. 0) ALLOCATE(CHARACTER(2) :: acf)
    CALL mpi_bcast(nu, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(hurst, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(acf, 2, mpi_character, 0, mpi_comm_world, ierr)
  ENDIF

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ----------------------------------------------------- parameters inversion -----------------------------------------------------

  CALL setup_interpolation('linear', 'copy', ok)

  IF (elastic) THEN
    lrange = [etass(1), etass2pp(1), etaps2pp(1)]
    hrange = [etass(2), etass2pp(2), etaps2pp(2)]
  ELSE
    lrange = [etass(1), nu(1), hurst(1)]
    hrange = [etass(2), nu(2), hurst(2)]
  ENDIF

  ALLOCATE(bestmodel(SIZE(lrange)))

  DO k = 1, SIZE(fbands, 2)                         !< loop over frequency bands

    drespl = 0.25_r32 / fbands(2, k)                !< Nyquist frequency after filtering is twice low-pass frequency

    current = ' '

    DO                                              !< loop needed only when "mode = 2"

      DO l = 1, SIZE(recvr)                           !< loop over receivers

        DO p = 1, SIZE(recvr(l)%event)                   !< loop over events for current receiver

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! --------------------------- check if event was recorded by other stations ("mode=2" only) ------------------------------

          IF (mode .eq. 2) THEN

            n = 0

            ! find how many receivers recorded current event and store this value in "n"
            DO j = 1, SIZE(recvr)
              DO i = 1, SIZE(recvr(j)%event)
                IF (recvr(j)%event(i) .eq. recvr(l)%event(p)) n = n + 1
              ENDDO
            ENDDO

            IF (n .lt. threshold) CYCLE                !< skip if event was recorded by too few receivers

            IF (LEN_TRIM(current) .ne. 0) THEN
              IF (recvr(l)%event(p) .ne. current) CYCLE               !< skip if event is not the one we are working on
            ELSE
              IF (ALLOCATED(inverted)) THEN
                IF (ANY(inverted .eq. recvr(l)%event(p))) CYCLE       !< skip if event has been already inverted
              ENDIF
              IF (ALLOCATED(blacklisted)) THEN
                IF (ANY(blacklisted .eq. recvr(l)%event(p))) CYCLE    !< skip if event has too few receivers with long enough windows
              ENDIF
            ENDIF

          ENDIF

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ------------------------------------- read recorded ground-motion, then broadcast  -------------------------------------

          IF (ALLOCATED(timeseries)) DEALLOCATE(timeseries, time)

          IF (world_rank .eq. 0) THEN

            fo = TRIM(recvr(l)%folder) + '/' + TRIM(recvr(l)%event(p)) + '_' + TRIM(recvr(l)%network(p)) + '_' +   &
                 TRIM(recvr(l)%code(l)) + '_' + TRIM(recvr(l)%channel(p)) + '[ENZ].mseed.ascii'

            CALL read_miniseed(ok, fo, dt, timeseries, recvr(l)%ts(p) * tlim)

            ! this is a hack to deal with fancy SED file naming convention
            IF (ok .ne. 0) THEN
              DO i = 1, 4
                fo = TRIM(recvr(l)%folder) + '/' + TRIM(recvr(l)%event(p)) + '_' + TRIM(recvr(l)%network(p)) + '_' +   &
                     TRIM(recvr(l)%code(l)) + '_' + TRIM(recvr(l)%channel(p)) + '[ENZ]_' + num2char(i) + '.mseed.ascii'
                CALL read_miniseed(ok, fo, dt, timeseries, recvr(l)%ts(p) * tlim)
                IF (ok .eq. 0) EXIT
              ENDDO
            ENDIF

          ENDIF

          CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

          IF (ok .ne. 0) THEN
            IF (world_rank .eq. 0) CALL report_error('File ' + fo + ' not found in folder ' + recvr(l)%folder)
            CALL mpi_barrier(mpi_comm_world, ierr)
            CALL mpi_abort(mpi_comm_world, ok, ierr)
          ENDIF

          IF (world_rank .eq. 0) v = SHAPE(timeseries)

          CALL mpi_bcast(v, 2, mpi_int, 0, mpi_comm_world, ierr)

          IF (world_rank .ne. 0) ALLOCATE(timeseries(v(1), v(2)))

          CALL mpi_bcast(timeseries, PRODUCT(v), mpi_real, 0, mpi_comm_world, ierr)
          CALL mpi_bcast(dt, 1, mpi_real, 0, mpi_comm_world, ierr)

          ALLOCATE(time(SIZE(timeseries, 1)))

          DO i = 1, SIZE(timeseries, 1)
            time(i) = (i - 1) * dt                                      !< original time vector (from recordings)
          ENDDO

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ----------------------------------------- check if windows are long enough  --------------------------------------------

          ! skip if recordings cannot contain all the desired windows
          IF ((recvr(l)%ts(p) + sdwindow*fwin + scwindow) .ge. time(SIZE(time))) CYCLE
          IF (elastic .and. ((recvr(l)%tp(p) + pdwindow*fwin + pcwindow) .ge. recvr(l)%ts(p))) CYCLE

          ! if we arrived at this point it means that the event can be inverted (another check will occur below for "mode=2")
          current = recvr(l)%event(p)
          IF (ALLOCATED(code)) THEN
            code = [code, recvr(l)%code(p)]
          ELSE
            code = [recvr(l)%code(p)]
          ENDIF

          ! add current event to list of already inverted events, just for making a summary later on
          IF (mode .le. 1) THEN
            IF (ALLOCATED(inverted)) THEN
              inverted = [inverted, current]
            ELSE
              inverted = [current]
            ENDIF
          ENDIF

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ---------------------------------------- filter, resample, compute envelope  -------------------------------------------

          n = NINT((SIZE(timeseries, 1) - 1) * dt / drespl) + 1         !< number of samples in decimated signal

          ALLOCATE(trespl(n), h(n), envlp(n), respl(n), analytic(n), spectrum(n))

          DO i = 1, n
            trespl(i) = (i - 1) * drespl                                !< time vector for resampling
          ENDDO

          ! "h" is weighting function to compute analytic signal
          h(:) = 0._r32
          h(1) = 1._r32

          IF (MOD(n, 2) .eq. 0) THEN
            h(n/2 + 1) = 1._r32
            h(2:n/2)   = 2._r32
          ELSE
            h(2:(n + 1)/2) = 2._r32
          ENDIF

          CALL make_fftw_plan([n])
          CALL make_iir_plan(ok, 'butter', dt, fbands(:, k), 'pass', 2, zphase = .true.)         !< two-pass, two-poles passband

          envlp(:) = 0._r32

          ! process components of motion and stack them
          DO j = 1, SIZE(timeseries, 2)

            timeseries(:, j) = iir(timeseries(:, j), ok)                 !< filter

            CALL interpolate(time, timeseries(:, j), trespl, respl)      !< resample

            CALL fft(respl, spectrum)

            DO i = 1, n
              spectrum(i) = spectrum(i) * h(i)
            ENDDO

            CALL ifft(analytic, spectrum)                             !< return analytic signal

            ! compute energy envelope
            DO i = 1, n
              envlp(i) = envlp(i) + 0.5_r32 * ABS(analytic(i))**2     !< compute envelope (eq. 2.27 Sato&Fehler) and stack
            ENDDO

          ENDDO

          CALL destroy_fftw_plan([n])
          CALL destroy_iir_plan()

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ------------------------------------------ add observables for P-wave windows  -----------------------------------------

          pts = 0                 !< "pts" store how many points are available for inversion

          IF (elastic) THEN

            IF (ALLOCATED(tpobs)) THEN
              tpobs = [tpobs, recvr(l)%tp(p)]                        !< append P-wave arrival time
            ELSE
              tpobs = [recvr(l)%tp(p)]
            ENDIF

            is = NINT((recvr(l)%tp(p) - pdwindow * (1._r32 - fwin)) / drespl) + 1               !< beginning direct P-window
            is = MAX(1, is)                                                                     !< make sure we don't go below 1
            ie = is + NINT(fwin*pdwindow / drespl)                                              !< end direct P-window
            is = NINT((recvr(l)%ts(p) - sdwindow * (1._r32 - fwin)) / drespl) + 1               !< beginning direct S-window

            ! we take one sample from direct P-window, then all the way up to just before direct S-window
            pts = 1 + (is - 1) - (ie + 1) + 1

          ENDIF

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ------------------------------------------ add observables for S-wave windows  -----------------------------------------

          is = NINT((recvr(l)%ts(p) - sdwindow * (1._r32 - fwin)) /  drespl) + 1                !< beginning direct S-window
          is = MAX(1, is)                                                                       !< make sure we don't go below 1
          ie = is + NINT(fwin*sdwindow / drespl)                                                !< end direct S-window

          ! we take one sample from direct S-window, then all the way up until the end
          pts = pts + 1 + n - (ie + 1) + 1

          IF (ALLOCATED(tsobs)) THEN
            tsobs  = [tsobs, recvr(l)%ts(p)]                         !< append S-wave arrival time
            envobs = [envobs, envlp]                                 !< append envelope
            iobs   = [iobs, n]                                       !< append total number of points for current envelope
            nobs   = [nobs, pts]                                     !< append actual number of points used during inversion
          ELSE
            tsobs  = [recvr(l)%ts(p)]
            envobs = [envlp]
            iobs   = [n]
            nobs   = [pts]
          ENDIF

          DEALLOCATE(trespl, h, envlp, respl, analytic, spectrum)

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! -------------------------------- inversion of each single coda recorded by a receiver ----------------------------------

          IF (mode .eq. 0) THEN
            IF (ALLOCATED(code)) THEN             !< obviously, do something only if receiver has recorded current event

              IF (world_rank .eq. 0) THEN
                CALL update_log('Inverting event ' + TRIM(current) + ' for receiver ' + TRIM(code(1)) + ' in the frequency band '+ &
                                num2char(fbands(1,k), notation='f', precision=1) + '-' +      &
                                num2char(fbands(2,k), notation='f', precision=1) + 'Hz')
              ENDIF

              CALL watch_start(tictoc(2), mpi_comm_self)

              ! arrange available processes into logic grid
              CALL build_comms()

              errors = 0

              IF (comm1 .ne. mpi_comm_null) THEN
                CALL na(comm1, misfit, lrange, hrange, seed, itermax, nsi, ns, nr, 0, bestmodel, sampled, fit, ok)
              ENDIF

#ifdef ERROR_TRAP
              IF (comm1 .ne. mpi_comm_null) CALL share_error(comm1)

              ! broadcast error code to tasks left outside cpus grid
              CALL mpi_bcast(errors, 1, mpi_int, 0, mpi_comm_world, ierr)

              IF (world_rank .eq. 0) CALL explain_error()

              IF (errors .gt. 1) CALL mpi_abort(mpi_comm_world, ok, ierr)

              CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

              IF (ok .ne. 0) THEN
                IF (world_rank .eq. 0) CALL report_error('NA exited with error code: ' + num2char(ok))
                CALL mpi_abort(mpi_comm_world, ok, ierr)
              ENDIF
#endif

              IF (comm1 .ne. mpi_comm_null) THEN
                CALL bestfit(k, [recvr(l)%code(1)], [current], bestmodel)
                CALL destroy_comms()
              ENDIF

              IF (world_rank .eq. 0) THEN

                ! write explored parameters space to disk
                CALL write_search(k, [recvr(l)%code(1)], [current], sampled, fit)

                CALL update_log(num2char('Models explored', width=29, fill='.') + num2char(SIZE(fit),width=17, justify='r') + '|')

                CALL watch_stop(tictoc(2), mpi_comm_self)

                tictoc(2) = tictoc(2) / 60._r64

                CALL update_log(num2char('Elapsed time (min)', width=29, fill='.') +    &
                                num2char(tictoc(2), notation='f', width=17, precision=1, justify='r') + '|', blankline = .false.)

                ! update log file
                IF (elastic) THEN
                  gss = bestmodel(1)
                  gpp = gss / bestmodel(2)
                  gps = gpp * bestmodel(3)
                  gsp = gps / 6._r32

                  CALL update_log(num2char('Best model parameters', width=29, fill='.') +    &
                                  num2char('EtaPP', width=17, justify='r') + '|' +   &
                                  num2char('EtaPS', width=13, justify='r') + '|' +   &
                                  num2char('EtaSP', width=13, justify='r') + '|' +   &
                                  num2char('EtaSS', width=13, justify='r') + '|', blankline = .false.)

                  CALL update_log(num2char('', width=29) + num2char(gpp, notation='s', width=17, precision=3, justify='r') + '|' + &
                                  num2char(gps, notation='s', width=13, precision=3, justify='r') + '|' + &
                                  num2char(gsp, notation='s', width=13, precision=3, justify='r') + '|' + &
                                  num2char(gss, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

                ELSE
                  gss   = bestmodel(1)
                  aksq  = bestmodel(2)
                  kappa = bestmodel(3)

                  CALL update_log(num2char('Best model parameters', width=29, fill='.') +    &
                                  num2char('EtaSS', width=17, justify='r') + '|' +    &
                                  num2char('Nu', width=13, justify='r')    + '|' +    &
                                  num2char('Hurst', width=13, justify='r') + '|', blankline = .false.)

                  CALL update_log(num2char('', width=29) + num2char(gss, notation='s', width=17, precision=3, justify='r') + '|' + &
                                  num2char(aksq, notation='s', width=13, precision=3, justify='r') + '|' +   &
                                  num2char(kappa, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

                ENDIF
              ENDIF

              DEALLOCATE(envobs, tsobs, nobs, iobs)
              DEALLOCATE(sampled, fit)
              IF (elastic) DEALLOCATE(tpobs)
              DEALLOCATE(code)

            ENDIF
          ENDIF

        ENDDO                                      !< end loop over events for current receiver

        IF (mode .eq. 0) THEN
          IF (ALLOCATED(inverted)) THEN
            IF (world_rank .eq. 0) CALL update_log('For receiver ' + TRIM(recvr(l)%code(1)) + ' inverted ' +    &
                                                    num2char(SIZE(inverted)) + ' event(s)')
            DEALLOCATE(inverted)
          ELSE
            IF (world_rank .eq. 0) CALL update_log('For receiver ' + TRIM(recvr(l)%code(1)) + ' no events inverted')
          ENDIF
        ENDIF

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ---------------------------------- single inversion of all coda recorded by a receiver -----------------------------------

        IF (mode .eq. 1) THEN
          IF (ALLOCATED(code)) THEN                  !< obviously, do something only if receiver has recorded some events

            IF (SIZE(code) .ge. threshold) THEN

              IF (world_rank .eq. 0) THEN
                CALL update_log('Inverting ' + num2char(SIZE(inverted)) + ' event(s) for receiver ' + TRIM(recvr(l)%code(1)) +  &
                                ' in the frequency band ' + num2char(fbands(1,k), notation='f', precision=1) + '-' +            &
                                num2char(fbands(2,k), notation='f', precision=1) + 'Hz')
              ENDIF

              CALL watch_start(tictoc(2), mpi_comm_self)

              ! arrange available processes into logic grid
              CALL build_comms()

              ! keep track of any error occurred during envelope calculation and NA search
              errors = 0

              IF (comm1 .ne. mpi_comm_null) THEN
                CALL na(comm1, misfit, lrange, hrange, seed, itermax, nsi, ns, nr, 0, bestmodel, sampled, fit, ok)
              ENDIF

#ifdef ERROR_TRAP
              IF (comm1 .ne. mpi_comm_null) CALL share_error(comm1)

              ! broadcast error code to tasks left outside cpus grid
              CALL mpi_bcast(errors, 1, mpi_int, 0, mpi_comm_world, ierr)

              IF (world_rank .eq. 0) CALL explain_error()

              ! do stop only for important errors
              IF (errors .gt. 1) CALL mpi_abort(mpi_comm_world, ok, ierr)

              CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

              IF (ok .ne. 0) THEN
                IF (world_rank .eq. 0) CALL report_error('NA exited with error code: ' + num2char(ok))
                CALL mpi_abort(mpi_comm_world, ok, ierr)
              ENDIF
#endif

              IF (comm1 .ne. mpi_comm_null) THEN
                ! write best-fitting model to disk
                CALL bestfit(k, [recvr(l)%code(1)], inverted, bestmodel)
                CALL destroy_comms()
              ENDIF

              IF (world_rank .eq. 0) THEN

                ! write explored parameters space to disk
                CALL write_search(k, [recvr(l)%code(1)], inverted, sampled, fit)

                CALL update_log(num2char('Models explored', width=29, fill='.') + num2char(SIZE(fit), width=17, justify='r') + '|')

                CALL watch_stop(tictoc(2), mpi_comm_self)

                tictoc(2) = tictoc(2) / 60._r64

                CALL update_log(num2char('Elapsed time (min)', width=29, fill='.') +    &
                                num2char(tictoc(2), notation='f', width=17, precision=1, justify='r') + '|', blankline = .false.)

                IF (elastic) THEN
                  gss = bestmodel(1)
                  gpp = gss / bestmodel(2)
                  gps = gpp * bestmodel(3)
                  gsp = gps / 6._r32

                  CALL update_log(num2char('Best model parameters', width=29, fill='.') +    &
                                  num2char('EtaPP', width=17, justify='r') + '|' +   &
                                  num2char('EtaPS', width=13, justify='r') + '|' +   &
                                  num2char('EtaSP', width=13, justify='r') + '|' +   &
                                  num2char('EtaSS', width=13, justify='r') + '|', blankline = .false.)

                  CALL update_log(num2char('', width=29) + num2char(gpp, notation='s', width=17, precision=3, justify='r') + '|' + &
                                  num2char(gps, notation='s', width=13, precision=3, justify='r') + '|' + &
                                  num2char(gsp, notation='s', width=13, precision=3, justify='r') + '|' + &
                                  num2char(gss, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

                ELSE
                  gss   = bestmodel(1)
                  aksq  = bestmodel(2)
                  kappa = bestmodel(3)

                  CALL update_log(num2char('Best model parameters', width=29, fill='.') +    &
                                  num2char('EtaSS', width=17, justify='r') + '|' +    &
                                  num2char('Nu', width=13, justify='r')    + '|' +    &
                                  num2char('Hurst', width=13, justify='r') + '|', blankline = .false.)

                  CALL update_log(num2char('', width=29) + num2char(gss, notation='s', width=17, precision=3, justify='r') + '|' + &
                                  num2char(aksq, notation='s', width=13, precision=3, justify='r') + '|' +   &
                                  num2char(kappa, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

                ENDIF
              ENDIF

            ELSE

              IF (world_rank .eq. 0) CALL update_log('For receiver ' + TRIM(recvr(l)%code(1)) + ' not enough events were found')

            ENDIF

            DEALLOCATE(envobs, tsobs, nobs, iobs)
            DEALLOCATE(sampled, fit)
            IF (elastic) DEALLOCATE(tpobs)
            DEALLOCATE(code, inverted)

          ELSE
            IF (world_rank .eq. 0) CALL update_log('For receiver ' + TRIM(recvr(l)%code(1)) + ' no events inverted')

          ENDIF
        ENDIF

      ENDDO                                       !< end loop over receivers

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------- single inversion of coda generated by same event ---------------------------------------

      IF (mode .eq. 2) THEN

        IF (LEN_TRIM(current) .eq. 0) EXIT         !< no further events were found, stop searching (exit second loop)

        IF (SIZE(code) .ge. threshold) THEN        !< do something only if enough receivers recorded the same event

          IF (world_rank .eq. 0) THEN
            CALL update_log('Inverting ' + num2char(SIZE(code)) + ' recordings for event ' + TRIM(current)  +  &
                            ' in the frequency band ' + num2char(fbands(1,k), precision=1) + '-' +             &
                            num2char(fbands(2,k), precision=1) + 'Hz')
            CALL update_log(num2char('List of receivers', width=29, fill='.') + show(code), blankline = .false.)
          ENDIF

          CALL watch_start(tictoc(2), mpi_comm_self)

          ! arrange available processes into logic grid
          CALL build_comms()

          ! keep track of any error occurred during envelope calculation and NA search
          errors = 0

          IF (comm1 .ne. mpi_comm_null) THEN
            CALL na(comm1, misfit, lrange, hrange, seed, itermax, nsi, ns, nr, 0, bestmodel, sampled, fit, ok)
          ENDIF

#ifdef ERROR_TRAP
          IF (comm1 .ne. mpi_comm_null) CALL share_error(comm1)

          ! broadcast error code to tasks left outside cpus grid
          CALL mpi_bcast(errors, 1, mpi_int, 0, mpi_comm_world, ierr)

          IF (world_rank .eq. 0) CALL explain_error()

          ! do stop only for important errors
          IF (errors .gt. 1) CALL mpi_abort(mpi_comm_world, ok, ierr)

          CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

          IF (ok .ne. 0) THEN
            IF (world_rank .eq. 0) CALL report_error('NA exited with error code: ' + num2char(ok))
            CALL mpi_abort(mpi_comm_world, ok, ierr)
          ENDIF
#endif

          IF (comm1 .ne. mpi_comm_null) THEN
            ! write best-fitting model to disk and explored parameters space to disk
            CALL bestfit(k, code, [current], bestmodel)
            CALL destroy_comms()
          ENDIF

          IF (world_rank .eq. 0) THEN

            ! write explored parameters space to disk
            CALL write_search(k, code, [current], sampled, fit)

            CALL update_log(num2char('Models explored', width=29, fill='.') + num2char(SIZE(fit),width=17, justify='r') + '|')

            CALL watch_stop(tictoc(2), mpi_comm_self)

            tictoc(2) = tictoc(2) / 60._r64

            CALL update_log(num2char('Elapsed time (min)', width=29, fill='.') + num2char(tictoc(2), notation = 'f', width = 17,  &
                            precision = 1, justify='r') + '|')

            IF (elastic) THEN
              gss = bestmodel(1)
              gpp = gss / bestmodel(2)
              gps = gpp * bestmodel(3)
              gsp = gps / 6._r32

              CALL update_log(num2char('Best model parameters', width=29, fill='.') + num2char('EtaPP', width=17, justify='r') + &
                              '|' + num2char('EtaPS', width=13, justify='r') + '|' + num2char('EtaSP', width=13, justify='r') +  &
                              '|' + num2char('EtaSS', width=13, justify='r') + '|', blankline = .false.)

              CALL update_log(num2char('', width=29) + num2char(gpp, notation='s', width=17, precision=3, justify='r') + '|' + &
                              num2char(gps, notation='s', width=13, precision=3, justify='r') + '|' + &
                              num2char(gsp, notation='s', width=13, precision=3, justify='r') + '|' + &
                              num2char(gss, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

            ELSE
              gss   = bestmodel(1)
              aksq  = bestmodel(2)
              kappa = bestmodel(3)

              CALL update_log(num2char('Best model parameters', width=29, fill='.') +    &
                              num2char('EtaSS', width=17, justify='r') + '|' +    &
                              num2char('Nu', width=13, justify='r')    + '|' +    &
                              num2char('Hurst', width=13, justify='r') + '|', blankline = .false.)

              CALL update_log(num2char('', width=29) + num2char(gss, notation='s', width=17, precision=3, justify='r') + '|' + &
                              num2char(aksq, notation='s', width=13, precision=3, justify='r') + '|' +   &
                              num2char(kappa, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

            ENDIF

            DEALLOCATE(sampled, fit)

          ENDIF

          ! add current event to list of already inverted events
          IF (ALLOCATED(inverted)) THEN
            inverted = [inverted, current]
          ELSE
            inverted = [current]
          ENDIF

        ELSE

          ! add event to blacklist if there are less receivers than expected due to too large windows: blacklisted events are
          ! skipped during search
          IF (ALLOCATED(blacklisted)) THEN
            blacklisted = [blacklisted, current]
          ELSE
            blacklisted = [current]
          ENDIF

        ENDIF

        DEALLOCATE(envobs, tsobs, nobs, iobs)
        IF (elastic) DEALLOCATE(tpobs)

        current = ' '

        DEALLOCATE(code)

      ELSE

        EXIT                                      !< no need to search for events if "mode" not equal 2

      ENDIF

    ENDDO                                  !< end of loop needed only if "mode = 2"

  ENDDO                                 !<  end loop over frequency bands

  CALL watch_stop(tictoc(1), mpi_comm_world)
  tictoc(1) = tictoc(1) / 60._r64

  IF (world_rank .eq. 0) CALL update_log('Program completed in' + num2char(tictoc(1), notation='f', width=10, precision=1) +  &
                                         ' minutes')

  CALL mpi_finalize(ierr)

END PROGRAM main

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

SUBROUTINE read_input_file(ok)

  ! Purpose:
  !   to read and check the INVETA input file. Name of input file is retrieved from command line. On exit, "ok" is not zero if an
  !   error occurred.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_inveta

  IMPLICIT none

  INTEGER(i32),             INTENT(OUT) :: ok
  CHARACTER(64)                         :: fo
  CHARACTER(:), ALLOCATABLE             :: str
  INTEGER(i32)                          :: lu, n, i

  !---------------------------------------------------------------------------------------------------------------------------------

  ok = 0

  CALL get_command_argument(1, fo, status = ok)                        !< get input file name from command line

  IF (ok .ne. 0) THEN
    CALL report_error('Could not read from command line')
    RETURN
  ENDIF

  OPEN(newunit = lu, file = TRIM(fo), status = 'old', form = 'formatted', access = 'sequential', action = 'read', iostat = ok)

  IF (ok .ne. 0) THEN
    CALL report_error('Error while opening file' + TRIM(fo))
    RETURN
  ENDIF

  n = recurrences(ok, lu, 'REC', com = '#')                        !< find number of receivers

  IF (ok .ne. 0) THEN
    CALL report_error(parser_error(ok))
    RETURN
  ENDIF

  IF (is_empty(n) .or. (n .le. 0)) THEN
    CALL report_error('No stations found in input file')
    ok = 1
    RETURN
  ENDIF

  ALLOCATE(recvr(n))

  ! loop over stations to store folders and read P/S-picking
  DO i = 1, n

    CALL parse(ok, str, lu, 'File', ["'", "'"], 'REC', nkey = i, com = '#')

    IF (ok .ne. 0) THEN
      CALL report_error(parser_error(ok))
      RETURN
    ENDIF

    IF (is_empty(str)) THEN
      CALL report_error('Argument "File" for receiver #' + num2char(i) + ' not found')
      ok = 1
      RETURN
    ENDIF

    ! go through file with events for current station and store its content
    CALL read_event_file(str, i, ok)

    IF (ok .ne. 0) RETURN

    CALL parse(ok, str, lu, 'Folder', ["'", "'"], 'REC', nkey = i, com = '#')

    IF (ok .ne. 0) THEN
      CALL report_error(parser_error(ok))
      RETURN
    ENDIF

    IF (is_empty(str)) THEN
      CALL report_error('Argument "Folder" for receiver #' + num2char(i) + ' not found')
      ok = 1
      RETURN
    ENDIF

    ! store argument "folder" for current receiver
    recvr(i)%folder = str

  ENDDO

  ! scan for frequency bands
  CALL parse(ok, fbands, lu, 'Bands', ['{', '}'], ',', ';', com = '#')                       !< fbands`
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (ALL(is_empty(fbands))) THEN
    CALL report_error('Argument "Bands" not found')
    ok = 1
    RETURN
  ENDIF

  ! scan for inversion parameters
  CALL parse(ok, etass, lu, 'EtaSS', ['{', '}'], ',', com = '#')                             !< etass
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (ALL(is_empty(etass))) THEN
    CALL report_error('Argument "EtaS" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, nu, lu, 'Nu', ['{', '}'], ',', com = '#')                                   !< nu
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  CALL parse(ok, mode, lu, 'Mode', ['=', ','], com = '#')                                    !< mode
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(mode)) THEN
    CALL report_error('Argument "Mode" not found')
    ok = 1
    RETURN
  ENDIF

  IF (mode .eq. 2) THEN

    CALL parse(ok, threshold, lu, 'Threshold', ['=', ','], com = '#')                       !< threshold
    IF (ok .ne. 0) CALL report_error(parser_error(ok))

    IF (is_empty(threshold)) THEN
      CALL report_error('Argument "Threshold" not found')
      ok = 1
      RETURN
    ENDIF

  ENDIF

  CALL parse(ok, sdwindow, lu, 'Swin', ['=', ','], 'DIRECT', com = '#')               !< sdwindow (for direct)
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(sdwindow)) THEN
    CALL report_error('Argument "Swin" for direct waves not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, fwin, lu, 'Factor', ['=', '%'], 'DIRECT', com = '#')                !< fwin
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  ! by default, window starts at Tp and Ts
  IF (is_empty(fwin)) fwin = 100._r32

  fwin = fwin / 100._r32

  CALL parse(ok, scwindow, lu, 'Swin', ['=', ','], 'CODA', com = '#')                 !< scwindow (for coda)
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(scwindow)) THEN
    CALL report_error('Argument "Swin" for coda waves not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, tlim, lu, 'Tlim', ['=', ','], 'CODA', com = '#')                 !< tlim (for coda)
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  ! by default, take whole time-series
  IF (is_empty(tlim)) tlim = 0._r32

  CALL parse(ok, nsi, lu, 'InitialModels', ['=', ','], com = '#')                           !< nsi
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(nsi)) THEN
    CALL report_error('Argument "InitialModels" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, ns, lu, 'Models', ['=', ','], com = '#')                                   !< ns
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(ns)) THEN
    CALL report_error('Argument "Models" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, nr, lu, 'Resampled', ['=', ','], com = '#')                                !< nr
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(nr)) THEN
    CALL report_error('Argument "Resampled" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, itermax, lu, 'Iterations', ['=', ','], com = '#')                          !< itermax
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(itermax)) THEN
    CALL report_error('Argument "Iterations" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, seed, lu, 'Seed', ['=', ','], com = '#')                                  !< seed
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(seed)) THEN
    CALL report_error('Argument "Seed" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, beta, lu, 'Beta', ['=', ','], com = '#')                                  !< beta
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(beta)) THEN
    CALL report_error('Argument "Beta" not found')
    ok = 1
    RETURN
  ENDIF

  CALL parse(ok, str, lu, 'Weight', ["'", "'"], com = '#')                                 !< noweight
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(str)) THEN
    CALL report_error('Argument "Weight" not found')
    ok = 1
    RETURN
  ENDIF

  noweight = .false.
  noweight = lowercase(str) .eq. 'n'

  ! if "nu" is not present, we are dealing with elastic RTT
  elastic = ALL(is_empty(nu))

  IF (elastic) THEN

    CALL parse(ok, etass2pp, lu, 'EtaSS/PP', ['{', '}'], ',', com = '#')                !< etass2pp
    IF (ok .ne. 0) CALL report_error(parser_error(ok))

    IF (ALL(is_empty(etass2pp))) THEN
      CALL report_error('Argument "EtaSS/PP" not found')
      ok = 1
      RETURN
    ENDIF

    CALL parse(ok, etaps2pp, lu, 'EtaPS/PP', ['{', '}'], ',', com = '#')                !< etaps2pp
    IF (ok .ne. 0) CALL report_error(parser_error(ok))

    IF (ALL(is_empty(etaps2pp))) THEN
      CALL report_error('Argument "EtaPS" not found')
      ok = 1
      RETURN
    ENDIF

    CALL parse(ok, pdwindow, lu, 'Pwin', ['=', ','], 'DIRECT', com = '#')               !< pdwindow (for direct)
    IF (ok .ne. 0) CALL report_error(parser_error(ok))

    IF (is_empty(pdwindow)) THEN
      CALL report_error('Argument "Pwin" for direct waves not found')
      ok = 1
      RETURN
    ENDIF

    CALL parse(ok, pcwindow, lu, 'Pwin', ['=', ','], 'CODA', com = '#')                 !< pcwindow (for coda)
    IF (ok .ne. 0) CALL report_error(parser_error(ok))

    IF (is_empty(pcwindow)) THEN
      CALL report_error('Argument "Pwin" for coda waves not found')
      ok = 1
      RETURN
    ENDIF

  ELSE

    ! assign some values to "acf" and "hurst" for acoustic isotropic RTT. Note: in such case, these are irrelevant
    IF (ALL(nu .eq. 0._r32)) THEN

      acf   = 'vk'
      hurst = [0.5_r32, 0.5_r32]

    ELSE

      CALL parse(ok, acf, lu, 'acf', ["'", "'"], com = '#')                                 !< acf
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(acf)) THEN
        CALL report_error('Argument "acf" not found')
        ok = 1
        RETURN
      ENDIF

      IF (acf .eq. 'vk') THEN
        CALL parse(ok, hurst, lu, 'Hurst', ['{', '}'], ',', com = '#')                     !< hurst
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (ALL(is_empty(hurst))) THEN
          CALL report_error('Argument "hurst" not found')
          ok = 1
          RETURN
        ENDIF
      ELSE
        ! for Gaussian media, set Hurst exponent to any value to avoid triggering errors
        hurst = [0.5_r32, 0.5_r32]
      ENDIF

    ENDIF

  ENDIF

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! -------------------------------------- check whether all input parameters make sense -------------------------------------------

  DO i = 1, SIZE(fbands, 2)
    IF (wrong_par(fbands(:, i), 0._r32, strict = .true.)) THEN
      CALL report_error('All frequency bands must be positive')
      ok = 1
      RETURN
    ENDIF
    IF (wrong_par(fbands(:, i), strict = .true.)) THEN
      CALL report_error('Lower bound for frequency bands must be smaller than upper bound')
      ok = 1
      RETURN
    ENDIF
  ENDDO

  IF (wrong_par([beta], 0._r32, strict = .true.)) THEN
    CALL report_error('Shear wave velocity ("beta") must be positive')
    ok = 1
    RETURN
  ENDIF

  IF (wrong_par([nsi, ns, nr, itermax], 0, strict = .true.)) THEN
    CALL report_error('NA parameters must be positive')
    ok = 1
    RETURN
  ENDIF

  IF (wrong_par([mode], 0, 2)) THEN
    CALL report_error('Mode parameter must lie in the range [0, 2]')
    ok = 1
    RETURN
  ENDIF

  IF (mode .eq. 2) THEN
    IF (wrong_par([threshold], 1)) THEN
      CALL report_error('Parameter "Threshold" must be larger than or equal to 1')
      ok = 1
      RETURN
    ENDIF
  ENDIF

  IF (wrong_par([sdwindow], 0._r32)) THEN
    CALL report_error('Parameter "Swin" for direct wave must be positive')
    ok = 1
    RETURN
  ENDIF
  IF (wrong_par([scwindow], 0._r32)) THEN
    CALL report_error('Parameter "Swin" for coda wave must be positive')
    ok = 1
    RETURN
  ENDIF
  IF (wrong_par([tlim], 0._r32)) THEN
    CALL report_error('Parameter "Tlim" for coda wave must be positive')
    ok = 1
    RETURN
  ENDIF

  IF (wrong_par([fwin], 0._r32)) THEN
    CALL report_error('Parameter "Factor" for direct wave must be positive')
    ok = 1
    RETURN
  ENDIF

  IF (wrong_par(etass, 0._r32, strict = .true.)) THEN
    CALL report_error('Parameter "EtaSS" must be positive')
    ok = 1
    RETURN
  ENDIF
  IF (wrong_par(etass)) THEN
    CALL report_error('Lower bound for parameter "EtaSS" must be smaller than or equal to upper bound')
    ok = 1
    RETURN
  ENDIF

  IF (elastic) THEN
    IF (wrong_par(etaps2pp, 0._r32, strict = .true.)) THEN
      CALL report_error('Parameter "EtaPS/PP" must be positive')
      ok = 1
      RETURN
    ENDIF
    IF (wrong_par(etaps2pp)) THEN
      CALL report_error('Lower bound for parameter "EtaPS/PP" must be smaller than or equal to upper bound')
      ok = 1
      RETURN
    ENDIF

    IF (wrong_par(etass2pp, 0._r32, strict = .true.)) THEN
      CALL report_error('Parameter "EtaSS/PP" must be positive')
      ok = 1
      RETURN
    ENDIF
    IF (wrong_par(etass2pp)) THEN
      CALL report_error('Lower bound for parameter "EtaSS/PP" must be smaller than or equal to upper bound')
      ok = 1
      RETURN
    ENDIF

    IF (wrong_par([pdwindow], 0._r32)) THEN
      CALL report_error('Parameter "Pwin" for direct wave must be positive')
      ok = 1
      RETURN
    ENDIF
    IF (wrong_par([pcwindow], 0._r32)) THEN
      CALL report_error('Parameter "Pwin" for coda wave must be positive')
      ok = 1
      RETURN
    ENDIF

  ELSE
    IF (wrong_par(nu, 0._r32)) THEN
      CALL report_error('"Nu" must be positive')
      ok = 1
      RETURN
    ENDIF
    IF (wrong_par(nu)) THEN
      CALL report_error('Lower bound for "Nu" must be smaller than or equal to upper bound')
      ok = 1
      RETURN
    ENDIF
    IF (wrong_par(hurst, 0._r32, 1._r32, .true.)) THEN
      CALL report_error('Hurst exponent must be in the interval (0, 1)')
      ok = 1
      RETURN
    ENDIF
    IF ( (acf .ne. 'vk') .and. (acf .ne. 'gs') ) THEN
      CALL report_error('Unknown value for "acf": ' + acf)
      ok = 1
      RETURN
    ENDIF
  ENDIF

END SUBROUTINE read_input_file

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

SUBROUTINE read_event_file(fo, i, ok)

  ! Purpose:
  !   to read info from file "fo" for i-th receiver containing information on each single event recorded. The file is expected to
  !   have the following format:
  !
  !   Line 1        -> header (skipped)
  !   Line 2 to EOF -> event ID, station ID, network ID, channel, direct P-wave arrival time, direct S-wave arrival time
  !
  !   Values are stored in global (module) variables. On exit, "ok" is not equal to zero if an error occurred.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_inveta

  IMPLICIT none

  CHARACTER(*),                            INTENT(IN)  :: fo
  INTEGER(i32),                            INTENT(IN)  :: i
  INTEGER(i32),                            INTENT(OUT) :: ok
  CHARACTER(10)                                        :: label
  CHARACTER(30)                                        :: id
  CHARACTER(200)                                       :: buffer
  CHARACTER(10), ALLOCATABLE, DIMENSION(:)             :: code, network, channel
  CHARACTER(30), ALLOCATABLE, DIMENSION(:)             :: event
  INTEGER(i32)                                         :: lu, l
  REAL(r32)                                            :: p, s, lat, long, z, alpha
  REAL(r32),     ALLOCATABLE, DIMENSION(:)             :: tp, ts

  !---------------------------------------------------------------------------------------------------------------------------------

  OPEN(newunit = lu, file = fo, status = 'old', form = 'formatted', access = 'sequential', action = 'read', iostat = ok)

  IF (ok .ne. 0) THEN
    CALL report_error('Error while opening file ' + fo)
    RETURN
  ENDIF

  READ(lu, *, iostat = ok)                        !< skip header

  l = 1                                           !< lines counter

  DO

    READ(lu, '(A200)', iostat = ok) buffer

    IF (ok .eq. -1) EXIT

    IF ( (ok .lt. -1) .and. (ok .gt. 0) ) THEN
      CALL report_error('Error while reading line #' + num2char(l) + ' of file ' + fo )
      RETURN
    ENDIF

    l = l + 1

    id = split(buffer, ACHAR(9), 1)

    id = split(split(id, ' ', 1), '-', 1) + split(split(id, ' ', 1), '-', 2) + split(split(id, ' ', 1), '-', 3) +   &
         split(split(id, ' ', 2), ':', 1) + split(split(id, ' ', 2), ':', 2) + split(split(id, ' ', 2), ':', 3)
    id = split(id, '.', 1)

    IF (ALLOCATED(tp)) THEN
      event = [event, id]
    ELSE
      event = [id]
    ENDIF

    label = split(buffer, ACHAR(9), 2)

    IF (ALLOCATED(tp)) THEN
      code = [code, label]
    ELSE
      code = [label]
    ENDIF

    label = split(buffer, ACHAR(9), 3)

    IF (ALLOCATED(tp)) THEN
      network = [network, label]
    ELSE
      network = [label]
    ENDIF

    label = split(buffer, ACHAR(9), 4)

    IF (ALLOCATED(tp)) THEN
      channel = [channel, label]
    ELSE
      channel = [label]
    ENDIF

    CALL char2num(split(buffer, ACHAR(9), 5), p)
    CALL char2num(split(buffer, ACHAR(9), 6), s)

    IF (ALLOCATED(tp)) THEN
      tp = [tp, p]
      ts = [ts, s]
    ELSE
      tp = [p]
      ts = [s]
    ENDIF

  ENDDO

  IF (.not.ALLOCATED(ts)) THEN
    CALL report_error('File ' + fo + ' seems to be empty or uncorrectly formatted')
    ok = 1
    RETURN
  ENDIF

  ! ALLOCATE(recvr(i)%event, source = event)
  ! ALLOCATE(recvr(i)%code,  source = code)
  ! ALLOCATE(recvr(i)%tp, source = tp)
  ! ALLOCATE(recvr(i)%ts, source = ts)
  recvr(i)%event   = event
  recvr(i)%code    = code
  recvr(i)%network = network
  recvr(i)%channel = channel
  recvr(i)%tp      = tp
  recvr(i)%ts      = ts

  CLOSE(lu, iostat = ok)

  IF (ok .ne. 0) THEN
    CALL report_error('Error while closing file ' + fo)
    RETURN
  ENDIF

END SUBROUTINE read_event_file

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

SUBROUTINE watch_start(tictoc, comm)

  ! Purpose:
  !   To start the MPI stopwatch. Timing is in double-precision. If specific communicator handle not given, mpi_comm_world is
  !   used.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions
  USE                :: mpi

  REAL(r64),    INTENT(OUT) :: tictoc                            !< initial time
  INTEGER(i32), INTENT(IN)  :: comm                              !< communicator handle
  INTEGER(i32)              :: ierr

  !-----------------------------------------------------------------------------------------------------------------------------

  CALL mpi_barrier(comm, ierr)

  tictoc = mpi_wtime()

END SUBROUTINE watch_start

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

SUBROUTINE watch_stop(tictoc, comm)

  ! Purpose:
  !   To stop the MPI stopwatch and return elapsed time. Timing is in double-precision.  If specific communicator handle not given,
  !   mpi_comm_world is used.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions
  USE                :: mpi

  REAL(r64),    INTENT(INOUT) :: tictoc                          !< elapsed time
  INTEGER(i32), INTENT(IN)    :: comm                            !< communicator handle
  INTEGER(i32)                :: ierr

  !-----------------------------------------------------------------------------------------------------------------------------

  CALL mpi_barrier(comm, ierr)

  tictoc = mpi_wtime() - tictoc

END SUBROUTINE watch_stop

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

SUBROUTINE coords2index(npts, dims, coords, fs, fe)

  ! Purpose:
  ! To return first ("fs") and last ("fe") grid index for a calling process having "coords" coordinates in a given cartesian topology,
  ! the latter characterized by "npts" points and "dims" processes.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   11/01/21                  original version
  !

  USE, NON_INTRINSIC :: m_precisions

  INTEGER(i32), DIMENSION(3), INTENT(IN)  :: npts, dims, coords
  INTEGER(i32), DIMENSION(3), INTENT(OUT) :: fs, fe
  INTEGER(i32)                            :: i

  !---------------------------------------------------------------------------------------------------------------------------------

  DO i = 1, 3
    fs(i) = 1 + INT( REAL(npts(i)) / REAL(dims(i)) * REAL(coords(i)) )
    fe(i) = INT( REAL(npts(i)) / REAL(dims(i)) * REAL(coords(i) + 1) )
  ENDDO

END SUBROUTINE coords2index

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
