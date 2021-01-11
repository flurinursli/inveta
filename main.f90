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
  USE                :: mpi

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
  INTEGER(i32)                              :: nsi, ns, nr, itermax, seed                 !< NA parameters
  INTEGER(i32), ALLOCATABLE, DIMENSION(:)   :: iobs, nobs                                 !< observations per recording
  INTEGER(i32), ALLOCATABLE, DIMENSION(:)   :: pprank1, pprank2, pprank3
  INTEGER(i32), ALLOCATABLE, DIMENSION(:,:) :: gs, ge                                     !< global indices for cartesian topology
  LOGICAL                                   :: elastic                                    !< scalar or elastic RTT
  LOGICAL                                   :: noweight                                   !< disable weighted linear regression
  REAL(r32)                                 :: hurst                                      !< Hurst exponent when acf=Von Karman
  REAL(r32)                                 :: drespl
  REAL(r32)                                 :: beta
  REAL(r32)                                 :: fwin
  REAL(r32)                                 :: pdwindow, sdwindow, pcwindow, scwindow     !< windows for direct and coda P-/S-waves
  REAL(r32),    ALLOCATABLE, DIMENSION(:)   :: etass, etass2pp, etaps2pp, nu              !< scattering parameters
  REAL(r32),    ALLOCATABLE, DIMENSION(:)   :: tobs, envobs, tpobs, tsobs                 !< observables used during inversion
  REAL(r32),    ALLOCATABLE, DIMENSION(:,:) :: fbands                                     !< frequency bands for inversion

  TYPE :: info
    CHARACTER(200)                           :: folder
    CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: code
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
      !   to write to disk the parameters space explored by NA, including average misfit values.
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
      INTEGER(i32)                                                             :: rank, ierr
      INTEGER(i32),               DIMENSION(0:SIZE(pprank2)-1)                 :: displs
      REAL(r32)                                                                :: gpp, gps, gsp, gss, gi, bnu, t, const
      REAL(r32),                                                    PARAMETER  :: tau = 0.25_r32, wp = 1._r32, ws = 23.4_r32
      REAL(r32),     ALLOCATABLE, DIMENSION(:)                                 :: time, envelope
      REAL(r32),                  DIMENSION(SUM(nobs))                         :: delta, weight, b, tobs
      REAL(r32),                  DIMENSION(SUM(nobs),SIZE(nobs)+1)            :: a

      !-----------------------------------------------------------------------------------------------------------------------------

      lf = NINT(fbands(1, f))
      hf = NINT(fbands(2, f))

      IF (elastic) THEN
        gss = bestmodel(1)
        gpp = gss / bestmodel(2)
        gps = gpp * bestmodel(3)
        gsp = gps / 6._r32
      ELSE
        gss = bestmodel(1)
        bnu = bestmodel(2)
      ENDIF

      gi = 0._r32

#include "linsys_incl.f90"

      n = SIZE(nobs)

      gi = b(n + 1) / beta

      IF (world_rank .eq. 0) THEN

        ! report best fitting parameters to a file
        IF (mode .le. 1) THEN

          fo = 'bestpar_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(1))  + '.txt'

          OPEN(newunit = lu, file = fo, status = 'unknown', form = 'formatted', access = 'sequential', position = 'append',    &
               action = 'write', iostat = ok)

          IF (elastic) THEN
            WRITE(lu, *) gpp, gps, gsp, gss, gi
          ELSE
            WRITE(lu, *) gss, bnu, gi
          ENDIF

          CLOSE(lu)

        ELSE

          DO j = 1, SIZE(nobs)

            fo = 'bestpar_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(j)) + '.txt'

            OPEN(newunit = lu, file = fo, status = 'unknown', form = 'formatted', access = 'sequential', position = 'append',    &
                 action = 'write', iostat = ok)

            IF (elastic) THEN
              WRITE(lu, *) gpp, gps, gsp, gss, gi
            ELSE
              WRITE(lu, *) gss, bnu, gi
            ENDIF

            CLOSE(lu)

          ENDDO

        ENDIF

      ENDIF

      ! write also envelopes associated to best fitting parameters
      DO j = 1, SIZE(nobs)

        n = iobs(j)

        ALLOCATE(time(n), envelope(n))

        DO i = 1, n
          time(i) = (i - 1) * drespl
        ENDDO

        IF (elastic) THEN
          CALL rtt(comm3, pprank3, time, tpobs(j) + tau, tsobs(j) + tau, gpp, gps, gsp, gss, gi, beta, wp, ws, tau, envelope, ok)
        ELSE
          CALL rtt(comm3, pprank3, time, tsobs(j) + tau, gss, gi, beta, acf, hurst, bnu, tau, envelope, ok)
        ENDIF

        IF (mode .le. 1) THEN
          fo = 'bestfit_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(1)) + '_' + TRIM(inverted(j)) + '.txt'
        ELSE
          fo = 'bestfit_' + num2char(lf) + '-' + num2char(hf) + '_' + TRIM(code(j)) + '_' + TRIM(inverted(1)) + '.txt'
        ENDIF

        IF (world_rank .eq. 0) THEN

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

    END SUBROUTINE bestfit

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION misfit(nd, mpar)

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
      INTEGER(i32)                                                            :: i, j, j0, j1, k, l, n, p, is, ie, ok, rank, ierr
      INTEGER(i32),              DIMENSION(0:SIZE(pprank2)-1)                 :: displs
      REAL(r32)                                                               :: gss, gsp, gpp, gps, bnu, t
      REAL(r32),                                                   PARAMETER  :: gi = 0._r32, tau = 0.25_r32
      REAL(r32),                                                   PARAMETER  :: wp = 1._r32, ws = 23.4_r32
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                                 :: time, envelope
      REAL(r32),                 DIMENSION(SUM(nobs))                         :: delta, weight, b, tobs
      REAL(r32),                 DIMENSION(SUM(nobs),SIZE(nobs)+1)            :: a

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (elastic) THEN
        gss = mpar(1)
        gpp = gss / mpar(2)
        gps = gpp * mpar(3)
        gsp = gps / 6._r32
      ELSE
        gss = mpar(1)
        bnu = mpar(2)
      ENDIF

#include "linsys_incl.f90"

      n = SIZE(nobs)

      misfit = 0._r32

      ! misfit is defined as "LN(obs) - (LN(syn) + LN(a) -b*t)"
      ! remember that "b = Qsi*beta = Qpi*alpha"
      DO j = 1, n
        ie = SUM(nobs(1:j))
        is = ie - nobs(j) + 1
        DO i = is, ie
          misfit = misfit + (delta(i) - b(j) + b(n + 1) * tobs(i))**2
        ENDDO
      ENDDO

    END FUNCTION misfit

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

    SUBROUTINE build_comms(nx, ny, nz)

      INTEGER(i32),                             INTENT(IN) :: nx, ny, nz
      INTEGER(i32)                                         :: cartopo, ierr
      INTEGER(i32),                   PARAMETER            :: ndims = 3
      INTEGER(i32), DIMENSION(ndims)                       :: dims, npts, coords
      LOGICAL,                        PARAMETER            :: reorder = .true.
      LOGICAL,      DIMENSION(ndims), PARAMETER            :: isperiodic = [.false., .false., .false.]

      !-----------------------------------------------------------------------------------------------------------------------------

      ! logic grid has "ns" points along x, number of observations along y and FFT points (for RTT) along z
      npts = [nx, ny, nz]

      ! define number of processes along each direction
      dims(1) = MIN(world_size, nx)
      dims(2) = MIN(world_size / dims(1), ny)                       !< equal to 1 when "mode=0"
      dims(3) = world_size / (dims(1) * dims(2))

if (world_rank == 0) print*, 'proc grid ', dims

      ! create topology
      CALL mpi_cart_create(mpi_comm_world, ndims, dims, isperiodic, reorder, cartopo, ierr)

      ! return process coordinates in current topology
      CALL mpi_cart_coords(cartopo, world_rank, 3, coords, ierr)

      ! return first/last index along each direction for the calling process. Note: first point has first index equal to 1.
      CALL coords2index(npts, dims, coords, gs(:, world_rank), ge(:, world_rank))

      ! make all processes aware of global indices
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 3, mpi_integer, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 3, mpi_integer, mpi_comm_world, ierr)

      ! release cartesian communicator
      CALL mpi_comm_free(cartopo, ierr)

      ! build all three communicators
      CALL build_pencil(0, comm1, pprank1)
      CALL build_pencil(1, comm2, pprank2)
      CALL build_pencil(2, comm3, pprank3)

    END SUBROUTINE build_comms

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

    SUBROUTINE build_pencil(dir, newcomm, n)

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

      INTEGER(i32),                                        INTENT(IN)  :: dir
      INTEGER(i32),                                        INTENT(OUT) :: newcomm
      INTEGER(i32), ALLOCATABLE, DIMENSION(:),             INTENT(OUT) :: n
      INTEGER(i32)                                                     :: i, ierr, rank, ntasks
      INTEGER(i32),              DIMENSION(0:world_size-1)             :: color
      LOGICAL,                   DIMENSION(2)                          :: bool

      !-----------------------------------------------------------------------------------------------------------------------------

      ! group processes into pencils
      DO i = 0, world_size - 1

        color(i) = 0

        ! pencil oriented along x-axis
        IF (dir .eq. 0) THEN
          bool(1) = (gs(2, i) .eq. gs(2, world_rank)) .and. (ge(2, i) .eq. ge(2, world_rank))
          bool(2) = (gs(3, i) .eq. gs(3, world_rank)) .and. (ge(3, i) .eq. ge(3, world_rank))
        ! pencil oriented along y-axis
        ELSEIF (dir .eq. 1) THEN
          bool(1) = (gs(1, i) .eq. gs(1, world_rank)) .and. (ge(1, i) .eq. ge(1, world_rank))
          bool(2) = (gs(3, i) .eq. gs(3, world_rank)) .and. (ge(3, i) .eq. ge(3, world_rank))
        ! pencil oriented along z-axis
        ELSEIF (dir .eq. 2) THEN
          bool(1) = (gs(1, i) .eq. gs(1, world_rank)) .and. (ge(1, i) .eq. ge(1, world_rank))
          bool(2) = (gs(2, i) .eq. gs(2, world_rank)) .and. (ge(2, i) .eq. ge(2, world_rank))
        ENDIF

        IF (ALL(bool .eqv. .true.)) color(i) = i + 1

      ENDDO

      ! process belonging to the same pencil have same color
      color(world_rank) = MAXVAL(color, dim = 1)

      ! create communicator subgroup
      CALL mpi_comm_split(mpi_comm_world, color(world_rank), world_rank, newcomm, ierr)

      ! process id and communicator size
      CALL mpi_comm_rank(newcomm, rank, ierr)
      CALL mpi_comm_size(newcomm, ntasks, ierr)

      ALLOCATE(n(0:ntasks - 1))

      ! number of points along pencil direction for calling process
      n(rank) = ge(dir + 1, world_rank) - gs(dir + 1, world_rank) + 1

      ! make whole communicator aware of points for each process
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, n, 1, mpi_integer, newcomm, ierr)

    END SUBROUTINE build_pencil

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
  INTEGER(i32)                               :: ierr, ok, lu, n, i, j, k, p, l, is, ie, pts, comm
  INTEGER(i32),               DIMENSION(2)   :: v
  REAL(r32)                                  :: dt, t, gss, gpp, gsp, gps, bnu
  REAL(r64)                                  :: tic
  REAL(r32),     ALLOCATABLE, DIMENSION(:)   :: time, trespl, h, envlp, respl, bestmodel, lrange, hrange, sampled, fit
  REAL(r32),     ALLOCATABLE, DIMENSION(:,:) :: timeseries

  !---------------------------------------------------------------------------------------------------------------------------------

  CALL mpi_init(ierr)

  CALL mpi_comm_size(mpi_comm_world, world_size, ierr)
  CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

  ! "gs"/"ge" are specified in "m_inveta" and store first/last global index along each direction for all processes
  ALLOCATE(gs(3, 0:world_size-1), ge(3, 0:world_size-1))

  IF (world_rank .eq. 0) CALL set_log_module(ok, screen = .true.)

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
                      num2char('S', width=13, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(pcwindow, notation = 'f', width=17, precision=3, justify='r') + '|' +  &
                      num2char(scwindow, notation='f', width=13, precision=3, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('Inversion parameters', width=29, fill='.') + num2char('EtaSS', width=17, justify='r') + '|' + &
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
      CALL update_log(num2char('Window width S-coda waves', width=29, fill='.') + num2char('S', width=17, justify='r') + '|',  &
                      blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(scwindow, notation = 'f', width=17, precision=3, justify='r') + '|',   &
                      blankline = .false.)
      CALL update_log(num2char('Inversion parameters', width=29, fill='.') + num2char('EtaSS', width=17, justify='r') + '|' +  &
                      num2char('Nu', width=13, justify='r') + '|', blankline = .false.)
      CALL update_log(num2char('', width=29) + num2char(etass(1), notation='s', width=8, precision=1, justify='r') + ',' +  &
                      num2char(etass(2),    notation='s', width=8, precision=1) + '|' +    &
                      num2char(nu(1), notation='f', width=6, precision=1, justify='r') + ',' +      &
                      num2char(nu(2), notation='f', width=6, precision=1) + '|' , blankline = .false.)
      IF (ANY(nu .ne. 0._r32)) THEN
        CALL update_log(num2char('acf of choice', width=29, fill='.') + num2char(acf, width=17, justify='r') + '|',    &
                        blankline = .false.)
        IF (acf .eq. 'vk')  &
        CALL update_log(num2char('Hurst exponent', width=29, fill='.') + num2char(hurst, notation='f', width=17, precision=3,   &
                        justify='r') + '|', blankline = .false.)
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

  IF (elastic) THEN
    IF (world_rank .ne. 0) ALLOCATE(etass2pp(2), etaps2pp(2))
    CALL mpi_bcast(etass2pp, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(etaps2pp, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(pdwindow, 1, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(pcwindow, 1, mpi_real, 0, mpi_comm_world, ierr)
  ELSE
    IF (world_rank .ne. 0) ALLOCATE(nu(2))
    IF (world_rank .ne. 0) ALLOCATE(CHARACTER(2) :: acf)
    CALL mpi_bcast(nu, 2, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(hurst, 1, mpi_real, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(acf, 2, mpi_character, 0, mpi_comm_world, ierr)
  ENDIF

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ----------------------------------------------------- parameters inversion -----------------------------------------------------

  CALL setup_interpolation('linear', 'copy', ok)

  IF (elastic) THEN
    lrange = [etass(1), etass2pp(1), etaps2pp(1)]
    hrange = [etass(2), etass2pp(2), etaps2pp(2)]
  ELSE
    lrange = [etass(1), nu(1)]
    hrange = [etass(2), nu(2)]
  ENDIF

  ALLOCATE(bestmodel(SIZE(lrange)))

  ! arrange available processes into logic grid
  IF (mode .eq. 0) CALL build_comms(ns, 1, query_fft_size())

  DO k = 1, SIZE(fbands, 2)                         !< loop over frequency bands

    drespl = 0.25_r32 / fbands(2, k)                !< Nyquist frequency after filtering is twice low-pass frequency

    current = ' '

    DO                                              !< loop needed only if "mode = 2"

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
            ! this is a hack to deal with irregular file naming convention
            DO i = 1, 4
              fo = TRIM(recvr(l)%folder) + '/' + TRIM(recvr(l)%event(p)) + '_CH_' + TRIM(recvr(l)%code(p)) + '_HH[ENZ]_' +  &
                   num2char(i) + '.mseed.ascii'
              CALL read_miniseed(ok, fo, dt, timeseries)
              IF (ok .eq. 0) EXIT
              fo = TRIM(recvr(l)%folder) + '/' + TRIM(recvr(l)%event(p)) + '_CH_' + TRIM(recvr(l)%code(p)) + '_HG[ENZ]_' +  &
                   num2char(i) + '.mseed.ascii'
              CALL read_miniseed(ok, fo, dt, timeseries)
              IF (ok .eq. 0) EXIT
              fo = TRIM(recvr(l)%folder) + '/' + TRIM(recvr(l)%event(p)) + '_CH_' + TRIM(recvr(l)%code(p)) + '_EH[ENZ]_' +  &
                   num2char(i) + '.mseed.ascii'
              CALL read_miniseed(ok, fo, dt, timeseries)
              IF (ok .eq. 0) EXIT
            ENDDO
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

          CALL watch_start(tic, mpi_comm_self)

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! -------------------------------- inversion of each single coda recorded by a receiver ----------------------------------

          IF (mode .eq. 0) THEN
            IF (ALLOCATED(code)) THEN             !< obviously, do something only if receiver has recorded current event

              IF (world_rank .eq. 0) CALL update_log('Inverting event ' + TRIM(current) + ' for receiver ' + TRIM(code(1)) +  &
                                               ' in the frequency band ' + num2char(fbands(1,k), notation='f', precision=1) +  &
                                               '-' + num2char(fbands(2,k), notation='f', precision=1) + 'Hz')

              CALL na(comm1, misfit, lrange, hrange, seed, itermax, nsi, ns, nr, 0, bestmodel, sampled, fit, ok)

              ! write best-fitting model to disk
              CALL bestfit(k, [recvr(l)%code(1)], [current], bestmodel)

              IF (world_rank .eq. 0) THEN

                ! write explored parameters space to disk
                CALL write_search(k, [recvr(l)%code(1)], [current], sampled, fit)

                CALL update_log(num2char('Models explored', width=29, fill='.') + num2char(SIZE(fit),width=17, justify='r') + '|')

                CALL watch_stop(tic, mpi_comm_self)

                CALL update_log(num2char('Elapsed time (s)', width=29, fill='.') + num2char(tic, notation = 'f', width = 17,  &
                                precision = 1, justify='r') + '|', blankline = .false.)

                ! update log file
                IF (elastic) THEN
                  gss = bestmodel(1)
                  gpp = gss / bestmodel(2)
                  gps = gpp * bestmodel(3)
                  gsp = gps / 6._r32

                  CALL update_log(num2char('Best model parameters', width=29, fill='.') + num2char('EtaPP', width=17, justify='r' &
                                  ) + '|' + num2char('EtaPS', width=13, justify='r') + '|' + num2char('EtaSP', width=13,    &
                                  justify='r') + '|' + num2char('EtaSS', width=13, justify='r') + '|', blankline = .false.)

                  CALL update_log(num2char('', width=29) + num2char(gpp, notation='s', width=17, precision=3, justify='r') + '|' + &
                                  num2char(gps, notation='s', width=13, precision=3, justify='r') + '|' + &
                                  num2char(gsp, notation='s', width=13, precision=3, justify='r') + '|' + &
                                  num2char(gss, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

                ELSE
                  gss = bestmodel(1)
                  bnu = bestmodel(2)

                  CALL update_log(num2char('Best model parameters', width=29, fill='.') + num2char('EtaSS', width=17, justify='r' &
                                  ) + '|' + num2char('Nu', width=13, justify='r') + '|', blankline = .false.)

                  CALL update_log(num2char('', width=29) + num2char(gss, notation='s', width=17, precision=3, justify='r') + '|' + &
                                  num2char(bnu, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

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

            IF (world_rank .eq. 0) CALL update_log('Inverting ' + num2char(SIZE(inverted)) + ' event(s) for receiver ' +  &
                                              TRIM(recvr(l)%code(1)) + ' in the frequency band ' + num2char(fbands(1,k),  &
                                              notation='f', precision=1) + '-' + num2char(fbands(2,k), notation='f', precision=1)+ &
                                              'Hz')

            ! arrange available processes into logic grid
            CALL build_comms(ns, SIZE(nobs), query_fft_size())

            CALL na(comm1, misfit, lrange, hrange, seed, itermax, nsi, ns, nr, 0, bestmodel, sampled, fit, ok)

            ! write best-fitting model to disk
            CALL bestfit(k, [recvr(l)%code(1)], inverted, bestmodel)

            CALL destroy_comms()

            IF (world_rank .eq. 0) THEN

              ! write explored parameters space to disk
              CALL write_search(k, [recvr(l)%code(1)], inverted, sampled, fit)

              CALL update_log(num2char('Models explored', width=29, fill='.') + num2char(SIZE(fit),width=17, justify='r') + '|')

              CALL watch_stop(tic, mpi_comm_self)

              CALL update_log(num2char('Elapsed time (s)', width=29, fill='.') + num2char(tic, notation = 'f', width = 17,  &
                              precision = 1, justify='r') + '|', blankline = .false.)

              IF (elastic) THEN
                gss = bestmodel(1)
                gpp = gss / bestmodel(2)
                gps = gpp * bestmodel(3)
                gsp = gps / 6._r32

                CALL update_log(num2char('Best model parameters', width=29, fill='.') + num2char('EtaPP', width=17, justify='r'  &
                                ) + '|' + num2char('EtaPS', width=13, justify='r') + '|' + num2char('EtaSP', width=13,    &
                                justify='r') + '|' + num2char('EtaSS', width=13, justify='r') + '|', blankline = .false.)

                CALL update_log(num2char('', width=29) + num2char(gpp, notation='s', width=17, precision=3, justify='r') + '|' + &
                                num2char(gps, notation='s', width=13, precision=3, justify='r') + '|' + &
                                num2char(gsp, notation='s', width=13, precision=3, justify='r') + '|' + &
                                num2char(gss, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

              ELSE
                gss = bestmodel(1)
                bnu = bestmodel(2)

                CALL update_log(num2char('Best model parameters', width=29, fill='.') + num2char('EtaSS', width=17, justify='r') + &
                                '|' + num2char('Nu', width=13, justify='r') + '|', blankline = .false.)

                CALL update_log(num2char('', width=29) + num2char(gss, notation='s', width=17, precision=3, justify='r') + '|' +  &
                                num2char(bnu, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

              ENDIF
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

          ! arrange available processes into logic grid
          CALL build_comms(ns, SIZE(nobs), query_fft_size())

          CALL na(comm1, misfit, lrange, hrange, seed, itermax, nsi, ns, nr, 0, bestmodel, sampled, fit, ok)

          ! write best-fitting model to disk and explored parameters space to disk
          CALL bestfit(k, code, [current], bestmodel)

          CALL destroy_comms()

          IF (world_rank .eq. 0) THEN

            ! write explored parameters space to disk
            CALL write_search(k, code, [current], sampled, fit)

            CALL update_log(num2char('Models explored', width=29, fill='.') + num2char(SIZE(fit),width=17, justify='r') + '|')

            CALL watch_stop(tic, mpi_comm_self)

            CALL update_log(num2char('Elapsed time (s)', width=29, fill='.') + num2char(tic, notation = 'f', width = 17,  &
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
              gss = bestmodel(1)
              bnu = bestmodel(2)

              CALL update_log(num2char('Best model parameters', width=29, fill='.') + num2char('EtaSS', width=17, justify='r') + &
                              '|' + num2char('Nu', width=13, justify='r') + '|', blankline = .false.)

              CALL update_log(num2char('', width=29) + num2char(gss, notation='s', width=17, precision=3, justify='r') + '|' + &
                              num2char(bnu, notation='s', width=13, precision=3, justify='r') + '|', blankline = .false.)

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

  IF (mode .eq. 0) CALL destroy_comms()

  DEALLOCATE(gs, ge)

  IF (world_rank .eq. 0) CALL update_log('Program completed')

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

  IF (is_empty(fwin)) fwin = 0._r32

  fwin = fwin / 100._r32

  CALL parse(ok, scwindow, lu, 'Swin', ['=', ','], 'CODA', com = '#')                 !< scwindow (for coda)
  IF (ok .ne. 0) CALL report_error(parser_error(ok))

  IF (is_empty(scwindow)) THEN
    CALL report_error('Argument "Swin" for coda waves not found')
    ok = 1
    RETURN
  ENDIF

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

    ! assign some values to "acf" and "hurst" for acoustic isotropic RTT
    IF (ALL(nu .eq. 0._r32)) THEN

      acf   = 'vk'
      hurst = 0.5_r32

    ELSE

      CALL parse(ok, acf, lu, 'acf', ["'", "'"], com = '#')                                 !< acf
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(acf)) THEN
        CALL report_error('Argument "acf" not found')
        ok = 1
        RETURN
      ENDIF

      IF (acf .eq. 'vk') THEN
        CALL parse(ok, hurst, lu, 'Hurst', ['=', ','], com = '#')                          !< hurst
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (is_empty(hurst)) THEN
          CALL report_error('Argument "hurst" not found')
          ok = 1
          RETURN
        ENDIF
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
    IF (wrong_par([hurst], 0._r32, 1._r32, .true.)) THEN
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
  !   Line 2 to EOF -> event ID, station ID, direct P-wave arrival time, direct S-wave arrival time
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
  CHARACTER(10), ALLOCATABLE, DIMENSION(:)             :: code
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

    label = split(buffer, ACHAR(9), 2)

    CALL char2num(split(buffer, ACHAR(9), 3), p)
    CALL char2num(split(buffer, ACHAR(9), 4), s)

    IF (ALLOCATED(tp)) THEN
      event = [event, id]
      code  = [code, label]
      tp = [tp, p]
      ts = [ts, s]
    ELSE
      event = [id]
      code  = [label]
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
  recvr(i)%event = event
  recvr(i)%code  = code
  recvr(i)%tp    = tp
  recvr(i)%ts    = ts

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
