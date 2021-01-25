! Purpose:
!   to solve the (weighted) linear least square problem "ln(y-x) = a - b*t" to retrieve "a" (amplitude) and "b" (attenuation) terms,
!   where "y" and "x" are observed and simulated coda and "b = Qi*velocity".
!
! Revisions:
!     Date                    Description of change
!     ====                    =====================
!   18/12/20                  original version
!   11/01/21                  added multiple communicators
!

weight(:) = 0._r32

CALL mpi_comm_rank(comm2, rank, ierr)

displs(0) = 0

DO i = 1, SIZE(pprank2) - 1
  displs(i) = displs(i - 1) + pprank2(i - 1)
ENDDO

j0 = displs(rank) + 1
j1 = j0 + pprank2(rank) - 1

! loop over observations and compute difference between observed and synthetic envelopes
DO j = j0, j1

  l = SUM(nobs(1:j))                   !< sum up number of points used for inversion
  n = iobs(j)                          !< number of points in current envelope

  ALLOCATE(time(n), envelope(n))

  DO i = 1, n
    time(i) = (i - 1) * drespl
  ENDDO

  ! solve forward problem
  IF (elastic) THEN
    CALL rtt(comm1, pprank1, time, tpobs(j) + tau, tsobs(j) + tau, gpp, gps, gsp, gss, gi, beta, wp, ws, tau, envelope, ok)
  ELSE
    CALL rtt(comm1, pprank1, time, tsobs(j) + tau, gss, gi, beta, acf, hurst, bnu, tau, envelope, ok)
  ENDIF

  k = l - nobs(j)                                        !< index of previous last point used for inversion
  p = SUM(iobs(1:j)) - iobs(j)                           !< index of last point for previous envelope

  IF (elastic) THEN

    is = NINT((tpobs(j) - pdwindow * (1._r32 - fwin)) / drespl) + 1                     !< beginning direct P-window
    is = MAX(1, is)                                                                     !< make sure we don't go below 1
    ie = is + NINT(fwin*pdwindow / drespl)                                              !< end direct P-window

    k        = k + 1
    delta(k) = LOG(mean(envobs(is + p:ie + p))) - LOG(mean(envelope(is:ie)))

    is = NINT((tsobs(j) - sdwindow * (1._r32 - fwin)) / drespl) + 1                     !< beginning direct S-window

    ! add points from just after direct P-window to just before direct S-window (e.g. P-coda)
    DO i = ie + 1, is - 1
      k        = k + 1
      delta(k) = LOG(envobs(i + p)) - LOG(envelope(i))
    ENDDO

  ENDIF

  is = NINT((tsobs(j) - sdwindow * (1._r32 - fwin)) / drespl) + 1                      !< beginning direct S-window
  is = MAX(1, is)                                                                      !< make sure we don't go below 1
  ie = is + NINT(fwin*sdwindow / drespl)                                               !< end direct S-window

  k        = k + 1
  delta(k) = LOG(mean(envobs(is + p:ie + p))) - LOG(mean(envelope(is:ie)))

  ! add points from just after direct S-window
  DO i = ie + 1, n
    k        = k + 1
    delta(k) = LOG(envobs(i + p)) - LOG(envelope(i))
  ENDDO

  DEALLOCATE(time, envelope)

ENDDO

! determine number of points (time points used for inversion) owned by each single rank in "comm2"
DO i = 0, SIZE(pprank) - 1
  j0 = displs(i) + 1
  j1 = j0 + pprank2(i) - 1
  pprank(i) = SUM(nobs(j0:j1))
ENDDO

! update displacements list
DO i = 1, SIZE(displs) - 1
  displs(i) = displs(i - 1) + pprank(i - 1)
ENDDO

! exchange data inside communicator
!CALL mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, delta, pprank, displs, mpi_real, comm2, ierr)
CALL mpi_iallgatherv(mpi_in_place, 0, mpi_datatype_null, delta, pprank, displs, mpi_real, comm2, req, ierr)

! setup weights and remaning parameters
DO j = 1, SIZE(nobs)

  l = SUM(nobs(1:j))                   !< sum up number of points used for inversion
  n = iobs(j)                          !< number of points in current envelope

  ALLOCATE(time(n))

  DO i = 1, n
    time(i) = (i - 1) * drespl
  ENDDO

  k = l - nobs(j)                                        !< index of previous last point used for inversion
  p = SUM(iobs(1:j)) - iobs(j)                           !< index of last point for previous envelope

  IF (elastic) THEN

    is = NINT((tpobs(j) - pdwindow * (1._r32 - fwin)) / drespl) + 1                     !< beginning direct P-window
    is = MAX(1, is)                                                                     !< make sure we don't go below 1
    ie = is + NINT(fwin*pdwindow / drespl)                                              !< end direct P-window

    t = 0._r32
    DO i = is, ie
      t = t + time(i) * envobs(i + p)
    ENDDO
    t = t / SUM(envobs(is + p:ie + p))                   !< balanced time point in direct P-window

    k         = k + 1
    tobs(k)   = t
    weight(k) = SQRT(REAL(ie - is + 1, r32))

    is = NINT((tsobs(j) - sdwindow * (1._r32 - fwin)) / drespl) + 1                     !< beginning direct S-window

    ! add points from just after direct P-window to just before direct S-window (e.g. P-coda)
    DO i = ie + 1, is - 1
      k         = k + 1
      tobs(k)   = time(i)
      weight(k) = 1._r32
    ENDDO

  ENDIF

  is = NINT((tsobs(j) - sdwindow * (1._r32 - fwin)) / drespl) + 1                      !< beginning direct S-window
  is = MAX(1, is)                                                                      !< make sure we don't go below 1
  ie = is + NINT(fwin*sdwindow / drespl)                                               !< end direct S-window

  t = 0._r32
  DO i = is, ie
    t = t + time(i) * envobs(i + p)
  ENDDO
  t = t / SUM(envobs(is + p:ie + p))                  !< balanced time point in S-window

  k         = k + 1
  tobs(k)   = t
  weight(k) = SQRT(REAL(ie - is + 1, r32))

  ! add points from just after direct S-window
  DO i = ie + 1, n
    k         = k + 1
    tobs(k)   = time(i)
    weight(k) = 1._r32
  ENDDO

  DEALLOCATE(time)

ENDDO

! set all weight to unity if we want ordinary LLSQ
IF (noweight) weight(:) = 1._r32

a(:,:) = 0._r32

DO j = 1, SIZE(nobs)
  ie = SUM(nobs(1:j))                 !< sum up number of points used for inversion
  is = ie - nobs(j) + 1
  DO i = is, ie
    a(i, j) = weight(i)
  ENDDO
ENDDO

j = SIZE(nobs) + 1

DO i = 1, SUM(nobs)
  a(i, j) = -tobs(i) * weight(i)
ENDDO

CALL mpi_wait(req, mpi_status_ignore, ierr)

DO i = 1, SUM(nobs)
  b(i) = delta(i) * weight(i)
ENDDO

! in output, first [SIZE(nobs) + 1] points contain solution: [ln(c1), ln(c2), ..., ln(cn), b]
CALL llsq_solver(a, b, ok)
