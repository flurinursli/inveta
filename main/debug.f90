PROGRAM driver

  ! gfortran -O2 -std=f2008 m_precisions.f90 m_strings.f90 m_fft_real.f90 m_interpolation.f90 m_llsq.f90 m_rtt.f90 hyp_2F1.f90 quadpack.f90 driver.f90 -I$FFTW_PATH/include -L$FFTW_PATH/lib -L$LAPACK_PATH -llapack -lblas -lfftw3 -lgsl -lgslcblas -cpp -DDEBUG -DDOUBLE_PREC

  USE, INTRINSIC     :: iso_fortran_env, only: stdout => output_unit
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_rtt
  USE, NON_INTRINSIC :: mpi

  IMPLICIT none

  CHARACTER(2)                                      :: acf
  CHARACTER(64)                                     :: fo
  INTEGER(i32)                                      :: i, world_size, rank, ierr, lu, ok, npts, flag
  INTEGER(i32), ALLOCATABLE, DIMENSION(:)           :: pprank
  REAL(r32)                                         :: beta, tmax, dt, tp, ts, gss, gss2pp, gps2pp, nu, hurst, alpha, gi
  REAL(r32)                                         :: gpp, gps, gsp
  REAL(r64)                                         :: tictoc
  REAL(r32),                              PARAMETER :: tau = 0.25_r32, wp = 1._r32, ws = 23.4_r32
  REAL(r32),    ALLOCATABLE, DIMENSION(:)           :: time, envelope

  !---------------------------------------------------------------------------------------------------------------------------------

  CALL mpi_init(ierr)

  CALL mpi_comm_size(mpi_comm_world, world_size, ierr)
  CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

  ALLOCATE(pprank(0:world_size-1))

  IF (rank .eq. 0) THEN

    CALL get_command_argument(1, fo, status = ok)                        !< get input file name from command line

    IF (ok .ne. 0) THEN

      WRITE(stdout, *) 'Could not read from command line'

    ELSE

      OPEN(newunit = lu, file = TRIM(fo), status = 'unknown', form = 'formatted', access = 'sequential', action = 'read',   &
           iostat = ok)

      READ(lu, *); READ(lu, *) flag            !< 0=acoustic, 1=elastic
      READ(lu, *); READ(lu, *) beta
      READ(lu, *); READ(lu, *) acf
      READ(lu, *); READ(lu, *) tmax
      READ(lu, *); READ(lu, *) dt
      READ(lu, *); READ(lu, *) tp
      READ(lu, *); READ(lu, *) ts
      READ(lu, *); READ(lu, *) gss
      READ(lu, *); READ(lu, *) gss2pp
      READ(lu, *); READ(lu, *) gps2pp
      READ(lu, *); READ(lu, *) nu
      READ(lu, *); READ(lu, *) hurst

    ENDIF

  ENDIF

  CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

  IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)

  CALL mpi_bcast(flag, 1, mpi_int, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(beta, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(acf, 2, mpi_character, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(tmax, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(dt, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(tp, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(ts, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(gss, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(gss2pp, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(gps2pp, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nu, 1, mpi_real, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(hurst, 1, mpi_real, 0, mpi_comm_world, ierr)

  npts = NINT(tmax / dt) + 1

  ALLOCATE(time(npts), envelope(npts))

  DO i = 1, npts
    time(i) = (i - 1) * dt
  ENDDO

  alpha = ts / tp * beta

  ! gi = b / beta
  gi = 0._r32

  CALL split_task(query_fft_size(), world_size, pprank)

  CALL mpi_barrier(mpi_comm_world, ierr)
  tictoc = mpi_wtime()

  IF (flag .eq. 0) THEN
    CALL rtt(mpi_comm_world, pprank, time, ts + tau, gss, gi, beta, acf, hurst, nu, tau, envelope, ok)
  ELSE
    gpp = gss / gss2pp
    gps = gps2pp * gpp
    gsp = gps / 6._r32
    CALL rtt(mpi_comm_world, pprank, time, tp + tau, ts + tau, gpp, gps, gsp, gss, gi, beta, wp, ws, tau, envelope, ok)
  ENDIF

  CALL mpi_barrier(mpi_comm_world, ierr)
  tictoc = mpi_wtime() - tictoc

  !CALL rtt(REAL(time, r64), REAL(ts, r64), REAL(beta, r64), REAL(gss, r64), REAL(gi, r64), denv)

  PRINT*, 'Elapsed time: ', REAL(tictoc, r32)

  IF (rank .eq. 0) THEN
    OPEN(newunit = lu, file = 'debug_envelope.txt', status = 'unknown', form = 'formatted', access = 'sequential', &
         action = 'write', iostat = ok)
    DO i = 1, npts
      WRITE(lu, '(G0, 3X, G0)') time(i), envelope(i)
    ENDDO
    CLOSE(lu)
  ENDIF

  CALL mpi_finalize(ierr)

  CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE split_task(npts, ntasks, pptask)

      ! Purpose:
      !   To evenly distribute elements of vector amongst processes, returning first and last index for each process.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(i32),                        INTENT(IN)  :: npts                       !< number of elements to be distributed
      INTEGER(i32),                        INTENT(IN)  :: ntasks                     !< number of MPI processes
      INTEGER(i32), DIMENSION(0:ntasks-1), INTENT(OUT) :: pptask                     !< points per task
      INTEGER(i32)                                     :: p, i0, i1

      !-----------------------------------------------------------------------------------------------------------------------------

      DO p = 0, ntasks - 1
        i0 = 1 + INT( REAL(npts, r32) / REAL(ntasks, r32) * REAL(p, r32) )
        i1 = INT( REAL(npts, r32) / REAL(ntasks, r32) * REAL(p + 1, r32) )
        pptask(p) = i1 - i0 + 1
      ENDDO

    END SUBROUTINE split_task


END PROGRAM driver
