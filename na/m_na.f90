MODULE m_na

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: mpi

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: na

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32) :: iseed, ic, nupd, ncald
  LOGICAL      :: timing
  REAL(r32)    :: tup, tcd, tres, tdev, taxis, tna

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --


  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE na(comm, misfun, lrange, hrange, seed, itermax, nsamplei, nsample, ncells, scaling, bestmodel, na_models, misfit, ok)

      INTEGER(i32),                                       INTENT(IN)  :: comm                               !< mpi communicator
      REAL(r32),                                          EXTERNAL    :: misfun                             !< misfit function
      REAL(r32),                 DIMENSION(:),            INTENT(IN)  :: lrange
      REAL(r32),                 DIMENSION(SIZE(lrange)), INTENT(IN)  :: hrange
      INTEGER(i32),                                       INTENT(IN)  :: seed, itermax, nsamplei, nsample, ncells, scaling
      REAL(r32),                 DIMENSION(SIZE(lrange)), INTENT(OUT) :: bestmodel
      REAL(r32),    ALLOCATABLE, DIMENSION(:),            INTENT(OUT) :: na_models
      REAL(r32),    ALLOCATABLE, DIMENSION(:),            INTENT(OUT) :: misfit
      INTEGER(i32),                                       INTENT(OUT) :: ok

      INTEGER(i32)                                                    :: nmodels, ntasks, rank, ierr, i, j, i0, i1, nclean
      INTEGER(i32)                                                    :: k, color, intercomm, mycomm
      INTEGER(i32)                                                    :: nd, ntot, ns, nsleep, mopt
      INTEGER(i32),              PARAMETER                            :: nh_max = 1000, nsleep_max = 1, maxseq = 50, nd_max = 1024
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)                         :: mfitord, iwork_NA1, iwork_NA2

      LOGICAL                                                         :: restartNA
      LOGICAL                                                         :: lroot

      REAL(r32)                                                       :: misfitval, mfitmin, mfitmean, mfitminc
      REAL(r64)                                                       :: tictoc, tic, time_fwd
      REAL(r32),                 DIMENSION(nd_max+1)                  :: scales
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                         :: sum, xcur, work_NA1, work_NA2, na_model
      REAL(r32),                 DIMENSION(2,nd_max)                  :: range, ranget



      REAL(r32) :: tmis, ttfwdelapsed, ttfor, minfwd, maxfwd, t1, t2, t3, tcdt, tdevt, telap, tnat, trest, tupt, taxist

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      CALL mpi_comm_size(comm, ntasks, ierr)
      CALL mpi_comm_rank(comm, rank, ierr)

      ! CALL watch_start(tictoc)

      !CALL user_init(nd,nd_max,range,scales)
      ! assign input parameters to NA arrays

      SELECT CASE(scaling)
      CASE(0)
        scales(1) = 0._r32
      CASE(1)
        scales(1) = -1._r32
      END SELECT

      nd = SIZE(lrange)

      DO i = 1, nd
        range(1, i) = lrange(i)
        range(2, i) = hrange(i)
      ENDDO

      IF (nd .le. 0) THEN
        ok = 1
        RETURN
      ENDIF

      IF (nd .gt. nd_max) THEN
        ok = 2
        RETURN
      ENDIF

      IF ( (ncells .gt. nsample) .or. (ncells .gt. nsamplei) ) THEN
        ok = 3
        RETURN
      ENDIF

      ! set parameters manually to replace CALL to "na_options"

      nsleep = 1

      iseed = ABS(seed)           !< ran3 assume negative seed numbers for initialisation

      nclean = 500

      timing = .false.

      nmodels = nsamplei + itermax*nsample

      ALLOCATE(misfit(nmodels), na_models(nd*nmodels), na_model(nd), xcur(nd), sum(max(nsample, nsamplei)))
      ALLOCATE(mfitord(nmodels), work_NA1(nmodels), work_NA2(nmodels))
      ALLOCATE(iwork_NA1(nmodels), iwork_NA2(nmodels))

      CALL initialize(range, ranget, scales, nmodels, misfit, na_models, nd, xcur, nsample, ncells, restartNA)

      CALL initial_sample(na_models, nd, ranget, nsamplei)

      ntot  = 0
      ncald = 0
      nupd  = 0
      time_fwd = 0._r64
      tic = 0._r64
      ns = nsamplei

      DO j = 1, itermax + 1

        DO i = 1, ns

          i0 = 1 + (i - 1 + ntot) * nd
          i1 = i0 + nd - 1

          CALL transform2raw(na_models(i0:i1), nd, range, scales, na_model)

          i0 = ntot + i
          misfit(i0) = misfun(nd, na_model)

        ENDDO

        ! CALL watch_stop(tic, mpi_comm_self)

        time_fwd = time_fwd + tic

        ! CALL mpi_allreduce(mpi_in_place, misfit(ntot + 1), ns, mpi_real, mpi_sum, comm, ierr)

        CALL statistics(misfit, ns, j, ntot, mfitmin, mfitminc, mfitmean, mopt, ncells, work_NA2, iwork_NA1, iwork_NA2, mfitord, &
                        ierr)

        ok = MAX(ok, ierr)

        i0 = 1 + (mopt - 1) * nd
        i1 = i0 + nd - 1

        CALL transform2raw(na_models(i0:i1), nd, range, scales, na_model)

        DO i = 1, nd
          bestmodel(i) = na_model(i)
        ENDDO

        ntot = ntot + ns
        ns = nsample

        ! t1 = cputime(t2,t3)
        ! tmis = tmis + t2
        ! tnat = tnat + t2

        IF (j .eq. itermax + 1) EXIT

        CALL sampling(na_models, ntot, nsample, nd, nsleep, ncells, misfit, mfitord, ranget, xcur, restartNA, nclean, work_NA1)

      ENDDO

      DO i = 1, ntot
        i0 = 1 + (i - 1) * nd
        i1 = i0 + nd - 1
        CALL transform2raw(na_models(i0:i1), nd, range, scales, na_models(i0:i1))
      ENDDO

      ! at this point all processes have same "na_models", "misfit" and "bestmodel"

      !IF (lroot) CALL writefun(nd, ntot, na_models, misfit)

      ! CALL watch_start(tictoc)

      ! IF (rank .eq. 0) THEN
      !
      !   WRITE(*,*)' '
      !   WRITE(*,*)' '
      !   WRITE(*,*)' '
      !   WRITE(*,*)'Performance statistics'
      !   WRITE(*,*)' '
      !   WRITE(*,*)'Total number of full dlist evaluations',ncald
      !   WRITE(*,*)'Total number of partial dlist updates ',nupd
      !   WRITE(*,*)'Lowest misfit found                   ',mfitmin
      !   WRITE(*,*)'Average misfit over all models        ',mfitmean
      !   WRITE(*,*)'Index of lowest misfit model          ',mopt
      !   WRITE(*,*)
      !
      !   WRITE(*,*) 'Total elapsed time for entire NA                   ', REAL(tictoc, r32)
      !   WRITE(*,*) 'Total cpu time for forward modeling only           ', REAL(time_fwd, r32)
      !
      !   WRITE(*,*)
      !
      ! ENDIF

      CALL mpi_barrier(comm, ierr)


    END SUBROUTINE na

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE initialize(range, ranget, scales, nmodels, misfit, na_models, nd, x, nsample, ncells, restartNA)

      REAL(r32),     DIMENSION(:,:), INTENT(IN)    :: range
      REAL(r32),     DIMENSION(:,:), INTENT(OUT)   :: ranget
      REAL(r32),     DIMENSION(:),   INTENT(INOUT) :: scales
      INTEGER(i32),                  INTENT(IN)    :: nmodels
      REAL(r32),     DIMENSION(:),   INTENT(OUT)   :: misfit, na_models, x
      INTEGER(i32),                  INTENT(IN)    :: nd, nsample, ncells
      LOGICAL,                       INTENT(OUT)   :: restartNA
      INTEGER(i32)                                 :: i
      REAL(r32)                                    :: tmp

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, nd*nmodels
        na_models(i) = 0._r32
      ENDDO
      DO i = 1, nmodels
        misfit(i) = 0._r32
      ENDDO

      restartNA = .true.

      ic = 1

      ! initialize rng (requires negative seed)
      iseed = -iseed
      tmp = ran3(iseed)

      IF (scales(1) .eq. 0._r32) THEN

        DO i = 1, nd
          ranget(1,i) = range(1,i)
          ranget(2,i) = range(2,i)
          scales(i+1) = 1._r32
        ENDDO

      ELSEIF (scales(1) .eq. -1._r32) THEN

        DO i = 1, nd
          ranget(1,i) = 0._r32
          ranget(2,i) = 1._r32
          scales(i+1) = range(2,i)-range(1,i)
        ENDDO

      ELSE

        DO i = 1, nd
          ranget(1,i)  = 0._r32
          ranget(2,i)  = (range(2,i) - range(1,i)) / scales(i+1)
        ENDDO

      ENDIF

      DO i = 1, nd
        x(i) = (ranget(2,i) + ranget(1,i)) / 2._r32
      ENDDO

    END SUBROUTINE initialize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE initial_sample(na_models, nd, range, nsample)

      REAL(r32),    DIMENSION(:),   INTENT(INOUT) :: na_models
      INTEGER(i32),                 INTENT(IN)    :: nd, nsample
      REAL(r32),    DIMENSION(:,:), INTENT(IN)    :: range
      INTEGER(i32)                                :: i, j
      REAL(r32)                                   :: a, b

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, nsample
        DO j = 1, nd
          a = ran3(iseed)
          b = 1._r32 - a
          na_models((i - 1)*nd + j) = b*range(1,j) + a*range(2,j)
        ENDDO
      ENDDO

    END SUBROUTINE initial_sample

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE transform2raw(model_sca, nd, range, scales, model_raw)

      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: model_sca
      INTEGER(i32),                 INTENT(IN)  :: nd
      REAL(r32),    DIMENSION(:,:), INTENT(IN)  :: range
      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: scales
      REAL(r32),    DIMENSION(:),   INTENT(OUT) :: model_raw
      INTEGER(i32)                              :: i
      REAL(r32)                                 :: a,  b

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (scales(1) .eq. 0._r32) THEN

        DO i = 1, nd
          model_raw(i) = model_sca(i)
        ENDDO

      ELSEIF (scales(1) .eq. -1._r32) THEN

        DO i = 1,nd
          b = model_sca(i)
          a = 1-b
          model_raw(i) = a*range(1,i) + b*range(2,i)
        ENDDO

      ELSE

        DO i = 1,nd
          model_raw(i) = range(1,i) + scales(i+1)*model_sca(i)
        ENDDO

      ENDIF

    END SUBROUTINE transform2raw

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE statistics(misfit, nsample, it, ntot, mfitmin, mfitminc, mfitmean, mopt, ncells, work, ind, iwork, mfitord, ok)

      REAL(r32),    DIMENSION(:), INTENT(IN)    :: misfit
      INTEGER(i32),               INTENT(IN)    :: nsample, it, ntot
      REAL(r32),                  INTENT(INOUT) :: mfitmin
      REAL(r32),                  INTENT(OUT)   :: mfitminc, mfitmean
      INTEGER(i32),               INTENT(INOUT) :: mopt
      INTEGER(i32),               INTENT(IN)    :: ncells
      REAL(r32),    DIMENSION(:), INTENT(INOUT) :: work
      INTEGER(i32), DIMENSION(:), INTENT(INOUT) :: ind, iwork, mfitord
      INTEGER(i32),               INTENT(OUT)   :: ok
      INTEGER(i32)                              :: i, j, iopt, flow, iselect, ntotal

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      mfitminc = misfit(ntot+1)
      mfitmean = mfitminc
      iopt     = ntot+1

      DO i = ntot + 2, ntot + nsample
        mfitmean = mfitmean + misfit(i)
        IF (misfit(i) .lt. mfitminc) THEN
          mfitminc = misfit(i)
          iopt = i
        ENDIF
      ENDDO

      mfitmean = mfitmean / REAL(nsample, r32)

      IF ( (mfitminc .lt. mfitmin) .or. (it .eq. 1) ) THEN
        mopt = iopt
        mfitmin = mfitminc
      ENDIF

      IF (ncells .eq. 1) THEN

        mfitord(1) = mopt

      ELSE

        ntotal = ntot + nsample

        DO i = 1,ntotal
          ind(i)  = i
          work(i) = misfit(i)
        ENDDO

        CALL jumble(ind, work, ntotal, ok)

        IF (ok .ne. 0) RETURN

        flow = select(ncells, ntotal, work, ind, iselect)

        DO j = 1, ncells
          iwork(j) = ind(j)
        ENDDO

        CALL indexx(ncells, work, ind, ok)

        IF (ok .ne. 0) RETURN

        DO j = 1, ncells
          mfitord(j) = iwork(ind(j))
        ENDDO

      ENDIF

    END SUBROUTINE statistics

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE jumble(iarr,arr,n,ok)

      INTEGER(i32),               INTENT(IN)    :: n
      INTEGER(i32), DIMENSION(n), INTENT(INOUT) :: iarr
      REAL(r32),    DIMENSION(n), INTENT(INOUT) :: arr
      INTEGER(i32),               INTENT(OUT)   :: ok
      INTEGER(i32)                              :: j, k, ival
      REAL(r32)                                 :: rn, val

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      rn = n

      DO j = 1, n

        val = ran3(iseed)
        k = 1 + INT(val * rn)

        IF (k .eq. n+1) THEN
          k = n
        ELSEIF (k .gt. n) THEN
          ok = 4
          RETURN
        ENDIF

        ival    = iarr(j)
        iarr(j) = iarr(k)
        iarr(k) = ival
        val     = arr(j)
        arr(j)  = arr(k)
        arr(k)  = val

      ENDDO

    END SUBROUTINE jumble

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION select(k, n, arr, ind, iselect)

      INTEGER(i32),               INTENT(IN)    :: k, n
      REAL(r32),    DIMENSION(n), INTENT(INOUT) :: arr
      INTEGER(i32), DIMENSION(n), INTENT(INOUT) :: ind
      INTEGER(i32),               INTENT(OUT)   :: iselect
      INTEGER(i32)                              :: l, ir, itemp, mid, ia, i, j
      REAL(r32)                                 :: temp, a

      !-----------------------------------------------------------------------------------------------------------------------------

      l = 1
      ir = n
      1     IF (ir-l .le. 1) THEN
        IF (ir-l .eq. 1) THEN
          IF (arr(ir) .lt. arr(l)) THEN
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
            itemp=ind(l)
            ind(l)=ind(ir)
            ind(ir)=itemp
          ENDIF
        ENDIF
        select=arr(k)
        iselect=ind(k)
        RETURN
      ELSE
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        itemp=ind(mid)
        ind(mid)=ind(l+1)
        ind(l+1)=itemp
        IF (arr(l+1) .gt. arr(ir)) THEN
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          itemp=ind(l+1)
          ind(l+1)=ind(ir)
          ind(ir)=itemp
        ENDIF
        IF (arr(l) .gt. arr(ir)) THEN
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          itemp=ind(l)
          ind(l)=ind(ir)
          ind(ir)=itemp
        ENDIF
        IF (arr(l+1) .gt. arr(l)) THEN
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          itemp=ind(l+1)
          ind(l+1)=ind(l)
          ind(l)=itemp
        ENDIF
        i=l+1
        j=ir
        a=arr(l)
        ia=ind(l)
        3       CONTINUE
        i=i+1
        IF (arr(i) .lt. a) GOTO 3
        4       CONTINUE
        j=j-1
        IF (arr(j) .gt. a) GOTO 4
        IF (j .lt. i) GOTO 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        itemp=ind(i)
        ind(i)=ind(j)
        ind(j)=itemp
        GOTO 3
        5       arr(l)=arr(j)
        arr(j)=a
        ind(l)=ind(j)
        ind(j)=ia
        IF (j .ge. k) ir=j-1
        IF (j .le. k) l=i
      ENDIF
      GOTO 1

    END FUNCTION select

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE indexx(n,arr,indx, ok)

      INTEGER(i32),                             INTENT(IN)  :: n
      REAL(r32),    DIMENSION(n),               INTENT(IN)  :: arr
      INTEGER(i32), DIMENSION(n),               INTENT(OUT) :: indx
      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32)                                          :: i, indxt, ir, itemp, j, jstack, k, l
      INTEGER(i32),                   PARAMETER             :: m = 7, nstack = 50
      INTEGER(i32), DIMENSION(nstack)                       :: istack
      REAL(r32)                                             :: a

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      DO j = 1,n
        indx(j)=j
      ENDDO

      jstack=0
      l=1
      ir=n

      1     IF (ir-l .lt. m) THEN

        DO j = l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          DO i = j-1,1,-1
            IF (arr(indx(i)) .le. a) GOTO 2
            indx(i+1)=indx(i)
          ENDDO
          i=0
      2         indx(i+1)=indxt
        ENDDO

        IF (jstack .eq. 0) RETURN
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2

      ELSE

        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp

        IF (arr(indx(l+1)) .gt. arr(indx(ir))) THEN
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        ENDIF

        IF (arr(indx(l)) .gt. arr(indx(ir))) THEN
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        ENDIF

        IF (arr(indx(l+1)) .gt. arr(indx(l))) THEN
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        ENDIF

        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)

        3       CONTINUE

        i=i+1
        IF (arr(indx(i)) .lt. a) GOTO 3

        4       CONTINUE

        j=j-1
        IF (arr(indx(j)) .gt. a) GOTO 4
        IF (j .lt. i) GOTO 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        GOTO 3
        5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2

        IF (jstack .gt. nstack) THEN
          ok = 5
          RETURN
        ENDIF

        IF (ir-i+1 .ge. j-l) THEN
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        ELSE
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        ENDIF

      ENDIF
      GOTO 1

    END SUBROUTINE indexx

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE sampling(na_models, ntot, nsample, nd, nsleep, ncells, misfit, mfitord, range, xcur, restartNA, nclean, dlist)

      REAL(r32),    DIMENSION(:),   INTENT(INOUT) :: na_models
      INTEGER(i32),                 INTENT(IN)    :: ntot, nsample, nd, nsleep, ncells
      REAL(r32),    DIMENSION(:),   INTENT(IN)    :: misfit                                  !< not used
      INTEGER(i32), DIMENSION(:),   INTENT(IN)    :: mfitord
      REAL(r32),    DIMENSION(:,:), INTENT(IN)    :: range
      REAL(r32),    DIMENSION(:),   INTENT(INOUT) :: xcur
      LOGICAL,                      INTENT(INOUT) :: restartNA
      INTEGER(i32),                 INTENT(IN)    :: nclean
      REAL(r32),    DIMENSION(:),   INTENT(INOUT) :: dlist
      INTEGER(i32)                                :: id, idnext, nsample0, cell, mopt, ind_cellnext, ind_celllast, nrem, is, il, iw
      INTEGER(i32)                                :: nsampercell, nsleep0, icount, j, kd, i, ind_cell, nodex
      LOGICAL                                     :: resetlist
      REAL(r32)                                   :: dsum, dcount, t1, t2, t3, dminx, x1, x2
      SAVE                                        :: id

      !-----------------------------------------------------------------------------------------------------------------------------

      idnext = 1 + INT(ran3(iseed)*nd)

      nsample0 = nsample

      ic = ic + 1

      IF (MOD(ic,nclean) .eq.0) resetlist = .true.

      ! global timing variables
      taxis = 0._r32
      tup   = 0._r32
      tna   = 0._r32
      tcd   = 0._r32
      tdev  = 0._r32
      tres  = 0._r32

      cell = 1
      mopt = mfitord(cell)
      ind_cellnext = mopt
      ind_celllast = 0
      dsum = 0._r32
      dcount = 0._r32

      nrem = MOD(nsample, ncells)

      IF (nrem .eq. 0) THEN
        nsampercell = nsample/ncells
      ELSE
        nsampercell = 1 + nsample/ncells
      ENDIF
      nsleep0 = nsleep

      ! IF (timing) THEN
      !   t1 = cputime(t2,t3)
      !   tna = tna + t2
      ! ENDIF

      icount = 0

      DO is = 1, nsample0

        ind_cell = ind_cellnext
        icount   = icount + 1

        IF (ind_cell .ne. ind_celllast) THEN

          ! IF (timing) THEN
          !   t1 = cputime(t2,t3)
          !   tna = tna + t2
          ! ENDIF

          CALL restart(na_models, nd, ind_cell, xcur, restartNA)

          ! IF (timing) THEN
          !   t1 = cputime(t2,t3)
          !   tres = tres + t2
          !   tna = tna + t2
          ! ENDIF

        ENDIF

        IF (restartNA) then
          resetlist = .true.
          restartNA = .false.
        ENDIF

        DO il = 1,nsleep0
          DO iw = 1,nd

            IF (.not.resetlist) THEN

              ! IF (timing) THEN
              !   t1 = cputime(t2,t3)
              !   tna = tna + t2
              ! ENDIF

              CALL update_dlist(idnext, id, dlist, na_models, nd, ntot, xcur, nodex, dminx)

              ! IF (timing) THEN
              !   t1 = cputime(t2,t3)
              !   tup = tup + t2
              !   tna = tna + t2
              ! ENDIF
              nupd = nupd + 1

            ELSE

              ! IF (timing) THEN
              !   t1 = cputime(t2,t3)
              !   tna = tna + t2
              ! ENDIF

              CALL calc_dlist(idnext, dlist, na_models, nd, ntot, xcur, nodex, dminx)

              ! IF (timing) THEN
              !   t1 = cputime(t2,t3)
              !   tcd = tcd + t2
              !   tna = tna + t2
              ! ENDIF

              ncald = ncald + 1
              resetlist = .false.

            ENDIF

            id = idnext

            ! IF (timing) THEN
            !   t1 = cputime(t2,t3)
            !   tna = tna + t2
            ! ENDIF

            CALL axis_intersect(xcur, id, dlist, na_models, nd, ntot, nodex, range(1,id), range(2,id), x1, x2)

            ! IF (timing) THEN
            !   t1 = cputime(t2,t3)
            !   tna = tna + t2
            !   taxis = taxis + t2
            ! ENDIF

            kd = id + (cell-1)*nd

            xcur(id) = x1 + (x2 - x1) * ran3(iseed)

            ! IF (timing) THEN
            !   t1 = cputime(t2,t3)
            !   tdev = tdev + t2
            !   tna = tna + t2
            ! ENDIF

            idnext = idnext + 1
            IF (idnext .gt. nd) idnext=1

          ENDDO
        ENDDO

        j = ntot + is
        DO i = 1,nd
          na_models(i + (j-1)*nd) = xcur(i)
        ENDDO

        ind_celllast = ind_cell

        IF (icount .eq. nsampercell) THEN
          icount = 0
          cell = cell + 1
          ind_cellnext = mfitord(cell)
          IF (cell .eq. nrem+1) nsampercell = nsampercell - 1
        ENDIF

      ENDDO

      ! IF (timing) THEN
      !   t1 = cputime(t2,t3)
      !   tna = tna + t2
      ! ENDIF

    END SUBROUTINE sampling

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE restart(na_models, nd, mreset, x, restartNA)

      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: na_models
      INTEGER(i32),                 INTENT(IN)  :: nd, mreset
      REAL(r32),    DIMENSION(:),   INTENT(OUT) :: x
      LOGICAL,                      INTENT(OUT) :: restartNA
      INTEGER(i32)                              :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1,nd
        x(i) = na_models((mreset - 1)*nd + i)
      ENDDO

      restartNA = .true.

    END SUBROUTINE restart

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE update_dlist(dim, dimlast, dlist, bp, nd, nb, x, node, dmin)

      INTEGER(i32),                 INTENT(IN)    :: dim, dimlast
      REAL(r32),    DIMENSION(:),   INTENT(INOUT) :: dlist
      REAL(r32),    DIMENSION(:),   INTENT(IN)    :: bp
      INTEGER(i32),                 INTENT(IN)    :: nb, nd
      REAL(r32),    DIMENSION(:),   INTENT(IN)    :: x
      INTEGER(i32),                 INTENT(OUT)   :: node
      REAL(r32),                    INTENT(OUT)   :: dmin
      INTEGER(i32)                                :: i
      REAL(r32)                                   :: d1, d2, ds

      !-----------------------------------------------------------------------------------------------------------------------------

      d1       = (x(dimlast) - bp(dimlast*1))
      d1       = d1*d1
      dmin     = dlist(1)+d1
      node     = 1
      d2       = (x(dim)-bp(dim*1))
      d2       = d2*d2
      dlist(1) = dmin-d2

      DO i = 2,nb
        d1 = (x(dimlast)-bp(dimlast + (i-1)*nd))
        ds = d1
        d1 = dlist(i) + d1*d1
        IF(dmin .gt. d1) THEN
          dmin = d1
          node = i
        ENDIF
        d2 = (x(dim)-bp(dim + (i-1)*nd))
        d2 = d2*d2
        dlist(i) = d1-d2
      ENDDO

    END SUBROUTINE update_dlist

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE calc_dlist(dim, dlist, bp, nd, nb, x, nodex, dminx)

        ! dminx is not used anywhere

      INTEGER(i32),                 INTENT(IN)  :: dim
      REAL(r32),    DIMENSION(:),   INTENT(OUT) :: dlist
      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: bp
      INTEGER(i32),                 INTENT(IN)  :: nb, nd
      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: x
      INTEGER(i32),                 INTENT(OUT) :: nodex
      REAL(r32),                    INTENT(IN)  :: dminx
      INTEGER(i32)                              :: j, i
      REAL(r32)                                 :: d, dmin, dsum, dnodex

      !-----------------------------------------------------------------------------------------------------------------------------

      dmin = 0._r32
      DO j = 1,dim-1
        d = (x(j)-bp(j*1))
        d = d*d
        dmin = dmin + d
      ENDDO

      DO j = dim+1,nd
        d = (x(j)-bp(j*1))
        d = d*d
        dmin = dmin + d
      ENDDO

      dlist(1) = dmin
      d = (x(dim)-bp(dim*1))
      d = d*d
      dmin = dmin + d
      nodex = 1

      DO i = 2,nb
        dsum = 0._r32
        DO j = 1,dim-1
          d = (x(j)-bp(j + (i-1)*nd))
          d = d*d
          dsum = dsum + d
        ENDDO
        DO j = dim+1,nd
          d = (x(j)-bp(j + (i-1)*nd))
          d = d*d
          dsum = dsum + d
        ENDDO
        dlist(i) = dsum
        d = (x(dim)-bp(dim + (i-1)*nd))
        d = d*d
        dsum = dsum + d
        IF (dmin .gt. dsum) THEN
          dmin = dsum
          nodex = i
        ENDIF
        dnodex = dmin
      ENDDO

    END SUBROUTINE calc_dlist

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE axis_intersect(x, dim, dlist, bp, nd, nb, nodex, xmin, xmax, x1, x2)

      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: x     !not used
      INTEGER(i32),                 INTENT(IN)  :: dim
      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: dlist
      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: bp
      INTEGER(i32),                 INTENT(IN)  :: nd, nb, nodex
      REAL(r32),                    INTENT(IN)  :: xmin, xmax
      REAL(r32),                    INTENT(OUT) :: x1, x2
      INTEGER(i32)                              :: j
      REAL(r32)                                 :: dp0, x0, xc, dpc, dx, xi

      !-----------------------------------------------------------------------------------------------------------------------------
      x1 = xmin
      x2 = xmax

      dp0 = dlist(nodex)
      x0  = bp(dim + (nodex-1)*nd)

      DO j = 1,nodex-1
        xc  = bp(dim + (j-1)*nd)
        dpc = dlist(j)
        dx  = x0 - xc

        IF (dx .ne. 0._r32) THEN
          xi = 0.5*(x0 + xc + (dp0 - dpc)/dx)
          IF( (xi .gt. xmin) .and. (xi .lt. xmax) ) THEN
            IF ( (xi .gt. x1) .and. (x0 .gt. xc) ) THEN
              x1 = xi
            ELSEIF ( (xi .lt. x2) .and. (x0 .lt. xc) ) THEN
              x2 = xi
            ENDIF
          ENDIF
        ENDIF

      ENDDO

      DO j = nodex+1,nb
        xc  = bp(dim + (j-1)*nd)
        dpc = dlist(j)
        dx  = x0 - xc

        IF (dx .ne. 0._r32) THEN
          xi = 0.5*(x0 + xc + (dp0 - dpc)/dx)
          IF ( (xi .gt. xmin) .and. (xi .lt. xmax) ) THEN
            IF ( (xi .gt. x1) .and. (x0 .gt. xc) ) THEN
              x1 = xi
            ELSEIF ( (xi .lt. x2) .and. (x0 .lt. xc) ) THEN
              x2 = xi
            ENDIF
          ENDIF
        ENDIF

      ENDDO

    END SUBROUTINE axis_intersect

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION ran3(idum)

      INTEGER(i32),                         INTENT(INOUT) :: idum
      INTEGER(i32)                                        :: i, iff = 0, ii, inext, inextp, k, mj, mk
      INTEGER(i32),               PARAMETER               :: mbig = 1000000000, mseed = 161803398, mz = 0
      INTEGER(i32), DIMENSION(55)                         :: ma
      REAL(r32),                  PARAMETER               :: fac = 1._r32/mbig
      SAVE                                                :: iff, inext, inextp, ma

      !-----------------------------------------------------------------------------------------------------------------------------

      IF ( (idum .lt. 0) .or. (iff .eq. 0 )) THEN
        iff = 1
        mj = mseed - IABS(idum)
        mj = MOD(mj, mbig)
        ma(55) = mj
        mk = 1
        DO i = 1,54
          ii = MOD(21*i, 55)
          ma(ii) = mk
          mk = mj-mk
          IF (mk .lt. mz) mk = mk + mbig
          mj=ma(ii)
        ENDDO
        DO k = 1,4
          DO i = 1,55
            ma(i) = ma(i)-ma(1 + MOD(i+30,55))
            IF (ma(i) .lt. mz) ma(i) = ma(i) + mbig
          ENDDO
        ENDDO
        inext = 0
        inextp = 31
        idum = 1
      ENDIF

      inext = inext + 1
      IF (inext .eq. 56) inext = 1
      inextp = inextp + 1
      IF (inextp .eq. 56) inextp = 1
      mj = ma(inext) - ma(inextp)
      IF (mj .lt. mz) mj = mj + mbig
      ma(inext) = mj
      ran3 = mj*fac

    END FUNCTION ran3

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! REAL(r32) FUNCTION cputime(t1, t2)
    !
    !   REAL(r32),              INTENT(OUT) :: t1, t2
    !   REAL(r32), DIMENSION(2)             :: tarray
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   cputime = dtime(tarray)
    !   t1 = tarray(1)
    !   t2 = tarray(2)
    !
    !
    !   ! CALL system_clock(count, rate)
    !   !
    !   ! tictoc = REAL(count, r64) / REAL(rate, r64)
    !
    ! END FUNCTION cputime

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE watch_start(tictoc, comm)
    !
    !   ! Purpose:
    !   !   To start the MPI stopwatch. Timing is in double-precision. If specific communicator handle not given, mpi_comm_world is
    !   !   used.
    !   !
    !   ! Revisions:
    !   !     Date                    Description of change
    !   !     ====                    =====================
    !   !   04/05/20                  original version
    !   !
    !
    !   REAL(r64),              INTENT(OUT) :: tictoc                            !< initial time
    !   INTEGER(i32), OPTIONAL, INTENT(IN)  :: comm                              !< communicator handle
    !   INTEGER(i32)                        :: ierr
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   IF (.not.PRESENT(comm)) THEN
    !     CALL mpi_barrier(mpi_comm_world, ierr)
    !   ELSE
    !     CALL mpi_barrier(comm, ierr)
    !   ENDIF
    !
    !   tictoc = mpi_wtime()
    !
    ! END SUBROUTINE watch_start
    !
    ! ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! !===============================================================================================================================
    ! ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !
    ! SUBROUTINE watch_stop(tictoc, comm)
    !
    !   ! Purpose:
    !   !   To stop the MPI stopwatch and return elapsed time. Timing is in double-precision.  If specific communicator handle not given,
    !   !   mpi_comm_world is used.
    !   !
    !   ! Revisions:
    !   !     Date                    Description of change
    !   !     ====                    =====================
    !   !   04/05/20                  original version
    !   !
    !
    !   REAL(r64),              INTENT(INOUT) :: tictoc                          !< elapsed time
    !   INTEGER(i32), OPTIONAL, INTENT(IN)    :: comm                            !< communicator handle
    !   INTEGER(i32)                          :: ierr
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   IF (.not.PRESENT(comm)) THEN
    !     CALL mpi_barrier(mpi_comm_world, ierr)
    !   ELSE
    !     CALL mpi_barrier(comm, ierr)
    !   ENDIF
    !
    !   tictoc = mpi_wtime() - tictoc
    !
    ! END SUBROUTINE watch_stop

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION na_error(ierr) RESULT(msg)

      ! Purpose:
      !   to translate an error code into a text for all routines in current module.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),  INTENT(IN) :: ierr
      CHARACTER(100)            :: msg

      !-----------------------------------------------------------------------------------------------------------------------------

      SELECT CASE(ierr)
        CASE(1)
          msg = 'wrong dimension for argument "lrange" in na: number of dimensions in parameter space cannot be zero'

        CASE(2)
          msg = 'wrong dimension for argument "lrange" in na: number of dimensions in parameter space cannot be larger than 1024'

        CASE(3)
          msg = 'wrong value for argument "ncells" in na: cannot resample more cells than the sample size'

        CASE(4)
          msg = 'internal error in routine jumble: index "k" larger than "n"'

        CASE(5)
          msg = 'internal error in routine indexx: parameter "nstack" too small'


      END SELECT

    END FUNCTION na_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *











END MODULE m_na
