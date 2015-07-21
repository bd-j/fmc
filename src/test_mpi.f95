      subroutine emcee_lnprob (ndim, p, lp)

        ! This subroutine demonstrates the syntax for implementing a
        ! log-probability function. This particular example implements
        ! an isotropic `ndim`-dimensional Gaussian with unit variance.
        !
        ! Inputs
        ! ------
        !
        ! ndim [integer]:
        !   The dimension of the parameter space.
        !
        ! pos [double precision (ndim)]:
        !   The position in parameter space where the probability should
        !   be computed.
        !
        ! Outputs
        ! -------
        !
        ! lp [double precision]:
        !   The log-probability computed at `pos`.

        implicit none

        integer, intent(in) :: ndim
        double precision, intent(in), dimension(ndim) :: p
        double precision, intent(out) :: lp

        integer :: i

        lp = 1.d0
        do i=1,ndim
          lp = lp + p(i)*p(i)
        enddo
        lp = -0.5d0 * lp

      end subroutine

      program main

        use mpi
        implicit none

        integer, parameter :: nwalkers=20, ndim=2
        double precision, dimension(ndim,nwalkers) :: pos
        double precision, dimension(nwalkers) :: lp
        integer, dimension(nwalkers) :: accept
        integer :: i, j

        INTEGER, parameter :: masterid=0
        DOUBLE PRECISION :: one_lnp
        INTEGER :: ierr, taskid, ntasks, rqst, received_tag, status(MPI_STATUS_SIZE)
        INTEGER :: k, npos, KILL=99, BEGIN=0
        LOGICAL :: wait=.TRUE.

        ! Initialize MPI, and get the total number of processes and
        ! your process number
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )

        ! The worker's only job is to calculate the value of a function
        ! after receiving a parameter vector.
        if (taskid.ne.masterid) then
           
           ! Start event loop
           do while (wait)
              ! Get the number of parameter positions that were sent
              call MPI_RECV(npos, 1, MPI_INTEGER, &
                   masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
              ! figure out what tag it was sent with.  This call does not return until
              ! until a parameter vector is received
              received_tag = status(MPI_TAG)
              ! Check if this is the kill tag
              if ((received_tag .EQ. KILL) .OR. (npos.EQ.0)) EXIT
              ! Otherwise look for data from the master
              call MPI_RECV(pos(1,1), npos*ndim, MPI_DOUBLE_PRECISION, &
                   masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

              ! Calculate the probability for these parameter positions.
              do k=1,npos
                 call emcee_lnprob(ndim, pos(:,k), one_lnp)
                 lp(k) = one_lnp
              enddo
              
              ! Send that back to the master
              call MPI_ISEND(lp(1), npos, MPI_DOUBLE_PRECISION, &
                   masterid, BEGIN, MPI_COMM_WORLD, rqst, ierr)
           enddo
           
        endif

        ! The master process must generate intial positions and then
        ! call emcee_advance_mpi
        if (taskid.eq.masterid) then
           
           ! First seed the random number generator... don't forget this!
           call init_random_seed ()

           ! Loop over the walkers and initialize their positions. In
           ! practice, you should initialize them more intelligently than
           ! this! Also, don't forget to compute the initial
           ! log-probabilities.
           do j=1,nwalkers
              ! Loop over the number of dimensions and initialize each one
              ! in the range `(0.5, 0.5)`.
              do i=1,ndim
                 call random_number(pos(i,j))
                 pos(i,j) = pos(i,j) - 0.5d0
              enddo
           enddo
           
           ! Compute the initial log-probability, in parallel
           write(*,*) "initialization"
           call function_parallel_map (ndim, nwalkers, ntasks-1, pos, lp)

           ! Start by running a burn-in of 200 steps.
           write(*,*) "burn-in"
           do i=1,200
              ! You'll notice that I'm overwriting the position and
              ! log-probability of the ensemble at each step. This works but
              ! you also have the option of saving the samples by giving
              ! different input and output arguments.
              call emcee_advance_mpi (ndim,nwalkers,2.d0,pos,lp,pos,lp,accept,ntasks-1)
           enddo

           ! Run a production chain of 500 steps and print to `stdout`.
           write(*,*) "production"
           do i=1,500
              call emcee_advance_mpi (ndim,nwalkers,2.d0,pos,lp,pos,lp,accept,ntasks-1)
              do j=1,nwalkers
                 write(*,*) pos(:, j), lp(j)
              enddo
           enddo
           write(*,*) "free comrades"
           ! Break the workers out of their event loops so they can
           ! close
           !call MPI_Barrier( MPI_COMM_WORLD )
           call free_workers(ntasks-1)
           
        endif
        
        call MPI_FINALIZE(ierr)
           
      end program
