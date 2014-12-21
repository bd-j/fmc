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
        DOUBLE PRECISION, dimension(ndim) :: one_pos
        DOUBLE PRECISION :: one_lnp
        INTEGER :: ierr, taskid, ntasks, received_tag, rqst, status(MPI_STATUS_SIZE)
        LOGICAL :: wait=.TRUE.

        ! Initialize MPI, and get the total number of processes and
        ! your process number
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )

        ! The slave's only job is calculate the value of a function
        ! after receiving a parameter vector.
        if (taskid.ne.masterid) then
           !Event loop
           do while (wait)
              ! Get the parameter value from the master, and figure out
              ! what tag it was sent with
              call MPI_RECV(one_pos, ndim, MPI_DOUBLE_PRECISION, &
                   masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
              received_tag = status(MPI_TAG)
           
              ! Calculate the probability for this parameter position.
              call emcee_lnprob(ndim, one_pos, one_lnp)
              !write(*,*) taskid, received_tag, one_pos, one_lnp
              
              ! Send that back to the master, with the correct tag
              call MPI_ISEND(one_lnp, 1, MPI_DOUBLE_PRECISION, &
                   masterid, received_tag, MPI_COMM_WORLD, rqst, ierr)
           enddo
        endif

        ! The master process must get intial positions and then
        ! call emcee_advance
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
           call function_parallel_map(ndim, nwalkers, ntasks-1, pos, lp)
           !do j=1,nwalkers
           !   write(*,*) j, pos(:, j), lp(j)
           !enddo

           write(*,*) '#initial lnp calculated'
           ! Start by running a burn-in of 200 steps.
           do i=1,200
              !write(*,*) 'iteration=', i
              ! You'll notice that I'm overwriting the position and
              ! log-probability of the ensemble at each step. This works but
              ! you also have the option of saving the samples by giving
              ! different input and output arguments.
              call emcee_advance_mpi (ndim,nwalkers,2.d0,pos,lp,pos,lp,accept,ntasks-1)
              !do j=1,nwalkers
              !   write(*,*) j, pos(:, j), lp(j)
              !enddo

           enddo

           ! Run a production chain of 500 steps and print to `stdout`.
           write(*,*) '# Production'
           do i=1,500
              call emcee_advance_mpi (ndim,nwalkers,2.d0,pos,lp,pos,lp,accept,ntasks-1)
              do j=1,nwalkers
                 write(*,*) pos(:, j), lp(j)
              enddo
           enddo
        endif
        
        call MPI_FINALIZE(ierr)
           
      end program
