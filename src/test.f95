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

        implicit none

        integer, parameter :: nwalkers=20, ndim=2
        double precision, dimension(ndim,nwalkers) :: pos
        double precision, dimension(nwalkers) :: lp
        integer, dimension(nwalkers) :: accept

        integer :: i, j

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

          ! Compute the initial log-probability.
          call emcee_lnprob(ndim, pos(:, j), lp(j))
        enddo

        ! Start by running a burn-in of 200 steps.
        do i=1,200
          ! You'll notice that I'm overwriting the position and
          ! log-probability of the ensemble at each step. This works but
          ! you also have the option of saving the samples by giving
          ! different input and output arguments.
          call emcee_advance (ndim,nwalkers,2.d0,pos,lp,pos,lp,accept)
        enddo

        ! Run a production chain of 500 steps and print to `stdout`.
        do i=1,500
          call emcee_advance (ndim,nwalkers,2.d0,pos,lp,pos,lp,accept)
          do j=1,nwalkers
            write(*,*) pos(:, j), lp(j)
          enddo
        enddo

      end program
