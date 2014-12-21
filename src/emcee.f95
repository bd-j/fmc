      subroutine emcee_advance (ndim, nwalkers, a, pin, lpin, &
                                pout, lpout, accept)

        ! This subroutine advances an ensemble of walkers using the
        ! Goodman & Weare stretch move.
        !
        ! Inputs
        ! ------
        !
        ! ndim [integer]:
        !   The dimension of the parameter space.
        !
        ! nwalkers [integer]:
        !   The number of walkers.
        !
        ! a [double precision]:
        !   The proposal scale (a tuning parameter). Using `a=2` is almost
        !   always the right move.
        !
        ! pin [double precision (ndim, nwalkers)]:
        !   The starting positions of the walkers in parameter space.
        !
        ! lpin [double precision (nwalkers)]:
        !   The value of the log-probability function at positions `pin`.
        !
        ! Outputs
        ! -------
        !
        ! pout [double precision (ndim, nwalkers)]:
        !   The final positions of the walkers in parameter space.
        !
        ! lpout [double precision (nwalkers)]:
        !   The value of the log-probability function at positions `pout`.
        !
        ! accept [integer (nwalkers)]:
        !   A binary list indicating whether or not each proposal was
        !   accepted.

        implicit none

        integer, intent(in) :: ndim, nwalkers
        double precision, intent(in) :: a
        double precision, intent(in), dimension(ndim,nwalkers) :: pin
        double precision, intent(in), dimension(nwalkers) :: lpin

        double precision, intent(out), dimension(ndim,nwalkers) :: pout
        double precision, intent(out), dimension(nwalkers) :: lpout
        integer, intent(out), dimension(nwalkers) :: accept

        integer :: k, ri
        double precision :: r, z, lp, diff
        double precision, dimension(ndim) :: q

        do k=1,nwalkers

          ! Compute a random stretch factor.
          call random_number(r)
          z = (a - 1.d0) * r + 1.d0
          z = z * z / a

          ! Select the helper walker.
          call random_number(r)
          ri = ceiling((nwalkers-1) * r)
          if (ri .ge. k) then
            ri = ri + 2
          !  q = pin(:, ri+1)
          !else
          !  q = pout(:, ri)
          endif
          q = pin(:, ri)
          ! Compute the proposal position.
          q = (1.d0 - z) * q + z * pin(:, k)

          ! Compute the new ln-probability.
          call emcee_lnprob (ndim, q, lp)
          diff = (ndim - 1.d0) * log(z) + lp - lpin(k)

          ! Accept or reject.
          if (diff .ge. 0.d0) then
            accept(k) = 1
          else
            call random_number(r)
            if (diff .ge. log(r)) then
              accept(k) = 1
            else
              accept(k) = 0
            endif
          endif

          ! Do the update.
          if (accept(k) .eq. 1) then
            pout(:, k) = q
            lpout(k) = lp
          else
            pout(:, k) = pin(:, k)
            lpout(k) = lpin(k)
          endif

        enddo

      end subroutine

      ! See: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
      subroutine init_random_seed ()

        implicit none

        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid, t(2), s
        integer(8) :: count, tms

        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
          form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
          read(un) seed
          close(un)
        else
          ! Fallback to XOR:ing the current time and pid. The PID is
          ! useful in case one launches multiple instances of the same
          ! program in parallel.
          call system_clock(count)
          if (count /= 0) then
            t = transfer(count, t)
          else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24 * 60 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
            t = transfer(tms, t)
          end if
          s = ieor(t(1), t(2))
          pid = getpid() + 1099279 ! Add a prime
          s = ieor(s, pid)
          if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
              seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
          else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
          end if
        end if
        call random_seed(put=seed)
      end subroutine init_random_seed
