        subroutine emcee_advance_mpi (ndim, nwalkers, a, pin, lpin, &
                                      pout, lpout, accept, nworkers)

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
        ! nworkers [integer]:
        !   The number of available workers, not including the master
        !   process, which is assumed to be workerid=0
        
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

        include 'mpif.h'
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

        INTEGER, intent(in) :: nworkers
        DOUBLE PRECISION, dimension(nwalkers) :: zarr
        INTEGER :: workerid
        
        ! Loop over the walkers to propose new positions and send them
        ! to workers
        do k=1,nwalkers
           ! Which worker does this go to?
           workerid = mod(k,nworkers) + 1
           
           ! Compute a random stretch factor and store it
           call random_number(r)
           z = (a - 1.d0) * r + 1.d0
           z = z * z / a
           zarr(k) = z
           
           ! Select the helper walker.
           call random_number(r)
           ri = ceiling((nwalkers-1) * r)
           if (ri .ge. k) then
              ri = ri + 1
           endif
           q = pin(:, ri)
           
           ! Compute the proposal position.
           q = (1.d0 - z) * q + z * pin(:, k)
           ! Dispatch proposal to a worker to figure out lnp
           MPI_ISEND(q, ndim, MPI_DOUBLE_PRECISION, &
                workerid, k, MPI_COMM_WORLD, ierr)   
        enddo

        ! Loop over the walkers to get the proposal lnp,
        ! accept/reject, and update
        do k=1,nwalkers
           ! Which worker had this proposal?
           workerid = mod(k,nworkers) + 1
           ! Get the answer from that worker
           MPI_RECV(lp, 1, MPI_DOUBLE_PRECISION, &
                workerid, k, MPI_COMM_WORLD, status, ierr)
           diff = (ndim - 1.d0) * log(zarr(k)) + lp - lpin(k)

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


      subroutine function_parallel_map(ndim, nk, nworkers, pos, lnp)
        !
        ! This subroutine sends rows of the pos array to whatever
        ! function is set up to receive them in a different process,
        ! and collects the results.  This could probably be done with
        ! scatter/gather, but I don't know how to use those yet.
        !
        include "mpif.h"
        implicit none

        integer, intent(in) :: ndim, nk, nworkers
        double precision, intent(in), dimension(ndim,nk) :: pos
        double precision, intent(out), dimension(nk) :: lnp

        integer :: k, workerid
        double precision :: one_lnp
        
        ! Send parameter positions to processes
        do k=1,nk
           workerid = mod(k,nworkers) + 1
           MPI_ISEND(pos(:,k), ndim, MPI_DOUBLE_PRECISION, &
                workerid, k, MPI_COMM_WORLD, ierr)   
        enddo

        ! Collect results
        do k=1,nk
           workerid = mod(k,nworkers) + 1
           MPI_RECV(one_lnp, 1, MPI_DOUBLE_PRECISION, &
                workerid, k, MPI_COMM_WORLD, status, ierr)
           lnp(k) = one_lnp
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
