      subroutine emcee_lnprob (ndim,p,lp)

        integer, intent(in) :: ndim
        double precision, intent(in), dimension(ndim) :: p
        double precision, intent(out) :: lp

        integer :: i

        lp = 1.d0
        do i=1,n
          lp = lp + p(i)*p(i)
        enddo

      end subroutine

      program main

        integer, parameter :: nwalkers=100, ndim=5
        double precision, dimension(ndim,nwalkers) :: pos
        double precision, dimension(nwalkers) :: lp

        integer :: i, j
        do j=1,nwalkers
          do i=1,ndim
            pos(i,j) = 1.d0
          enddo
        enddo

        call emcee_advance (ndim,nwalkers,pos,lp)

      end program
