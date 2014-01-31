      subroutine emcee_advance (ndim,nwalkers,pars,lnprob)

        integer, intent(in) :: ndim, nwalkers
        double precision, intent(inout), dimension(ndim,nwalkers) :: pars
        double precision, intent(inout), dimension(nwalkers) :: lnprob

        double precision :: lp

        call emcee_lnprob (ndim,pars(:,1),lp)
        write(*,*) lp

      end subroutine
