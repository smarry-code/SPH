program SPH
c--------------------------------------------------------------
c This is a three dimensional SPH code. The followings are the
c basic parameters needed in this code or calculated by this code
c
c   mass   : mass of particles                    [in]
c   ntotal : total particle number used           [in]
c   dt     : Time step used in the time integration [in]
c   itype  : types of particles                   [in]
c   x      : coordinates of particles             [in/out]
c   vx     : velocities of particles              [in/out]
c   rho    : densities of particles               [in/out]
c   p      : pressure of particles                [in/out]
c   u      : internal energy of particles         [in/out]
c   hsml   : smoothing lengths of particles       [in/out]
c   c      : sound velocity of particles          [out]
c   s      : entropy of particles (not used here) [out]
c   e      : total energy of particles            [out]
c--------------------------------------------------------------

    implicit none
    include 'param.inc'
    
    integer ntotal, itype(maxn), maxtimestep, d, m, i, yesorno
    double precision x(dim,maxn), vx(dim, maxn), mass(maxn), rho(maxn), &
                    p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt
    double precision s1, s2
    
    call time_print
    call time_elapsed(s1)
    
    if (shocktube) dt = 0.005
    if (shearcavity) dt = 5.e-5
    
1   call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
    
    write(*,*)
    write(*,*) '***************************************************'
    write(*,*) 'Please input the maximal time steps '
    write(*,*) '***************************************************'
    read(*,*) maxtimestep
    
    call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, &
                          hsml, ntotal, maxtimestep, dt)
    
    call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
    
    write(*,*)
    write(*,*) '***************************************************'
    write(*,*) 'Are you going to run more time steps? (0=No, 1=yes)'
    write(*,*) '***************************************************'
    read(*,*) yesorno
    
    if (yesorno.ne.0) go to 1
    
    call time_print
    call time_elapsed(s2)
    write(*,*) 'Elapsed CPU time = ', s2-s1

end