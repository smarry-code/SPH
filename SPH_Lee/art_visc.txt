subroutine art_visc(ntotal, hsml, mass, x, vx, niac, rho, c, &
                    pair_i, pair_j, w, dwdx, dvxdt, dedt)
c--------------------------------------------------------------
c Subroutine to calculate the artificial viscosity (Monaghan, 1992)
c
c ntotal : Number of particles (including virtual particles) [in]
c hsml   : Smoothing length [in]
c mass   : Particle masses [in]
c x      : Coordinates of all particles [in]
c vx     : Velocities of all particles [in]
c niac   : Number of interaction pairs [in]
c rho    : Density [in]
c c      : Sound velocity [in]
c pair_i : List of first partner of interaction pair [in]
c pair_j : List of second partner of interaction pair [in]
c w      : Kernel for all interaction pairs [in]
c dwdx   : Derivative of kernel with respect to x, y and z [in]
c dvxdt  : Acceleration with respect to x, y and z [out]
c dedt   : Change of specific internal energy [out]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer ntotal, niac, pair_i(max_interaction), &
        pair_j(max_interaction)
double precision hsml(maxn), mass(maxn), x(dim,maxn), vx(dim,maxn), &
                 rho(maxn), c(maxn), w(max_interaction), &
                 dwdx(dim,max_interaction), dvxdt(dim,maxn), dedt(maxn)

integer i, j, k, d
double precision dx, dvx(dim), alpha, beta, etq, piv, &
                 muv, vr, rr, h, mc, mrho, mhsml

! Parameters for artificial viscosity
parameter(alpha = 1.0e0)  ! Shear viscosity coefficient
parameter(beta  = 1.0e0)  ! Bulk viscosity coefficient
parameter(etq   = 0.1e0)  ! Parameter to avoid singularities

! Initialize arrays
do i = 1, ntotal
    do d = 1, dim
        dvxdt(d,i) = 0.0e0
    enddo
    dedt(i) = 0.0e0
enddo

! Main loop over interaction pairs
do k = 1, niac
    i = pair_i(k)
    j = pair_j(k)
    mhsml = (hsml(i) + hsml(j)) / 2.0  ! Average smoothing length
    
    vr = 0.0e0  ! Relative velocity projection
    rr = 0.0e0  ! Squared distance
    
    ! Calculate relative velocity and distance
    do d = 1, dim
        dvx(d) = vx(d,i) - vx(d,j)
        dx     = x(d,i) - x(d,j)
        vr     = vr + dvx(d) * dx  ! v_ij · r_ij
        rr     = rr + dx * dx      ! |r_ij|^2
    enddo
    
    ! Apply artificial viscosity only for approaching particles
    if (vr < 0.0e0) then
        ! Calculate viscosity coefficient μ_ij
        muv = mhsml * vr / (rr + mhsml*mhsml*etq*etq)
        
        ! Calculate artificial pressure Π_ij
        mc   = 0.5e0 * (c(i) + c(j))       ! Average sound speed
        mrho = 0.5e0 * (rho(i) + rho(j))   ! Average density
        piv  = (beta * muv - alpha * mc) * muv / mrho
        
        ! Calculate force contribution
        do d = 1, dim
            h = -piv * dwdx(d,k)  ! Force per unit mass
            
            ! Apply forces symmetrically
            dvxdt(d,i) = dvxdt(d,i) + mass(j) * h
            dvxdt(d,j) = dvxdt(d,j) - mass(i) * h
            
            ! Energy change from viscosity work
            dedt(i) = dedt(i) - mass(j) * dvx(d) * h
            dedt(j) = dedt(j) - mass(i) * dvx(d) * h
        enddo
    endif
enddo

! Finalize energy change
do i = 1, ntotal
    dedt(i) = 0.5e0 * dedt(i)
enddo

end