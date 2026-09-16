subroutine h_upgrade(dt, ntotal, mass, vx, rho, niac, &
                    pair_i, pair_j, dwdx, hsml)
c--------------------------------------------------------------
c Subroutine to evolve smoothing length
c
c dt     : Time step [in]
c ntotal : Number of particles [in]
c mass   : Particle masses [in]
c vx     : Velocities of all particles [in]
c rho    : Density [in]
c niac   : Number of interaction pairs [in]
c pair_i : List of first partner of interaction pair [in]
c pair_j : List of second partner of interaction pair [in]
c dwdx   : Derivative of kernel with respect to x, y and z [in]
c hsml   : Smoothing length [in/out]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer ntotal, niac, pair_i(max_interaction), &
        pair_j(max_interaction)
double precision mass(maxn), vx(dim, maxn), rho(maxn), &
                 dwdx(dim, max_interaction), hsml(maxn), dt

integer i, j, k, d
double precision dt, fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)

! Method 0: Keep smoothing length unchanged
if (sle .eq. 0) then
    return

! Method 1: h proportional to (m/rho)^(1/dim)
else if (sle .eq. 1) then
    fac = 2.0  ! Scaling factor
    do i = 1, ntotal
        hsml(i) = fac * (mass(i)/rho(i))**(1./dim)
    enddo

! Method 2: Continuity equation approach
else if (sle .eq. 2) then
    do i = 1, ntotal
        vcc(i) = 0.0
    enddo
    
    ! Calculate velocity divergence
    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        hvcc = 0.0
        
        do d = 1, dim
            dvx(d) = vx(d,j) - vx(d,i)
            hvcc = hvcc + dvx(d) * dwdx(d,k)
        enddo
        
        vcc(i) = vcc(i) + mass(j) * hvcc / rho(j)
        vcc(j) = vcc(j) + mass(i) * hvcc / rho(i)
    enddo
    
    ! Update smoothing length
    do i = 1, ntotal
        dhsml(i) = (hsml(i)/dim) * vcc(i)
        hsml(i) = hsml(i) + dt * dhsml(i)
        
        ! Prevent negative smoothing length
        if (hsml(i) .le. 0) hsml(i) = hsml(i) - dt * dhsml(i)
    enddo
endif

end