subroutine av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)
c--------------------------------------------------------------
c Subroutine to calculate the average velocity to correct velocity
c for preventing penetration (Monaghan, 1992)
c
c ntotal : Number of particles [in]
c mass   : Particle masses [in]
c niac   : Number of interaction pairs [in]
c pair_i : List of first partner of interaction pair [in]
c pair_j : List of second partner of interaction pair [in]
c w      : Kernel for all interaction pairs [in]
c vx     : Velocity of each particle [in]
c rho    : Density of each particle [in]
c av     : Average velocity of each particle [out]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer ntotal, niac, pair_i(max_interaction), &
        pair_j(max_interaction)
double precision mass(maxn), w(max_interaction), vx(dim, maxn), &
                 rho(maxn), av(dim, maxn)

integer i, j, k, d
double precision vcc, dvx(dim), epsilon

! Parameter for velocity averaging
epsilon = 0.3  ! Tuning constant (0.1-0.5)

! Initialize average velocity array
do i = 1, ntotal
    do d = 1, dim
        av(d,i) = 0.0
    enddo
enddo

! Calculate SPH sum for average velocity
do k = 1, niac
    i = pair_i(k)
    j = pair_j(k)
    
    ! Compute relative velocity
    do d = 1, dim
        dvx(d) = vx(d,i) - vx(d,j)
        av(d,i) = av(d,i) - 2.0 * mass(j) * dvx(d) / (rho(i) + rho(j)) * w(k)
        av(d,j) = av(d,j) + 2.0 * mass(i) * dvx(d) / (rho(i) + rho(j)) * w(k)
    enddo
enddo

! Scale by epsilon parameter
do i = 1, ntotal
    do d = 1, dim
        av(d,i) = epsilon * av(d,i)
    enddo
enddo

end