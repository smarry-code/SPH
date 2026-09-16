subroutine ext_force(ntotal, mass, x, niac, pair_i, pair_j, &
                    itype, hsml, dvxdt)
c--------------------------------------------------------------
c Subroutine to calculate the external forces, e.g. gravitational forces.
c The forces from the interactions with boundary virtual particles
c are also calculated here as external forces.
c
c ntotal : Number of particles [in]
c mass   : Particle masses [in]
c x      : Coordinates of all particles [in]
c niac   : Number of interaction pairs [in]
c pair_i : List of first partner of interaction pair [in]
c pair_j : List of second partner of interaction pair [in]
c itype  : Type of particles [in]
c hsml   : Smoothing Length [in]
c dvxdt  : Acceleration with respect to x, y and z [out]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer ntotal, itype(maxn), niac, &
        pair_i(max_interaction), pair_j(max_interaction)
double precision mass(maxn), x(dim,maxn), hsml(maxn), &
                 dvxdt(dim,maxn)

integer i, j, k, d
double precision dx(dim), rr, f, rr0, dd, p1, p2

! Initialize acceleration array
do i = 1, ntotal
    do d = 1, dim
        dvxdt(d,i) = 0.0
    enddo
enddo

! Consider self-gravity (if enabled in param.inc)
if (self_gravity) then
    do i = 1, ntotal
        dvxdt(dim,i) = -9.8  ! Gravity in last dimension (y in 2D, z in 3D)
    enddo
endif

! Boundary particle force and penalty anti-penetration force
rr0 = 1.25e-5  ! Cutoff distance
dd = 1.e-2     ! Repulsion strength
p1 = 12        ! Power for repulsive force
p2 = 4         ! Power for attractive force

do k = 1, niac
    i = pair_i(k)
    j = pair_j(k)
    
    ! Only interactions between real and virtual particles
    if (itype(i) > 0 .and. itype(j) < 0) then
        rr = 0.0
        
        ! Calculate distance between particles
        do d = 1, dim
            dx(d) = x(d,i) - x(d,j)
            rr = rr + dx(d)*dx(d)
        enddo
        rr = sqrt(rr)
        
        ! Apply repulsive force if within cutoff
        if (rr < rr0) then
            f = ((rr0/rr)**p1 - (rr0/rr)**p2) / rr**2
            
            ! Apply force to real particle
            do d = 1, dim
                dvxdt(d,i) = dvxdt(d,i) + dd * dx(d) * f
            enddo
        endif
    endif
enddo

end