subroutine sum_density(ntotal, hsml, mass, niac, pair_i, pair_j, w, &
                      itype, rho)
c--------------------------------------------------------------
c Subroutine to calculate the density with SPH summation algorithm
c
c ntotal   : Number of particles [in]
c hsml     : Smoothing Length [in]
c mass     : Particle masses [in]
c niac     : Number of interaction pairs [in]
c pair_i   : List of first partner of interaction pair [in]
c pair_j   : List of second partner of interaction pair [in]
c w        : Kernel for all interaction pairs [in]
c itype    : type of particles [in]
c rho      : Density [out]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer ntotal, niac, pair_i(max_interaction), &
        pair_j(max_interaction), itype(maxn)
double precision hsml(maxn), mass(maxn), w(max_interaction), rho(maxn)

integer i, j, k, d
double precision selfdens, hv(dim), r, wi(maxn)

! wi(maxn) - integration of the kernel itself
do d = 1, dim
    hv(d) = 0.e0
enddo

! Self density of each particle: Wii (Kernel for distance 0)
r = 0.0  ! Distance = 0 for self-interaction

do i = 1, ntotal    
    call kernel(r, hv, hsml(i), selfdens, hv)
    wi(i) = selfdens * mass(i) / rho(i)  ! Initial approximation
enddo


! Calculate normalization denominator
do k = 1, niac
    i = pair_i(k)
    j = pair_j(k)
    wi(i) = wi(i) + mass(j) / rho(j) * w(k)
    wi(j) = wi(j) + mass(i) / rho(i) * w(k)
enddo

! Calculate SPH sum for density
do i = 1, ntotal
    rho(i) = 0.0
enddo
! Add self-contribution
do i = 1, ntotal
    call kernel(r, hv, hsml(i), selfdens, hv)
    rho(i) = rho(i) + selfdens * mass(i)
enddo

do k = 1, niac
    i = pair_i(k)
    j = pair_j(k)
    rho(i) = rho(i) + mass(j) * w(k)
    rho(j) = rho(j) + mass(i) * w(k)
enddo

! Normalize density
if (nor_density) then
    do i = 1, ntotal
        rho(i) = rho(i) / wi(i)
    enddo
endif

end

subroutine con_density(ntotal, mass, niac, pair_i, pair_j, dwdx, vx, &
                      itype, x, rho, drhodt)
c--------------------------------------------------------------
c Subroutine to calculate the density with SPH continuity approach
c
c ntotal   : Number of particles [in]
c mass     : Particle masses [in]
c niac     : Number of interaction pairs [in]
c pair_i   : List of first partner of interaction pair [in]
c pair_j   : List of second partner of interaction pair [in]
c dwdx     : Derivative of kernel [in]
c vx       : Velocities of all particles [in]
c itype    : Type of particles [in]
c x        : Coordinates of all particles [in]
c rho      : Density [in]
c drhodt   : Density change rate [out]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer ntotal, niac, pair_i(max_interaction), &
        pair_j(max_interaction), itype(maxn)
double precision mass(maxn), dwdx(dim,max_interaction), vx(dim,maxn), &
                 x(dim,maxn), rho(maxn), drhodt(maxn)

integer i, j, k, d
double precision vcc, dvx(dim)

do i = 1, ntotal
    drhodt(i) = 0.0
enddo

do k = 1, niac
    i = pair_i(k)
    j = pair_j(k)
    
    do d = 1, dim
        dvx(d) = vx(d,i) - vx(d,j)
    enddo
    
    vcc = 0.0
    do d = 1, dim
        vcc = vcc + dvx(d) * dwdx(d,k)
    enddo
    
    drhodt(i) = drhodt(i) + mass(j) * vcc
    drhodt(j) = drhodt(j) + mass(i) * vcc
enddo

! For axisymmetric case (1D with symmetry)
!if (dim .eq. 1 .and. nsym .eq. 1) then
!    do i = 1, ntotal
!         if (x(1,i) .gt. 1e-8) then
!            drhodt(i) = drhodt(i) - rho(i) * vx(1,i) / x(1,i)
!        endif
!    enddo
!endif

end