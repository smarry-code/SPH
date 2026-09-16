subroutine output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
c--------------------------------------------------------------
c Subroutine for saving particle information to external disk file
c
c x     : coordinates of particles [in]
c vx    : velocities of particles [in]
c mass  : mass of particles [in]
c rho   : densities of particles [in]
c p     : pressure of particles [in]
c u     : internal energy of particles [in]
c c     : sound velocity of particles [in]
c itype : types of particles [in]
c hsml  : smoothing lengths of particles [in]
c ntotal: total particle number [in]
c--------------------------------------------------------------

implicit none
include 'param.inc'
integer itype(maxn), ntotal
double precision x(dim, maxn), vx(dim, maxn), mass(maxn), &
                 rho(maxn), p(maxn), u(maxn), c(maxn), hsml(maxn)

integer i, d
character(len=100) :: filename

! Create timestamped filename
write(filename, "(A,I10.10,A)") "../sph_output_", itimestep, ".dat"
open(1, file=filename)

! Write header with simulation parameters
write(1, "(A)") "SPH Simulation Output"
write(1, "(A, I10)") "Timestep: ", itimestep
write(1, "(A, F12.6)") "Time: ", time
write(1, "(A, I8)") "Total particles: ", ntotal
write(1, "(A)") "-----------------------------------------"

! Write column headers
write(1, "(A)") "ID Type Mass X Y"
if (dim > 2) write(1, "(A)") " Z"
write(1, "(A)") " VX VY"
if (dim > 2) write(1, "(A)") " VZ"
write(1, "(A)") " Density Pressure Energy SoundSpeed HSML"

! Write particle data
do i = 1, ntotal
    write(1, 1001) i, itype(i), mass(i), (x(d,i), d=1,dim), &
                  (vx(d,i), d=1,dim), rho(i), p(i), u(i), c(i), hsml(i)
enddo

1001 format(1x, I8, 2x, I3, 2x, E14.8, 6(2x, E14.8), 4(2x, E14.8))

close(1)

end