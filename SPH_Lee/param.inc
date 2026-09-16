c--------------------------------------------------------------
c    Including file for parameters and constants used
c    in the entire SPH software packages.
c--------------------------------------------------------------

c    dim : Dimension of the problem (1, 2 or 3)
    integer dim
    parameter (dim = 1)
    
c    maxn : Maximum number of particles
c    max_interaction : Maximum number of interaction pairs
    integer maxn, max_interaction
    parameter (maxn = 12000, max_interaction = 100 * maxn)
    
c    Parameters for the computational geometry
    double precision x_maxgeom, x_mingeom, y_maxgeom, 
    &                y_mingeom, z_maxgeom, z_mingeom
    parameter (x_maxgeom = 10.e0, x_mingeom = -10.e0,
    &          y_maxgeom = 10.e0, y_mingeom = -10.e0,
    &          z_maxgeom = 10.e0, z_mingeom = -10.e0)

c    SPH algorithm for particle approximation (pa_sph)
c      1: (e.g. (p(i)+p(j))/rho(i)/rho(j))
c      2: (e.g. p(i)/rho(i)**2 + p(j)/rho(j)**2)
    integer pa_sph
    parameter(pa_sph = 2)

c    Nearest neighboring particle searching (nnps) method
c      1: Simplest and direct searching
c      2: Sorting grid linked list
c      3: Tree algorithm
    integer nnps
    parameter(nnps = 1)

c    Smoothing length evolution (sle) algorithm
c      0: Keep unchanged
c      1: h = fac * (m/rho)**(1/dim)
c      2: dh/dt = (-1/dim)*(h/rho)*(drho/dt)
c      3: Other approaches
    integer sle
    parameter(sle = 0)

c    Smoothing kernel function
c      1: cubic spline kernel by W4-Spline (Monaghan 1985)
c      2: Gauss kernel (Gingold and Monaghan 1981)
c      3: Quintic kernel (Morris 1997)
    integer skf
    parameter(skf = 1)

c    Switches for different scenarios [logical flags]
    logical summation_density, average_velocity, config_input,
    &       virtual_part, vp_input, visc, ex_force, heat_artificial,
    &       visc_artificial, self_gravity, nor_density
    
c    Density approximation
    parameter(summation_density = .true.)   ! .TRUE.: density summation
    
c    Velocity treatment
    parameter(average_velocity = .false.)   ! .TRUE.: Monaghan average velocity
    
c    Initial configuration
    parameter(config_input = .false.)       ! .TRUE.: load from file
    
c    Virtual particles
    parameter(virtual_part = .false.)       ! .TRUE.: use virtual particles
    parameter(vp_input = .false.)           ! .TRUE.: load virtual particles
    
c    Physical effects
    parameter(visc = .false.)               ! .TRUE.: physical viscosity
    parameter(ex_force = .false.)           ! .TRUE.: external forces
    parameter(visc_artificial = .true.)     ! .TRUE.: artificial viscosity
    parameter(heat_artificial = .false.)    ! .TRUE.: artificial heat
    parameter(self_gravity = .false.)       ! .TRUE.: self-gravity
    
c    Density normalization
    parameter(nor_density = .false.)        ! .TRUE.: density normalization

c    Symmetry of the problem
c      0: no symmetry
c      1: axis symmetry
c      2: center symmetry
    integer nsym
    parameter(nsym = 0)

c    Control parameters for output
    logical int_stat
    parameter(int_stat = .true.)            ! Print interaction statistics
    
    integer print_step, save_step, moni_particle
    parameter(print_step = 100,             ! Timesteps for screen output
    &         save_step = 500,              ! Timesteps for disk output
    &         moni_particle = 1600)         ! Monitored particle number

c    Mathematical constants
    double precision pi
    parameter(pi = 3.14159265358979323846)

c    Simulation cases [logical flags]
    logical shocktube, shearcavity
    parameter(shocktube = .true.)           ! 1D shock tube simulation
    parameter(shearcavity = .false.)        ! 2D shear cavity simulation