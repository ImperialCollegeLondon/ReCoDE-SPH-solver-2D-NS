# Smooth Particle Hydrodynamics (SPH)
Smooth particle hydrodynamics (SPH) is a branch of computational fluid dynamics, which belongs in the category of particle-based and mesh-free methods. It was first developed for astrophysical flows related problems and has gained an increasing popularity over the last few years due to its very good applicability and easy extensibility in problems describing free surface and multiphase flows in both simple and complex geometries.

In the SPH formulation, the fluid is being discretised by fictitious particles whose interpolated properties can approximate any function f(x) that is of interest over a domain $\Omega$.


$$f(x) \sim \int_{\Omega} f(x')W(x-x',h)dx'$$


In this  W is a kernel function and h is defined as the smoothing length that characterizes the size of the support domain of the kernel. The discrete equivalent of the above expression for the $i^{th}$ particle can be written as:

$$f_i = \sum_{j}^{N} \frac{m_j}{\rho_{j}} f_j W_{ij} $$

The kernel function provides the weight which is assigned to each particle based on their distance from the point of interest (i.e. the point where the function f(x) is to be evaluated). The movement of the particles obeys Newton's second law of motion and the forces applied on the particles are being calculated as explained above.

More information on SPH and its applications can be found in the following resources

- http://dx.doi.org/10.1098/rspa.2019.0801
- https://doi.org/10.3389/fams.2021.797455
- https://www.ansys.com/en-gb/blog/what-is-sph-simulation
- https://www.dive-solutions.de/blog/sph-basics
- https://doi.org/10.5194/gmd-11-2691-2018
- https://arxiv.org/pdf/2104.00537.pdf

# The algorithm
In this exemplar the following algorithm which describes the solution steps of a 2D formulation of the Navier-Stokes equation is implemented:

## Density 

The density of the fluid associated with each particle i is approximated as

$$\rho_i = \sum_{j} m \phi_d(r_{ij} ,h)  $$

where $r_{ij} = x_{i} −x_{j}$ and $m$ is the mass of a particle and the kernel density function for density, $\phi_{d}(r_{ij},h)$ is given by 

$$\phi_{d}(r_{ij},h) = \begin{cases}
\frac{4}{\pi h^2{(1 −q_{ij}^2)^3}} & \text{if $q_{ij} < 1$}\\
0 & \text{otherwise} 
\end{cases}$$

where $q_{ij}$ is the distance between particle $i$ and particle $j$, normalised by the interaction radius, given by

$$q_{ij} = \frac{||r_{ij}||} {h} $$
 
## Pressure

The pressure is calculated based on the ideal gas law

$$p_i = k(\rho_{i} −\rho_{0})$$

where $\rho_{0}$ is a resting density and $k$ is a gas constant.


## Pressure force

The force exerted on the particle due to pressure from neighboring fluid particles is calculated as

$$F_{pi} = −\sum_{j} \frac{m}{\rho_{j}} \frac{(p_i + p_j)}{2} \nabla(\phi_{p})(r_{ij} ,h)$$ 

where $\phi_d$ is the the kernel density function for pressure :

$$\nabla(\phi_{p})(r_{ij} ,h) = \begin{cases}
− 30 \pi h^3 r_{ij} \frac{(1−q)^2}{q} & \text{for } q_{ij} < 1 \text{ and } i \neq j\\
0 & \text{otherwise}
\end{cases}
$$

while otherwise it is set to 0.

## Viscous force

The force acting on each particle due to viscous effects is calculated as

$$F_{vi} = −\mu \sum_{j} m\rho_{j} v_{ij} \nabla^{2} \phi v(r_i,h)$$

where $v_{ij} = v_i −v_j$, $v_{i}$ is the velocity of particle i and $\mu$ is the dynamic viscosity.

If $q < 1 $ and $i \neq j$

$$\nabla^{2} \phi v(r_i,h) = 40 \pi h^4 (1 −q)$$

while otherwise it is set to 0.

## Gravity force: 

Finally, the force due to gravity is calculated as

$$F_{gi} = (0, −\rho_{i}g)$$

## Acceleration

The acceleration of each particle is calculated as:

$$a_i =\frac{F_{pi} + F_{vi} + F_{gi}} {\rho_{i}}$$

## Time integration

We solve the equation as a function of time by finding the velocity and position of each particle at each of a number of time steps. We denote a property $x$ of particle $i$ at time step $t$ as $x^{t}_i$. The state of the property half way between time steps $t$ and $t + 1$ is denoted as $x^{t + \frac{1}{2}}_i$.

We begin with the initial conditions of the system, which are the positions and velocities of the particles at time $t = 0$. We iteratively use the state of the system at time step $t$ to find the state of the system at time step $t + 1$ using a leap-frog scheme which provides improved stability characteristics.

$$v^{t + \frac{1}{2}}_i = v^{t−\frac{1}{2}}_i + a_{i}^t \delta{t}$$

$$x^{(t+1)}_i = x^{t}_i + v^{t+ \frac{1}{2}}_i \delta{t}$$

However, because the velocity is calculated at half-steps, we need to initialise the scheme on the first time step using:

$$v^{\frac{1}{2}}_i = v^{0}_i + a_i^0 \frac{∆t}{2}$$

where $∆t$ is the time step size. To ensure convergence, a small time-step is required. A value of $∆t = 10^{−4}$s is suggested.


## Initial Condition

To initialise the simulation, one or more particles must be specified. These should be evenly distributed. 

Once the particles are placed, the particle densities should be evaluated with an assumed mass of m = 1 and the mass subsequently scaled so that the density is equal to the reference density:

$$ m = \frac{N \rho_{0}}{\sum_{i} \rho_{i}} $$

##  Boundary Conditions
All the boundaries are solid walls with damped reflecting boundary conditions. If particles are within a distance h of the boundary, their velocity and position should be updated to reverse the motion of the particle and keep it within the domain. For example, on the right boundary, the x-component of position and velocity would be modified as

$$ x  = 1 - h $$
$$ u  = - eu $$ 