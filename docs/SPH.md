# Smooth particle hydrodynamics (SPH)
Smooth particle hydrodynamics (SPH) is a branch of computational fluid dynamics, which belongs in the category of particle-based and mesh-free methods. It was first developed for astrophyscial flows related problems and has gained an increasing popularity over the last few years due to its very good applicabilty and easy extensibility in problems describing free surface and multiplhase flows in both simple and complex geometries.

In the SPH formulation, the fluid is being discretised by ficticious particles whose interpolated properties can approximate any function f(x) that is of interest over a domain $\Omega$.

$f(x) \sim \int_{\Omega} f(x')W(x-x',h)dx'$

In this  W is a kernel function and h is defined as the smoothing length that characterizes the size of the support domain of the kernel. The discrete equivalent of the above expression can be written as:

$f(x_i) = \sum_{j}^{N} \frac{m_j}{\rho_{j}} f_j W_{ij} $

The kernel function provides the weight which is assigned to each particle based on their distance from the point of interest (i.e. the point where the function f(x) is to be evaluated). The movement of the particles obeys Newton's second law of motion and the forces applied on the particles are being calculated as explained above.

More information on SPH and its applications can be found in the following resources

• http://dx.doi.org/10.1098/rspa.2019.0801
• https://doi.org/10.3389/fams.2021.797455
• https://www.ansys.com/en-gb/blog/what-is-sph-simulation
• https://www.dive-solutions.de/blog/sph-basics
• https://doi.org/10.5194/gmd-11-2691-2018
• https://arxiv.org/pdf/2104.00537.pdf

# The algorithm
In this exemplar the following algorithm which describes the solution steps of a 2D formulation of the Navier-Stokes equation is implemented:

### Calculate density 

The density of the fluid associated with each particle i is approximated as

$\rho_i = \sum_{j} m \phi d(r_{ij} ,h) $

where $r_{ij} = x_{i} −x_{j}$. Let the distance between particles, normalised by the interaction radius, be 

$q = ||rij|| h $

then the kernel density function for density, $\phi_{d}(r_{i},h)$ for q<1, is given by 

$\phi_{d}(r_{ij},h) = \frac{4}{\pi h^2{(1 −q^2)^3}}$ 

while otherwise it is set to 0.
 
### Calculate pressure

The pressure is being calculated based on the ideal gas law

$p_i = k(\rho_{i} −\rho_{0})$

where $\rho_{0}$ is a resting density and k is a gas constant.


## Calculate pressure force: 

The force exerted on the particle due to pressure from neighbouring fluid particles is calculated as

$Fp_{i} = −\sum_{j} \frac{m}{\rho_{j}} \frac{(p_i + p_j)}{2} \nabla(\phi_{p})(r_{ij} ,h)$

where if q < 1 and $i \neq j$ 

$\nabla(\phi_{p})(r_{ij} ,h) = − 30 \pi h^3 r_{ij} \frac{(1−q)^2}{q}$ 

while otherwise it is set to 0.

### Calculate viscous force

The force acting on each particle due to viscous effects is calculated as

$Fv_i = −\mu \sum_{j} m\rho_{j} v_{ij} \nabla^{2} \phi v(r_i,h)$

where $v_{ij} = v_i −v_j$. 

If $q < 1 $ and $i \neq j$

$\nabla^{2} \phi v(r_i,h) = 40 \pi h^4 (1 −q)$

while otherwise it is set to 0.

### Calculate gravity force: 

Finally, the force due to gravity is calculated as

$Fg_i = (0, −\rho_{i}g).$


### Time integration

Advancing the particles in time is performed explicitly via the use of a leap-frog scheme which provides improved stability characteristics.  Here we use the superscript in brackets to denote the time-level.

$a_i = Fp_i + Fv_i + Fg_i \rho_{i}$

$v^{t + \frac{1}{2}}_i = v^{t−\frac{1}{2}}_i + a_{i} \delta{t}$


$x^{(t+1)}_i = x^{t}_i + v^{t+ \frac{1}{2}}_i \delta{t}$

However, because the velocity is calculated at half-steps, we need to initialise the scheme on the first time step by instead doing:

$a_i = Fp_i + Fv_i + Fg_i ρ_i$

$v^{\frac{1}{2}}_i = v^{0}_i + a_i \frac{∆t}{2} x^{1}_i = x^{0}_i + v^{\frac{1}{2}}_i \delta{t}$

To ensure convergence in time a time-step of $\delta{t}$ = 10−4 is suggested.