include("../src/dlra.jl")
using Plots


# Dimension of the problem
d = 3
# Rank r to be projected on
r = 2

# --
# Utils function to get the affine plane from three (non-colinear) points 
function plane_from_points(p1, p2, p3)
    # Extract coordinates
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    x3, y3, z3 = p3

    # Construct matrix and right-hand side vector
    A = [x1 y1 1; x2 y2 1; x3 y3 1]
    b = [z1, z2, z3]

    # Solve for a, b, c
    abc = A \ b  # Equivalent to inv(A) * b
    return abc[1], abc[2], abc[3]  # Return (a, b, c)
end
# -- 


# Number of particules to use (i.e. number of trajectories)
N = 1000

# Starting point 
x0 = rand(d, N); x0 = projection_matrix_rank_r(x0, r)

# Diffusion matrix of size d \times m (m can be whatever...)
σ(x) = Diagonal(ones(d)) # Let's just take σ to be the identity

# Drift coefficient of size d (again, can be whatsoever)
b(x) = - x # Let's just take the identity function...


# Time horizon
T = 10
# Number of steps in EM scheme (and the subsequent timestep h)
n = 1000; h = T/n 

# Simulation of the N trajectories
Xs = generate_dlra_traj(x0, n, N, σ, b, h, r);

# Let's determine the r-dim affine space in which the particle are 
x1 = Xs[1, :, 1]
x2 = Xs[1, :, 2]
x3 = Xs[1, :, 3]
v, w, c = plane_from_points(x1, x2, x3)


# Anim plot 
particules = 1:30 # -- too much = long to compute
every_n = 10 # -- do not display all all time 
# -- Display within fixed boxed
x_min, x_max = minimum(Xs[:, 1, particules]), maximum(Xs[:, 1, particules])
y_min, y_max = minimum(Xs[:, 2, particules]), maximum(Xs[:, 2, particules])
z_min, z_max = minimum(Xs[:, 3, particules]), maximum(Xs[:, 3, particules])
# -- To display the two dimensional plane
xs = range(x_min, x_max, length=50)  # X range
ys = range(y_min, y_max, length=50)  # Y range
zs = [v * xi + w * yi + c for xi in xs, yi in ys'][:, 1, :]  
anim = @animate for i = 1:n
    if i % every_n == 0
        println(i)
        p = plot(xlims = (x_min, x_max), ylims = (y_min, y_max),
        zlims = (z_min, z_max), camera=(70, 15), dpi = 200)
        surface!(p, xs, ys, zs, alpha=0.1, color=:blue, legend=false)
        for j = 1:d
            plot!(p, Xs[1:i, 1, particules], Xs[1:i, 2, particules], Xs[1:i, 3, particules], label = "")
        end 
    end
end every every_n
 
gif(anim)
            





# -- Let's take b = - ∇V where V is the Coulomb potential in d=1
function b(x)
    nb_particles = length(x)
    force = zeros(nb_particles)
    for i = 1:nb_particles
       xi = x[i]; res = 0
       for j = union(1:i-1, i+1:nb_particles)
            res += -sign(xi - x[j])
       end 
       force[i] = res 
    end 
    # -- + harmonic potential 
    for i = 1:nb_particles
        force[i] += d * x[i]
    end 
    return - force 
end
