include("../src/dlra.jl")
using Plots


# Dimension of the problem
d = 3
# Rank r to be projected on
r = 2

# Number of particules to use (i.e. number of trajectories)
N = 2500

# Starting point 
x0 = rand(d, N); x0 = projection_matrix_rank_r(x0, r)

# Diffusion matrix of size d \times m (m can be whatever...)
σ(x) = Diagonal(ones(d)) # Let's just take σ to be the identity

# Drift coefficient of size d (again, can be whatsoever)
b(x) = - x # Let's just take the identity function...


# Time horizon
T = 2
# Number of steps in EM scheme (and the subsequent timestep h)
n = 1000; h = T/n 

# Simulation of the N trajectories
Xs = generate_dlra_traj(x0, n, N, σ, b, h, r);

# Anim plot 
gr()
particules = 1:10 # -- too much = long to compute
every_n = 10 # -- do not display all all time 
# -- Display within fixed boxed
x_min, x_max = minimum(Xs[:, 1, particules]), maximum(Xs[:, 1, particules])
y_min, y_max = minimum(Xs[:, 2, particules]), maximum(Xs[:, 2, particules])
z_min, z_max = minimum(Xs[:, 3, particules]), maximum(Xs[:, 3, particules])
anim = @animate for i = 1:n
    if i % every_n == 0
        println(i)
        p = plot(xlims = (x_min, x_max), ylims = (y_min, y_max),
        zlims = (z_min, z_max), camera= (-5, 15), dpi = 250)
        #surface!(p, xs, ys, zs, alpha=0.1, color=:blue, legend=false, colorbar = false)
        for j = 1:d
            plot!(p, Xs[1:i, 1, particules], Xs[1:i, 2, particules], Xs[1:i, 3, particules], label = "")
        end 
    end
end every every_n
 
gif(anim, "DLRA/media/Dimension3Rank2.gif")
            

