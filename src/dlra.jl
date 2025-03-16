using TSVD, Optim, ProgressBars, LinearAlgebra



# -- Compute the projection onto the set of rank r matrices
# -- this is just the truncated SVD (using TSVD.jl)
function projection_matrix_rank_r(A, r)
    U, s, V = tsvd(A, r)
    return U*Diagonal(s)*transpose(V)
end

# -- Compute the projection of the matrix Z onto the tangent 
# of the rank r matrices submanifold at the point (i.e matrix) A
# --> see Koch & Lubich '07 (https://epubs.siam.org/doi/10.1137/050639703)
function projection_tangent_matrix_rank_r(Z, A, r)
    U, s, V = tsvd(A, r)
    P_U = I - U * U'
    P_V = I - V * V' 
    return Z - P_U' * Z * P_V' 
end 


# -- Generate N trajectories starting from x_0 = (x_0^1, \dots, x_0^N)
# with : n -- number of Euler--Maruyama steps 
#        σ -- the diffusion matrix function 
#        b -- the drift vector function 
#        h -- the timestep 
#        r -- the rank r
# --> This function is certainly not optimised e.g. I do the tSVD twice...

function generate_dlra_traj(x0, n, N, σ, b, h, r)
    d, m = size(σ(x0[:, 1]))
    bms = randn(n, N, m)
    res = zeros(n, d, N)
    for i = 1:N
        res[1, :, i] .= x0[:, i]
    end 
    Δb = zeros(d, N)
    Δσ = zeros(d, N)
    for k = ProgressBar(2:n)
        x = copy(res[k - 1, :, :])
        # -- computing the increment 
        for i = 1:N
            Δb[:, i] .= b(x[:, i])
            Δσ[:, i] .= σ(x[:,i])*bms[k, i, :]
        end 
        Δ = h * Δb + √h * Δσ
        # -- projecting the increment on the tangent of rank r matrices 
        x .+= projection_tangent_matrix_rank_r(Δ, x, r)
        # -- projecting the point on the set of rank r matrices 
        x .= projection_matrix_rank_r(x, r)
        res[k, :, :] .= x
    end
    return res 
end 

