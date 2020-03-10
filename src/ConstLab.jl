module ConstLab

using Tensors, Parameters
import NLsolve, ForwardDiff, DiffResults

# Stress and strain control
abstract type ControlType end
struct StrainControl <: ControlType
    f # f(t) -> strain
end
struct StressControl <: ControlType
    # TODO: Make this just a MixedControl instead,
    # no gain in having a separate type
    f # f(t) -> stress
end
struct MixedControl <: ControlType
    # TODO: improve this masking by e.g. making masks 4th order tensors
    strain
    strain_mask::Vector{Bool}
    stress
    stress_mask::Vector{Bool}
    function MixedControl(strain, strain_mask, stress, stress_mask)
        strain_mask = collect(strain_mask.data)
        stress_mask = collect(stress_mask.data)
        for (a, b) in zip(strain_mask, stress_mask)
            a ⊻ b || error("no can do")
        end
        return new(strain, strain_mask, stress, stress_mask)
    end
end

###############################################
###############################################
###############################################
###############################################
abstract type Material end

abstract type MaterialState end

####################
## Linear Elastic ##
####################
struct LinearElastic <: Material
    G::Float64
    K::Float64
end
struct LinearElasticState{D,N} <: MaterialState
    ε::SymmetricTensor{2,D,Float64,N}
    σ::SymmetricTensor{2,D,Float64,N}
end
function LinearElasticState{D}() where D
    z = zero(SymmetricTensor{2,D})
    return LinearElasticState(z, z)
end

function constitutive_driver(model::LinearElastic, ε::SymmetricTensor, state::LinearElasticState)
    G, K = model.G, model.K
    σ = 2G * dev(ε) + 3K * vol(ε)
    return LinearElasticState(ε, σ)
end

#############################
## Plastic Prototype Model ##
#############################
struct Plastic <: Material
    G::Float64
    K::Float64
    σ_y::Float64 # Yield stress
    H::Float64
    κ∞::Float64
    α∞::Float64
    r::Float64
end
struct PlasticState{D,N} <: MaterialState
    ε::SymmetricTensor{2,D,Float64,N}
    σ::SymmetricTensor{2,D,Float64,N}
    # εₚ::SymmetricTensor{2,D,Float64,N} # ??
    κ::Float64
    α::SymmetricTensor{2,D,Float64,N}
    μ::Float64
end
# Initial state
function PlasticState{D}() where D
    z = zero(SymmetricTensor{2,D})
    return PlasticState(z, z, #=z,=# 0.0, z, 0.0)
end

function elastic_tangent(model::Plastic)
    G, K = model.G, model.K
    # Use Lamé parameters
    λ = K - 2G/3
    μ = G
    δ(i,j) = i == j ? 1.0 : 0.0
    E = SymmetricTensor{4,3}() do i, j, k, l
        λ * δ(i,j) * δ(k,l) + μ * (δ(i,k) * δ(j,l) + δ(i,l) * δ(j,k))
    end
    return E
end

mise(s) = √(3/2) * norm(dev(s))

function constitutive_driver(model::Plastic, ε::SymmetricTensor, state::PlasticState;
        solver_params=Dict{Symbol,Any}())
    G, K, σ_y = model.G, model.K, model.σ_y

    # Assume elastic step
    Δε = ε - state.ε
    σ_tr = state.σ + 2G * dev(Δε) + K * tr(Δε) * one(Δε)
    Φ = mise(σ_tr - state.α) - σ_y - state.κ
    if Φ < 0
        # elastic
        σ = σ_tr
        κ = state.κ
        α = state.α
        μ = state.μ
        dσdε = elastic_tangent(model)
    else
        # plastic; need to iterate
        σ, κ, α, μ, dσdε = solve_local_problem(model, Δε, state; solver_params=solver_params)
    end
    return σ, dσdε, PlasticState(ε, σ, κ, α, μ)
end

function pack_variables(σ, κ, α, μ::T) where T
    n = length(σ.data) + length(κ) + length(α.data) + length(μ)
    @assert n == 14 # hehe
    X = Vector{T}(undef, n)
    return pack_variables!(X, σ, κ, α, μ)
end

function pack_variables!(X, σ, κ, α, μ)
    tovoigt!(X, σ)
    X[7] = κ
    tovoigt!(X, α; offset=6+1)
    X[end] = μ
    return X
end

function extract_variables(X)
    N = 6
    σ = fromvoigt(SymmetricTensor{2,3}, X)
    κ = X[N+1]
    α = fromvoigt(SymmetricTensor{2,3}, X; offset=N+1)
    μ = X[end]
    return σ, κ, α, μ
end


function solve_local_problem(model::Plastic, Δε::SymmetricTensor, state::PlasticState{D,N}; solver_params=Dict{Symbol,Any}()) where {D,N}
    # extract material parameters
    @unpack G, K, σ_y, κ∞, α∞, r, H = model

    σₙ = state.σ
    κₙ = state.κ
    αₙ = state.α
    μₙ = state.μ * 0 # ??
    σ_tr = σₙ + 2G * dev(Δε) + K * tr(Δε) * one(Δε)

    # Initial guess
    X0 = pack_variables(σ_tr, κₙ, αₙ, μₙ)

    # Residual function
    function residual!(R, X)
        # fill!(R, 0) # pack_variables! should overwrite
        # Unpack variables
        σ, κ, α, μ = extract_variables(X)
        σ_red_dev = dev(σ - α)

        # Compute residuals
        Rσ = σ - σ_tr + 3G*μ/mise(σ_red_dev) * σ_red_dev
        Rκ = κ - κₙ - r * H * μ * (1 - κ/κ∞)
        Rα = α - αₙ - (1-r) * H * μ * (σ_red_dev / mise(σ_red_dev) - α/α∞)
        RΦ = mise(σ_red_dev) - σ_y - κ

        # Scale stresses with σ_y to make the residual dimensionless
        Rσ /= σ_y
        Rκ /= σ_y
        Rα /= σ_y
        RΦ /= σ_y

        # Populate residual vector R
        pack_variables!(R, Rσ, Rκ, Rα, RΦ)
        return R
    end

    # Use NLsolve to solve non-linear system
    params = Dict{Symbol,Any}(:ftol=>1e-14, :show_trace=>false)
    merge!(params, solver_params)
    R = similar(X0) # residual cache for constructing OnceDifferentiable (not used :( )
    cache = NLsolve.OnceDifferentiable(residual!, X0, R; autodiff=:forward)
    res = NLsolve.nlsolve(cache, X0; params...)
    NLsolve.converged(res) || error("local problem did not converge: $res")
    X = res.zero # Solution to the problem
    J = cache.DF # Extract Jacobian from last step

    # # Custom implementation of Newton's algorithm
    # X = X0 # zeros(14)
    # R = zeros(14)
    # cfg = ForwardDiff.JacobianConfig(residual!, R, X)
    # out = DiffResults.JacobianResult(X)
    # iters = 0
    # local J
    # while true; iters += 1
    #     # residual!(R, X)
    #     ForwardDiff.jacobian!(out, residual!, R, X, cfg)
    #     g = DiffResults.value(out)
    #     J = DiffResults.jacobian(out)
    #     if norm(g) < 1e-6
    #         @info "local problem converged" iters norm(g)
    #         break
    #     elseif iters > 20
    #         error("did not converge")
    #     end
    #     ΔX = J \ g
    #     X .-= ΔX
    # end

    # Extract variables from the solution X
    σ, κ, α, μ = extract_variables(X)
    # Extract ATS tensor from final Jacobian
    dRdε = elastic_tangent(model) # For fixed X
    Jinv = inv(J)
    Jinv_t = fromvoigt(SymmetricTensor{4,3}, Jinv; offset_i=0, offset_j=0)
    dσdε = (Jinv_t / σ_y) ⊡ dRdε # scale Jacobian with σ_y since residual scaled with σ_y

    return σ, κ, α, μ, dσdε
end


###############################################
###############################################
###############################################
###############################################

"""
    integrate(model, time, s0)

Integrate the material response.
"""
function integrate(model::Material, ctrl::StrainControl, time, state::T;
        solver_params=Dict{Symbol,Any}()) where {T <: MaterialState}
    states = T[]
    for (i, t) in enumerate(time)
        ε = ctrl.f(t)
        σ, dσdε, state′ = constitutive_driver(model, ε, state; solver_params=solver_params)
        push!(states, state′)
        state = state′
    end
    return states
end

function integrate(model::Material, ctrl::MixedControl, time, state::T;
        solver_params=Dict{Symbol,Any}()) where {T <: MaterialState}
    states = T[]
    for (i, t) in enumerate(time)
        # Compute controlled strain
        ε_c = ctrl.strain(t)
        # Compute controlled stress
        σ_c = ctrl.stress(t)

        # Iterate to find the correct strain
        local state′
        iters = 0
        Δε = zero(state.ε) # initial guess
        # println("----- strain iterations -----")
        while true; iters += 1
            # Current guess
            ε = SymmetricTensor{2,3}(ε_c.data .* ctrl.strain_mask .+ (state.ε.data .+ Δε.data) .* ctrl.stress_mask)

            σ′, dσdε, state′ = constitutive_driver(model, ε, state; solver_params=solver_params)

            # Compute stress residual
            R = tovoigt(σ′ - σ_c)[ctrl.stress_mask]
            # println("norm(R) = ", norm(R))
            if norm(R)/model.σ_y < 1e-6 # TOL
                # println("strain iterations converged in $iters iterations")
                break
            elseif iters > 20
                error("strain iterations did not converge")
            end
            # Take a step in strain
            J = tovoigt(dσdε)[ctrl.stress_mask, ctrl.stress_mask]
            ΔΔε = - J \ R
            ΔΔε_a = zeros(6)
            ΔΔε_a[ctrl.stress_mask] .= ΔΔε
            ΔΔε_t = fromvoigt(SymmetricTensor{2,3}, ΔΔε_a)
            Δε += ΔΔε_t
        end
        # Done, store stuff
        push!(states, state′)
        state = state′
    end
    return states
end


end # module
