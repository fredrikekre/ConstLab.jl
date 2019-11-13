module ConstLab

using Tensors, Parameters, ForwardDiff

# Stress and strain control
abstract type ControlType end
struct StrainControl <: ControlType
    f # f(t) -> strain
end
struct StressControl <: ControlType
    f # f(t) -> stress
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

mise(s) = √(3/2) * norm(dev(s))

function constitutive_driver(model::Plastic, ε::SymmetricTensor, state::PlasticState)
    G, K, σ_y = model.G, model.K, model.σ_y

    # Assume elastic step
    Δε = ε - state.ε
    σ_tr_dev = dev(state.σ + 2G * Δε)
    @show mise(σ_tr_dev - state.α)
    Φ = mise(σ_tr_dev - state.α) - σ_y - state.κ
    # @show Φ
    if Φ < 0 # elastic
        @info "elastic"
        σ_dev = σ_tr_dev # + K * tr(ε) * one(σ_tr_dev)
        κ = state.κ
        α = state.α
        μ = state.μ
    else
        @info "plastic"
        # need to iterate
        σ_dev, κ, α, μ = solve_local_problem(model, Δε, state)
        # σ = σ_dev + K * tr(ε) * one(σ_tr_dev)
        # @show mise(σ_dev - α) - σ_y - κ
    end
    σ = σ_dev + K * tr(ε) * one(σ_dev)
    return PlasticState(ε, σ, κ, α, μ)
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

function solve_local_problem(model::Plastic, Δε::SymmetricTensor, state::PlasticState{D,N}) where {D,N}
    # extract material paremeters
    @unpack G, K, σ_y, κ∞, α∞, r, H = model

    σₙ = state.σ
    κₙ = state.κ
    αₙ = state.α
    μₙ = state.μ * 0 # ??

    # σₙ_dev = dev(σₙ - αₙ)
    σ_tr_dev = dev(σₙ) + 2G * dev(Δε)

    # pack variables
    X = zeros(14)
    # ΔX = zeros(14)
    R = zeros(14)
    pack_variables!(X, σ_tr_dev, κₙ, αₙ, μₙ)

    function residual!(R, X)
        fill!(R, 0)
        # unpack variables
        σ, κ, α, μ = extract_variables(X)
        σ_red_dev = dev(σ - α)

        # Compute residuals
        Rσ = dev(σ) - σ_tr_dev + 3G*μ/mise(σ_red_dev) * σ_red_dev
        Rκ = κ - κₙ - r * H * μ * (1 - κ/κ∞)
        Rα = dev(α) - dev(αₙ) - (1-r) * H * μ * (σ_red_dev / mise(σ_red_dev) - 1/α∞ * dev(α))
        RΦ = mise(σ_red_dev) - σ_y - κ

        # populate R
        pack_variables!(R, Rσ, Rκ, Rα, RΦ)
        return R
    end

    cfg = ForwardDiff.JacobianConfig(residual!, R, X)
    out = DiffResults.JacobianResult(X)

    iters = 0
    while true; iters += 1
        # residual!(R, X)
        ForwardDiff.jacobian!(out, residual!, R, X, cfg)
        g = DiffResults.value(out)
        J = DiffResults.jacobian(out)
        if norm(g, Inf) < 1e-6
            @info "local problem converged" iters norm(g)
            break
        elseif iters > 20
            error("did not converge")
        end
        ΔX = J \ g
        X .-= ΔX
    end
    σ, κ, α, μ = extract_variables(X)
    @info "local problem converged" dev(σ) κ dev(α) μ
    return dev(σ), κ, dev(α), μ
end



###############################################
###############################################
###############################################
###############################################

"""
    integrate(model, time, s0)

Integrate the material response.
"""
function integrate(model::Material, ctrl::ControlType, time, state::T) where {T <: MaterialState}
    states = T[]
    for (i, t) in enumerate(time)
        ε = ctrl.f(t)
        @show ε
        state′ = constitutive_driver(model, ε, state)
        push!(states, state′)
        state = state′
    end
    return states
end

end # module
