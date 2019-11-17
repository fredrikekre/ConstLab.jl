# # Von Mises Plasticity with Mixed Hardening
#
#
# Load packages
using Tensors
import ConstLab
using ConstLab: LinearElastic, LinearElasticState, StrainControl,
      Plastic, PlasticState

# ## Define the strain control
σ_amp = 0.005
ctrl = StrainControl((t) -> SymmetricTensor{2,3}((i,j) -> i == j == 1 ? σ_amp*sin(t) : 0.0))
time = range(0.0, 2pi; length=200)

# ## Elastic model
## model = LinearElastic(2.0, 3.0)
## s0 = LinearElasticState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}))
## res = ConstLab.integrate(model, ctrl, time, s0)

# ## Set material parameters
E = 200.0e9
ν =  0.3
G = E/(2*(1+ν))
K = E/(3*(1-2*ν))
σ_y = 200.0e6
H = 0.1*E
κ∞ = 0.2*σ_y
α∞ = 0.1*σ_y
r = 0.5
## σ_c = 0.1*E
## t_star = 100.0
## n = 1.5
model = Plastic(G, K, σ_y, H, κ∞, α∞, r)

# ## Set initial conditions
s0 = let z = zero(SymmetricTensor{2,3})
    PlasticState(z, z, 0.0, z, 0.0)
end

# ## Do integration
res = ConstLab.integrate(model, ctrl, time, s0);


# ## Plot results
using PGFPlotsX, LaTeXStrings

push!(empty!(PGFPlotsX.CUSTOM_PREAMBLE), raw"\usepackage{amsmath}")

mise(s) = √(3/2) * norm(dev(s))

# ### deviatoric 1,1
p = @pgf Axis({xlabel=L"\varepsilon_{11}", ylabel=L"(\boldsymbol{\sigma}_\mathrm{dev})_{11}"},
    PlotInc({mark="none"},
        Coordinates([r.ε[1,1] for r in res],
                    [dev(r.σ)[1,1] for r in res])
    ),
)

# ### vonMises
p = @pgf Axis({xlabel=L"\varepsilon_{11}", ylabel=L"\sigma_\mathrm{e}^\mathrm{vM}"},
    PlotInc({mark="none"},
        Coordinates([r.ε[1,1] for r in res],
                    [mise(r.σ) for r in res])
    ),
)

