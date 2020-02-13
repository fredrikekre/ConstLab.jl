using Revise
using ConstLab
using ConstLab: StrainControl, StressControl, Plastic, PlasticState
using Debugger
using Tensors

# Values from Kim
ε_amp = 0.05
# ctrl = StrainControl((t) -> SymmetricTensor{2,3}((i,j) -> i == j == 1 ? ε_amp*t : 0.0))
time = range(0.0, 1.0; length=101)

# Set material parameters
E = 200e9
ν = 0.3
G = E/(2*(1+ν))
K = E/(3*(1-2*ν))
σ_y = 543e6 # 200.0e6
H = 6.5e9 # 0.1*E
κ∞ = 132e6 # 0.2*σ_y
α∞ = 22.8e6 # 0.1*σ_y
r = 0.5


model = Plastic(G, K, σ_y, H, κ∞, α∞, r)
# ## Set initial conditions
s0 = let z = zero(SymmetricTensor{2,3})
    PlasticState(z, z, 0.0, z, 0.0)
end

# main() = ConstLab.integrate(model, ctrl, time, s0; solver_params=Dict(:method => :newton, :ftol=>1e-4));

# @run main()

# Mixed


ctrl = ConstLab.MixedControl(
    t -> SymmetricTensor{2,3}((i,j) -> i == j == 1 ? ε_amp*t : 0.0),
    SymmetricTensor{2,3}((true, false, false, false, false, false)),
    t -> zero(SymmetricTensor{2,3}),
    SymmetricTensor{2,3}((false, true, true, true, true, true))
)

main() = ConstLab.integrate(model, ctrl, time, s0; solver_params=Dict(:method => :newton, :ftol=>1e-4));

@enter main()
