using Revise
using ConstLab
using ConstLab: StrainControl, StressControl, MixedControl, Plastic, PlasticState
using Debugger
using Tensors

ncycles = 2
σ_amp = 801e6
σ_0 = 201e6 * 0
ctrl = MixedControl(
    (t) -> SymmetricTensor{2,3}((i,j) -> 0.0), # strain
    SymmetricTensor{2,3}((false, false, false, false, false, false)), # controlled strain components
    t -> SymmetricTensor{2,3}((i,j) -> i == j ? σ_0 + σ_amp*sin(t) : 0.0), # stress
    SymmetricTensor{2,3}((true, true, true, true, true, true)) # controlled stress components
)
time = range(0.0, ncycles * 2pi; length=ncycles*500)

# ## Set initial conditions
s0 = let z = zero(SymmetricTensor{2,3})
    PlasticState(z, z, 0.0, z, 0.0)
end

# Material parameters
E = 200.0e9
ν = 0.3
G = E/(2*(1+ν))
K = E/(3*(1-2*ν))
σ_y = 543e6 # 200.0e6
H = 6.5e9 # 0.1*E
κ∞ = 228e6 # 0.2*σ_y
α∞ = 132e6 # 0.1*σ_y
r = 0.5
model = Plastic(G, K, σ_y, H, κ∞, α∞, r);

main() = ConstLab.integrate(model, ctrl, time, s0; solver_params=Dict(:method => :newton));

@enter main()
