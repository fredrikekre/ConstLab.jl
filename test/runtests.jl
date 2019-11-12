

using Revise, Tensors
import ConstLab
using ConstLab: LinearElastic, LinearElasticState, StrainControl,
      Plastic, PlasticState


ctrl = StrainControl((t) -> SymmetricTensor{2,3}((i,j) -> i == j == 1 ? t : 0.0))
time = range(0.0, 2.0; length=50)

# # LinearElastic
# model = LinearElastic(2.0, 3.0)
# s0 = LinearElasticState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}))

# res = ConstLab.integrate(model, ctrl, time, s0)

# Plastic
model = Plastic(2.0, 3.0, 4.0, 5.0, 10, 10, 0.5)
s0 = let z = zero(SymmetricTensor{2,3})
    PlasticState(z, z, 0.0, z, 0.0)
end

res = ConstLab.integrate(model, ctrl, time, s0)



using PGFPlotsX, LaTeXStrings

push!(empty!(PGFPlotsX.CUSTOM_PREAMBLE), raw"\usepackage{amsmath}")

mise(s) = √(3/2) * norm(s)

p = @pgf Axis({xlabel=L"\varepsilon_{11}", ylabel=L"\sigma_\mathrm{e}^\mathrm{vM}"},
    PlotInc({},
        Coordinates([r.ε[1,1] for r in res],
                    [mise(r.σ) for r in res])
    ),
)



p = @pgf Axis({xlabel=L"\varepsilon_{11}"},
    PlotInc({},
        Coordinates([r.ε[1,1] for r in res],
                    [r.σ[1,1] for r in res])
    ),
    LegendEntry(L"\sigma_{11}"),
    PlotInc({},
        Coordinates([r.ε[1,1] for r in res],
                    [r.σ[2,2] for r in res])
    ),
    LegendEntry(L"\sigma_{22}"),
    PlotInc({},
        Coordinates([r.ε[1,1] for r in res],
                    [r.σ[1,2] for r in res])
    ),
    LegendEntry(L"\sigma_{12}"),
)

