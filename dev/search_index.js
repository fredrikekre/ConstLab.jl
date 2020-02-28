var documenterSearchIndex = {"docs":
[{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"EditURL = \"https://github.com/fredrikekre/ConstLab.jl/blob/master/docs/src/literate/von_mises_plasticity_mixed_hardening.jl\"","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Von-Mises-Plasticity-with-Mixed-Hardening-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"Load packages","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"using Tensors\nimport ConstLab\nusing ConstLab: LinearElastic, LinearElasticState, StrainControl,\n      Plastic, PlasticState","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Define-the-strain-control-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Define the strain control","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"σ_amp = 0.005\nctrl = StrainControl((t) -> SymmetricTensor{2,3}((i,j) -> i == j == 1 ? σ_amp*sin(t) : 0.0))\ntime = range(0.0, 2pi; length=200)","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Elastic-model-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Elastic model","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"# model = LinearElastic(2.0, 3.0)\n# s0 = LinearElasticState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}))\n# res = ConstLab.integrate(model, ctrl, time, s0)","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Set-material-parameters-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Set material parameters","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"E = 200.0e9\nν =  0.3\nG = E/(2*(1+ν))\nK = E/(3*(1-2*ν))\nσ_y = 200.0e6\nH = 0.1*E\nκ∞ = 0.2*σ_y\nα∞ = 0.1*σ_y\nr = 0.5\n# σ_c = 0.1*E\n# t_star = 100.0\n# n = 1.5\nmodel = Plastic(G, K, σ_y, H, κ∞, α∞, r)","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Set-initial-conditions-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Set initial conditions","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"s0 = let z = zero(SymmetricTensor{2,3})\n    PlasticState(z, z, 0.0, z, 0.0)\nend","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Do-integration-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Do integration","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"res = ConstLab.integrate(model, ctrl, time, s0; solver_params=Dict(:method => :newton));\nnothing #hide","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#Plot-results-1","page":"Von Mises Plasticity with Mixed Hardening","title":"Plot results","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"using PGFPlotsX, LaTeXStrings\n\npush!(empty!(PGFPlotsX.CUSTOM_PREAMBLE), raw\"\\usepackage{amsmath}\")\n\nmise(s) = √(3/2) * norm(dev(s))","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#deviatoric-1,1-1","page":"Von Mises Plasticity with Mixed Hardening","title":"deviatoric 1,1","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"p = @pgf Axis({xlabel=L\"\\varepsilon_{11}\", ylabel=L\"(\\boldsymbol{\\sigma}_\\mathrm{dev})_{11}\"},\n    PlotInc({mark=\"none\"},\n        Coordinates([r.ε[1,1] for r in res],\n                    [dev(r.σ)[1,1] for r in res])\n    ),\n)","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#vonMises-1","page":"Von Mises Plasticity with Mixed Hardening","title":"vonMises","text":"","category":"section"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"p = @pgf Axis({xlabel=L\"\\varepsilon_{11}\", ylabel=L\"\\sigma_\\mathrm{e}^\\mathrm{vM}\"},\n    PlotInc({mark=\"none\"},\n        Coordinates([r.ε[1,1] for r in res],\n                    [mise(r.σ) for r in res])\n    ),\n)","category":"page"},{"location":"examples/von_mises_plasticity_mixed_hardening/#","page":"Von Mises Plasticity with Mixed Hardening","title":"Von Mises Plasticity with Mixed Hardening","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"DocTestSetup = :(using JuAFEM)","category":"page"},{"location":"#ConstLab.jl-1","page":"Home","title":"ConstLab.jl","text":"","category":"section"}]
}