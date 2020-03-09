using Documenter, ConstLab
ENV["JULIA_DEBUG"] = "Documenter,LiveServer"
# Load remaining dependencies to avoid leaking
# precompiling ... into the docs.
import PGFPlotsX, LaTeXStrings, Literate, Tensors

if haskey(ENV, "GITHUB_ACTIONS")
    PGFPlotsX.latexengine!(PGFPlotsX.PDFLATEX)
    ENV["GITHUB_EVENT_PATH"] = joinpath(@__DIR__, "actions-event.json")
end

# Generate examples
include("generate.jl")

GENERATEDEXAMPLES = [joinpath("examples", f) for f in (
    "von_mises_plasticity_mixed_hardening.md",
    )]

# Build documentation.
makedocs(
    format = Documenter.HTML(
        prettyurls = haskey(ENV, "GITHUB_ACTIONS"), # disable for local builds
        assets = ["assets/custom.css"], # white background for images in documenter-dark
    ),
    sitename = "ConstLab.jl",
    doctest = false,
    strict = false,
    clean = !isdefined(Main, :LiveServer),
    pages = Any[
        "Home" => "index.md",
        "Examples" => GENERATEDEXAMPLES,
        ],
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/fredrikekre/ConstLab.jl.git",
    push_preview = true,
)
