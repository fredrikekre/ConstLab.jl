using Documenter, ConstLab

# Generate examples
include("generate.jl")

GENERATEDEXAMPLES = [joinpath("examples", f) for f in (
    "von_mises_plasticity_mixed_hardening.md",
    )]

# Build documentation.
makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "HAS_JOSH_K_SEAL_OF_APPROVAL")), # disable for local builds
    sitename = "ConstLab.jl",
    doctest = false,
    strict = false,
    pages = Any[
        "Home" => "index.md",
        "Examples" => GENERATEDEXAMPLES,
        ],
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/fredrikekre/ConstLab.jl.git",
)
