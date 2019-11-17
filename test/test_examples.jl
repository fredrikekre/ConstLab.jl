module VonMisesMixedPlasticity
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/von_mises_plasticity_mixed_hardening.jl"))
        end
    end
end

