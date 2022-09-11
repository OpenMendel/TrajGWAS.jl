module TrajGWAS
    using DataFrames, Tables, LinearAlgebra
    using Printf, Reexport, Statistics, SnpArrays, VCFTools, VariantCallFormat
    using Distributions, CSV, BGEN
    import LinearAlgebra: BlasReal, copytri!
    import DataFrames: DataFrame

    @reexport using WiSER
    @reexport using StatsModels
    @reexport using Distributions
    export WSVarScoreTest, WSVarScoreTestInvariant, test!
    export trajgwas
    export matchindices
    export Adjustor
    using MathOptInterface
    const MOI = MathOptInterface

    """
        ◺(n::Integer)

    Triangular number `n * (n + 1) / 2`.
    """
    @inline ◺(n::Integer) = (n * (n + 1)) >> 1
    include("scoretest.jl")
    include("scoretest_invariant.jl")
    include("pvalues.jl")
    include("gwas.jl")
    include("spa.jl")
end
