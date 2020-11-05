module vGWAS
    using DataFrames, Tables, LinearAlgebra
    using Printf, Reexport, Statistics
    using Distributions
    import LinearAlgebra: BlasReal, copytri!
    import DataFrames: DataFrame

    @reexport using WiSER
    @reexport using StatsModels
    @reexport using Distributions

    """
        ◺(n::Integer)

    Triangular number `n * (n + 1) / 2`.
    """
    @inline ◺(n::Integer) = (n * (n + 1)) >> 1
    include("scoretest.jl")
    include("scoretest_invariant.jl")
    include("pvalues.jl")
end
