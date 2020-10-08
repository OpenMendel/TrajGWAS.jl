module vGWAS

using DataFrames, Tables, LinearAlgebra
using Printf, Reexport, Statistics
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

"""
    WSVarScoreTestBaseObs

A base per-observation object for the score test of within-subject variance linear mixed model data instance without information on X1 or W1. 
Contains base variables for testing
H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0, 
for the full model of WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ), 
We make use of the fitted null model.
"""
struct WSVarScoreTestBaseObs{T <: BlasReal}
    nullObs       :: WSVarLmmObs{T} # is it redundant?
    n :: Int
    p :: Int
    q :: Int
    l :: Int

    # working arrays
    A_21_β2β1_pre   :: Matrix{T} # becomes A_21_β2β1 when right-multiplied by X1
    A_21_τ2τ1_pre   :: Matrix{T} # becomes A_21_τ2τ1 when right-multiplied by W1
    A_21_Lγτ1_pre   :: Matrix{T} # becomes A_21_Lγτ1 when right-multiplied by W1
end

function WSVarScoreTestBaseObs(nullObs::WSVarLmmObs{T}) where T <: BlasReal
    n, p, q, l = size(nullObs.Xt, 2), size(nullObs.Xt, 1), size(nullObs.Zt, 1), size(nullObs.Wt, 1)
    q◺ = ◺(q)
    A_21_β2β1_pre = Matrix{T}(undef, p, n)
    A_21_τ2τ1_pre = Matrix{T}(undef, l, n)
    A_21_Lγτ1_pre = Matrix{T}(undef, q◺, n)

    # construct A_21_β2β1_pre
    # X2t * Vinv, Vinv = Dinv - UUt
    X2t_U = nullObs.Xt * transpose(nullObs.Ut)
    mul!(A_21_β2β1_pre, -X2t_U, nullObs.Ut)
    @inbounds @simd for j in 1:n
        for i in 1:p
            A_21_β2β1_pre[i, j] += nullObs.Xt[i, j] * nullObs.Dinv[j] # first term
        end
    end    

    # construct A_21_τ2τ1_pre
    # W2t * D * Vinv .* Vinv * D, this is no longer symmetric.
    # compute W2t * D * Vinv .* Vinv 
    mul!(A_21_τ2τ1_pre, nullObs.Wt_D_Ut_kr_Utt, nullObs.Ut_kr_Ut) # third term
    @inbounds @simd for j in 1:n
        for i in 1:l
            A_21_τ2τ1_pre[i, j] += nullObs.Wt_D_Dinv[i, j] * nullObs.Dinv[j] # first term
            A_21_τ2τ1_pre[i, j] += -2 * nullObs.Wt_D_sqrtdiagDinv_UUt[i, j] * nullObs.sqrtDinv_UUt[j] # second term
        end
    end
    # right-multiply by D. 
    @inbounds @simd for j in 1:n
        for i in 1:l
            A_21_τ2τ1_pre[i, j] = A_21_τ2τ1_pre[i, j] * nullObs.expwτ[j]
        end
    end    

    # construct A_21_Lγτ1_pre
    # 2 * Cq' * (L'Z'(V^-1) ⊙ Z'(V^-1)) * Diagonal(expwτ)
    # nullObs.storage_q◺n is always Cq' * (L'Z'(V^-1) ⊙ Z'(V^-1)).
    A_21_Lγτ1_pre .=  2 * nullObs.storage_q◺n
    @inbounds @simd for j in 1:n
        for i in 1:l
            A_21_Lγτ1_pre[i, j] = A_21_Lγτ1_pre[i, j] * nullObs.expwτ[j]
        end
    end

    WSVarScoreTestBaseObs{T}(nullObs, n, p, q, l, A_21_β2β1_pre, A_21_τ2τ1_pre, A_21_Lγτ1_pre)
end

"""
    WSVarScoreTestObs

A base per-observation object for the score test of within-subject variance linear mixed model data instance.
H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0, 
for the full model of WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ).
We make use of the fitted null model.
"""
struct WSVarScoreTestObs{T <: BlasReal}
    # data
    testbaseobs         :: WSVarScoreTestBaseObs{T}   # the base object without information on X1 and W1. is it redundant?
    r_X1                :: Int              # number of test variables in X1
    r_W1                :: Int              # number of test variables in W1
    r                   :: Int              # number of total test variables, = r_X1 + r_W1.
    X1t                 :: Matrix{T}        # test variables for X
    W1t                 :: Matrix{T}        # test variables for W

    # working arrays
    ψ_1                 :: Vector{T}
    # ψ_β2 = testbaseobs.nullObs.∇β         
    # ψ_τ2 = testbaseobs.nullObs.∇τ
    # ψ_Lγ = vech(testbaseobs.nullObs.∇Lγ), be consistent with order of variables in WiSER.sandwich!()

    A_21_β2β1   :: Matrix{T} # p x r_X1, testbaseobs.A_21_β2β1_pre * X1
    A_21_τ2τ1   :: Matrix{T} # l x r_W1, testbaseobs.A_21_τ2τ1_pre * W1
    A_21_Lγτ1   :: Matrix{T} # q◺ x r_W1, testbaseobs.A_21_Lγτ1_pre * W1
end

function WSVarScoreTestObs(testbaseobs::WSVarScoreTestBaseObs{T}, X1obs::AbstractMatrix{T}, W1obs::AbstractMatrix{T}) where T <: BlasReal
    n, p, q, l = testbaseobs.n, testbaseobs.p, testbaseobs.q, testbaseobs.l
    q◺ = ◺(q)
    r_X1 = size(X1obs, 2)
    r_W1 = size(W1obs, 2)
    r = r_X1 + r_W1
    X1t = transpose(X1obs)
    W1t = transpose(W1obs)

    ψ_1 = Vector{T}(undef, r)
    ψ_β1 = @view(ψ_1[1:r_X1])
    ψ_τ1 = @view(ψ_1[r_X1+1:end])

    A_21_β2β1 = Matrix{T}(undef, p, r_X1)
    A_21_τ2τ1 = Matrix{T}(undef, l, r_W1)
    A_21_Lγτ1 = Matrix{T}(undef, q◺, r_W1)

    mul!(ψ_β1, X1t, reshape(testbaseobs.nullObs.Dinv_r - transpose(testbaseobs.nullObs.rt_UUt), :))
    mul!(ψ_τ1, -W1t, testbaseobs.nullObs.diagDVRV)

    mul!(A_21_β2β1, testbaseobs.A_21_β2β1_pre, X1obs)
    mul!(A_21_τ2τ1, testbaseobs.A_21_τ2τ1_pre, W1obs)
    mul!(A_21_Lγτ1, testbaseobs.A_21_Lγτ1_pre, W1obs)

    WSVarScoreTestObs{T}(testbaseobs, r_X1, r_W1, r, X1t, W1t, ψ_1, A_21_β2β1, A_21_τ2τ1, A_21_Lγτ1)
end

"""
    WSVarScoreTest

A base object for the score test of within-subject variance linear mixed model data instance.
H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0, 
for the full model of WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ).
We make use of the fitted null model.
"""
struct WSVarScoreTest{T <: BlasReal}
    nullmodel   :: WSVarLmmModel{T}
    testobs :: Vector{WSVarScoreTestObs{T}}
    
    # dimensions
    p          :: Int       # number of mean parameters in linear regression
    q          :: Int       # number of random effects
    l          :: Int       # number of parameters for modeling WS variability
    m          :: Int       # number of individuals/clusters
    nsum       :: Int       # number of observations (summed across individuals)
    r_X1                :: Int              # number of test variables in X1
    r_W1                :: Int              # number of test variables in W1
    r                   :: Int              # number of total test variables, = r_X1 + r_W1.

    # working arrays 
    ψ_1             :: Vector{T}        # length-r vector, sum_i testobs.ψ_1
    B_11            :: Matrix{T}        # r x r matrix.
    B_21            :: Matrix{T}        # (p + l + q◺) x r matrix. 
    A_21            :: Matrix{T}        # (p + l + q◺) x r matrix.
    # B_22 = nullmodel.B
    # inv(A_22) = nullmodel.Ainv
end

function WSVarScoreTest(nullmodel::WSVarLmmModel{T}, X1vec::Vector{<:AbstractMatrix{T}}, W1vec::Vector{<:AbstractMatrix{T}}; 
    testbaseobsvec::Union{Vector{WSVarScoreTestBaseObs{T}}, Nothing} = nothing) where T <: BlasReal
    if testbaseobsvec === nothing
        testbaseobsvec = [WSVarScoreTestBaseObs(obs) for obs in nullmodel.data]
    end
    testobsvec = [WSVarScoreTestObs(testbaseobs, X1obs, W1obs) for (testbaseobs, X1obs, W1obs) in zip(testbaseobsvec, X1vec, W1vec)]

    @assert nullmodel.isfitted[1] "Please fit the model first."

    p, q, l, m, nsum = nullmodel.p, nullmodel.q, nullmodel.l, nullmodel.m, nullmodel.nsum
    q◺ = ◺(q)
    r_X1, r_W1, r = testobsvec[1].r_X1, testobsvec[1].r_W1, testobsvec[1].r

    ψ_1 = Vector{T}(undef, r)
    B_11 = Matrix{T}(undef, r, r)
    B_21 = Matrix{T}(undef, p + l + q◺, r)
    A_21 = Matrix{T}(undef, p + l + q◺, r)

    # build ψ_1: sum_i testobs.ψ_1
    fill!(ψ_1, zero(T))
    for testobs in testobsvec
        ψ_1 .+= testobs.ψ_1
    end

    # build B_11: using BLAS.syr!()
    fill!(B_11, zero(T))
    for testobs in testobsvec 
        BLAS.syr!('U', T(1), testobs.ψ_1, B_11)
    end
    copytri!(B_11, 'U')
    lmul!(1 / m, B_11)

    # build B_21: using BLAS.ger!()
    fill!(B_21, zero(T))
    ψ_2 = Vector{T}(undef, p + l + q◺)
    ψ_β2 = @view ψ_2[1 : p]
    ψ_τ2 = @view ψ_2[p + 1 : p + l]
    ψ_Lγ = @view ψ_2[p + l + 1: end]
    for (obs, testobs) in zip(nullmodel.data, testobsvec)
        ψ_β2 .= obs.∇β
        ψ_τ2 .= obs.∇τ
        offset = 1
        @inbounds for j in 1:q, i in j:q
            ψ_Lγ[offset] = obs.∇Lγ[i, j]
            offset += 1
        end
        BLAS.ger!(T(1), ψ_2, testobs.ψ_1, B_21)
    end
    lmul!(1 / m, B_21)

    # build A_21: 1/m sum_i Ai_21.
    fill!(A_21, zero(T))
    A_21_β2β1 = @view A_21[1         : p    , 1        : r_X1]  
    A_21_τ2τ1 = @view A_21[p + 1     : p + l, r_X1 + 1 : r] 
    A_21_Lγτ1 = @view A_21[p + l + 1 : end  , r_X1 + 1 : r]
    for testobs in testobsvec
        A_21_β2β1 .+= testobs.A_21_β2β1
        A_21_τ2τ1 .+= testobs.A_21_τ2τ1
        A_21_Lγτ1 .+= testobs.A_21_Lγτ1
    end
    lmul!(1 / m, A_21)

    WSVarScoreTest{T}(nullmodel, testobsvec, p, q, l, m, nsum, r_X1, r_W1, r, ψ_1, B_11, B_21, A_21)
end

end