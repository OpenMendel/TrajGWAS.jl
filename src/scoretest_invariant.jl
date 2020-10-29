# """
#     WSVarScoreTestInvariant

# A base object for the score test of within-subject variance linear mixed model
# data instance, independent of the input data.
# Numbers of test variables should be predetermined.

# H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0, for the full model of
# WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ).
# We make use of the fitted null model.
# """
# struct WSVarScoreTestInvariant{T <: BlasReal}
#     nullmodel   :: WSVarLmmModel{T}

#     # dimensions
#     p          :: Int       # #mean parameters in linear regression
#     q          :: Int       # #random effects
#     q◺         :: Int       # triangle number of q
#     l          :: Int       # #parameters for modeling WS variability
#     m          :: Int       # #individuals/clusters
#     r_X1                :: Int      # #test variables in X1
#     r_W1                :: Int      # #test variables in W1
#     r                   :: Int      # #total test variables, = r_X1 + r_W1.

#     # working arrays
#     A_21_β2β1_rowsums :: Matrix{T} # for time-invariant X1, p x m.
#     A_21_τ2τ1_rowsums :: Matrix{T} # for time-invariant W1, l x m.
#     A_21_Lγτ1_rowsums :: Matrix{T} # for time-invariant W1, q◺ x m.

#     ψ_1pre :: Matrix{T} # r x m,
#     ψ_1 :: Vector{T} # length-r, sum_i
#     # ψ_β2 = testbaseobs.nullObs.∇β
#     # ψ_τ2 = testbaseobs.nullObs.∇τ
#     # ψ_Lγ = vech(testbaseobs.nullObs.∇Lγ), be consistent with order of variables in WiSER.sandwich!()

#     ψ_β1 :: AbstractMatrix{T} # = @view(ψ_1[1:r_X1, :])
#     ψ_τ1 :: AbstractMatrix{T} # = @view(ψ_1[r_X1+1:end, :])

#     A_21_β2β1   :: Matrix{T} # p x r_X1, testbaseobs.A_21_β2β1_pre * X1
#     A_21_τ2τ1   :: Matrix{T} # l x r_W1, testbaseobs.A_21_τ2τ1_pre * W1
#     A_21_Lγτ1   :: Matrix{T} # q◺ x r_W1, testbaseobs.A_21_Lγτ1_pre * W1

# end

# function WSVarScoreTestInvariant(
#     nullmodel::WSVarLmmModel{T}) where T <: BlasReal

#     @assert nullmodel.isfitted[1] "Please fit the model first."
#     testbaseobsvec = [WSVarScoreTestBaseObs(obs) for obs in nullmodel.data]

#     p, q, l, m = nullmodel.p, nullmodel.q, nullmodel.l, nullmodel.m
#     q◺ = ◺(q)

#     for (i, baseobs) in enumerate(testbaseobsvec)
#         A_21_β2β1_rowsums[:, i] .= baseobs.A_21_β2β1_rowsum
#         A_21_τ2τ1_rowsums[:, i] .= baseobs.A_21_τ2τ1_rowsum
#         A_21_Lγτ1_rowsums[:, i] .= baseobs.A_21_Lγτ1_rowsum
#     end
#     WSVarScoreTestInvariantPre{T}(nullmodel, p, q, q◺, l, m,
#         A_21_β2β1_rowsums, A_21_τ2τ1_rowsums, A_21_Lγτ1_rowsums)
# end
