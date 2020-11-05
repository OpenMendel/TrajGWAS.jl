"""
    WSVarScoreTestBaseObs

A base per-observation object for the score test of within-subject variance
linear mixed model data instance without information on X1 or W1.
Contains base variables for testing
H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0,
for the full model of WiSER (with parameters β = [β1, β2],
     τ = [τ1, τ2], and Lγ),
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

    A_21_β2β1_rowsum :: Vector{T} # for time-invariant X1, length p
    A_21_τ2τ1_rowsum :: Vector{T} # for time-invariant W1, length l
    A_21_Lγτ1_rowsum :: Vector{T} # for time-invariant W1, length q◺
end

function WSVarScoreTestBaseObs(nullObs::WSVarLmmObs{T}) where T <: BlasReal
    n, p, q, l = size(nullObs.Xt, 2), size(nullObs.Xt, 1),
        size(nullObs.Zt, 1), size(nullObs.Wt, 1)
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
            A_21_β2β1_pre[i, j] += nullObs.Xt[i, j] *
                nullObs.Dinv[j] # first term
        end
    end

    # construct A_21_τ2τ1_pre
    # W2t * D * Vinv .* Vinv * D, this is no longer symmetric.
    # compute W2t * D * Vinv .* Vinv
    mul!(A_21_τ2τ1_pre, nullObs.Wt_D_Ut_kr_Utt, nullObs.Ut_kr_Ut) # third term
    @inbounds @simd for j in 1:n
        for i in 1:l
            # first term
            A_21_τ2τ1_pre[i, j] += nullObs.Wt_D_Dinv[i, j] * nullObs.Dinv[j]
            # second term
            A_21_τ2τ1_pre[i, j] += -2 * nullObs.Wt_D_sqrtdiagDinv_UUt[i, j] *
                nullObs.sqrtDinv_UUt[j]
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

    A_21_β2β1_rowsum = reshape(sum(A_21_β2β1_pre; dims=2), :)
    A_21_τ2τ1_rowsum = reshape(sum(A_21_τ2τ1_pre; dims=2), :)
    A_21_Lγτ1_rowsum = reshape(sum(A_21_Lγτ1_pre; dims=2), :)

    WSVarScoreTestBaseObs{T}(nullObs, n, p, q, l, A_21_β2β1_pre, A_21_τ2τ1_pre,
        A_21_Lγτ1_pre, A_21_β2β1_rowsum, A_21_τ2τ1_rowsum, A_21_Lγτ1_rowsum)
end

"""
    WSVarScoreTestObs

A base per-observation object for the score test of
within-subject variance linear mixed model data instance.
H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0,
for the full model of WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ).
We make use of the fitted null model.
"""
struct WSVarScoreTestObs{T <: BlasReal}
    # data
    testbaseobs         :: WSVarScoreTestBaseObs{T}
    r_X1                :: Int  # number of test variables in X1
    r_W1                :: Int  # number of test variables in W1
    r                   :: Int  # number of total test variables, = r_X1 + r_W1.
    #X1t                 :: Matrix{T}        # test variables for X
    #W1t                 :: Matrix{T}        # test variables for W

    # working arrays
    ψ_1                 :: Vector{T}
    ψ_β1                :: AbstractVector{T}
    ψ_τ1                :: AbstractVector{T}
    # ψ_β2 = testbaseobs.nullObs.∇β
    # ψ_τ2 = testbaseobs.nullObs.∇τ
    # ψ_Lγ = vech(testbaseobs.nullObs.∇Lγ),
    # be consistent with order of variables in WiSER.sandwich!()

    A_21_β2β1   :: Matrix{T} # p x r_X1, testbaseobs.A_21_β2β1_pre * X1
    A_21_τ2τ1   :: Matrix{T} # l x r_W1, testbaseobs.A_21_τ2τ1_pre * W1
    A_21_Lγτ1   :: Matrix{T} # q◺ x r_W1, testbaseobs.A_21_Lγτ1_pre * W1
end

function WSVarScoreTestObs(testbaseobs::WSVarScoreTestBaseObs{T},
    r_X1, r_W1
    #X1obs::AbstractMatrix{T}, W1obs::AbstractMatrix{T}
) where {T <: BlasReal}
    n, p, q, l = testbaseobs.n, testbaseobs.p, testbaseobs.q, testbaseobs.l
    q◺ = ◺(q)
    #r_X1 = size(X1obs, 2)
    #r_W1 = size(W1obs, 2)
    r = r_X1 + r_W1
    #X1t = Matrix{T}(undef, r_X1, n)
    #W1t = Matrix{T}(undef, r_W1, n)

    ψ_1 = Vector{T}(undef, r)
    ψ_β1 = @view(ψ_1[1:r_X1])
    ψ_τ1 = @view(ψ_1[r_X1+1:end])

    A_21_β2β1 = Matrix{T}(undef, p, r_X1)
    A_21_τ2τ1 = Matrix{T}(undef, l, r_W1)
    A_21_Lγτ1 = Matrix{T}(undef, q◺, r_W1)

    WSVarScoreTestObs{T}(testbaseobs, r_X1, r_W1, r, #X1t, W1t,
        ψ_1, ψ_β1, ψ_τ1, A_21_β2β1, A_21_τ2τ1, A_21_Lγτ1)
end

"""
    update!(testobs::WSVarScoreTestObs, X1obs, W1obs)

Updates WSVarScoreTestObs based on the test data.
"""
function update!(testobs::WSVarScoreTestObs,
    X1obs::Union{Nothing, AbstractMatrix{T}},
    W1obs::Union{Nothing, AbstractMatrix{T}}
) where {T <: BlasReal}
    testbaseobs = testobs.testbaseobs
    if testobs.r_X1 > 0
        mul!(testobs.ψ_β1,
            transpose(X1obs),
            reshape(testbaseobs.nullObs.Dinv_r -
                transpose(testbaseobs.nullObs.rt_UUt), :
            )
        )
        mul!(testobs.A_21_β2β1, testbaseobs.A_21_β2β1_pre, X1obs)
    end
    if testobs.r_W1 > 0
        mul!(testobs.ψ_τ1, -transpose(W1obs), testbaseobs.nullObs.diagDVRV)
        mul!(testobs.A_21_τ2τ1, testbaseobs.A_21_τ2τ1_pre, W1obs)
        mul!(testobs.A_21_Lγτ1, testbaseobs.A_21_Lγτ1_pre, W1obs)
    end
end

"""
    WSVarScoreTest

A base object for the score test of within-subject vLMM data instance.
H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0,
for the full model of WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ).
We make use of the fitted null model.
"""
struct WSVarScoreTest{T <: BlasReal}
    nullmodel   :: WSVarLmmModel{T}
    testobsvec :: Vector{WSVarScoreTestObs{T}}

    # dimensions
    p          :: Int       # #mean parameters in linear regression
    q          :: Int       # #random effects
    l          :: Int       # #parameters for modeling WS variability
    m          :: Int       # #individuals/clusters
    nsum       :: Int       # #observations (summed across individuals)
    r_X1                :: Int      # #test variables in X1
    r_W1                :: Int      # #test variables in W1
    r                   :: Int      # #total test variables, = r_X1 + r_W1.

    # working arrays
    ψ_1             :: Vector{T}        # length-r vector, sum_i testobs.ψ_1
    ψ_1obs          :: Matrix{T}        # r x m
    ψ_2obs          :: Matrix{T}        # (p + l + q◺) x m
    B_11            :: Matrix{T}        # r x r matrix.
    B_21            :: Matrix{T}        # (p + l + q◺) x r matrix.
    A_21            :: Matrix{T}        # (p + l + q◺) x r matrix.
    # B_22 = nullmodel.B
    # inv(A_22) = nullmodel.Ainv
end

function WSVarScoreTest(nullmodel::WSVarLmmModel{T},
    r_X1::Int, r_W1::Int;
    testbaseobsvec::Union{Vector{WSVarScoreTestBaseObs{T}}, Nothing} = nothing
) where {T <: BlasReal}

    @assert r_X1 >= 0 && r_W1 >= 0
    @assert r_X1 > 0 || r_W1 > 0
    if testbaseobsvec === nothing
        testbaseobsvec = [WSVarScoreTestBaseObs(obs) for obs in nullmodel.data]
    end
    testobsvec = [WSVarScoreTestObs(testbaseobs, r_X1, r_W1) for
        testbaseobs in testbaseobsvec]

    @assert nullmodel.isfitted[1] "Please fit the model first."

    p, q, l, m, nsum = nullmodel.p, nullmodel.q,
                        nullmodel.l, nullmodel.m, nullmodel.nsum
    q◺ = ◺(q)
    r = r_X1 + r_W1

    ψ_1 = Vector{T}(undef, r)
    ψ_1obs = Matrix{T}(undef, r, m)
    ψ_2obs = Matrix{T}(undef, p + l + q◺, m)
    fill!(ψ_2obs, zero(T))
    for (i, obs) in enumerate(nullmodel.data)
        ψ_β2 = @view ψ_2obs[1:p, i]
        ψ_τ2 = @view ψ_2obs[(p + 1):(p + l), i]
        ψ_Lγ = @view ψ_2obs[(p + l + 1):end, i]
        ψ_β2 .= obs.∇β
        ψ_τ2 .= obs.∇τ
        offset = 1
        @inbounds for j in 1:q, i in j:q
            ψ_Lγ[offset] = obs.∇Lγ[i, j]
            offset += 1
        end
    end

    B_11 = Matrix{T}(undef, r, r)
    B_21 = Matrix{T}(undef, p + l + q◺, r)
    A_21 = Matrix{T}(undef, p + l + q◺, r)

    WSVarScoreTest{T}(nullmodel, testobsvec, p, q, l, m, nsum, r_X1, r_W1, r,
        ψ_1, ψ_1obs, ψ_2obs, B_11, B_21, A_21)
end

"""
    test!(st::WSVarScoreTest, X1Vec, W1Vec)

Performs the score test, and returns the three p-values on the
time-variant test data
"""
function test!(st::WSVarScoreTest,
    X1vec::Union{Nothing, Vector{<:AbstractMatrix{T}}},
    W1vec::Union{Nothing, Vector{<:AbstractMatrix{T}}}) where {T <: BlasReal}

    if X1vec === nothing
        @assert st.r_X1 == 0
        X1vec = repeat([nothing], length(st.testobsvec))
    elseif W1vec === nothing
        @assert st.r_W1 == 0
        W1vec = repeat([nothing], length(st.testobsvec))
    end

    nm = st.nullmodel
    p, q, l, m, nsum = nm.p, nm.q, nm.l, nm.m, nm.nsum
    r_X1, r_W1, r = st.r_X1, st.r_W1, st.r

    # update testobsvec first
    for (testobs, X1obs, W1obs) in zip(st.testobsvec, X1vec, W1vec)
        if r_X1 > 0
            @assert size(X1obs, 1) == testobs.testbaseobs.n
            @assert size(X1obs, 2) == st.r_X1
        end
        if r_W1 > 0
            @assert size(W1obs, 1) == testobs.testbaseobs.n
            @assert size(W1obs, 2) == st.r_W1
        end
        update!(testobs, X1obs, W1obs)
    end

    # build ψ_1: sum_i testobs.ψ_1
    fill!(st.ψ_1, zero(T))
    for (i, testobs) in enumerate(st.testobsvec)
        st.ψ_1obs[:, i] .= testobs.ψ_1
        st.ψ_1 .+= testobs.ψ_1
    end

    # build B_11: using BLAS.syr!()
    fill!(st.B_11, zero(T))
    for testobs in st.testobsvec
        BLAS.syr!('U', T(1), testobs.ψ_1, st.B_11)
    end
    copytri!(st.B_11, 'U')
    lmul!(one(T) / m, st.B_11)

    # build B_21
    mul!(st.B_21, st.ψ_2obs, transpose(st.ψ_1obs), one(T) / m, zero(T))

    # build A_21: 1/m sum_i Ai_21.
    fill!(st.A_21, zero(T))
    A_21_β2β1 = @view st.A_21[1         : p    , 1        : r_X1]
    A_21_τ2τ1 = @view st.A_21[p + 1     : p + l, r_X1 + 1 : r]
    A_21_Lγτ1 = @view st.A_21[p + l + 1 : end  , r_X1 + 1 : r]
    if r_X1 > 0
        for testobs in st.testobsvec
            A_21_β2β1 .+= testobs.A_21_β2β1
        end
    end
    if r_W1 > 0
        for testobs in st.testobsvec
            A_21_τ2τ1 .+= testobs.A_21_τ2τ1
            A_21_Lγτ1 .+= testobs.A_21_Lγτ1
        end
    end
    lmul!(1 / m, st.A_21)

    pvalues(st)
end
