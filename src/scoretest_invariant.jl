"""
    WSVarScoreTestInvariant

A base object for the score test of within-subject variance linear mixed model
data instance, independent of the input data.
Numbers of test variables should be predetermined.

H0: β1 = 0 and τ1 = 0, H1: β1 ≠ 0 or τ1 ≠ 0, for the full model of
WiSER (with parameters β = [β1, β2], τ = [τ1, τ2], and Lγ).
We make use of the fitted null model.
"""
struct WSVarScoreTestInvariant{T <: BlasReal}
    nullmodel   :: WSVarLmmModel{T}

    # dimensions
    p          :: Int       # #mean parameters in linear regression
    q          :: Int       # #random effects
    q◺         :: Int       # triangle number of q
    l          :: Int       # #parameters for modeling WS variability
    m          :: Int       # #individuals/clusters
    r_X1                :: Int      # #test variables in X1
    r_W1                :: Int      # #test variables in W1
    r                   :: Int      # #total test variables, = r_X1 + r_W1.

    # working arrays
    A_21_β2β1_rowsums :: Matrix{T} # for time-invariant X1, p x m.
    A_21_τ2τ1_rowsums :: Matrix{T} # for time-invariant W1, l x m.
    A_21_Lγτ1_rowsums :: Matrix{T} # for time-invariant W1, q◺ x m.

    #ψ_1pre :: Matrix{T} # r x m,
    ψ_1  :: Vector{T} # length-r, sum_i
    ψ_1obs :: Matrix{T} # r x m
    ψ_β1 :: AbstractVector{T} # = @view(ψ_1[1:r_X1])
    ψ_τ1 :: AbstractVector{T} # = @view(ψ_1[r_X1+1:end])

    # ψ_β1 = transpose(X1) * ψ_β1_pre.
    ψ_β1_pre :: Vector{T} # length-m.
    # ψ_τ1 = transpose(W1) * ψ_τ1_pre.
    ψ_τ1_pre :: Vector{T} # length-m.

    ψ_2obs  :: Matrix{T} # (p + l + q◺) x m.

    B_11 :: Matrix{T} # size r x r
    B_21 :: Matrix{T} # p + l + q◺ x r
    A_21 :: Matrix{T} # p + l + q◺ x r
    A_21_β2β1   :: AbstractMatrix{T} # p x r_X1, testbaseobs.A_21_β2β1_pre * X1
    A_21_τ2τ1   :: AbstractMatrix{T} # l x r_W1, testbaseobs.A_21_τ2τ1_pre * W1
    A_21_Lγτ1   :: AbstractMatrix{T} # q◺ x r_W1, testbaseobs.A_21_Lγτ1_pre * W1
    AinvBAinv   :: AbstractMatrix{T} # r x r.
    tmp_sr      :: AbstractMatrix{T} # p + l + q◺ x r.
    tmp_srx1    :: AbstractMatrix{T} # p + l + q◺ x r_X1.
    tmp_srw1    :: AbstractMatrix{T} # p + l + q◺ x r_W1.
    tmp_rr      :: AbstractMatrix{T} # r x r.
    tmp_rx1rx1  :: AbstractMatrix{T} # r_X1 x r_X1.
    tmp_rw1rw1  :: AbstractMatrix{T} # r_W1 x r_W1.
    tmp_r       :: AbstractVector{T}
    tmp_rx1     :: AbstractVector{T}
    tmp_rw1     :: AbstractVector{T}
end

function WSVarScoreTestInvariant(
    nullmodel::WSVarLmmModel{T},
    r_X1::Int, r_W1::Int;
    testbaseobsvec::Union{Vector{WSVarScoreTestBaseObs{T}}, Nothing} = nothing
) where {T <: BlasReal}
    @assert r_X1 >= 0 && r_W1 >= 0
    @assert r_X1 > 0 || r_W1 > 0
    if testbaseobsvec === nothing
        testbaseobsvec = [WSVarScoreTestBaseObs(obs) for obs in nullmodel.data]
    end
    @assert nullmodel.isfitted[1] "Please fit the model first."

    p, q, l, m = nullmodel.p, nullmodel.q, nullmodel.l, nullmodel.m
    q◺ = ◺(q)
    r = r_X1 + r_W1

    A_21_β2β1_rowsums = Matrix{T}(undef, p,  m)
    A_21_τ2τ1_rowsums = Matrix{T}(undef, l,  m)
    A_21_Lγτ1_rowsums = Matrix{T}(undef, q◺, m)
    for (i, baseobs) in enumerate(testbaseobsvec)
        A_21_β2β1_rowsums[:, i] .= baseobs.A_21_β2β1_rowsum
        A_21_τ2τ1_rowsums[:, i] .= baseobs.A_21_τ2τ1_rowsum
        A_21_Lγτ1_rowsums[:, i] .= baseobs.A_21_Lγτ1_rowsum
    end

    ψ_1  = Vector{T}(undef, r)
    ψ_1obs = Matrix{T}(undef, r, m)
    ψ_β1 = @view(ψ_1[1:r_X1])
    ψ_τ1 = @view(ψ_1[(r_X1 + 1):end])

    # ψ_β1 = transpose(X1) * ψ_β1_pre.
    ψ_β1_pre = Vector{T}(undef, m)
    # ψ_τ1 = transpose(W1) * ψ_τ1_pre.
    ψ_τ1_pre = Vector{T}(undef, m)
    for (i, obs) in enumerate(nullmodel.data)
        ψ_β1_pre[i] = sum(obs.Dinv_r - transpose(obs.rt_UUt))
        ψ_τ1_pre[i] = -sum(obs.diagDVRV)
    end

    # build ψ_2obs
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
    A_21_β2β1 = @view A_21[1:p, 1:r_X1]
    A_21_τ2τ1 = @view A_21[(p + 1):(p + l), (r_X1 + 1):r]
    A_21_Lγτ1 = @view A_21[(p + l + 1):end, (r_X1 + 1):r]

    AinvBAinv = nullmodel.Ainv * nullmodel.B * nullmodel.Ainv
    tmp_sr      = Matrix{T}(undef, p + l + q◺, r)
    tmp_srx1    = Matrix{T}(undef, p + l + q◺, r_X1)
    tmp_srw1    = Matrix{T}(undef, p + l + q◺, r_W1)
    tmp_rr      = Matrix{T}(undef, r, r)
    tmp_rx1rx1  = Matrix{T}(undef, r_X1, r_X1)
    tmp_rw1rw1  = Matrix{T}(undef, r_W1, r_W1)
    tmp_r       = Vector{T}(undef, r)
    tmp_rx1     = Vector{T}(undef, r_X1)
    tmp_rw1     = Vector{T}(undef, r_W1)

    WSVarScoreTestInvariant{T}(nullmodel, p, q, q◺, l, m, r_X1, r_W1, r,
        A_21_β2β1_rowsums, A_21_τ2τ1_rowsums, A_21_Lγτ1_rowsums,
        ψ_1, ψ_1obs, ψ_β1, ψ_τ1, ψ_β1_pre, ψ_τ1_pre, ψ_2obs,
        B_11, B_21, A_21, A_21_β2β1, A_21_τ2τ1, A_21_Lγτ1,
        AinvBAinv, tmp_sr, tmp_srx1, tmp_srw1, tmp_rr, tmp_rx1rx1, tmp_rw1rw1,
        tmp_r, tmp_rx1, tmp_rw1
    )
end

"""
    test!(st::WSVarScoreTestInvariant, X1, W1)

Performs the score test, and returns the three p-values on the
time-invariant test data.

# Arguments

- `st::WSVarScoreTestInvariant`: data structure for the test
- `X1::Union{Nothing, AbstractMatrix{T}}`: m x r_X1 data, nothing if r_X1 == 0
- `W1::Union{Nothing, AbstractMatrix{T}}`: m x r_W1 data, nothing if r_W1 == 0
"""
function test!(st::WSVarScoreTestInvariant,
    X1::Union{Nothing, AbstractVecOrMat{T}},
    W1::Union{Nothing, AbstractVecOrMat{T}}) where {T <: BlasReal}

    if st.r_X1 > 0
        @assert size(X1, 1) == st.m
        @assert size(X1, 2) == st.r_X1
    end
    if st.r_W1 > 0
        @assert size(W1, 1) == st.m
        @assert size(W1, 2) == st.r_W1
    end

    nm = st.nullmodel
    p, q, l, m = nm.p, nm.q, nm.l, nm.m
    r_X1, r_W1, r = st.r_X1, st.r_W1, st.r

    # build ψ_1: sum_i testobs.ψ_1
    if r_X1 > 0
        st.ψ_1obs[1:r_X1, :] .= transpose(X1) .* transpose(st.ψ_β1_pre)
        mul!(st.ψ_β1, transpose(X1), st.ψ_β1_pre)
    end
    if r_W1 > 0
        st.ψ_1obs[(r_X1 + 1):end, :] .= transpose(W1) .* transpose(st.ψ_τ1_pre)
        mul!(st.ψ_τ1, transpose(W1), st.ψ_τ1_pre)
    end

    # build B_11: using BLAS.syrk!()
    BLAS.syrk!('U', 'N', one(T), st.ψ_1obs, zero(T), st.B_11)
    copytri!(st.B_11, 'U')
    lmul!(one(T) / m, st.B_11)

    # build B_21
    mul!(st.B_21, st.ψ_2obs, transpose(st.ψ_1obs), one(T) / m, zero(T))

    # build A_21: 1/m sum_i Ai_21.
    fill!(st.A_21, zero(T))
    if r_X1 > 0
        mul!(st.A_21_β2β1, st.A_21_β2β1_rowsums, X1, one(T) / m, zero(T))
    end
    if r_W1 > 0
        mul!(st.A_21_τ2τ1, st.A_21_τ2τ1_rowsums, W1, one(T) / m, zero(T))
        mul!(st.A_21_Lγτ1, st.A_21_Lγτ1_rowsums, W1, one(T) / m, zero(T))
    end

    pvalues!(st)
end
