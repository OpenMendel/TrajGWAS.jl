"""
    test_statistic(ψ_1, B_11, A_21, B_21)

Computes the test statistic for the score test
"""
function test_statistic(
    st::Union{WSVarScoreTest{T},WSVarScoreTestInvariant{T}},
    ψ_1::AbstractVector{T},
    B_11::AbstractMatrix{T},
    A_21::AbstractMatrix{T},
    B_21::AbstractMatrix{T},
    nm::WSVarLmmModel{T},
    m::Integer
) where {T <: BlasReal}
    Vψ_1 = B_11 - transpose(A_21) * nm.Ainv * B_21 -
        transpose(B_21) * nm.Ainv * A_21 +
        transpose(A_21) * nm.Ainv * nm.B * nm.Ainv * A_21
    eigfact = eigen!(Symmetric(Vψ_1))
    v1 = one(T) / sqrt(m) * transpose(eigfact.vectors) * ψ_1
    atol = 1e-8 # tolerance for determining rank
    rk = 0 # rank
    ts = zero(T)
    for j in 1:length(ψ_1)
        if eigfact.values[j] > atol
            ts += abs2(v1[j]) / eigfact.values[j]
            rk += 1
        end
    end
    return ts, rk
    # TODO: devise something to avoid reallocation in st.
end

"""
    pvalues!(st::WSVarScoreTest)

Returns three p-values for the score test:
for mean (`p1`), within-sample variability (`p2`), and both (`p3`).
-1 is returned for any skipped tests.
"""
function pvalues!(st::Union{WSVarScoreTest{T},WSVarScoreTestInvariant{T}}
) where {T <: BlasReal}
    nm = st.nullmodel
    r_X1, r_W1, r = st.r_X1, st.r_W1, st.r
    if st.r_X1 > 0
        ψ_1p  = @view st.ψ_1[1:r_X1]
        B_11p = @view st.B_11[1:r_X1, 1:r_X1]
        A_21p = @view st.A_21[:, 1:r_X1]
        B_21p = @view st.B_21[:, 1:r_X1]
        v1, r1 = test_statistic(st, ψ_1p, B_11p, A_21p, B_21p, nm, st.m)
        p1 = v1 ≤ 0 ? 1.0 : ccdf(Chisq(r1), v1)
    else
        p1 = -one(T)
    end

    if st.r_W1 > 0
        ψ_1p  = @view st.ψ_1[(r_X1 + 1):end]
        B_11p = @view st.B_11[(r_X1 + 1):end, (r_X1 + 1):end]
        A_21p = @view st.A_21[:, (r_X1 + 1):end]
        B_21p = @view st.B_21[:, (r_X1 + 1):end]
        v2, r2 = test_statistic(st, ψ_1p, B_11p, A_21p, B_21p, nm, st.m)
        p2 = v2 ≤ 0 ? 1.0 : ccdf(Chisq(r2), v2)
    else
        p2 = -one(T)
    end

    # for both
    if st.r_X1 > 0 && st.r_W1 > 0
        v3, r3 = test_statistic(st, st.ψ_1, st.B_11, st.A_21, st.B_21, nm, st.m)
        p3 = v3 ≤ 0 ? 1.0 : ccdf(Chisq(r3), v3)
    else
        p3 = -one(T)
    end
    p1, p2, p3
end
