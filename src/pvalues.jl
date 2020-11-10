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
    one(T) / m * (transpose(ψ_1) * inv(B_11 -
        transpose(A_21) * nm.Ainv * B_21 -
        transpose(B_21) * nm.Ainv * A_21 +
        transpose(A_21) * nm.Ainv * nm.B * nm.Ainv * A_21
        ) * ψ_1
    )
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
        v1 = test_statistic(st, ψ_1p, B_11p, A_21p, B_21p, nm, st.m)
        p1 = Distributions.ccdf.(Chisq(st.r_X1), v1)
    else
        p1 = -one(T)
    end

    if st.r_W1 > 0
        ψ_1p  = @view st.ψ_1[(r_X1 + 1):end]
        B_11p = @view st.B_11[(r_X1 + 1):end, (r_X1 + 1):end]
        A_21p = @view st.A_21[:, (r_X1 + 1):end]
        B_21p = @view st.B_21[:, (r_X1 + 1):end]
        v2 = test_statistic(st, ψ_1p, B_11p, A_21p, B_21p, nm, st.m)
        p2 = Distributions.ccdf.(Chisq(st.r_W1), v2)
    else
        p2 = -one(T)
    end

    # for both
    if st.r_X1 > 0 && st.r_W1 > 0
        v3 = test_statistic(st, st.ψ_1, st.B_11, st.A_21, st.B_21, nm, st.m)
        p3 = Distributions.ccdf.(Chisq(st.r), v3)
    else
        p3 = -one(T)
    end
    p1, p2, p3
end
