#Code to do type I error simulation for vGWAS.jl
# after reworking SPA

using Interpolations
import Optim: minimizer, optimize, LBFGS, NelderMead

function ecgf(z::AbstractVector; knots::Integer=10000)

    range = (-200.0, 200.0)
    idx0 = quantile.(Cauchy(), (1:knots)/(knots+1))
    idx1 = idx0 .* maximum(range) / maximum(idx0)

    K0_int, K1_int, K2_int = ecgf(idx1, z)

    return K0_int, K1_int, K2_int
end

function ecgf(vs::AbstractVector, zs::AbstractVector)
    exp_zv = Vector{Float64}(undef, length(zs))
    zexp_zv = Vector{Float64}(undef, length(zs))
    z2exp_zv = Vector{Float64}(undef, length(zs))
    K0 = Vector{Float64}(undef, length(vs))
    K1 = similar(K0)
    K2 = similar(K0)
    for (i, v) in enumerate(vs)
        exp_zv .= zs .* v
        alpha = maximum(exp_zv)
        exp_zv .= exp.(exp_zv .- alpha)

        zexp_zv .= zs .* exp_zv
        z2exp_zv .= zs .^ 2 .* exp_zv
        M0 = mean(exp_zv)
        M1 = mean(zexp_zv)
        M2 = mean(z2exp_zv)
        K0[i] = log(M0) + alpha
        K1[i] = M1 / M0
        K2[i] = M2 / M0 - (M1 / M0) ^ 2
    end
    K0_int = LinearInterpolation(vs, convert(Vector{Float64}, K0);
        extrapolation_bc=Linear())
    K1_int = LinearInterpolation(vs, convert(Vector{Float64}, K1);
        extrapolation_bc=Linear())
    K2_int = LinearInterpolation(vs, convert(Vector{Float64}, K2);
        extrapolation_bc=Linear())
    return K0_int, K1_int, K2_int
end

struct vGWASEcgfCollection
    K0_β
    K1_β
    K2_β
    K0_τ
    K1_τ
    K2_τ
    K0_βτ
    K1_βτ
    K2_βτ
end
function ecgf(st::WSVarScoreTestInvariant)
    K0_β, K1_β, K2_β = ecgf(st.ψ_β1_pre)
    K0_τ, K1_τ, K2_τ = ecgf(st.ψ_τ1_pre)
    K0_βτ, K1_βτ, K2_βτ = ecgf(st.ψ_βτ_pre)
    vGWASEcgfCollection(K0_β, K1_β, K2_β, K0_τ, K1_τ, K2_τ, K0_βτ, K1_βτ, K2_βτ)
end

function ecgf_g(cnts, vals_norm, K0, K1, K2)
    function K0_(z; tmp=Vector(undef, 3))
        for i in 1:3
            @inbounds tmp[i] = K0(vals_norm[i] * z)
        end
        cnts[1] * tmp[1] + cnts[2] * tmp[2] + cnts[3] * tmp[3]
    end
    function K1_(z; tmp=Vector(undef, 3))
        for i in 1:3
            @inbounds tmp[i] = vals_norm[i] * K1(vals_norm[i] * z)
        end
        cnts[1] * tmp[1] + cnts[2] * tmp[2] + cnts[3] * tmp[3]
    end
    function K2_(z; tmp=Vector(undef, 3))
        for i in 1:3
            @inbounds tmp[i] = vals_norm[i] ^ 2 * K2(vals_norm[i] * z)
        end
        cnts[1] * tmp[1] + cnts[2] * tmp[2] + cnts[3] * tmp[3]
    end
    K0_, K1_, K2_
end

function ecgf_g(vals_norm, K0, K1, K2)
    function K0_(z; tmp=Vector(undef, n))
        for i in 1:length(vals)
            @inbounds tmp[i] = K0(vals_norm[i] * z)
        end
        sum(tmp)
    end
    function K1_(z; tmp=Vector(undef, n))
        for i in 1:length(vals)
            @inbounds tmp[i] = vals[i] * K1(vals_norm[i] * z)
        end
        sum(tmp)
    end
    function K2_(z; tmp=Vector(undef, n))
        for i in 1:length(vals)
            @inbounds tmp[i] = vals[i] ^ 2 * K2(vals_norm[i] * z)
        end
        sum(tmp)
    end
    K0_, K1_, K2_
end

function spa(g::AbstractVector, pre_vec::AbstractVector, p_alt, K0, K1, K2;
    tmp_g=similar(g), tmp_g2=similar(g, 3), r=2.0, r_var=var(pre_vec))
    # involves internal normalization
    m = mean(skipmissing(g))
    s = std(skipmissing(g))
    tmp_g .= (g .- m) ./ s
    s_tau = dot(tmp_g, pre_vec)
    cutoff = r * sqrt(r_var * sum(x -> x ^ 2, tmp_g))
    return _get_pval(s_tau, cutoff, p_alt, cnts, vals_norm, K0, K1, K2, tmp_g2)
end

function spa(g::AbstractVector,
    st::WSVarScoreTestInvariant,
    p_alt,
    Ks::vGWASEcgfCollection;
    genotypes=true,
    tmp_g=similar(g),
    tmp_g2=genotypes ? similar(g, 3) : similar(g),
    r=2.0)

    # involves internal normalization
    m = mean(skipmissing(g))
    s = std(skipmissing(g))
    tmp_g .= (g .- m) ./ s

    s_β = dot(tmp_g, st.ψ_β1_pre)
    s_τ = dot(tmp_g, st.ψ_τ1_pre)
    s_βτ = dot(tmp_g, st.ψ_βτ_pre)

    cutoff_factor = sum(x -> x^2, tmp_g)
    cutoff_β = r * sqrt(st.var_β1_pre * cutoff_factor)
    cutoff_τ = r * sqrt(st.var_τ1_pre * cutoff_factor)
    cutoff_βτ = r * sqrt(st.var_βτ_pre * cutoff_factor)
    if genotypes
        cnts = map(x -> count(x .== g), [0.0, 1.0, 2.0])
        @assert sum(cnts) == length(g) "With genotypes == true, the values in g must be 0, 1, or 2."
        vals_norm = ([0.0, 1.0, 2.0] .- m) ./ s
        p_β = _get_pval(s_β, cutoff_β, p_alt[1], cnts, vals_norm,
            Ks.K0_β, Ks.K1_β, Ks.K2_β, tmp_g2)
        p_τ = _get_pval(s_τ, cutoff_τ, p_alt[2], cnts, vals_norm,
            Ks.K0_τ, Ks.K1_τ, Ks.K2_τ, tmp_g2)
        p_βτ = _get_pval(s_βτ, cutoff_βτ, p_alt[3], cnts, vals_norm,
            Ks.K0_βτ, Ks.K1_βτ, Ks.K2_βτ, tmp_g2)
    else
        vals_norm = (g .- m) ./ s
        p_β = _get_pval(s_β, cutoff_β, p_alt[1], vals_norm,
            Ks.K0_β, Ks.K1_β, Ks.K2_β, tmp_g2)
        p_τ = _get_pval(s_τ, cutoff_τ, p_alt[2], vals_norm,
            Ks.K0_τ, Ks.K1_τ, Ks.K2_τ, tmp_g2)
        p_βτ = _get_pval(s_βτ, cutoff_βτ, p_alt[3], vals_norm,
            Ks.K0_βτ, Ks.K1_βτ, Ks.K2_βτ, tmp_g2)  
    end
    p_β, p_τ, p_βτ
end

function _get_pval(s, cutoff, p_alt, cnts, vals_norm, K0, K1, K2, tmp3)
    if abs(s) < cutoff
        return p_alt
    else
        K0_, K1_, K2_ = ecgf_g(cnts, vals_norm, K0, K1, K2)
        K0__(x) = K0_(x; tmp=tmp3)
        K1__(x) = K1_(x; tmp=tmp3)
        K2__(x) = K2_(x; tmp=tmp3)
        return _spa_pval(s, K0__, K1__, K2__)
    end
end

function _get_pval(s, cutoff, p_alt, vals_norm, K0, K1, K2, tmp3)
    if abs(s) < cutoff
        return p_alt
    else
        K0_, K1_, K2_ = ecgf_g(vals_norm, K0, K1, K2)
        K0__(x) = K0_(x; tmp=tmp3)
        K1__(x) = K1_(x; tmp=tmp3)
        K2__(x) = K2_(x; tmp=tmp3)
        return _spa_pval(s, K0__, K1__, K2__)
    end
end

function _spa_pval(s::AbstractFloat, K0, K1, K2)
    f = x -> K0(x[1]) - x[1] * s
    function g!(storage, x)
        storage[1] = K1(x[1]) - s
    end
    function h!(storage, x)
        storage[1] = K2(x[1])
    end
    r = optimize(f, g!, h!, [0.0], LBFGS())
    zeta = minimizer(r)[1]

    omega = sign(zeta) * sqrt(max(0, 2 * (zeta * s - K0(zeta))))
    nu = zeta * sqrt(K2(zeta))
    z2 = omega + 1.0/omega * log(nu / omega)

    return ccdf(Normal(), abs(z2)) + cdf(Normal(), -abs(z2))
end
