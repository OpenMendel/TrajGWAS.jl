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
        extrapolation_bc=Line())
    K1_int = LinearInterpolation(vs, convert(Vector{Float64}, K1);
        extrapolation_bc=Line())
    K2_int = LinearInterpolation(vs, convert(Vector{Float64}, K2);
        extrapolation_bc=Line())
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

"""
    spa(g, st, p_alt, Ks; mode, g_norm, tmp_ecgf, r)

Perform saddlepoint approximation for a single variant. 

## arguments 
- `g`::AbstractVector:  genotype vector
- `st::WSVarScoreTestInvariant`
- `p_alt` : alternative p-value if saddlepoint approximation is skipped due to cutoff `r`
- `Ks::vGWASEcgfCollection`
- `mode`: `:bed` for hard calls of 0, 1, 2; `:ukbbgen` for bgen format provided by UK Biobank; anything else for others
"""
function spa(g::AbstractVector,
    st::WSVarScoreTestInvariant,
    p_alt,
    Ks::vGWASEcgfCollection;
    g_norm=similar(g),
    cnts=nothing,
    ref_vals= cnts !== nothing ? similar(g, length(cnts)) : nothing,
    vals_norm= cnts !== nothing ? similar(g, length(cnts)) : g_norm,
    tmp_ecgf = cnts !== nothing ? similar(g, length(cnts)) : similar(g),
    r=2.0
    )
    # involves internal normalization
    m = mean(skipmissing(g))
    s = std(skipmissing(g))
    g_norm .= (g .- m) ./ s

    if ref_vals !== nothing
        vals_norm .= (ref_vals .- m) ./ s
        @inbounds for i in 1:length(vals_norm)
            if isnan(vals_norm[i])
                vals_norm[i] = zero(eltype(vals_norm))
            end
        end
    end
    cutoff_factor = sum(x -> x^2, g_norm)
    p_β = spa(g_norm, st.ψ_β1_pre, vals_norm,
        cutoff_factor, r, p_alt[1], Ks.K0_β, Ks.K1_β, Ks.K2_β;
        cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_β1_pre)
    p_τ = spa(g_norm, st.ψ_τ1_pre, vals_norm,
        cutoff_factor, r, p_alt[2], Ks.K0_τ, Ks.K1_τ, Ks.K2_τ;
        cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_τ1_pre)
    p_βτ = spa(g_norm, st.ψ_βτ_pre, vals_norm,
        cutoff_factor, r, p_alt[3], Ks.K0_βτ, Ks.K1_βτ, Ks.K2_βτ;
        cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_βτ_pre)
    # s_β = dot(g_norm, st.ψ_β1_pre)
    # s_τ = dot(g_norm, st.ψ_τ1_pre)
    # s_βτ = dot(g_norm, st.ψ_βτ_pre)


    # cutoff_β = r * sqrt(st.var_β1_pre * cutoff_factor)
    # cutoff_τ = r * sqrt(st.var_τ1_pre * cutoff_factor)
    # cutoff_βτ = r * sqrt(st.var_βτ_pre * cutoff_factor)
    # if mode == :bed
    #     cnts = begin
    #         r = zeros(Int, 4)
    #         for v in g
    #             try
    #                 r[convert(Int, round(v)) + 1] += 1
    #             catch e
    #                 r[4] += 1
    #             end
    #         end
    #         r
    #     end
    #     #@assert sum(cnts) == length(g) "With genotypes == true, the values in g must be 0, 1, or 2."
    #     vals_norm = ([0.0, 1.0, 2.0, m] .- m) ./ s
    #     p_β = _get_pval(s_β, cutoff_β, p_alt[1], vals_norm,
    #         Ks.K0_β, Ks.K1_β, Ks.K2_β, tmp_ecgf; cnts=cnts)
    #     p_τ = _get_pval(s_τ, cutoff_τ, p_alt[2], vals_norm,
    #         Ks.K0_τ, Ks.K1_τ, Ks.K2_τ, tmp_ecgf; cnts=cnts)
    #     p_βτ = _get_pval(s_βτ, cutoff_βτ, p_alt[3], vals_norm,
    #         Ks.K0_βτ, Ks.K1_βτ, Ks.K2_βτ, tmp_ecgf; cnts=cnts)
    # elseif mode == :ukbbgen
    #     cnts = begin
    #         r = zeros(Int, 512)
    #         for v in g
    #             try
    #                 r[convert(Int, round(v*255)) + 1] += 1
    #             catch e
    #                 r[512] += 1
    #             end
    #         end
    #         r
    #     end
    #     vals_norm = vcat((bgenlookup .- m) ./ s, [0])
    #     p_β = _get_pval(s_β, cutoff_β, p_alt[1], vals_norm,
    #         Ks.K0_β, Ks.K1_β, Ks.K2_β, tmp_ecgf; cnts=cnts)
    #     p_τ = _get_pval(s_τ, cutoff_τ, p_alt[2], vals_norm,
    #         Ks.K0_τ, Ks.K1_τ, Ks.K2_τ, tmp_ecgf; cnts=cnts)
    #     p_βτ = _get_pval(s_βτ, cutoff_βτ, p_alt[3], vals_norm,
    #         Ks.K0_βτ, Ks.K1_βτ, Ks.K2_βτ, tmp_ecgf; cnts=cnts)        
    # else
    #     vals_norm = (g .- m) ./ s
    #     p_β = _get_pval(s_β, cutoff_β, p_alt[1], vals_norm,
    #         Ks.K0_β, Ks.K1_β, Ks.K2_β, tmp_ecgf)
    #     p_τ = _get_pval(s_τ, cutoff_τ, p_alt[2], vals_norm,
    #         Ks.K0_τ, Ks.K1_τ, Ks.K2_τ, tmp_ecgf)
    #     p_βτ = _get_pval(s_βτ, cutoff_βτ, p_alt[3], vals_norm,
    #         Ks.K0_βτ, Ks.K1_βτ, Ks.K2_βτ, tmp_ecgf)  
    # end
    p_β, p_τ, p_βτ
end

function spa(g_norm::AbstractVector, pre_vec::AbstractVector,
    vals_norm::AbstractVector,
    cutoff_factor::Real, r::Real, p_alt::Real, K0, K1, K2;
    cnts=nothing, 
    tmp_ecgf = g === nothing ? similar(g) : similar(g, length(cnts)), 
    pre_vec_var::Real=var(pre_vec)
    )
    s = dot(g_norm, pre_vec)
    cutoff = r * sqrt(pre_vec_var * cutoff_factor)
    return _get_pval(s, cutoff, p_alt, vals_norm, K0, K1, K2, tmp_ecgf; cnts=cnts)
end

function _get_pval(s, cutoff, p_alt, vals_norm, K0, K1, K2, tmp_ecgf; cnts=nothing)
    function K0_(z)
        if cnts === nothing
            @inbounds for i in 1:length(vals_norm)
                tmp_ecgf[i] = K0(vals_norm[i] * z)
            end
            return sum(tmp_ecgf)
        else
            @inbounds for i in 1:length(tmp_ecgf)
                if cnts[i] != 0
                    tmp_ecgf[i] = K0(vals_norm[i] * z)
                end
            end
            return dot(cnts, tmp_ecgf)
        end
    end
    function K1_(z)
        if cnts === nothing
            @inbounds for i in 1:length(vals_norm)
                tmp_ecgf[i] = vals_norm[i] * K1(vals_norm[i] * z)
            end
            return sum(tmp_ecgf)
        else
            @inbounds for i in 1:length(tmp_ecgf)
                if cnts[i] != 0
                    tmp_ecgf[i] = vals_norm[i] * K1(vals_norm[i] * z)
                end
            end
            return dot(cnts, tmp_ecgf)
        end
    end
    function K2_(z)
        if cnts === nothing
            @inbounds for i in 1:length(vals_norm)
                tmp_ecgf[i] = vals_norm[i] ^ 2 * K2(vals_norm[i] * z)
            end
            return sum(tmp_ecgf)
        else
            @inbounds for i in 1:length(tmp_ecgf)
                if cnts[i] != 0
                    tmp_ecgf[i] = vals_norm[i] ^ 2 * K2(vals_norm[i] * z)
                end
            end
            return dot(cnts, tmp_ecgf)
        end
    end
    function f(x)
        r = K0_(x[1]) - x[1] * s
        return r
    end
    function g!(storage, x)
        storage[1] = K1_(x[1]) - s
    end
    function h!(storage, x)
        storage[1] = K2_(x[1])
    end
    if abs(s) < cutoff
        return p_alt
    else
        r = optimize(f, g!, h!, [0.0], LBFGS())
        zeta = minimizer(r)[1]
        omega = sign(zeta) * sqrt(max(0, 2 * (zeta * s - K0_(zeta))))
        nu = zeta * sqrt(max(0, K2_(zeta)))
        z2 = omega + 1.0/omega * log(nu / omega)
    
        return ccdf(Normal(), abs(z2)) + cdf(Normal(), -abs(z2))
    end
end
