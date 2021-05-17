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
end

function ecgf(st::WSVarScoreTestInvariant)
    K0_β, K1_β, K2_β = ecgf(st.ψ_β1_pre)
    K0_τ, K1_τ, K2_τ = ecgf(st.ψ_τ1_pre)
    vGWASEcgfCollection(K0_β, K1_β, K2_β, K0_τ, K1_τ, K2_τ)
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
    dir_alt,
    Ks::vGWASEcgfCollection;
    g_norm=similar(g),
    cnts=nothing,
    ref_vals= cnts !== nothing ? similar(g, length(cnts)) : nothing,
    vals_norm= cnts !== nothing ? similar(g, length(cnts)) : g_norm,
    tmp_ecgf = cnts !== nothing ? similar(g, length(cnts)) : similar(g),
    r=2.0,
    adjustor::Union{Adjustor, Nothing}=nothing,
    adj_cutoff::Real=5e-5
    )
    # involves internal normalization
    m = mean(skipmissing(g))
    #s = std(skipmissing(g))

    # screening for β
    g_norm .= g .- m
    @inbounds for i in 1:length(g_norm)
        g_norm[i] = isnan(g_norm[i]) ? 0.0 : g_norm[i]
    end
    cutoff_factor = sum(x -> x^2, g_norm)
    std_β = sqrt(st.var_β1_pre * cutoff_factor)
    score_β = dot(g_norm, st.ψ_β1_pre)
    z_β = score_β/std_β

    if abs(z_β) < r
        p_β = p_alt[1]
        dir_β = dir_alt[1]
    else
        g_norm .= (g .- m) ./ std_β
        if ref_vals !== nothing
            vals_norm .= (ref_vals .- m) ./ std_β
            @inbounds for i in 1:length(vals_norm)
                if isnan(vals_norm[i])
                    vals_norm[i] = zero(eltype(vals_norm))
                end
            end
        end
        p_β, dir_β = spa(g_norm, st.ψ_β1_pre, vals_norm,
            cutoff_factor, 0.0, p_alt[1], dir_alt[1], Ks.K0_β, Ks.K1_β, Ks.K2_β;
            cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_β1_pre, 
            adjustor=adjustor, adj_cutoff=adj_cutoff)
    end

    # screening for τ

    g_norm .= g .- m
    @inbounds for i in 1:length(g_norm)
        g_norm[i] = isnan(g_norm[i]) ? 0.0 : g_norm[i]
    end
    # cutoff_factor = sum(x -> x^2, g_norm)
    std_τ = sqrt(st.var_τ1_pre * cutoff_factor)
    score_τ = dot(g_norm, st.ψ_τ1_pre)
    z_τ = score_τ/std_τ

    if abs(z_τ) < r
        p_τ = p_alt[2]
        dir_τ = dir_alt[2]
    else
        g_norm .= (g .- m) ./ std_τ
        if ref_vals !== nothing
            vals_norm .= (ref_vals .- m) ./ std_τ
            @inbounds for i in 1:length(vals_norm)
                if isnan(vals_norm[i])
                    vals_norm[i] = zero(eltype(vals_norm))
                end
            end
        end
        p_τ, dir_τ = spa(g_norm, st.ψ_τ1_pre, vals_norm,
            cutoff_factor, 0.0, p_alt[2], dir_alt[2], Ks.K0_τ, Ks.K1_τ, Ks.K2_τ;
            cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_τ1_pre, 
            adjustor=adjustor, adj_cutoff=adj_cutoff)
    end
    # p_β, dir_β = spa(g_norm, st.ψ_β1_pre, vals_norm,
    #     cutoff_factor, r, p_alt[1], dir_alt[1], Ks.K0_β, Ks.K1_β, Ks.K2_β;
    #     cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_β1_pre, 
    #     adjustor=adjustor, adj_cutoff=adj_cutoff)
    # p_τ, dir_τ = spa(g_norm, st.ψ_τ1_pre, vals_norm,
    #     cutoff_factor, r, p_alt[2], dir_alt[2], Ks.K0_τ, Ks.K1_τ, Ks.K2_τ;
    #     cnts=cnts, tmp_ecgf=tmp_ecgf, pre_vec_var=st.var_τ1_pre,
    #     adjustor=adjustor, adj_cutoff=adj_cutoff)

    p_β, p_τ, dir_β, -dir_τ
end

function spa(g_norm::AbstractVector, pre_vec::AbstractVector,
    vals_norm::AbstractVector,
    cutoff_factor::Real, r::Real, p_alt::Real, dir_alt, K0, K1, K2;
    cnts=nothing, 
    tmp_ecgf = g_norm === nothing ? similar(g_norm) : similar(g_norm, length(cnts)), 
    pre_vec_var::Real=var(pre_vec), 
    adjustor::Union{Adjustor, Nothing}=nothing,
    adj_cutoff::Real=1e-3
    )
    s = dot(g_norm, pre_vec)
    cutoff = r * sqrt(pre_vec_var * cutoff_factor)
    pval, dir = _get_pval(s, cutoff, p_alt, dir_alt, vals_norm, K0, K1, K2, tmp_ecgf; cnts=cnts)
    if pval < adj_cutoff && adjustor !== nothing
        g_norm .= g_norm .- adjustor.X_XtXinv * adjustor.Xt * g_norm
        g_norm /= sqrt(sum(x -> x^2, g_norm) * pre_vec_var) # renormalization
        s = dot(g_norm, pre_vec)
        pval, dir = _get_pval(s, Inf, 0.0, 0, g_norm, K0, K1, K2, tmp_ecgf; cnts=nothing)
    end
    return pval, dir
end

const normal = Normal()
const lbfgs = LBFGS()
const singleton_zero = [0.0]
function _get_pval(s, cutoff, p_alt, dir_alt, vals_norm, K0, K1, K2, tmp_ecgf; cnts=nothing)
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
        return p_alt, dir_alt
    else
        r = optimize(f, g!, h!, singleton_zero, lbfgs)
        zeta = minimizer(r)[1]
        omega = sign(zeta) * sqrt(max(0.0, 2 * (zeta * s - K0_(zeta))))
        nu = zeta * sqrt(max(0, K2_(zeta)))
        z2 = omega + 1.0/omega * log(nu / omega)
    
        return ccdf(normal, abs(z2)) + cdf(normal, -abs(z2)), isnan(z2) ? 0 : convert(Int, sign(z2))
    end
end
