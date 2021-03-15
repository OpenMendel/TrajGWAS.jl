using vGWAS
using Test
using CSV, DataFrames, WiSER

filepath = normpath(joinpath(dirname(pathof(WiSER)), "../data/"))
df = DataFrame(CSV.File(filepath * "sbp.csv"), copycols = false)
vlmm = WSVarLmmModel(
    @formula(sbp ~ 1 + agegroup + bmi_std + meds), # removed obswt, gender
    @formula(sbp ~ 1 + bmi_std),
    @formula(sbp ~ 1 + agegroup + bmi_std), # removed obswt, gender, meds
    :id, df
);
WiSER.fit!(vlmm)
df_cut = df[!, [:obswt, :gender, :meds]]
df_cut[!, :gender_bool] = map(x -> x == "Male", df_cut[!, :gender])
df_cut[!, :meds_bool] = map(x -> x == "OnMeds", df_cut[!, :meds])
df_cut = df_cut[!, setdiff(names(df_cut), ["meds", "gender"])];
X1vec = Matrix{Float64}[]
W1vec = Matrix{Float64}[]
for i in unique(df[!, :id])
    push!(X1vec, convert(Matrix{Float64}, df_cut[df[!, :id] .== i, :]))
    push!(W1vec, convert(Matrix{Float64}, df_cut[df[!, :id] .== i, :]))
end

@testset "generic" begin
st = vGWAS.WSVarScoreTest(vlmm, 3, 3)
_, v1, _ = vGWAS.test!(st, X1vec, W1vec)
@test isapprox(v1, 6.236443105947029e-11)

st = vGWAS.WSVarScoreTest(vlmm, 1, 1);

# testing for obswt
X1vec_obswt = map(x -> x[:, 1], X1vec)
v1, v2, v3 = vGWAS.test!(st, X1vec_obswt, X1vec_obswt)
@test all(isapprox.((v1, v2, v3),
    (0.4409730648525017, 0.7156750808059059, 0.6854385016534825)))

# testing for gender
X1vec_gender = map(x -> x[:, 2], X1vec)
v1, v2, v3 = vGWAS.test!(st, X1vec_gender, X1vec_gender)
@test isapprox(v2, 0.8535798925691855)

# testing for meds. meds is already included in the null model.
st = vGWAS.WSVarScoreTest(vlmm, 0, 1)
X1vec_meds = map(x -> x[:, 3], X1vec)
v1, v2, v3 = vGWAS.test!(st, nothing, X1vec_meds)
@test all(isapprox.((v1, v2, v3),
    (-1.0, 1.3935485985991808e-12, -1.0)))

end

@testset "generic_invariant" begin
for i in 1:length(X1vec)
    for j in 1:size(X1vec[i], 1)
        @assert all(X1vec[i][j, :] .== X1vec[i][1, :]) # in fact, we are dealing with time-invariant features
    end
end
X1 = Matrix{Float64}(undef, length(X1vec), 3)
for i in 1:length(X1vec)
    X1[i, :] .= X1vec[i][1, :]
end
W1 = copy(X1);

st_i = vGWAS.WSVarScoreTestInvariant(vlmm, 3, 3)
_, v1, _ = vGWAS.test!(st_i, X1, W1)
@test isapprox(v1, 6.236443105947029e-11)

st_i = vGWAS.WSVarScoreTestInvariant(vlmm, 1, 1);
# obswt
v1, v2, v3 = vGWAS.test!(st_i, X1[:, 1], W1[:, 1])
@test all(isapprox.((v1, v2, v3),
    (0.4409730648525017, 0.7156750808059059, 0.6854385016534825)))

# testing for gender
v1, v2, v3 = vGWAS.test!(st_i, X1[:, 2], W1[:, 2])
@test isapprox(v2, 0.8535798925691855)

# testing for meds. meds is already included in the null model.
# meds
st_i = vGWAS.WSVarScoreTestInvariant(vlmm, 0, 1);
v1, v2, v3 = vGWAS.test!(st_i, nothing, W1[:, 3])
@test all(isapprox.((v1, v2, v3),
    (-1.0, 1.3935485985991808e-12, -1.0)))
end

filepath = joinpath(dirname(pathof(vGWAS)), "../data/")

pvalpath = "testpval"
@testset "vgwas_singlesnp_plink" begin
vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "vgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        usespa = false)
results = CSV.read(pvalpath, DataFrame)
@test all(isapprox.((mean(results.betapval),
    mean(results.taupval),
    mean(results.jointpval)),
    (0.20892290777343045, 0.1905578803635017, 0.15579304125062576)))

#spa 
vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "vgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        usespa = true)
results = CSV.read(pvalpath, DataFrame)

@test all(isapprox.((mean(results.betapval),
    mean(results.taupval),
    mean(results.jointpval)),
    (0.20841900621550116, 0.18913803849088423, 0.15402200977008706); rtol=2e-2))
end

@testset "vgwas_singlesnp_vcf" begin
    vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "vgwas_vcfex.csv",
        filepath * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS,
        usespa = false)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval),
        mean(results.jointpval)),
        (0.47580382024875084, 0.4752535666464088, 0.48936769742780345); rtol=2e-2))

    #spa version
    vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "vgwas_vcfex.csv",
        filepath * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS,
        usespa = true)
    results = CSV.read(pvalpath, DataFrame)
    println(mean(results.betapval))
    println(mean(results.taupval))
    println(mean(results.jointpval))
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval),
        mean(results.jointpval)),
        (0.4709957279235137, 0.4691934843663542, 0.4766095197667451); rtol=7e-2))
end

@testset "vgwas_singlesnp_bgen" begin
    vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        filepath * "vgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        usespa = false,
        geneticformat = "BGEN", 
        pvalfile = pvalpath)
    results = CSV.read(pvalpath, DataFrame)

    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval),
        mean(results.jointpval)),
        (0.46429181963000665, 0.4931759050587741, 0.5046081797409608); rtol=2e-2))

    #spa version
    vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        filepath * "vgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        usespa = true,
        geneticformat = "BGEN", 
        pvalfile = pvalpath)
    results = CSV.read(pvalpath, DataFrame)
    println(mean(results.betapval))
    println(mean(results.taupval))
    println(mean(results.jointpval))
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval),
        mean(results.jointpval)),
        (0.46428062017147764, 0.49230409589144597, 0.5024177397539207); rtol=2e-2))
end

@testset "vgwas_snpset_plink" begin
vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex),
        :id,
        filepath * "vgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        analysistype = "snpset",
        snpset = 2)
results = CSV.read(pvalpath, DataFrame)
@test all(isapprox.((mean(results.betapval),
    mean(results.taupval),
    mean(results.jointpval)),
    (0.10426027515117509, 0.08863103713602266, 0.059319604727780056); rtol=2e-2))
end

@testset "vgwas_snpset_vcf" begin
vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "vgwas_vcfex.csv",
        filepath * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS,
        analysistype = "snpset",
        snpset = 20)
results = CSV.read(pvalpath, DataFrame)
@test all(isapprox.((mean(results.betapval),
    mean(results.taupval),
    mean(results.jointpval)),
    (0.25347423101832717, 0.19024751040999272, 0.22966140947913333); rtol=2e-2))
end

@testset "vgwas_gxe_plink" begin
vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "vgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        analysistype = "gxe",
        e = :sex,
        snpinds = 1:10)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval),
        mean(results.jointpval)),
        (0.5924867510209822, 0.40587655481713425, 0.5209862444588567); rtol=2e-2))
end

@testset "vgwas_gxe_bgen" begin
    vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        filepath * "vgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        analysistype = "gxe",
        e = :sex,
        snpinds = 1:3)
    results = CSV.read(pvalpath, DataFrame)
    println(mean(results.betapval))
    println(mean(results.taupval))
    println(mean(results.jointpval))
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval),
        mean(results.jointpval)),
        (0.5161655222831536, 0.4049021864084856, 0.5439672480147976); rtol=2e-2))
    @test all(isapprox.((mean(results.snpeffectnullbeta),
        mean(results.snpeffectnulltau)),
        (0.6076516729747372, -0.3486224003739454); rtol=2e-2))
end



rm(pvalpath)
rm("vgwas.null.txt")
