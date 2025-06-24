using BinnedDistributionFit
using Test, FHist
using TestItems

@testitem "ExtendPdf" begin
    d = ExtendPdf((x, _)->x, (1,3))
    @test BinnedDistributionFit.scalar_eval(d, 3) == 3
    @test BinnedDistributionFit.vector_eval(d, [3, 5]) == [3, 5]

    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    @test BinnedDistributionFit.scalar_eval(d2, 3) == 9
    @test BinnedDistributionFit.vector_eval(d2, [3, 5]) == [9, 25]

    d3 = ExtendPdf((x, _) -> xor((x << 2), x) % (x^2 + 1), ())
    @test BinnedDistributionFit.scalar_eval(d3, 3) == 5
    @test BinnedDistributionFit.vector_eval(d3, [3, 5]) == [5, 17]
    end

@testitem "SumOfPdfs" begin
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    d3 = d1+d2
    @test BinnedDistributionFit.scalar_eval(d3, 3) == 12
    @test BinnedDistributionFit.vector_eval(d3, [3, 5]) == [12, 30]
    end

@testitem "RooFitNLL ExtendPdf" begin
    d = ExtendPdf((x, _)->x, (1,3))
    h = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF = RooFitNLL_functor(d, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF([2.0]) ≈ 1.6827899396467232
end


@testset "RooFit SumOfPdfs" begin
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->x^2, (1,3))
    sd = d1+d2
    h = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF = RooFitNLL_functor(sd, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF([[2.0, 0.5], [], []]) ≈ 0.8497340452248006
end
