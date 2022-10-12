using HousingLifeCycleWorkhorse
using Test


@testset "Fundamentals" begin
    par = benchpar()
    🍔 = 2.0
    🏰 = 2.0

    @test util(🍔,🏰,0.0,par) < util(🍔,🏰,1.0,par)

    ## LoM
    par=benchpar()
    pn = 2.0
    sc = .0
    hc = 1.0
    bc = 0.0
    @test LoM(par,bc,sc,hc,pn) == 0.0

    ## LoM
    if par.rm > 0
        @test r(par,1.0) < r(par,-1.0)
    else
        @test r(par,1.0) == r(par,-1.0)
    end

end

@testset "Scaling independence" begin
    par = benchpar(nx=21,ne=1,nh=1)
    V,pol = solve_decproblems(par);
    aggmom, pan = simulate(par,pol,V);

    factor = 0.5
    par2 = rescale_parameters(par,factor)
    V,pol = solve_decproblems(par2);
    aggmom_scale, pan_scale = simulate(par2,pol,V);
    @test isapprox(aggmom.x,aggmom_scale.x/factor,rtol=0.000001)
    @test isapprox(aggmom.c,aggmom_scale.c/factor,rtol=0.000001)
    @test isapprox(aggmom.b,aggmom_scale.b/factor,rtol=0.000001)
end
