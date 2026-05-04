using Test
using AutoTortoise

@testset "Schwarzschild tortoise" begin
    f = schwarzschild(M=1.0)
    h = horizons(f)
    I = static_region(h; prefer=:exterior)
    tm = tortoise_map(f, I; convention=:full_line)
    r = 10.0
    @test isapprox(tortoise(tm, r), r + 2log(r - 2), rtol=1e-10, atol=1e-10)
    x = tortoise(tm, 8.0)
    @test isapprox(inverse_tortoise(tm, x), 8.0, rtol=1e-9, atol=1e-9)
end

@testset "Schwarzschild--dS finite static interval" begin
    f = schwarzschild_ds(M=1.0, Lambda=0.08)
    h = horizons(f)
    @test length(h.positive_roots) >= 2
    I = static_region(h; prefer=:finite)
    @test isfinite(I.right)
    @test I.sign == 1
    tm = tortoise_map(f, I)
    r = (I.left + I.right) / 2
    x = tortoise(tm, r)
    @test isapprox(inverse_tortoise(tm, x), r, rtol=1e-9, atol=1e-9)
end

@testset "Schwarzschild--AdS finite boundary convention" begin
    f = schwarzschild_ads(M=1.0, L=10.0)
    h = horizons(f)
    I = static_region(h; prefer=:exterior)
    tm = tortoise_map(f, I; convention=:boundary_zero)
    @test tm.convention == :boundary_zero
    @test abs(tortoise(tm, 1.0e6)) < 2.0e-4
end
