using AutoTortoise

M = 1.0
Lambda = 0.08

f = schwarzschild_ds(M=M, Lambda=Lambda)
h = horizons(f)
I = static_region(h; prefer=:finite)
tm = tortoise_map(f, I)

println("positive roots = ", h.positive_roots)
println("selected interval = ", I)

r = (I.left + I.right) / 2
x = tortoise(tm, r)
r_back = inverse_tortoise(tm, x)

println("r      = ", r)
println("r_*    = ", x)
println("inverse= ", r_back)
