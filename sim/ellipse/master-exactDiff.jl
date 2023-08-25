include("../../startup.jl")

γ = Curve(t -> [cos(π * t), 0.5 * sin(π * t)], t -> [-π * sin(π * t), 0.5 * π * cos(π * t)], t -> [-π^2 * cos(π * t), -0.5π^2 * sin(π * t)])
g = Profile(x -> x, x -> 1)
#g = Profile(x -> sin(1.5x), x -> 1.5cos(1.5x))
f = Sleeve(γ, g)


fRec = sleeveApproximation(f, 0.5, 1e-3, 1e-4)

display(fRec.γ)

x = [fRec.γ[k][1] for k=1:length(fRec.γ)]
y = [fRec.γ[k][2] for k=1:length(fRec.γ)]

#plotlyjs()
pgfplotsx()

K = 100
tx = LinRange(-1.5, 1.5, K)
ty = LinRange(-0.5, 1., K)
F = [f([tx[l], ty[k]]) for k=1:K, l=1:K]
FRec = [fRec([tx[l], ty[k]]) for k=1:K, l=1:K]

display(heatmap(tx, ty, abs.(F - FRec), c=:viridis, tex_output_standalone=true))

savefig("sim/ellipse/error-exactDiff.tex")