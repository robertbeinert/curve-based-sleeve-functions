include("../../startup.jl")

γ = Curve(t -> [cos(π * t), 0.5 * sin(π * t)], t -> [-π * sin(π * t), 0.5 * π * cos(π * t)], t -> [-π^2 * cos(π * t), -0.5π^2 * sin(π * t)])
g = Profile(x -> x, x -> 1)
#g = Profile(x -> sin(1.5x), x -> 1.5cos(1.5x))
f = Sleeve(γ, g; h=1e-8)


fRec = sleeveApproximation(f, 0.5, 1e-3, 1e-4)

display(fRec.γ)

x = [fRec.γ[k][1] for k=1:length(fRec.γ)]
y = [fRec.γ[k][2] for k=1:length(fRec.γ)]

#plotlyjs()
pgfplotsx()

P1 = plot(x,y, xlim=(-1.05,1.05), ylim=(-0.05,1.05), label="", tex_output_standalone=true)
P2 = plot(fRec.g2.knots.^2, fRec.g2.samples, label="", tex_output_standalone=true)

K = 100
tx = LinRange(-1.5, 1.5, K)
ty = LinRange(-0.5, 1., K)
F = [f([tx[l], ty[k]]) for k=1:K, l=1:K]
FRec = [fRec([tx[l], ty[k]]) for k=1:K, l=1:K]

P3 = heatmap(tx, ty, FRec, c=:viridis, tex_output_standalone=true)
P4 = heatmap(tx, ty, abs.(F - FRec), c=:viridis, tex_output_standalone=true)

display(P1)
display(P2)
display(P3)
display(P4)

savefig(P1, "sim/ellipse/curve-finiteDiff.tex")
savefig(P2, "sim/ellipse/profile-finiteDiff.tex")
savefig(P3, "sim/ellipse/sleeve-finiteDiff.tex")
savefig(P4, "sim/ellipse/error-finiteDiff.tex")
