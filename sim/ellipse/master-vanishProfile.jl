include("../../startup.jl")

γ = Curve(t -> [cos(π * t), 0.5 * sin(π * t)], t -> [-π * sin(π * t), 0.5 * π * cos(π * t)], t -> [-π^2 * cos(π * t), -0.5π^2 * sin(π * t)])
g = Profile(x -> x^2, x -> 2x)
#g = Profile(x -> sin(1.5x), x -> 1.5cos(1.5x))
f = Sleeve(γ, g)


fRec = sleeveApproximation(f, 0.5, 1e-2, 1e-4)

display(fRec.γ)

x = [fRec.γ[k][1] for k=1:length(fRec.γ)]
y = [fRec.γ[k][2] for k=1:length(fRec.γ)]

#plotlyjs()
pgfplotsx()

t = range(0.0, 1.0, length=100)
Γ = γ.(t)
xTrue = [Γ[k][1] for k=1:100]
yTrue = [Γ[k][2] for k=1:100]

P1 = plot(xTrue, yTrue, label="GTRUE")
plot!(P1, x,y,  ls=:solid, label="GREC", tex_output_standalone=true)

P2 = plot(t, t.^2, label="GTRUE")
plot!(P2, fRec.g2.knots.^2, fRec.g2.samples,  ls=:dash, label="GREC", tex_output_standalone=true)
png("sim/ellipse/profile-vanishProfile")

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

savefig(P1, "sim/ellipse/curve-vanishProfile.tex")
savefig(P2, "sim/ellipse/profile-vanishProfile.tex")
savefig(P3, "sim/ellipse/sleeve-vanishProfile.tex")
savefig(P4, "sim/ellipse/error-vanishProfile.tex")