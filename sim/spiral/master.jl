include("../../startup.jl")

c(t) = (0.125 + 0.375t) * [cos(3π * t), sin(3π * t)]
cp(t) = 3π * (0.125 + 0.375t) * [-sin(3π * t), cos(3π * t)] +
        0.375 * [cos(3π * t), sin(3π * t)]
cpp(t) = 9π^2 * (0.125 + 0.375t) * [-cos(3π * t), -sin(3π * t)] +
         2.25π * [-sin(3π * t), cos(3π * t)]

γ = Curve(c, cp, cpp)
#g = Profile(x -> x, x -> 1)
g = Profile(x -> sin(0.5π * x), x -> 0.5π * cos(0.5π * x))
f = Sleeve(γ, g)


fRec = sleeveApproximation(f, 0.125, 1e-3, 1e-4)

display(fRec.γ)

x = [fRec.γ[k][1] for k=1:length(fRec.γ)]
y = [fRec.γ[k][2] for k=1:length(fRec.γ)]

# plotlyjs()
# display(plot(x,y, ls=:solid, shape=:auto, ms=2))
# png("sim/spiral/curve2d")
# display(plot(fRec.g2.knots.^2, fRec.g2.samples))
# png("sim/spiral/profile2d")

# K = 100
# tx = LinRange(-1., 1., K)
# ty = LinRange(-1., 1., K)
# F = [norm([tx[k], ty[l]]) ≤ 1 ? f([tx[k], ty[l]]) : NaN for k=1:K, l=1:K]
# FRec = [norm([tx[k], ty[l]]) ≤ 1 ? fRec([tx[k], ty[l]]) : NaN for k=1:K, l=1:K]

# display(heatmap(tx, ty, FRec, c=:viridis))
# png("sim/spiral/sleeve2d")

# display(heatmap(tx, ty, abs.(F - FRec), c=:viridis))
# png("sim/spiral/error2d")


pgfplotsx()
P = plot(x,y, ls=:solid, label="", xlim=(-0.55,0.55), ylim=(-0.55,0.55), tex_output_standalone=true)
display(P)
savefig("sim/spiral/curve2d.tex")

P = plot(fRec.g2.knots.^2, fRec.g2.samples, label="", tex_output_standalone=true)
display(P)
savefig("sim/spiral/profile2d.tex")

K = 100
tx = LinRange(-1., 1., K)
ty = LinRange(-1., 1., K)
F = [norm([tx[k], ty[l]]) ≤ 1 ? f([tx[k], ty[l]]) : NaN for l=1:K, k=1:K]
FRec = [norm([tx[k], ty[l]]) ≤ 1 ? fRec([tx[k], ty[l]]) : NaN for l=1:K, k=1:K]

P = heatmap(tx, ty, FRec, c=:viridis, tex_output_standalone=true)
display(P)
savefig("sim/spiral/sleeve2d.tex")

P = heatmap(tx, ty, abs.(F - FRec), c=:viridis, tex_output_standalone=true)
display(P)
savefig("sim/spiral/error2d.tex")