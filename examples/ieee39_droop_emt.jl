using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
import ForwardDiff as FD
import CairoMakie as CM

include(joinpath(@__DIR__, "..", "src", "PowerEMT.jl"))
import .PowerEMT as PE
using .PowerEMT: t

case_data = PE.ieee39_case()
const BASE_MVA = case_data.base_mva
bus = case_data.bus
gen = case_data.gen
branch = case_data.branch
gen_idx = Int.(gen[:, 1])

const G_MIN = 1e-5
const LG_MIN = 1e-10

pf = PE.newton_pf(
	bus, branch, gen;
	baseMVA = BASE_MVA,
	tol = 1e-12,
	g_min = G_MIN,
	lg_min = LG_MIN,
)

reference_name = :vg31
dev_specs = Dict(
	Symbol("vg", bus_id) => (
		bus=bus_id,
		model=bus_id == 31 ? PE.SG6Ref : PE.Droop,
		power=nothing,
		params=(; ratio=1.5),
	)
	for bus_id in Int.(gen[:, 1])
)

PE.align_reference!(pf, dev_specs, reference_name)
components = PE.generate_components(
	pf, branch, dev_specs;
	g_min = G_MIN,
	lg_min = LG_MIN,
)
(; devs, nodes, edges) = components

network_eqs = PE.composing(nodes, edges, devs, reference_name)
subsystems = [collect(values(nodes)); collect(values(edges)); collect(values(devs))]
emt = System(network_eqs, t; systems = subsystems, name = :ieee39_droop_emt)
emt_compiled = mtkcompile(emt)

tspan = (0.0, 1.0)
base_prob = ODEProblem(emt_compiled, [], tspan; build_initializeprob = false)
prob = PE.remake_operating_point(
	base_prob, devs, pf, dev_specs,
)
initial_rhs = prob.f(prob.u0, prob.p, first(tspan))
emt_rhs(x) = prob.f(x, prob.p, first(tspan))
A = FD.jacobian(emt_rhs, prob.u0)
λ = eigvals(A)

condition(u, t, integrator) = t - 0.1 == 0
affect(integrator) = integrator.ps[devs[:vg30].v_m] = 1.05
cb = DiscreteCallback(condition, affect)
sol = solve(
	prob,
	Rodas5P();
	reltol = 1e-7,
	abstol = 1e-9,
	saveat = 1e-3,
	callback = cb,
	tstops = [0.1],
)




SciMLBase.successful_retcode(sol) || error("EMT solve failed with return code $(sol.retcode)")
println("IEEE-39 all-droop EMT simulation completed: $(length(sol.t)) saved points")
for bus_id in gen_idx
	dev = devs[Symbol("vg", bus_id)]
	println("bus $bus_id final frequency = $(round(60 * sol[dev.ω][end], digits = 6)) Hz")
end


fig = CM.Figure()
ax = CM.Axis(
	fig[1, 1];
	xlabel = "Time (s)",
	ylabel = "|uc| (p.u.)",
	title = "Generator capacitor-voltage magnitudes",
)
for bus_id in gen_idx
	dev = devs[Symbol("vg", bus_id)]
	uc_magnitude = if hasproperty(dev, :ucd)
		hypot.(sol[dev.ucd], sol[dev.ucq])
	else
		# SG6Ref has no converter capacitor; use its machine-terminal voltage.
		hypot.(sol[dev.vd], sol[dev.vq])
	end
	CM.lines!(ax, sol.t, uc_magnitude; label = "bus $bus_id")
end

ax = CM.Axis(
	fig[2, 1];
	xlabel = "Time (s)",
	ylabel = "qg (p.u.)",
	title = "Generator reactive power",
)
for bus_id in gen_idx
	dev = devs[Symbol("vg", bus_id)]
	CM.lines!(ax, sol.t, sol[dev.qg]; label = "bus $bus_id")
end

ax = CM.Axis(
	fig[3, 1];
	xlabel = "Time (s)",
	ylabel = "pg (p.u.)",
	title = "Generator active power",
)
for bus_id in gen_idx
	dev = devs[Symbol("vg", bus_id)]
	CM.lines!(ax, sol.t, sol[dev.pg]; label = "bus $bus_id")
end
fig
