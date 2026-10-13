using ModelingToolkit
using OrdinaryDiffEq
import CairoMakie as CM
using OrdinaryDiffEqSDIRK
using LinearAlgebra
import ForwardDiff as FD

include(joinpath(@__DIR__, "..", "src", "DiffPower.jl"))
import .DiffPower as PE
using .DiffPower: t

case_data = PE.ieee14_case()
const BASE_MVA = case_data.base_mva
bus = case_data.bus
gen = case_data.gen
branch = case_data.branch
gen_idx = Int.(gen[:, 1])

const G_MIN = 1e-5
const LG_MIN = 1e-10

# AC power-flow operating point used to initialize the EMT states.
pf = PE.newton_pf(
	bus, branch, gen;
	baseMVA = BASE_MVA,
	tol = 1e-12,
	g_min = G_MIN,
	lg_min = LG_MIN,
)


reference_name = :vg1
dev_specs = Dict(
	:vg1 => (bus=1, model=PE.SG6Ref, power=nothing, params=(; ratio=1.5)),
	:vg2 => (bus=2, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
	:vg3 => (bus=3, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
	:vg6 => (bus=6, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
	:vg8 => (bus=8, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
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
emt = System(network_eqs, t; systems = subsystems, name = :ieee14_droop_emt)
emt_compiled = mtkcompile(emt)


tspan = (0.0, 1.0)
base_prob = ODEProblem(emt_compiled, [], tspan; build_initializeprob = false)
prob = PE.remake_operating_point(
	base_prob, devs, pf, dev_specs,
)

initial_rhs = prob.f(prob.u0, prob.p, first(tspan))
state_names = string.(unknowns(emt_compiled))
largest_initial_terms = sortperm(abs.(initial_rhs); rev = true)[1:min(10, length(initial_rhs))]
println("Largest initial state derivatives:")
for index in largest_initial_terms
	println("  ", rpad(state_names[index], 48), " ", initial_rhs[index])
end

# At 0.1 s, increase the active-power reference of the GFM at bus 2.
condition(u, t, integrator) = t - 0.1 == 0
affect!(integrator) = integrator.ps[devs[:vg2].pr] += 0.05
cb = DiscreteCallback(condition, affect!)

sol = solve(
	prob,
	KenCarp4();
	reltol = 1e-7,
	abstol = 1e-9,
	saveat = 1e-3,
	callback = cb,
	tstops = [0.1],
)

SciMLBase.successful_retcode(sol) || error("EMT solve failed with return code $(sol.retcode)")
println("IEEE-14 SG6Ref + OBS EMT simulation completed: $(length(sol.t)) saved points")
for bus_id in gen_idx
	dev = devs[Symbol("vg", bus_id)]
	frequency = hasproperty(dev, :ωr) ? sol[dev.ωr][end] : sol[dev.ω][end]
	println("bus $bus_id final frequency = $(round(60 * frequency, digits = 6)) Hz")
end

fig = CM.Figure(size = (900, 900))
ax_uc = CM.Axis(
	fig[1, 1];
	ylabel = "|u| (p.u.)",
	title = "Generator terminal/capacitor-voltage magnitudes",
)
ax_pg = CM.Axis(
	fig[2, 1];
	ylabel = "pg (p.u.)",
	title = "Generator active power",
)
ax_qg = CM.Axis(
	fig[3, 1];
	xlabel = "Time (s)",
	ylabel = "qg (p.u.)",
	title = "Generator reactive power",
)

for bus_id in gen_idx
	dev = devs[Symbol("vg", bus_id)]
	voltage_magnitude = if bus_id == pf.slack_bus
		hypot.(sol[dev.Vd], sol[dev.Vq])
	else
		hypot.(sol[dev.ucd], sol[dev.ucq])
	end
	label = "bus $bus_id"
	CM.lines!(ax_uc, sol.t, voltage_magnitude; label = label)
	CM.lines!(ax_pg, sol.t, sol[dev.pg]; label = label)
	CM.lines!(ax_qg, sol.t, sol[dev.qg]; label = label)
end

for ax in (ax_uc, ax_pg, ax_qg)
	CM.vlines!(ax, [0.1]; color = :black, linestyle = :dash, linewidth = 1)
end
CM.linkxaxes!(ax_uc, ax_pg, ax_qg)
CM.axislegend(ax_uc; position = :rb)
fig
