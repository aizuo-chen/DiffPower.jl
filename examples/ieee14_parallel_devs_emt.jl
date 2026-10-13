using ModelingToolkit
using OrdinaryDiffEq

include(joinpath(@__DIR__, "..", "src", "PowerEMT.jl"))
import .PowerEMT as PE
using .PowerEMT: t

data = PE.ieee14_case()
const G_MIN = 1e-5
const LG_MIN = 1e-10

pf = PE.newton_pf(data.bus, data.branch, data.gen;
	baseMVA=data.base_mva, tol=1e-12, g_min=G_MIN, lg_min=LG_MIN)

reference_name = :vg1_1
dev_specs = Dict(
	:vg1_1 => (bus=1, model=PE.SG6Ref, power=nothing, params=(; ratio=1.5)),
	:vg2_1 => (bus=2, model=PE.Droop, power=0.15+0.05im, params=(; ratio=1.5)),
	:vg2_2 => (bus=2, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
	:vg3 => (bus=3, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
	:vg6 => (bus=6, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
	:vg8 => (bus=8, model=PE.OBS, power=nothing, params=(; ratio=1.5)),
)

PE.align_reference!(pf, dev_specs, reference_name)
components = PE.generate_components(pf, data.branch, dev_specs;
	g_min=G_MIN, lg_min=LG_MIN)
(; devs, nodes, edges, dev_powers, dev_currents) = components

network_eqs = PE.composing(nodes, edges, devs, reference_name)
subsystems = [collect(values(nodes)); collect(values(edges)); collect(values(devs))]
emt = System(network_eqs, t; systems=subsystems, name=:ieee14_parallel_devs_emt)
compiled = mtkcompile(emt)
base_prob = ODEProblem(compiled, [], (0.0, 0.2); build_initializeprob=false)
prob = PE.remake_operating_point(base_prob, devs, pf, dev_specs)

bus2_names = [:vg2_1, :vg2_2]
power_error = sum(dev_powers[name] for name in bus2_names) -
	pf.V_bus[2]*conj(pf.I_gen[2])
current_error = sum(dev_currents[name] for name in bus2_names) - pf.I_gen[2]
initial_residual = maximum(abs, prob.f(prob.u0, prob.p, 0.0))

println("Resolved powers at bus 2: ", Dict(name => dev_powers[name] for name in bus2_names))
println("Bus-2 power allocation error: ", abs(power_error))
println("Bus-2 current allocation error: ", abs(current_error))
println("Maximum initial state derivative: ", initial_residual)

sol = solve(prob, Rodas5P(); reltol=1e-8, abstol=1e-10, saveat=1e-3)
SciMLBase.successful_retcode(sol) || error("parallel-device simulation failed")
println("IEEE-14 parallel-device EMT simulation completed: ", length(sol.t), " saved points")
