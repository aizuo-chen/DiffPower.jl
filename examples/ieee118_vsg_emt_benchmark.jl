using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using LinearAlgebra
using Printf
using SparseArrays
using LinearSolve

include(joinpath(@__DIR__, "..", "src", "DiffPower.jl"))
import .DiffPower as PE
using .DiffPower: t

function measured(f, label)
	GC.gc()
	result = @timed f()
	@printf("%-28s %10.3f s  %10.1f MiB\n", label, result.time, result.bytes / 2.0^20)
	flush(stdout)
	return result.value
end

case_data = PE.ieee118_case()
bus = case_data.bus
gen = case_data.gen
branch = case_data.branch

const G_MIN = 1e-5
const LG_MIN = 1e-10
const TEND = let option = findfirst(arg -> startswith(arg, "--tend="), ARGS)
	option === nothing ? 0.1 : parse(Float64, split(ARGS[option], "="; limit = 2)[2])
end
const TSPAN = (0.0, TEND)
const STEP_TIME = 0.02
const STEP_SIZE = 0.05
const SPECIALIZATION = "--autospecialize" in ARGS ?
	SciMLBase.AutoSpecialize : SciMLBase.NoSpecialize
TEND > STEP_TIME || error("simulation end time must exceed callback time $STEP_TIME")
const SOLVER_NAME = "--kencarp4" in ARGS ? :KenCarp4 : :TRBDF2

println("IEEE-118 all-VSG DiffPower benchmark")
println("Julia: ", VERSION, ", threads: ", Threads.nthreads())
println("buses=$(size(bus, 1)), generators=$(size(gen, 1)), branches=$(size(branch, 1))")

pf = measured("AC power flow") do
	PE.newton_pf(
		bus,
		branch,
		gen;
		baseMVA = case_data.base_mva,
		tol = 1e-12,
		maxiter = 30,
		g_min = G_MIN,
		lg_min = LG_MIN,
	)
end

reference_bus = pf.slack_bus
reference_name = Symbol("vg", reference_bus)
generator_buses = Int.(gen[:, 1])

dev_specs = Dict(
	Symbol("vg", bus_id) => (
		bus = bus_id,
		model = bus_id == reference_bus ? PE.VSGRef : PE.VSG,
		power = nothing,
		params = (;),
	)
	for bus_id in generator_buses
)

measured("reference-frame alignment") do
	PE.align_reference!(pf, dev_specs, reference_name)
end

components = measured("component generation") do
	PE.generate_components(
		pf,
		branch,
		dev_specs;
		g_min = G_MIN,
		lg_min = LG_MIN,
	)
end

(; devs, nodes, edges) = components

network_eqs = measured("network equation assembly") do
	PE.composing(nodes, edges, devs, reference_name)
end

emt = measured("symbolic System assembly") do
	subsystems = [collect(values(nodes)); collect(values(edges)); collect(values(devs))]
	System(network_eqs, t; systems = subsystems, name = :ieee118_all_vsg)
end

emt_compiled = measured("mtkcompile") do
	mtkcompile(emt; split = true)
end

base_prob = measured("ODEProblem construction") do
	ODEProblem{true}(
		emt_compiled,
		[],
		TSPAN;
		build_initializeprob = false,
	)
end

prob = measured("operating-point remake") do
	PE.remake_operating_point(base_prob, devs, pf, dev_specs)
end

equilibrium_rhs = similar(prob.u0)
prob.f(equilibrium_rhs, prob.u0, prob.p, first(prob.tspan))
equilibrium_residual = norm(equilibrium_rhs, Inf)

pattern = measured("topology Jacobian pattern") do
	PE.topology_jacobian_pattern(emt_compiled, components, dev_specs, reference_name)
end

pattern_check = measured("finite-difference check") do
	PE.check_jacobian_pattern(prob, pattern)
end
pattern_check.valid || error(
	"Jacobian pattern omitted $(length(pattern_check.missing)) significant entries; first entries: $(first(pattern_check.missing, min(10, length(pattern_check.missing))))",
)

target_bus = first(filter(!=(reference_bus), generator_buses))
target_device = devs[Symbol("vg", target_bus)]

jacobian_cache = measured("colored Jacobian prepare") do
	PE.prepare_colored_jacobian(prob, pattern)
end

J0 = copy(jacobian_cache.prototype)
measured("first colored Jacobian") do
	PE.colored_jacobian!(
		J0,
		jacobian_cache,
		prob.u0,
		prob.p,
		first(prob.tspan),
	)
end


measured("10 warm Jacobians") do
	for _ in 1:10
		PE.colored_jacobian!(
			J0,
			jacobian_cache,
			prob.u0,
			prob.p,
			first(prob.tspan),
		)
	end
end

numeric_prob = measured("numerical sparse problem") do
	PE.with_colored_jacobian(prob, jacobian_cache;
		specialization = SPECIALIZATION)
end

initial_rhs = similar(numeric_prob.u0)
measured("first RHS evaluation") do
	numeric_prob.f(initial_rhs, numeric_prob.u0, numeric_prob.p, first(TSPAN))
end
initial_residual = norm(initial_rhs, Inf)
println("initial RHS infinity norm: ", initial_residual)
state_names = string.(unknowns(emt_compiled))
largest_initial_terms = sortperm(abs.(initial_rhs); rev = true)[1:min(10, length(initial_rhs))]
println("largest initial derivatives:")
for index in largest_initial_terms
	println("  ", rpad(state_names[index], 48), " ", initial_rhs[index])
end
flush(stdout)

measured("100 warm RHS evaluations") do
	for _ in 1:100
		numeric_prob.f(initial_rhs, numeric_prob.u0, numeric_prob.p, first(TSPAN))
	end
end

println("dynamic states: ", length(numeric_prob.u0))
println("Jacobian structural nnz: ", nnz(pattern))
println("Jacobian density: ", nnz(pattern) / length(numeric_prob.u0)^2)
println("Jacobian colors: ", jacobian_cache.ncolors)
println("equilibrium RHS infinity norm: ", equilibrium_residual)

eigenpairs = measured("full equilibrium spectrum") do
	eigen(Matrix(J0))
end
lambda = eigenpairs.values
all(isfinite, real.(lambda)) && all(isfinite, imag.(lambda)) ||
	error("equilibrium Jacobian has non-finite eigenvalues")
rightmost = sortperm(real.(lambda); rev = true)
near_zero = count(value -> abs(value) <= 1e-6, lambda)
unstable = count(value -> real(value) > 1e-6, lambda)
println("equilibrium eigenvalue count: ", length(lambda))
println("spectral abscissa (1/s): ", maximum(real, lambda))
println("unstable roots (Re > 1e-6): ", unstable,
	", near-zero roots (abs <= 1e-6): ", near_zero)
println("12 rightmost eigenvalues (1/s):")
for index in first(rightmost, min(12, length(rightmost)))
	@printf("  % .8e %+.8ei\n", real(lambda[index]), imag(lambda[index]))
end
for index in first(rightmost, min(4, length(rightmost)))
	imag(lambda[index]) < 0 && continue
	mode = eigenpairs.vectors[:, index]
	residual = norm(J0*mode-lambda[index]*mode, Inf) /
		(norm(J0, Inf)*norm(mode, Inf))
	@printf("rightmost mode % .5e %+.5ei: relative eigen residual %.3e\n",
		real(lambda[index]), imag(lambda[index]), residual)
	for state_index in first(sortperm(abs.(mode); rev = true), 8)
		@printf("  %-42s %.6f\n", state_names[state_index], abs(mode[state_index]))
	end
end
for bus_index in findall(<(0.0), pf.c_bus)
	println("negative dynamic cg at bus ", bus_index, ": ", pf.c_bus[bus_index])
end
flush(stdout)

if "--stability-only" in ARGS
	println("stability-only diagnostic complete")
	exit()
end

W0 = copy(J0)
W0 .*= -1e-5
for index in axes(W0, 1)
	W0[index, index] += 1.0
end

F0 = measured("first sparse LU") do
	lu(W0)
end
measured("10 warm sparse LU") do
	for _ in 1:10
		lu(W0)
	end
end
println("sparse LU fill nnz: ", nnz(F0.L)+nnz(F0.U),
	", fill ratio: ", (nnz(F0.L)+nnz(F0.U))/nnz(W0))

if "--linear-only" in ARGS
	println("linear-only diagnostic complete")
	exit()
end

function run_simulation(problem = numeric_prob)
	# A callback mutates pr; give every run its own parameter storage so the
	# first and warm solves represent the same physical experiment.
	run_problem = remake(problem; p = deepcopy(problem.p))
	linsolve = KLUFactorization()
	algorithm = SOLVER_NAME === :TRBDF2 ?
		TRBDF2(; linsolve) :
		KenCarp4(; linsolve)
	event_count = Ref(0)
	function affect!(integrator)
		event_count[] += 1
		old_reference = integrator.ps[target_device.pr]
		integrator.ps[target_device.pr] = old_reference+STEP_SIZE
		println("callback at t=", integrator.t,
			": bus ", target_bus, " pr ", old_reference,
			" -> ", integrator.ps[target_device.pr])
		flush(stdout)
	end
	cb = DiscreteCallback((u, t, integrator) -> t == STEP_TIME,
		affect!; save_positions = (true, true))
	sol = solve(
		run_problem,
		algorithm;
		callback = cb,
		tstops = [STEP_TIME],
		saveat = 1e-3,
		reltol = 1e-7,
		abstol = 1e-9,
		dtmin = 1e-5,
		maxiters = 10^7,
	)
	event_count[] == 1 || error("expected one callback event, observed $(event_count[])")
	return sol
end

sol_first = measured("first solve") do
	run_simulation()
end

SciMLBase.successful_retcode(sol_first) || error("first solve failed: $(sol_first.retcode)")
println("first-solve statistics: ", sol_first.stats)
println("first-solve final magnitude: ", maximum(abs, sol_first.u[end]))
flush(stdout)

if "--single" in ARGS
	println("single-solve diagnostic complete at t=", sol_first.t[end])
	exit()
end

sol_warm = measured("warm solve") do
	run_simulation()
end

SciMLBase.successful_retcode(sol_warm) || error("warm solve failed: $(sol_warm.retcode)")

println("initial RHS infinity norm: ", initial_residual)
println("perturbation: bus $target_bus VSG pr += $STEP_SIZE at t=$STEP_TIME s via callback")
println("solver: ", SOLVER_NAME, " with topology sparsity and colored finite differences")
println("specialization: ", SPECIALIZATION)
println("saved points: ", length(sol_warm.t), ", final time: ", sol_warm.t[end])
println("maximum final state magnitude: ", maximum(abs, sol_warm.u[end]))
println("solver statistics: ", sol_warm.stats)

function print_stats(label, sol)
	println(label, " retcode=", sol.retcode,
		", final_t=", sol.t[end],
		", max_state=", maximum(abs, sol.u[end]))
	for name in propertynames(sol.stats)
		println("  ", name, " = ", getproperty(sol.stats, name))
	end
end

if "--scan" in ARGS
	print_stats("scan t=$(TEND)", sol_warm)
	for horizon in (1e-3, 1e-2, 5e-2)
		horizon <= TEND && continue
		scan_prob = remake(numeric_prob; tspan = (0.0, horizon))
		scan_result = measured("warm solve to $(horizon)") do
			run_simulation(scan_prob)
		end
		print_stats("scan t=$(horizon)", scan_result)
		flush(stdout)
	end
end
