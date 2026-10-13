using JLD2
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using Printf

include(raw"E:\papers\EMT_dev - 6\DiffPower.jl")
import .DiffPower as DP

function measured(f, label)
    GC.gc()
    result = @timed f()
    @printf("%-27s %9.3f s %9.1f MiB\n",
        label, result.time, result.bytes / 2.0^20)
    flush(stdout)
    return result.value
end

data = JLD2.load(raw"E:\papers\EMT_dev - 6\EMT-network\case118.jld2")
pf = measured("EMT_dev6 power flow") do
    DP.solve_power_flow(data["bus"], data["gen"], data["branch"])
end
op_raw = DP.operating_point_from_power_flow(pf)

# This matches the device mix in EMT_dev - 6/Untitled-3.jl.
device_models = (
    DP.GFMVSGRefDevice,
    ntuple(_ -> DP.GFLPLLDevice, 11)...,
    ntuple(_ -> DP.GFMDroopDevice, 9)...,
    ntuple(_ -> DP.GFMVSGDevice, 33)...,
)
source_buses = op_raw.types.source_buses
device_paras = (
    (source_buses[1], (ki1 = 0.5, kp1 = 25.0)),
    (source_buses[6], (ki_pll = 0.05, kp_pll = 0.1)),
    (source_buses[40], (ki1 = 0.5, kp1 = 25.0)),
)
device_params = DP.device_initialize_params(op_raw, device_models;
    device_paras)
op = DP.rotate_operating_point_to_reference_frame(op_raw, device_params)
components = measured("EMT_dev6 components") do
    DP.generate_components(op, device_models, device_params)
end
model, theta = measured("EMT_dev6 build_model") do
    DP.build_model(components.nodes, components.lines, components.devices;
        nonref_device_types = (
            DP.GFMVSGDevice, DP.GFMDroopDevice, DP.GFLPLLDevice))
end
x0 = measured("EMT_dev6 build_x0") do
    DP.build_x0(model, theta, components.init)
end
f! = (dx, x, p, t) -> DP.system_rhs!(dx, x, p, model, t)
dx0 = similar(x0)
f!(dx0, x0, theta, 0.0)
println("states=", length(x0), ", initial RHS infinity norm=", norm(dx0, Inf))

const STEP_TIME = 0.02
const STEP_SIZE = 0.05
target_block = model.device_blocks[28] # Device 29; device 1 is the reference.
pr_index = target_block.pmap.pr

function run_simulation()
    prob = ODEProblem(f!, copy(x0), (0.0, 0.1), copy(theta))
    event_count = Ref(0)
    function affect!(integrator)
        event_count[] += 1
        old_pr = integrator.p[pr_index]
        integrator.p[pr_index] = old_pr + STEP_SIZE
        println("callback at t=", integrator.t, ": pr ", old_pr,
            " -> ", integrator.p[pr_index])
        flush(stdout)
    end
    cb = DiscreteCallback((u, t, integrator) -> t == STEP_TIME,
        affect!; save_positions = (true, true))
    sol = solve(prob, Rodas5P();
        callback = cb, tstops = [STEP_TIME], saveat = 1e-3,
        reltol = 1e-7, abstol = 1e-9, dtmin = 1e-5)
    event_count[] == 1 || error("expected one callback; observed $(event_count[])")
    return sol
end

first_sol = measured("EMT_dev6 first solve") do
    run_simulation()
end
println("first retcode=", first_sol.retcode, ", final t=", first_sol.t[end],
    ", stats=", first_sol.stats)
SciMLBase.successful_retcode(first_sol) || exit(1)

warm_sol = measured("EMT_dev6 warm solve") do
    run_simulation()
end
println("warm retcode=", warm_sol.retcode, ", final t=", warm_sol.t[end],
    ", stats=", warm_sol.stats)
println("repeat endpoint difference=",
    norm(first_sol.u[end] - warm_sol.u[end], Inf))
