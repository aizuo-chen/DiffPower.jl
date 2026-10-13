using ModelingToolkit
using OrdinaryDiffEq
import CairoMakie as CM
using OrdinaryDiffEqSDIRK
using LinearAlgebra
import ForwardDiff as FD
using ModelingToolkit: t_nounits as t
include(joinpath(@__DIR__, "..", "src", "PowerEMT.jl"))
using .PowerEMT

# ============================================================
# single GFM test
# ============================================================
@named sg = PowerEMT.VSG(params = (; ratio = 5, lg = 0.2))
input_eqs = [
	sg.ud ~ 0.0,
	sg.uq ~ 1.0,
	sg.ωr ~ 1.0,
]
@named model = System(input_eqs, t; systems = [sg])
model_c = mtkcompile(model)
tspan = (0.0, 10.0)
prob = ODEProblem(model_c, [], tspan)
sol = solve(prob, KenCarp4();
	reltol = 1e-9,
	abstol = 1e-10,
	saveat = 0.001,
)

res = prob.f(prob.u0, prob.p, 0.0)  # Evaluate the ODE function at the initial condition
maximum(res)
prob.u0

f(x) = prob.f(x, prob.p, 0.0)

A = FD.jacobian(f, prob.u0)
eigen(A).values
# condition(u, t, integrator) = t - 0.5 == 0

# function affect!(integrator)
# 	integrator.u[2] += 0.02
# end

# cb = DiscreteCallback(condition, affect!)

# sol2 = solve(
# 	prob,
# 	KenCarp4();
# 	reltol = 1e-9,
# 	abstol = 1e-10,
# 	saveat = 0.001,
# 	callback = cb,
# 	tstops = [0.5],
# )

# CM.lines(sol2.t, sol2[sg.ucq], color = :blue)
