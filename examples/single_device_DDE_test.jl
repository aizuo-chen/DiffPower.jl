# 从 single_device_test.jl 复制并改成 ODE / DDE 对照。
# 按 ## 分段执行，或在已有 Julia 包环境里 include 本文件。
# 保留原来的单台 VSG、LCL 参数和固定电网输入；只延时最终桥电压。

## 1. 加载包和本地 DiffPower
using ModelingToolkit, OrdinaryDiffEq, DelayDiffEq, LinearAlgebra
using OrdinaryDiffEqSDIRK: KenCarp4
import CairoMakie as CM
import ForwardDiff as FD
using SciMLBase: SciMLBase
import SymbolicIndexingInterface as SII
using ModelingToolkit: t_nounits as t
include(joinpath(@__DIR__, "..", "src", "DiffPower.jl"))
using .DiffPower

## 2. 与原 single_device_test 相同的模型
@named sg = DiffPower.VSG(params = (; ratio = 5, lg = 0.2))
input_eqs = [sg.ud ~ 0.0, sg.uq ~ 1.0, sg.ωr ~ 1.0]
@named model = System(input_eqs, t; systems = [sg])
model_c = mtkcompile(model)
tspan = (0.0, 10.0)
prob = ODEProblem(model_c, [], tspan)
@assert prob.f.mass_matrix == I
x0 = copy(prob.u0)

Td = 200e-6                # 最终桥电压指令延时，秒
# 原 VSG 参数在 200 μs 下会出现振荡增长；保留这个结果供对照。
# 可自行改 Td 后重跑后续步骤，但延时大小不保证闭环稳定。
dtmax = 20e-6              # 积分最大步长，不是采样周期
event_time = 0.5
power_step = 0.02          # pr 阶跃；设为 0.0 就是原脚本的无扰动工况
@assert 0 < dtmax <= Td && tspan[1] < event_time < tspan[2]

## 3. 分开保存当前 / 过去可能用到的参数
p_before = copy(prob.p)
p_after = copy(p_before)
pr0 = SII.getp(model_c, sg.pr)(p_before)
SII.setp(model_c, sg.pr)(p_after, pr0 + power_step)

# 不用 callback 修改同一份 p，避免查过去时读到跳变后的 pr。
model_parameters(parameters, time) = time < parameters.event_time ?
									 parameters.p_before : parameters.p_after

# vds / vqs 是限幅后、变换到固定网络 dq 坐标的最终桥电压。
command_observed = SII.observed(model_c, [sg.vds, sg.vqs])
sample_observed = SII.observed(model_c, [sg.pg, sg.ω, sg.ucd, sg.ucq])
i1d_index = SII.variable_index(model_c, sg.i1d)
i1q_index = SII.variable_index(model_c, sg.i1q)
omega_b, lf = SII.getp(model_c, [sg.ωb, sg.lf])(p_before)
Bv = zeros(length(x0), 2)
Bv[i1d_index, 1] = omega_b / lf
Bv[i1q_index, 2] = omega_b / lf
parameters = (; p_before, p_after, event_time, Td, x0, Bv, command_observed)
@show unknowns(model_c) pr0 Td Bv

## 4. 原 ODE，以及只多一个 history 查询的 DDE
function ode_rhs!(dx, x, parameters, time)
	p_model = model_parameters(parameters, time)
	# 兼容 MTK 生成的原地 / 非原地数值函数。
	if SciMLBase.isinplace(prob)
		prob.f(dx, x, p_model, time)
	else
		dx .= prob.f(x, p_model, time)
	end
	nothing
end


function dde_rhs!(dx, x, history, parameters, time)
	ode_rhs!(dx, x, parameters, time)
	past_time = time - parameters.Td
	past_state = history(parameters, past_time)
	p_now = model_parameters(parameters, time)
	p_past = model_parameters(parameters, past_time)
	v_now = parameters.command_observed(x, p_now, time)
	v_past = parameters.command_observed(past_state, p_past, past_time)
	dx .+= parameters.Bv * (v_past - v_now)
	nothing
end

# 原 ODE 已含 Bv*v_now，所以减当前电压、加过去电压。
# 只改 i1d / i1q 两个导数，VSG 和电容、网侧电感仍使用当前状态。
# 代码直接延时固定标称 dq 指令，没有单独执行相位补偿算法。
# 若对应 abc 硬件延时，这等价于假设其标称相位滞后已由接口补偿。
# 未补偿的 abc 延时，还需把过去的复电压乘 exp(-im*omega_b*Td)。
# 用户给 t<=0 的历史，求解器负责 t>0 的历史保存和插值。
initial_history(parameters, time) = parameters.x0

res = similar(x0)
ode_rhs!(res, x0, parameters, 0.0)
@assert norm(res, Inf) < 1e-7   # 常值历史必须对应原模型稳态
@show norm(res, Inf)
# 零延时检查：DDE 必须完全退化为相同的 ODE。
probe_x = x0 .+ 1e-4
dx_ode = similar(x0);
dx_dde = similar(x0)
ode_rhs!(dx_ode, probe_x, parameters, event_time + 0.1)
zero_delay_parameters = merge(parameters, (; Td = 0.0))
dde_rhs!(dx_dde, probe_x, (p, time)->probe_x,
	zero_delay_parameters, event_time + 0.1)
@assert dx_ode == dx_dde

## 5. 同一个 x0 和电路，分别运行无延时 / 真实延时
prob_ode = ODEProblem(ode_rhs!, copy(x0), tspan, parameters)
prob_dde = DDEProblem(dde_rhs!, copy(x0), initial_history, tspan,
	parameters; constant_lags = (Td,))
sol_ode = solve(prob_ode, KenCarp4(); dtmax, reltol = 1e-9, abstol = 1e-10,
	saveat = 0.001, tstops = [event_time])
sol_dde = solve(prob_dde, MethodOfSteps(KenCarp4()); dtmax,
	reltol = 1e-9, abstol = 1e-10, saveat = 0.001,
	tstops = [event_time, event_time + Td])
@assert SciMLBase.successful_retcode(sol_ode) && SciMLBase.successful_retcode(sol_dde)
@assert sol_ode.t[end] == tspan[2] && sol_dde.t[end] == tspan[2]
@assert all(x->all(isfinite, x), sol_dde.u)
pre_step_drift = maximum(norm(x-x0, Inf) for (x, time) in zip(sol_dde.u, sol_dde.t)
											 if time < event_time)
if pre_step_drift >= 1e-7
	@warn "DDE 在阶跃前已偏离稳态；求解成功不表示延时闭环稳定" pre_step_drift Td
end
@show sol_ode.retcode sol_dde.retcode pre_step_drift

## 6. 保留原脚本的无延时 ODE 线性化；这些不是 DDE 特征根
f(x) = prob.f(x, p_before, 0.0)
A = FD.jacobian(f, x0)
ode_eigenvalues = eigen(A).values
@show ode_eigenvalues

## 7. 看两种仿真的功率和频率；变量均保留在 REPL
y_ode = [sample_observed(x, model_parameters(parameters, time), time)
		 for (x, time) in zip(sol_ode.u, sol_ode.t)]
y_dde = [sample_observed(x, model_parameters(parameters, time), time)
		 for (x, time) in zip(sol_dde.u, sol_dde.t)]
fig = CM.Figure(size = (950, 400))
axp = CM.Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "p (p.u.)")
axf = CM.Axis(fig[1, 2]; xlabel = "t (s)", ylabel = "f (Hz)")
for (label, solution, ys, style) in (("ODE", sol_ode, y_ode, :dash),
	("DDE", sol_dde, y_dde, :solid))
	CM.lines!(axp, solution.t, [y[1] for y in ys]; label, linestyle = style)
	CM.lines!(axf, solution.t, [omega_b/(2pi)*y[2] for y in ys]; label, linestyle = style)
end
CM.axislegend(axp)
CM.vlines!(axp, [event_time]; color = :gray, linestyle = :dot)
CM.vlines!(axf, [event_time]; color = :gray, linestyle = :dot)
output_dir = joinpath(@__DIR__, "single_device_DDE_results")
mkpath(output_dir)
CM.save(joinpath(output_dir, "ODE_DDE.png"), fig)
open(joinpath(output_dir, "checks.txt"), "w") do io
	println(io, "Original single-device VSG parameters: ratio=5, lg=0.2.")
	println(io, "Td_s=", Td, "; dtmax_s=", dtmax, "; event_s=", event_time, "; power_step=", power_step)
	println(io, "states=", length(x0), "; initial RHS=", norm(res, Inf), "; zero-delay RHS error=", norm(dx_ode-dx_dde, Inf))
	println(io, "ODE=", sol_ode.retcode, "; DDE=", sol_dde.retcode, "; final time=", sol_dde.t[end])
	println(io, "DDE pre-step drift=", pre_step_drift)
	println(io, "ODE eigenvalues above are not delayed-system characteristic roots.")
end
display(fig)

# 逐行查看：
# command_observed(x0, p_before, 0.0)
# command_observed(sol_dde(1.0-Td), p_after, 1.0-Td) # 1 s 时施加的过去指令
# sol_ode(1.0)-sol_dde(1.0)
# 在步骤 5 构造 prob_dde 后，也可以自己推进一步：
# integrator = init(prob_dde,MethodOfSteps(KenCarp4());dtmax,
#     reltol=1e-9,abstol=1e-10,tstops=[event_time,event_time+Td])
# step!(integrator)
# integrator.t
# integrator.u
