function InfBus end

const INF_BUS_DEFAULTS = (
	ωb = 2pi*60,
	rg = 1e-3,
	lg = 0.0,
	id_m = 0.0,
	iq_m = 0.1,
	ud_m = 0.0,
	uq_m = 1.0,
)

function validate_infbus_parameters(p)
	isfinite(p.ωb) && p.ωb > 0 || throw(ArgumentError("InfBus ωb must be positive"))
	isfinite(p.rg) && p.rg >= 0 || throw(ArgumentError("InfBus rg must be nonnegative"))
	isfinite(p.lg) && p.lg >= 0 || throw(ArgumentError("InfBus lg must be nonnegative"))
	p.rg > 0 || p.lg > 0 || throw(ArgumentError("InfBus requires rg > 0 or lg > 0"))
	return nothing
end

# Current is positive from the internal source to its terminal. At nominal
# frequency, ug = u_terminal + (rg + j*lg)*i in both construction and initialization.
function operating_point(::typeof(InfBus), voltage, current, params::NamedTuple = (;))
	p = merge(INF_BUS_DEFAULTS, params)
	validate_infbus_parameters(p)
	ug = voltage + (p.rg + 1im*p.lg)*current
	states = iszero(p.lg) ? (;) : (;id=real(current), iq=imag(current))
	parameters = (;ω=1.0, ωb=p.ωb, rg=p.rg, lg=p.lg,
		ugd=real(ug), ugq=imag(ug))
	return (;states, parameters)
end

"""
    InfBus(; name, params=(;))

Infinite bus behind a series impedance. `lg=0` (default) retains the algebraic
resistance interface. `lg>0` adds the two RL current differential states.
`rg` and `lg` are per-unit resistance/reactance; `ωb` converts to physical time.
Changing between zero/nonzero `lg` changes the equations, so rebuild the model.
"""
function InfBus(; name, params::NamedTuple = (;))
	p = merge(INF_BUS_DEFAULTS, params)
	validate_infbus_parameters(p)
	(; ωb, rg, lg, id_m, iq_m, ud_m, uq_m) = p
	ug = (id_m+1im*iq_m) * (rg+1im*lg)+ud_m+1im*uq_m
	ps = @parameters begin
		ω = 1.0
		ωb = ωb
		rg = rg
		lg = lg
		ugd = real(ug)
		ugq = imag(ug)
	end
	# 定义变量
	sts = @variables begin
		ωr(t)
		id(t) = id_m
		iq(t) = iq_m
		ud(t)
		uq(t)
	end
	# lg=0 不生成电流动态，也不人为夹到一个 lg_min。
	eq = if iszero(p.lg)
		[id ~ (ugd - ud) / rg, iq ~ (ugq - uq) / rg]
	else
		[
			D(id) / ωb ~ (ugd - ud - rg*id) / lg + ωr*iq,
			D(iq) / ωb ~ (ugq - uq - rg*iq) / lg - ωr*id,
		]
	end
	# 创建 System
	return System(eq, t, sts, ps; name = name)
end

jacobian_state_names(::typeof(InfBus)) = (:id, :iq)
