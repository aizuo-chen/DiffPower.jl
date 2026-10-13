function InfBusRef(; name, params::NamedTuple = (;))
	defaults = (
		rg = 1e-3,
		id_m = 0.0,
		iq_m = 0.1,
		ud_m = 0.0,
		uq_m = 1.0,
	);
	p = merge(defaults, params);
	# 解构初始条件向量
	(; rg, id_m, iq_m, ud_m, uq_m) = p
	ug = (id_m+1im*iq_m) * rg+ud_m+1im*uq_m;
	ps = @parameters begin
		ω = 1.0
		rg = rg
		ugd = real(ug)
		ugq = imag(ug)
	end
	# 定义变量
	sts = @variables begin
		ωr(t)
		id(t)
		iq(t)
		ud(t)
		uq(t)
	end
	# 定义方程
	eq = [
		id ~ (ugd - ud) / rg,
		iq ~ (ugq - uq) / rg,
	]
	# 创建 ODESystem
	return System(eq, t, sts, ps; name = name)
end
