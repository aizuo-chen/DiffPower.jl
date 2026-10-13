function Line_fault(; name, params::NamedTuple = (;))
	# default parameter values
	defaults = (rg = 0.015, lg = 0.100, lg_min = 1e-10, k = 1.0, ωb = 2pi*60, id_m = 0.0, iq_m = 0.0)
	p = merge(defaults, params)
	(; rg, lg, k, ωb, id_m, iq_m, tf, tb, tr) = p # 解构参数向量
	lg = max(lg, p.lg_min)
	# 定义参数
	ps = [@parameters begin
			rg = rg
			lg = lg
			k = k
			ωb = ωb
			tf = tf ##故障时间
			tb = tb ##线路切除
			tr = tr ##线路重合闸
		end
		@discretes begin
			rl(t) = 1e9
			S(t) = 1
		end]

	# 定义变量
	sts = @variables begin
		ωr(t)
		isd(t)
		isq(t)
		itd(t)
		itq(t)
		id1(t) = id_m
		iq1(t) = iq_m
		id2(t) = id_m
		iq2(t) = iq_m
		usd(t)
		usq(t)
		utd(t)
		utq(t)
	end
	# 定义方程
	eqs = [
		D(id1) / ωb ~ ((usd / k - rg * id1 - (id1 - id2) * rl) / lg + ωr * iq1)*S,
		D(iq1) / ωb ~ ((usq / k - rg * iq1 - (iq1 - iq2) * rl) / lg - ωr * id1)*S,
		D(id2) / ωb ~ (((id1 - id2) * rl - rg * id2 - utd) / lg + ωr * iq2)*S,
		D(iq2) / ωb ~ (((iq1 - iq2) * rl - rg * iq2 - utq) / lg - ωr * id2)*S,
		isd ~ id1 / k,
		isq ~ iq1 / k,
		itd ~ id2,
		itq ~ iq2,
	]
	# 定义离散事件
	discrete_events = [
		ModelingToolkit.SymbolicDiscreteCallback(
			(t == tf) => [rl ~ 1e-2]; discrete_parameters = rl, iv = t),
		ModelingToolkit.SymbolicDiscreteCallback(
			(t == tb) => [rl ~ 1e9, S~0]; discrete_parameters = [rl, S], iv = t),
		(t == tb) => [id1~0, iq1~0, id2~0, iq2~0],
		ModelingToolkit.SymbolicDiscreteCallback(
			(t == tr) => [S~1]; discrete_parameters = [S], iv = t),
	]
	# 创建 ODESystem
	return System(eqs, t, sts, ps; name = name, discrete_events = discrete_events)
end
