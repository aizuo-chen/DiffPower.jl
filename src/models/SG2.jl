function SG_2(; name, params::NamedTuple = (;))
	defaults = (
		ωb = 2pi*60,
		lg = 0.05,
		rg = 0.01,
		H = 500,
		Da = 50,
		ra = 0.0043,
		ld′ = 0.46,
		lq′ = 0.46,
		id_m = 0,
		iq_m = 0.1,
		ud_m = 0,
		uq_m = 1,
	);
	p = merge(defaults, params);
	# 解构初始条件向量
	(; ωb, lg, rg, H, Da, ra, ld′, lq′, id_m, iq_m, ud_m, uq_m) = p
	## 迭代初始值  
	v_m = (id_m+1im*iq_m) * (rg+1im*lg)+ud_m+1im*uq_m
	Ed′_m = real(v_m) - lq′ * iq_m + ra*id_m;
	Eq′_m = imag(v_m) + ld′ * id_m + ra*iq_m;
	pr_m = real(v_m*(id_m-1im*iq_m));
	ps = @parameters begin
		ωb = ωb
		ld′ = ld′
		lq′ = lq′
		ra = ra
		H = H
		Da = Da
		Ed′ = Ed′_m
		Eq′ = Eq′_m
		pr = pr_m
		lg = lg
		rg = rg
	end
	# 定义变量
	sts = @variables begin
		δ(t) = 0
		ω(t) = 1.0
		ωr(t)
		pg(t)
		qg(t)
		id(t) = id_m
		iq(t) = iq_m
		id_c(t)
		iq_c(t)
		vd(t)
		vq(t)
		vd_c(t)
		vq_c(t)
		ud(t)
		uq(t)
	end

	T1 = [cos(δ) sin(δ); -sin(δ) cos(δ)] # from grid to component
	T2 = [cos(δ) -sin(δ); sin(δ) cos(δ)] # from component to grid
	eq1 = [
		# 系统方程
		pg ~ vd_c * id_c + vq_c * iq_c,
		qg ~ -vd_c * iq_c + vq_c * id_c,
		2H*D(ω) ~ (pr - pg)/ω - Da * (ω - 1),
		D(δ) / ωb ~ ω - ωr,
		vd_c ~ Ed′ + ωr * lq′ * iq_c - ra*id_c,
		vq_c ~ Eq′ - ωr * ld′ * id_c - ra*iq_c,
		# 线路方程
		D(id) / ωb ~ -rg * id / lg + (vd - ud) / lg + ωr * iq,
		D(iq) / ωb ~ -rg * iq / lg + (vq - uq) / lg - ωr * id,
	]
	# 坐标变换
	eq2 = [id_c; iq_c] .~ T1 * [id; iq]# from grid to component
	eq3 = [vd; vq] .~ T2 * [vd_c; vq_c]# from component to grid
	eqs = vcat([eq1, eq2, eq3]...)
	# 创建 ODESystem
	return System(eqs, t, sts, ps; name = name)
end


