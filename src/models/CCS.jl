function CCS(; name, params::NamedTuple = (;))
	defaults = (
		ωb = 2pi*60,
		cf = 0.05,
		lf = 0.12,
		rf = 0.035,
		lg = 0.1,
		rg = 0.02,
		hp = 0.2,
		hc = 11,
		id_m = 0.0,
		iq_m = 0.1,
		ud_m = 0,
		uq_m = 1,
		h_sat = 1.0,
		imax = 5,
	);
	p = merge(defaults, params);# 解构参数向量
	# 解构初始条件向量
	(; ωb, cf, lf, rf, lg, rg, hp, hc, id_m, iq_m, ud_m, uq_m, h_sat, imax) = p
	# 控制器参数
	ki2 = hc*rf; ## ki = ω*rf(二阶系统, 零极点对消除)
	kp2 = hc*lf; ## kp = ω*Lf
	# 迭代初始值
	uc = (id_m+1im*iq_m) * (rg+1im*lg)+ud_m+1im*uq_m
	i1 = id_m+1im*iq_m+1im*cf*uc
	sr = uc*conj(i1)
	δ_m = angle(uc) - π/2
	idr_m, iqr_m = [cos(δ_m) sin(δ_m); -sin(δ_m) cos(δ_m)] * [real(i1); imag(i1)]
	iid_m = -rf*idr_m/ki2
	iiq_m = -rf*iqr_m/ki2
	ps = @parameters begin
		ωb = ωb
		lf = lf
		rf = rf
		cf = cf
		lg = lg
		rg = rg
		pr = real(sr)
		qr = imag(sr)
		hp = hp
		h_sat = h_sat
		imax = imax
		ki2 = ki2
		kp2 = kp2
	end
	# 定义变量
	sts = @variables begin
		δ(t) = δ_m
		ω(t) = 1.0
		ωr(t)
		pg(t)
		qg(t)
		vd_c(t)
		vq_c(t)
		vd(t)
		vq(t)
		i1d(t) = real(i1)
		i1q(t) = imag(i1)
		i1d_c(t)
		i1q_c(t)
		ucd(t) = real(uc)
		ucq(t) = imag(uc)
		ucd_c(t)
		ucq_c(t)
		idrr(t) = idr_m
		iqrr(t) = iqr_m
		idr(t)
		iqr(t)
		imr(t)
		iid(t) = iid_m
		iiq(t) = iiq_m
		id(t) = id_m
		iq(t) = iq_m
		ud(t)
		uq(t)
	end

	T1 = [cos(δ) sin(δ); -sin(δ) cos(δ)] # from grid to component
	T2 = [cos(δ) -sin(δ); sin(δ) cos(δ)] # from component to grid
	eq1 = [
		# 控制方程 
		pg ~ ucd_c * i1d_c + ucq_c * i1q_c,
		qg ~ -ucd_c * i1q_c + ucq_c * i1d_c,
		D(δ) / ωb ~ ω - ωr,
		D(ω) / ωb ~ 0,
		D(idrr) / ωb ~ hp * ((pr * ucd_c + qr * ucq_c) / clamp((ucd_c^2 + ucq_c^2), 0.64, 1.44) - i1d_c)+h_sat*(idr-idrr),
		D(iqrr) / ωb ~ hp * ((pr * ucq_c - qr * ucd_c) / clamp((ucd_c^2 + ucq_c^2), 0.64, 1.44) - i1q_c)+h_sat*(iqr-iqrr),
		imr ~ ifelse(sqrt(idrr^2 + iqrr^2) <= imax, 1.0, imax / sqrt(idrr^2 + iqrr^2)),
		idr ~ idrr * imr,
		iqr ~ iqrr * imr,
		D(iid) / ωb ~ i1d_c - idr,
		D(iiq) / ωb ~ i1q_c - iqr,
		vd_c ~ -kp2 * (i1d_c - idr) - ki2 * iid + ucd_c - ω * lf * i1q_c,
		vq_c ~ -kp2 * (i1q_c - iqr) - ki2 * iiq + ucq_c + ω * lf * i1d_c,
		# 系统方程
		D(i1d) / ωb ~ -rf * i1d / lf + (vd - ucd) / lf + ωr * i1q,
		D(i1q) / ωb ~ -rf * i1q / lf + (vq - ucq) / lf - ωr * i1d,
		D(ucd) / ωb ~ (i1d - id) / cf + ωr * ucq,
		D(ucq) / ωb ~ (i1q - iq) / cf - ωr * ucd,
		D(id) / ωb ~ -rg * id / lg + (ucd - ud) / lg + ωr * iq,
		D(iq) / ωb ~ -rg * iq / lg + (ucq - uq) / lg - ωr * id,
	]
	# 坐标变换
	eq2 = [i1d_c; i1q_c] .~ T1 * [i1d; i1q]# from grid to component
	eq3 = [ucd_c, ucq_c] .~ T1 * [ucd; ucq]
	eq4 = [vd, vq] .~ T2 * [vd_c; vq_c]# from component to grid
	eqs = vcat([eq1, eq2, eq3, eq4]...)
	# 定义离散事件
	discrete_events = [
		0.02 => [ω ~ clamp(Pre(ω) - (qr-qg)/pg*hp, 0.95, 1.05)],
	]
	# 创建 System
	return System(eqs, t, sts, ps; name = name, discrete_events = discrete_events)
end

