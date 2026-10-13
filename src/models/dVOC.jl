function dVOC(; name, params::NamedTuple = (;))
	defaults = (
		ωb = 2pi*60,
		cf = 0.05,
		lf = 0.17,
		rf = 0.04,
		lg = 0.2,
		rg = 0.01,
		mp = 0.015, # equal mp in the droop
		mq = 0.75, # 
		hv = 7,
		ζv = 1.0,
		hc = 11,
		id_m = 0.0,
		iq_m = 0.1,
		ud_m = 0,
		uq_m = 1,
		h_sat = -1,
		imax = 5,
		vmax = 1.4,
	);
	p = merge(defaults, params);
	# 解构初始条件向量
	(; ωb, cf, lf, rf, lg, rg, mp, mq, hv, ζv, hc, id_m, iq_m, ud_m, uq_m, h_sat, imax, vmax) = p
	# 控制器参数
	ki1 = hv^2*cf; ## ki1 = ωₙ²*cf(二阶系统，标准化)
	kp1 = 2*ζv*hv*cf; ## kp1 = 2ζωₙcf
	ki2 = hc*rf; ## ki = ω*rf(二阶系统, 零极点对消除)
	kp2 = hc*lf;  ## kp = ω*Lf
	## 迭代初始值
	uc = (id_m+1im*iq_m) * (rg+1im*lg)+ud_m+1im*uq_m
	i1 = id_m+1im*iq_m+1im*cf*uc
	sr = uc*conj(i1)
	δ_m = angle(uc) - π/2
	V_m = abs(uc)
	idr_m, iqr_m = [cos(δ_m) sin(δ_m); -sin(δ_m) cos(δ_m)] * [real(i1); imag(i1)]
	uud_m = -(idr_m+cf*V_m) / ki1
	uuq_m = -iqr_m/ki1
	iid_m = -rf*idr_m/ki2
	iiq_m = -rf*iqr_m/ki2

	ps = @parameters begin
		ωb = ωb
		lf = lf
		rf = rf
		cf = cf
		lg = lg
		rg = rg
		mp = mp
		mq = mq
		pr = real(sr)
		qr = imag(sr)
		v_m = V_m
		ki1 = ki1
		kp1 = kp1
		ki2 = ki2
		kp2 = kp2
		h_sat = h_sat
		imax = imax
		vmax = vmax
	end
	# 定义变量
	sts = @variables begin
		δ(t) = δ_m
		ω(t)
		ωr(t)
		pg(t)
		qg(t)
		vr(t) = V_m
		vd_c(t)
		vq_c(t)
		vd(t)
		vq(t)
		vds(t)
		vqs(t)
		vmr(t)
		i1d(t) = real(i1)
		i1q(t) = imag(i1)
		i1d_c(t)
		i1q_c(t)
		ucd(t) = real(uc)
		ucq(t) = imag(uc)
		ucd_c(t)
		ucq_c(t)
		idr(t)
		iqr(t)
		idrr(t)
		iqrr(t)
		imr(t)
		uud(t) = uud_m
		uuq(t) = uuq_m
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
		ω ~ clamp(1 + mp / vr^2 * (pr - pg), 0.95, 1.05),
		D(δ) / ωb ~ ω - ωr,
		D(vr) / ωb ~ mq * vr * (v_m^2 - vr^2) + mp / vr * (qr - qg),
		D(uud) / ωb ~ ucd_c - 0+h_sat*(idr-idrr),
		D(uuq) / ωb ~ ucq_c - vr+h_sat*(iqr-iqrr),
		idrr ~ -ki1 * uud - kp1 * (ucd_c - 0) - ω * cf * ucq_c,
		iqrr ~ -ki1 * uuq - kp1 * (ucq_c - vr) + ω * cf * ucd_c,
		imr ~ ifelse(sqrt(idrr^2 + iqrr^2) <= imax, 1.0, imax / sqrt(idrr^2 + iqrr^2)),
		idr ~ idrr * imr,
		iqr ~ iqrr * imr,
		D(iid) / ωb ~ i1d_c - idr,
		D(iiq) / ωb ~ i1q_c - iqr,
		vd_c ~ -kp2 * (i1d_c - idr) - ki2 * iid + ucd_c - ω * lf * i1q_c,
		vq_c ~ -kp2 * (i1q_c - iqr) - ki2 * iiq + ucq_c + ω * lf * i1d_c,
		vmr ~ ifelse(sqrt(vd^2 + vq^2) <= vmax, 1.0, vmax / sqrt(vd^2 + vq^2)),
		vds ~ vd * vmr,
		vqs ~ vq * vmr,
		# 系统方程
		D(i1d) / ωb ~ -rf * i1d / lf + (vds - ucd) / lf + ωr * i1q,
		D(i1q) / ωb ~ -rf * i1q / lf + (vqs - ucq) / lf - ωr * i1d,
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

	return System(eqs, t, sts, ps; name = name)
end

