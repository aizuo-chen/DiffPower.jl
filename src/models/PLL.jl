function PLL end

jacobian_state_names(::typeof(PLL)) =
	(:δ, :pll, :pp, :qq, :iid, :iiq, :i1d, :i1q, :ucd, :ucq, :id, :iq)

const PLL_DEFAULTS = (
	ωb = 2pi*60,
	cf = 0.05,
	lf = 0.12,
	rf = 0.035,
	lg = 0.05,
	rg = 0.02,
	hp = 0.1,
	hpll = 0.5,
	ξpll = 0.707,
	hc = 11,
	h_sat = 1,
	hsv = -1,
	imax = 5,
	vmax = 1.3,
	id_m = 0,
	iq_m = 0.1,
	ud_m = 0,
	uq_m = 1,
	ratio = 1.0,
)

function operating_point(::typeof(PLL), voltage, current, params::NamedTuple = (;))
	p = merge(PLL_DEFAULTS, params,
		(id_m = real(current), iq_m = imag(current), ud_m = real(voltage), uq_m = imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	cf = p.cf*scale
	lf, rf, lg, rg = p.lf/scale, p.rf/scale, p.lg/scale, p.rg/scale
	imax = p.imax*scale
	# ki2 scales as 1/scale, so scale hsv to retain the same voltage tracking rate.
	hsv = p.hsv*scale
	ki_pll = p.hpll^2
	kp_pll = 2*p.ξpll*p.hpll
	ki2 = p.hc*rf
	kp2 = p.hc*lf
	uc = current*(rg+1im*lg)+voltage
	i1 = current+1im*cf*uc
	sr = uc*conj(i1)
	δ = angle(uc)-pi/2
	idr, iqr = [cos(δ) sin(δ); -sin(δ) cos(δ)]*[real(i1); imag(i1)]
	states = (δ = δ, pll = 0.0, pp = iqr/p.hp, qq = idr/p.hp,
		iid = -rf*idr/ki2, iiq = -rf*iqr/ki2,
		i1d = real(i1), i1q = imag(i1), ucd = real(uc), ucq = imag(uc),
		id = real(current), iq = imag(current))
	parameters = (ωb = p.ωb, lf = lf, rf = rf, cf = cf, lg = lg, rg = rg,
		pr = real(sr), qr = imag(sr), hp = p.hp, h_sat = p.h_sat,
		hsv = hsv, imax = imax, vmax = p.vmax,
		ki_pll = ki_pll, kp_pll = kp_pll, ki2 = ki2, kp2 = kp2,
		ratio = p.ratio, scale = scale)
	return (; states, parameters)
end

function PLL(; name, params::NamedTuple = (;))
	p = merge(PLL_DEFAULTS, params)
	point = operating_point(PLL, p.ud_m+1im*p.uq_m, p.id_m+1im*p.iq_m, p)
	(; ωb, lf, rf, cf, lg, rg, pr, qr, hp, h_sat, hsv, imax, vmax,
		ki_pll, kp_pll, ki2, kp2, ratio, scale) = point.parameters
	δ_m, pll_m, pp_m, qq_m, iid_m, iiq_m =
		(point.states[k] for k in (:δ, :pll, :pp, :qq, :iid, :iiq))
	i1d_m, i1q_m, ucd_m, ucq_m, id_m, iq_m =
		(point.states[k] for k in (:i1d, :i1q, :ucd, :ucq, :id, :iq))

	ps = @parameters begin
		ωb = ωb
		lf = lf
		rf = rf
		cf = cf
		lg = lg
		rg = rg
		pr = pr
		qr = qr
		hp = hp
		h_sat = h_sat
		hsv = hsv
		imax = imax
		vmax = vmax
		ki_pll = ki_pll
		kp_pll = kp_pll
		ki2 = ki2
		kp2 = kp2
		ratio = ratio
		scale = scale
	end
	# 定义变量
	sts = @variables begin
		δ(t) = δ_m
		ω(t)
		ωr(t)
		pg(t)
		qg(t)
		vd_c(t)
		vq_c(t)
		vds_c(t)
		vqs_c(t)
		vd(t)
		vq(t)
		vds(t)
		vqs(t)
		vmr(t)
		i1d(t) = i1d_m
		i1q(t) = i1q_m
		i1d_c(t)
		i1q_c(t)
		ucd(t) = ucd_m
		ucq(t) = ucq_m
		ucd_c(t)
		ucq_c(t)
		idr(t)
		iqr(t)
		idrr(t)
		iqrr(t)
		imr(t)
		pll(t) = pll_m
		pp(t) = pp_m
		qq(t) = qq_m
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
		D(pll) / ωb ~ ucd_c,
		ω ~ clamp(-ki_pll * pll - kp_pll * ucd_c + 1, 0.95, 1.05),
		D(δ) / ωb ~ ω - ωr,
		D(pp) / ωb ~ pr - pg+h_sat*(iqr-iqrr),
		D(qq) / ωb ~ qr - qg+h_sat*(idr-idrr),
		idrr ~ hp*qq,
		iqrr ~ hp*pp,
		imr ~ ifelse(sqrt(idrr^2 + iqrr^2) <= imax, 1.0, imax / sqrt(idrr^2 + iqrr^2)),
		idr ~ idrr * imr,
		iqr ~ iqrr * imr,
		D(iid) / ωb ~ i1d_c - idr + hsv*(vds_c-vd_c),
		D(iiq) / ωb ~ i1q_c - iqr + hsv*(vqs_c-vq_c),
		vd_c ~ -kp2 * (i1d_c - idr) - ki2 * iid + ucd_c - ω * lf * i1q_c,
		vq_c ~ -kp2 * (i1q_c - iqr) - ki2 * iiq + ucq_c + ω * lf * i1d_c,
		vmr ~ ifelse(sqrt(vd^2 + vq^2) <= vmax, 1.0, vmax / sqrt(vd^2 + vq^2)),
		vds_c ~ vd_c * vmr,
		vqs_c ~ vq_c * vmr,
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
	# 创建 ODESystem
	return System(eqs, t, sts, ps; name = name)
end
