function OBS end

const OBS_DEFAULTS = (
		ωb = 2pi*60,
		cf = 0.05,
		lf = 0.12,
		rf = 0.035,
		lg = 0.1,
		rg = 0.02,
		kv = 1.5,
		ko = 0.5,
		kp = 0.03,
		lo = 0.05,
		kc = 10.0,
		id_m = 0,
		iq_m = 0.1,
		ud_m = 0,
		uq_m = 1,
		imax = 1.2,
		vmax = 1.3,
		ratio = 1.0,
)

function operating_point(::typeof(OBS), voltage, current, params::NamedTuple = (;))
	p = merge(OBS_DEFAULTS, params,
		(id_m = real(current), iq_m = imag(current), ud_m = real(voltage), uq_m = imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	cf = p.cf*scale
	lf, rf, lg, rg = p.lf/scale, p.rf/scale, p.lg/scale, p.rg/scale
	lo, kp, kc, imax = p.lo/scale, p.kp/scale, p.kc/scale, p.imax*scale
	uc = current*(rg+1im*lg)+voltage
	i1 = current+1im*cf*uc
	vr = i1*(rf+1im*lf)+uc
	sr = vr*conj(i1)
	uu = vr+(p.ko-1im)*lo*i1
	states = (uud=real(uu), uuq=imag(uu), i1d=real(i1), i1q=imag(i1),
		ucd=real(uc), ucq=imag(uc), id=real(current), iq=imag(current))
	parameters = (ωb=p.ωb, cf=cf, lf=lf, rf=rf, lg=lg, rg=rg, lo=lo, kp=kp,
		kv=p.kv, ko=p.ko, kc=kc, pr=real(sr), v_m=abs(vr), imax=imax,
		vmax=p.vmax, ratio=p.ratio, scale=scale)
	return (; states, parameters)
end

function OBS(; name, params::NamedTuple = (;))
	p = merge(OBS_DEFAULTS, params);
	# unpack parameters
	(; ωb, cf, lf, rf, lg, rg, lo, kv, ko, kp, kc, id_m, iq_m, ud_m, uq_m,
		imax, vmax) = p
	scale = current_capacity_scale(id_m, iq_m, p.ratio)
	cf *= scale
	lf /= scale
	rf /= scale
	lg /= scale
	rg /= scale
	lo /= scale
	kp /= scale
	kc /= scale
	imax *= scale
	# initial conditions
	u_m = ud_m+1im*uq_m
	i_m = id_m + 1im*iq_m
	uc_m = i_m * (rg+1im*lg) + u_m
	i1_m = i_m + 1im*cf*uc_m
	vr_m = i1_m * (rf+1im*lf) + uc_m
	sr = vr_m*conj(i1_m)
	uu = vr_m + (ko-1im)*lo*i1_m
	ps = @parameters begin
		ωb = ωb
		cf = cf
		lf = lf
		rf = rf
		lg = lg
		rg = rg
		lo = lo
		kp = kp
		kv = kv
		ko = ko
		kc = kc
		pr = real(sr)
		v_m = abs(vr_m)
		imax = imax
		vmax = vmax
		ratio = p.ratio
		scale = scale
	end
	# 定义变量
	sts = @variables begin
		ω(t)
		ωr(t)
		pg(t)
		qg(t)
		p̂g(t)
		q̂g(t)
		ecd(t)
		ecq(t)
		eecd(t)
		eecq(t)
		em(t)
		uud(t) = real(uu)
		uuq(t) = imag(uu)
		vr(t)
		vrd(t)
		vrq(t)
		vd(t)
		vq(t)
		vds(t)
		vqs(t)
		vmr(t)
		i1d(t) = real(i1_m)
		i1q(t) = imag(i1_m)
		ucd(t) = real(uc_m)
		ucq(t) = imag(uc_m)
		id(t) = id_m
		iq(t) = iq_m
		ud(t)
		uq(t)
	end

	eqs = [
		# Control equations
		pg ~ ucd * i1d + ucq * i1q,
		qg ~ -ucd * i1q + ucq * i1d,
		p̂g ~ vrd * i1d + vrq * i1q,
		q̂g ~ -vrd * i1q + vrq * i1d,
		ω ~ 1 - ko * kv * (v_m - vr) / vr -
			(ωr - 1) * lo * (ko * p̂g - q̂g) / vr^2,
		vr ~ sqrt(vrd^2 + vrq^2),
		eecd ~ i1d + 1/kc*(kp*(pr-p̂g)*vrd/vr + (vrd + kv*vrq)*(v_m-vr)/vr),
		eecq ~ i1q + 1/kc*(kp*(pr-p̂g)*vrq/vr + (vrq - kv*vrd)*(v_m-vr)/vr),
		em ~ min(sqrt(eecd^2 + eecq^2), imax),
		ecd ~ kc * (em * eecd / sqrt(eecd^2 + eecq^2) - i1d),
		ecq ~ kc * (em * eecq / sqrt(eecd^2 + eecq^2) - i1q),
		D(uud) / ωb ~ ko * ecd + (ωr - 1) * uuq,
		D(uuq) / ωb ~ ko * ecq - (ωr - 1) * uud,
		vrd ~ uud-ko*lo*i1d-i1q*lo,
		vrq ~ uuq-ko*lo*i1q+i1d*lo,
		vd ~ vrd+ecd,
		vq ~ vrq+ecq,
		vmr ~ ifelse(sqrt(vd^2 + vq^2) <= vmax, 1.0, vmax / sqrt(vd^2 + vq^2)),
		vds ~ vd * vmr,
		vqs ~ vq * vmr,
		# Network equations
		D(i1d) / ωb ~ -rf * i1d / lf + (vds - ucd) / lf + ωr * i1q,
		D(i1q) / ωb ~ -rf * i1q / lf + (vqs - ucq) / lf - ωr * i1d,
		D(ucd) / ωb ~ (i1d - id) / cf + ωr * ucq,
		D(ucq) / ωb ~ (i1q - iq) / cf - ωr * ucd,
		D(id) / ωb ~ -rg * id / lg + (ucd - ud) / lg + ωr * iq,
		D(iq) / ωb ~ -rg * iq / lg + (ucq - uq) / lg - ωr * id,
	]
	# System construction
	return System(eqs, t, sts, ps; name = name)
end

