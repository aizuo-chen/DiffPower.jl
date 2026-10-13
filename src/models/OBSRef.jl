function OBSRef end

const OBS_REF_DEFAULTS = OBS_DEFAULTS

function operating_point(::typeof(OBSRef), voltage, current, params::NamedTuple = (;))
	p = merge(OBS_REF_DEFAULTS, params,
		(id_m=real(current), iq_m=imag(current), ud_m=real(voltage), uq_m=imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	cf = p.cf*scale
	lf, rf, lg, rg = p.lf/scale, p.rf/scale, p.lg/scale, p.rg/scale
	lo, kp, kc, imax = p.lo/scale, p.kp/scale, p.kc/scale, p.imax*scale
	uc = current*(rg+1im*lg)+voltage
	i1 = current+1im*cf*uc
	vr = i1*(rf+1im*lf)+uc
	sr = vr*conj(i1)
	uu = vr+(p.ko-1im)*lo*i1
	states = (uuq=imag(uu), i1d=real(i1), i1q=imag(i1), ucd=real(uc),
		ucq=imag(uc), id=real(current), iq=imag(current))
	parameters = (ωb=p.ωb, cf=cf, lf=lf, rf=rf, lg=lg, rg=rg, lo=lo, kp=kp,
		kv=p.kv, ko=p.ko, kc=kc, pr=real(sr), v_m=abs(vr), imax=imax,
		vmax=p.vmax, ratio=p.ratio, scale=scale)
	return (; states, parameters)
end

function OBSRef(; name, params::NamedTuple = (;))
	p = merge(OBS_REF_DEFAULTS, params)
	(; ωb, cf, lf, rf, lg, rg, lo, kv, ko, kp, kc,
		id_m, iq_m, ud_m, uq_m, imax, vmax) = p
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

	u_m = ud_m+1im*uq_m
	i_m = id_m+1im*iq_m
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

	sts = @variables begin
		ω(t) = 1.0
		pg(t)
		qg(t)
		p̂g(t)
		ecd(t)
		ecq(t)
		eecd(t)
		eecq(t)
		em(t)
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
		pg ~ ucd * i1d + ucq * i1q,
		qg ~ -ucd * i1q + ucq * i1d,
		p̂g ~ vrd * i1d + vrq * i1q,
		vr ~ sqrt(vrd^2 + vrq^2),
		eecd ~ i1d + 1/kc*(kp*(pr-p̂g)*vrd/vr + (vrd + kv*vrq)*(v_m-vr)/vr),
		eecq ~ i1q + 1/kc*(kp*(pr-p̂g)*vrq/vr + (vrq - kv*vrd)*(v_m-vr)/vr),
		em ~ min(sqrt(eecd^2 + eecq^2), imax),
		ecd ~ kc * (em * eecd / sqrt(eecd^2 + eecq^2) - i1d),
		ecq ~ kc * (em * eecq / sqrt(eecd^2 + eecq^2) - i1q),
		ω ~ 1 - ko * ecd / uuq,
		D(uuq) / ωb ~ ko * ecq,
		vrd ~ -ko*lo*i1d-i1q*lo,
		vrq ~ uuq-ko*lo*i1q+i1d*lo,
		vd ~ vrd+ecd,
		vq ~ vrq+ecq,
		vmr ~ ifelse(sqrt(vd^2 + vq^2) <= vmax, 1.0, vmax / sqrt(vd^2 + vq^2)),
		vds ~ vd * vmr,
		vqs ~ vq * vmr,
		D(i1d) / ωb ~ -rf * i1d / lf + (vds - ucd) / lf + ω * i1q,
		D(i1q) / ωb ~ -rf * i1q / lf + (vqs - ucq) / lf - ω * i1d,
		D(ucd) / ωb ~ (i1d - id) / cf + ω * ucq,
		D(ucq) / ωb ~ (i1q - iq) / cf - ω * ucd,
		D(id) / ωb ~ -rg * id / lg + (ucd - ud) / lg + ω * iq,
		D(iq) / ωb ~ -rg * iq / lg + (ucq - uq) / lg - ω * id,
	]

	return System(eqs, t, sts, ps; name = name)
end

function reference_axis_phasor(::typeof(OBSRef), voltage, current, params::NamedTuple)
	ratio = get(params, :ratio, OBS_REF_DEFAULTS.ratio)
	scale = current_capacity_scale(real(current), imag(current), ratio)
	cf = get(params, :cf, OBS_REF_DEFAULTS.cf) * scale
	lf = get(params, :lf, OBS_REF_DEFAULTS.lf) / scale
	rf = get(params, :rf, OBS_REF_DEFAULTS.rf) / scale
	lg = get(params, :lg, OBS_REF_DEFAULTS.lg) / scale
	rg = get(params, :rg, OBS_REF_DEFAULTS.rg) / scale
	lo = get(params, :lo, OBS_REF_DEFAULTS.lo) / scale
	ko = get(params, :ko, OBS_REF_DEFAULTS.ko)
	uc = voltage + (rg + 1im*lg) * current
	i1 = current + 1im*cf*uc
	vr = uc + (rf + 1im*lf) * i1
	uu = vr + (ko - 1im) * lo * i1
	return uu
end

