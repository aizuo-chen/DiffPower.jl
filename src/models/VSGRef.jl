function VSGRef end

jacobian_state_names(::typeof(VSGRef)) =
	(:ω, :uud, :uuq, :iid, :iiq, :i1d, :i1q, :ucd, :ucq, :id, :iq)

const VSG_REF_DEFAULTS = (
	ωb = 2pi*60,
	cf = 0.05,
	lf = 0.17,
	rf = 0.04,
	lg = 0.2,
	rg = 0.01,
	Mv = 500.0,
	Dv = 70.7,
	mq = 0.0,
	kp1 = 2.0,
	ki1 = 5.0,
	kp2 = 1.7,
	ki2 = 0.4,
	hsc = -1,
	hsv = -1,
	imax = 1.2,
	vmax = 1.3,
	id_m = 0,
	iq_m = 0.1,
	ud_m = 0,
	uq_m = 1,
	ratio = 1.0,
)

function operating_point(::typeof(VSGRef), voltage, current, params::NamedTuple = (;))
	p = merge(VSG_REF_DEFAULTS, params,
		(id_m = real(current), iq_m = imag(current), ud_m = real(voltage), uq_m = imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	cf = p.cf*scale
	lf, rf, lg, rg = p.lf/scale, p.rf/scale, p.lg/scale, p.rg/scale
	mq, Mv, Dv, hsc = p.mq/scale, p.Mv*scale, p.Dv*scale, p.hsc/scale
	imax = p.imax*scale
	(; kp1, ki1, kp2, ki2) = p
	uc = current*(rg+1im*lg)+voltage
	i1 = current+1im*cf*uc
	sr = uc*conj(i1)
	V = abs(uc)
	idr, iqr = real(i1), imag(i1)
	states = (ω = 1.0, uud = -(idr+cf*V)/ki1, uuq = -iqr/ki1,
		iid = -rf*idr/ki2, iiq = -rf*iqr/ki2, i1d = real(i1), i1q = imag(i1),
		ucd = real(uc), ucq = imag(uc), id = real(current), iq = imag(current))
	parameters = (ωb = p.ωb, lf = lf, rf = rf, cf = cf, lg = lg, rg = rg, mq = mq,
		pr = real(sr), qr = imag(sr), v_m = V, Mv = Mv, Dv = Dv,
		ki1 = ki1, kp1 = kp1, ki2 = ki2, kp2 = kp2,
		hsc = hsc, hsv = p.hsv, imax = imax, vmax = p.vmax, ratio = p.ratio, scale = scale)
	return (; states, parameters)
end

function VSGRef(; name, params::NamedTuple = (;))
	p = merge(VSG_REF_DEFAULTS, params)
	(; ωb, cf, lf, rf, lg, rg, mq, Mv, Dv, kp1, ki1, kp2, ki2, hsc, hsv,
		id_m, iq_m, ud_m, uq_m, imax, vmax) = p
	scale = current_capacity_scale(id_m, iq_m, p.ratio)
	cf *= scale
	lf /= scale
	rf /= scale
	lg /= scale
	rg /= scale
	mq /= scale
	Mv *= scale
	Dv *= scale
	hsc /= scale
	imax *= scale

	uc = (id_m+1im*iq_m)*(rg+1im*lg)+ud_m+1im*uq_m
	i1 = id_m+1im*iq_m+1im*cf*uc
	sr = uc*conj(i1)
	V_m = abs(uc)
	idr_m, iqr_m = real(i1), imag(i1)
	uud_m = -(idr_m+cf*V_m)/ki1
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
		mq = mq
		pr = real(sr)
		qr = imag(sr)
		v_m = V_m
		Mv = Mv
		Dv = Dv
		ki1 = ki1
		kp1 = kp1
		ki2 = ki2
		kp2 = kp2
		hsc = hsc
		hsv = hsv
		imax = imax
		vmax = vmax
		ratio = p.ratio
		scale = scale
	end

	sts = @variables begin
		ω(t) = 1.0
		pg(t)
		qg(t)
		vr(t)
		vd(t)
		vq(t)
		vds(t)
		vqs(t)
		vmr(t)
		i1d(t) = real(i1)
		i1q(t) = imag(i1)
		ucd(t) = real(uc)
		ucq(t) = imag(uc)
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

	eqs = [
		pg ~ ucd*i1d+ucq*i1q,
		qg ~ -ucd*i1q+ucq*i1d,
		Mv*D(ω) ~ pr-pg-Dv*(ω-1),
		vr ~ v_m+mq*(qr-qg),
		D(uud)/ωb ~ ucd+hsc*(idr-idrr),
		D(uuq)/ωb ~ ucq-vr+hsc*(iqr-iqrr),
		idrr ~ -ki1*uud-kp1*ucd-ω*cf*ucq,
		iqrr ~ -ki1*uuq-kp1*(ucq-vr)+ω*cf*ucd,
		imr ~ ifelse(sqrt(idrr^2+iqrr^2) <= imax, 1.0, imax/sqrt(idrr^2+iqrr^2)),
		idr ~ idrr*imr,
		iqr ~ iqrr*imr,
		# Voltage saturation feedback uses the same dq frame as the current PI.
		D(iid)/ωb ~ i1d-idr+hsv*(vds-vd),
		D(iiq)/ωb ~ i1q-iqr+hsv*(vqs-vq),
		vd ~ -kp2*(i1d-idr)-ki2*iid+ucd-ω*lf*i1q,
		vq ~ -kp2*(i1q-iqr)-ki2*iiq+ucq+ω*lf*i1d,
		vmr ~ ifelse(sqrt(vd^2+vq^2) <= vmax, 1.0, vmax/sqrt(vd^2+vq^2)),
		vds ~ vd*vmr,
		vqs ~ vq*vmr,
		D(i1d)/ωb ~ -rf*i1d/lf+(vds-ucd)/lf+ω*i1q,
		D(i1q)/ωb ~ -rf*i1q/lf+(vqs-ucq)/lf-ω*i1d,
		D(ucd)/ωb ~ (i1d-id)/cf+ω*ucq,
		D(ucq)/ωb ~ (i1q-iq)/cf-ω*ucd,
		D(id)/ωb ~ -rg*id/lg+(ucd-ud)/lg+ω*iq,
		D(iq)/ωb ~ -rg*iq/lg+(ucq-uq)/lg-ω*id,
	]

	return System(eqs, t, sts, ps; name = name)
end

function reference_axis_phasor(::typeof(VSGRef), voltage, current, params::NamedTuple)
	ratio = get(params, :ratio, VSG_REF_DEFAULTS.ratio)
	scale = current_capacity_scale(real(current), imag(current), ratio)
	rg = get(params, :rg, VSG_REF_DEFAULTS.rg)/scale
	lg = get(params, :lg, VSG_REF_DEFAULTS.lg)/scale
	return voltage+(rg+1im*lg)*current
end
