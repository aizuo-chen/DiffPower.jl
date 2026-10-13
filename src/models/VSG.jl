function VSG end

jacobian_state_names(::typeof(VSG)) =
	(:δ, :ω, :uud, :uuq, :iid, :iiq, :i1d, :i1q, :ucd, :ucq, :id, :iq)

const VSG_DEFAULTS = (
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
	# Negative gains because both PI outputs contain minus ki * integrator.
	# hsc: current saturation -> voltage PI; hsv: voltage saturation -> current PI.
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

function operating_point(::typeof(VSG), voltage, current, params::NamedTuple = (;))
	p = merge(VSG_DEFAULTS, params,
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
	δ = angle(uc)-pi/2
	V = abs(uc)
	idr, iqr = [cos(δ) sin(δ); -sin(δ) cos(δ)]*[real(i1); imag(i1)]
	states = (δ = δ, ω = 1.0, uud = -(idr+cf*V)/ki1, uuq = -iqr/ki1,
		iid = -rf*idr/ki2, iiq = -rf*iqr/ki2, i1d = real(i1), i1q = imag(i1),
		ucd = real(uc), ucq = imag(uc), id = real(current), iq = imag(current))
	parameters = (ωb = p.ωb, lf = lf, rf = rf, cf = cf, lg = lg, rg = rg, mq = mq,
		pr = real(sr), qr = imag(sr), v_m = V, Mv = Mv, Dv = Dv,
		ki1 = ki1, kp1 = kp1, ki2 = ki2, kp2 = kp2,
		hsc = hsc, hsv = p.hsv, imax = imax, vmax = p.vmax, ratio = p.ratio, scale = scale)
	return (; states, parameters)
end

function VSG(; name, params::NamedTuple = (;))
	p = merge(VSG_DEFAULTS, params)
	(; ωb, cf, lf, rf, lg, rg, mq, Mv, Dv, kp1, ki1, kp2, ki2,
		id_m, iq_m, ud_m, uq_m, hsc, hsv, imax, vmax) = p
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
	δ_m = angle(uc)-pi/2
	V_m = abs(uc)
	idr_m, iqr_m = [cos(δ_m) sin(δ_m); -sin(δ_m) cos(δ_m)]*[real(i1); imag(i1)]
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
		δ(t) = δ_m
		ω(t) = 1.0
		ωr(t)
		pg(t)
		qg(t)
		vr(t)
		vd_c(t)
		vq_c(t)
		vds_c(t)
		vqs_c(t)
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

	T1 = [cos(δ) sin(δ); -sin(δ) cos(δ)]
	T2 = [cos(δ) -sin(δ); sin(δ) cos(δ)]
	eq1 = [
		pg ~ ucd_c*i1d_c+ucq_c*i1q_c,
		qg ~ -ucd_c*i1q_c+ucq_c*i1d_c,
		Mv*D(ω)/ωb ~ pr-pg-Dv*(ω-1),
		D(δ)/ωb ~ ω-ωr,
		vr ~ v_m+mq*(qr-qg),
		D(uud)/ωb ~ ucd_c+hsc*(idr-idrr),
		D(uuq)/ωb ~ ucq_c-vr+hsc*(iqr-iqrr),
		idrr ~ -ki1*uud-kp1*ucd_c-ω*cf*ucq_c,
		iqrr ~ -ki1*uuq-kp1*(ucq_c-vr)+ω*cf*ucd_c,
		imr ~ ifelse(sqrt(idrr^2+iqrr^2) <= imax, 1.0, imax/sqrt(idrr^2+iqrr^2)),
		idr ~ idrr*imr,
		iqr ~ iqrr*imr,
		# Feed back saturated minus unsaturated voltage in the controller dq frame.
		D(iid)/ωb ~ i1d_c-idr+hsv*(vds_c-vd_c),
		D(iiq)/ωb ~ i1q_c-iqr+hsv*(vqs_c-vq_c),
		vd_c ~ -kp2*(i1d_c-idr)-ki2*iid+ucd_c-ω*lf*i1q_c,
		vq_c ~ -kp2*(i1q_c-iqr)-ki2*iiq+ucq_c+ω*lf*i1d_c,
		vmr ~ ifelse(sqrt(vd^2+vq^2) <= vmax, 1.0, vmax/sqrt(vd^2+vq^2)),
		vds_c ~ vd_c*vmr,
		vqs_c ~ vq_c*vmr,
		vds ~ vd*vmr,
		vqs ~ vq*vmr,
		D(i1d)/ωb ~ -rf*i1d/lf+(vds-ucd)/lf+ωr*i1q,
		D(i1q)/ωb ~ -rf*i1q/lf+(vqs-ucq)/lf-ωr*i1d,
		D(ucd)/ωb ~ (i1d-id)/cf+ωr*ucq,
		D(ucq)/ωb ~ (i1q-iq)/cf-ωr*ucd,
		D(id)/ωb ~ -rg*id/lg+(ucd-ud)/lg+ωr*iq,
		D(iq)/ωb ~ -rg*iq/lg+(ucq-uq)/lg-ωr*id,
	]
	eq2 = [i1d_c; i1q_c] .~ T1*[i1d; i1q]
	eq3 = [ucd_c; ucq_c] .~ T1*[ucd; ucq]
	eq4 = [vd; vq] .~ T2*[vd_c; vq_c]
	eqs = vcat(eq1, eq2, eq3, eq4)

	return System(eqs, t, sts, ps; name = name)
end
