function COBS end

const COBS_DEFAULTS = (
	ωb = 2pi*60,
	ωhat = 1.0,
	cf = 0.05,
	lf = 0.12,
	rf = 0.035,
	lg = 0.1,
	rg = 0.02,
	kv = 1.5,
	ko = 0.5,
	kp = 0.03,
	kc = 10.0,
	id_m = 0,
	iq_m = 0.1,
	ud_m = 0,
	uq_m = 1,
	imax = 1.2,
	vmax = 1.3,
	ratio = 1.0,
	ra = 0.10,
	wa = 1.00,
)

function operating_point(::typeof(COBS), voltage, current, params::NamedTuple = (;))
	p = merge(COBS_DEFAULTS, params,
		(id_m=real(current), iq_m=imag(current), ud_m=real(voltage), uq_m=imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	cf = p.cf*scale
	lf, rf, lg, rg = p.lf/scale, p.rf/scale, p.lg/scale, p.rg/scale
	kp, kc, ra, imax = p.kp/scale, p.kc/scale, p.ra/scale, p.imax*scale
	uc = voltage+(rg+1im*lg)*current
	i1 = current+1im*cf*uc
	sc = uc*conj(i1)
	n = uc/abs(uc)
	ia = i1*conj(n)
	states = (vrd=real(uc), vrq=imag(uc), χpar=real(ia), χperp=imag(ia),
		i1d=real(i1), i1q=imag(i1), ucd=real(uc), ucq=imag(uc),
		id=real(current), iq=imag(current))
	parameters = (ωb=p.ωb, ωhat=p.ωhat, cf=cf, lf=lf, rf=rf, lg=lg, rg=rg, kp=kp,
		kv=p.kv, ko=p.ko, kc=kc, ra=ra, wa=p.wa, pr=real(sc), v_m=abs(uc),
		imax=imax, vmax=p.vmax, ratio=p.ratio, scale=scale)
	return (; states, parameters)
end

function COBS(; name, params::NamedTuple = (;))
	p = merge(COBS_DEFAULTS, params)
	(; ωb, ωhat, cf, lf, rf, lg, rg, kv, ko, kp, kc, ra, wa,
		id_m, iq_m, ud_m, uq_m, imax, vmax) = p
	scale = current_capacity_scale(id_m, iq_m, p.ratio)
	cf *= scale
	lf /= scale
	rf /= scale
	lg /= scale
	rg /= scale
	kp /= scale
	kc /= scale
	ra /= scale
	imax *= scale
	u_m = ud_m+1im*uq_m
	i_m = id_m+1im*iq_m
	uc_m = u_m+(rg+1im*lg)*i_m
	i1_m = i_m+1im*cf*uc_m
	sc = uc_m*conj(i1_m)
	ia_m = i1_m*conj(uc_m/abs(uc_m))

	ps = @parameters begin
		ωb=ωb; ωhat=ωhat; cf=cf; lf=lf; rf=rf; lg=lg; rg=rg
		kp=kp; kv=kv; ko=ko; kc=kc; ra=ra; wa=wa
		pr=real(sc); v_m=abs(uc_m); imax=imax; vmax=vmax
		ratio=p.ratio; scale=scale
	end
	sts = @variables begin
		ω(t); ωr(t); pg(t); qg(t); p̂g(t); q̂g(t); vr(t)
		nd(t); nq(t); ipar(t); iperp(t); ihpar(t); ihperp(t)
		χpar(t)=real(ia_m); χperp(t)=imag(ia_m)
		eobsd(t); eobsq(t); edampd(t); edampq(t); etotd(t); etotq(t)
		eecd(t); eecq(t); em(t); ecd(t); ecq(t)
		vrd(t)=real(uc_m); vrq(t)=imag(uc_m)
		vd(t); vq(t); vds(t); vqs(t); vmr(t)
		i1d(t)=real(i1_m); i1q(t)=imag(i1_m)
		ucd(t)=real(uc_m); ucq(t)=imag(uc_m)
		id(t)=id_m; iq(t)=iq_m; ud(t); uq(t)
	end
	eqs = [
		pg~ucd*i1d+ucq*i1q, qg~-ucd*i1q+ucq*i1d,
		p̂g~vrd*i1d+vrq*i1q, q̂g~-vrd*i1q+vrq*i1d,
		vr~sqrt(vrd^2+vrq^2), nd~vrd/vr, nq~vrq/vr,
		ipar~nd*i1d+nq*i1q, iperp~-nq*i1d+nd*i1q,
		D(χpar)/ωb~wa*(ipar-χpar), D(χperp)/ωb~wa*(iperp-χperp),
		ihpar~ipar-χpar, ihperp~iperp-χperp,
		eobsd~(kp*(pr-p̂g)*vrd+(vrd+kv*vrq)*(v_m-vr))/vr,
		eobsq~(kp*(pr-p̂g)*vrq+(vrq-kv*vrd)*(v_m-vr))/vr,
		edampd~-ra*(nd*ihpar-nq*ihperp),
		edampq~-ra*(nq*ihpar+nd*ihperp),
		etotd~eobsd+edampd, etotq~eobsq+edampq,
		eecd~i1d+etotd/kc, eecq~i1q+etotq/kc,
		em~min(sqrt(eecd^2+eecq^2),imax),
		ecd~kc*(em*eecd/sqrt(eecd^2+eecq^2)-i1d),
		ecq~kc*(em*eecq/sqrt(eecd^2+eecq^2)-i1q),
		# The transient-current damping acts only on the inverter-voltage command.
		# Keep the observer synchronization law identical to the original OBS law.
		ω~ωhat-ko*kv*(v_m-vr)/vr,
		D(vrd)/ωb~ko*eobsd+(ωr-ωhat)*vrq,
		D(vrq)/ωb~ko*eobsq-(ωr-ωhat)*vrd,
		vd~vrd+ecd+rf*i1d-lf*i1q,
		vq~vrq+ecq+rf*i1q+lf*i1d,
		vmr~ifelse(sqrt(vd^2+vq^2)<=vmax,1.0,vmax/sqrt(vd^2+vq^2)),
		vds~vd*vmr, vqs~vq*vmr,
		D(i1d)/ωb~-rf*i1d/lf+(vds-ucd)/lf+ωr*i1q,
		D(i1q)/ωb~-rf*i1q/lf+(vqs-ucq)/lf-ωr*i1d,
		D(ucd)/ωb~(i1d-id)/cf+ωr*ucq,
		D(ucq)/ωb~(i1q-iq)/cf-ωr*ucd,
		D(id)/ωb~-rg*id/lg+(ucd-ud)/lg+ωr*iq,
		D(iq)/ωb~-rg*iq/lg+(ucq-uq)/lg-ωr*id,
	]
	return System(eqs,t,sts,ps;name=name)
end
