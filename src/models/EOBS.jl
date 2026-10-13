function EOBS end

const EOBS_DEFAULTS = merge(COBS_DEFAULTS, (
	km = 1.0,
	ra = 0.5,
	wa = 10.0,
	kq = 0.10,
))

function operating_point(::typeof(EOBS), voltage, current, params::NamedTuple = (;))
	p = merge(EOBS_DEFAULTS, params,
		(id_m=real(current), iq_m=imag(current), ud_m=real(voltage), uq_m=imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	cf = p.cf*scale
	lf, rf, lg, rg = p.lf/scale, p.rf/scale, p.lg/scale, p.rg/scale
	kp, kq, kc, ra, imax = p.kp/scale, p.kq/scale, p.kc/scale, p.ra/scale, p.imax*scale
	uc = voltage+(rg+1im*lg)*current
	i1 = current+1im*cf*uc
	sc = uc*conj(i1)
	n = uc/abs(uc)
	ia = i1*conj(n)
	states = (vrd=real(uc), vrq=imag(uc), χpar=real(ia), χperp=imag(ia),
		i1d=real(i1), i1q=imag(i1), ucd=real(uc), ucq=imag(uc),
		id=real(current), iq=imag(current))
	parameters = (ωb=p.ωb, ωhat=p.ωhat, cf=cf, lf=lf, rf=rf, lg=lg, rg=rg, kp=kp,
		km=p.km, kv=p.kv, ko=p.ko, kc=kc, ra=ra, wa=p.wa, kq=kq,
		pr=real(sc), qr=imag(sc), v_m=abs(uc),
		imax=imax, vmax=p.vmax, ratio=p.ratio, scale=scale)
	return (; states, parameters)
end

function EOBS(; name, params::NamedTuple = (;))
	p = merge(EOBS_DEFAULTS, params)
	(; ωb, ωhat, cf, lf, rf, lg, rg, km, kv, ko, kp, kq, kc, ra, wa,
		id_m, iq_m, ud_m, uq_m, imax, vmax) = p
	scale = current_capacity_scale(id_m, iq_m, p.ratio)
	cf *= scale
	lf /= scale
	rf /= scale
	lg /= scale
	rg /= scale
	kp /= scale
	kq /= scale
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
		kp=kp; km=km; kv=kv; ko=ko; kc=kc; ra=ra; wa=wa; kq=kq
		pr=real(sc); qr=imag(sc); v_m=abs(uc_m); imax=imax; vmax=vmax
		ratio=p.ratio; scale=scale
	end
	sts = @variables begin
		ω(t); ωr(t); pg(t); qg(t); vr(t)
		nd(t); nq(t); ipar(t); iperp(t); ihpar(t); ihperp(t)
		χpar(t)=real(ia_m); χperp(t)=imag(ia_m)
		eobsd(t); eobsq(t); eqd(t); eqq(t)
		edampd(t); edampq(t); etotd(t); etotq(t)
		eecd(t); eecq(t); em(t); ecd(t); ecq(t)
		vrd(t)=real(uc_m); vrq(t)=imag(uc_m)
		vd(t); vq(t); vds(t); vqs(t); vmr(t)
		i1d(t)=real(i1_m); i1q(t)=imag(i1_m)
		ucd(t)=real(uc_m); ucq(t)=imag(uc_m)
		id(t)=id_m; iq(t)=iq_m; ud(t); uq(t)
	end
	eqs = [
		pg~ucd*i1d+ucq*i1q, qg~-ucd*i1q+ucq*i1d,
		vr~sqrt(vrd^2+vrq^2), nd~vrd/vr, nq~vrq/vr,
		ipar~nd*i1d+nq*i1q, iperp~-nq*i1d+nd*i1q,
		D(χpar)/ωb~wa*(ipar-χpar), D(χperp)/ωb~wa*(iperp-χperp),
		ihpar~ipar-χpar, ihperp~iperp-χperp,
		# The paper's observer uses measured capacitor-port active power.
		eobsd~(kp*(pr-pg)*vrd+(km*vrd+kv*vrq)*(v_m-vr))/vr,
		eobsq~(kp*(pr-pg)*vrq+(km*vrq-kv*vrd)*(v_m-vr))/vr,
		# Explicit Q-V droop is parallel to the observer voltage and acts only
		# through the inverter command, not through the observer dynamics.
		eqd~kq*(qr-qg)*nd, eqq~kq*(qr-qg)*nq,
		edampd~-ra*(nd*ihpar-nq*ihperp),
		edampq~-ra*(nq*ihpar+nd*ihperp),
		etotd~eobsd+eqd+edampd, etotq~eobsq+eqq+edampq,
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
