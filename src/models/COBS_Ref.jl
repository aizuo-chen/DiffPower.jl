function COBS_Ref end

const COBS_REF_DEFAULTS = COBS_DEFAULTS

function operating_point(::typeof(COBS_Ref),voltage,current,params::NamedTuple=(;))
	p=merge(COBS_REF_DEFAULTS,params,
		(id_m=real(current),iq_m=imag(current),ud_m=real(voltage),uq_m=imag(voltage)))
	scale=current_capacity_scale(p.id_m,p.iq_m,p.ratio)
	cf=p.cf*scale
	lf,rf,lg,rg=p.lf/scale,p.rf/scale,p.lg/scale,p.rg/scale
	kp,kc,ra,imax=p.kp/scale,p.kc/scale,p.ra/scale,p.imax*scale
	uc=voltage+(rg+1im*lg)*current
	i1=current+1im*cf*uc
	sc=uc*conj(i1)
	ia=i1*conj(uc/abs(uc))
	states=(vrq=abs(uc),χpar=real(ia),χperp=imag(ia),i1d=real(i1),
		i1q=imag(i1),ucd=real(uc),ucq=imag(uc),id=real(current),iq=imag(current))
	parameters=(ωb=p.ωb,ωhat=p.ωhat,cf=cf,lf=lf,rf=rf,lg=lg,rg=rg,
		kp=kp,kv=p.kv,ko=p.ko,kc=kc,ra=ra,wa=p.wa,pr=real(sc),
		v_m=abs(uc),imax=imax,vmax=p.vmax,ratio=p.ratio,scale=scale)
	return (;states,parameters)
end

function COBS_Ref(;name,params::NamedTuple=(;))
	p=merge(COBS_REF_DEFAULTS,params)
	(;ωb,ωhat,cf,lf,rf,lg,rg,kv,ko,kp,kc,ra,wa,
		id_m,iq_m,ud_m,uq_m,imax,vmax)=p
	scale=current_capacity_scale(id_m,iq_m,p.ratio)
	cf*=scale; lf/=scale; rf/=scale; lg/=scale; rg/=scale
	kp/=scale; kc/=scale; ra/=scale; imax*=scale
	u_m=ud_m+1im*uq_m; i_m=id_m+1im*iq_m
	uc_m=u_m+(rg+1im*lg)*i_m
	i1_m=i_m+1im*cf*uc_m
	sc=uc_m*conj(i1_m)
	ia_m=i1_m*conj(uc_m/abs(uc_m))

	ps=@parameters begin
		ωb=ωb; ωhat=ωhat; cf=cf; lf=lf; rf=rf; lg=lg; rg=rg
		kp=kp; kv=kv; ko=ko; kc=kc; ra=ra; wa=wa
		pr=real(sc); v_m=abs(uc_m); imax=imax; vmax=vmax
		ratio=p.ratio; scale=scale
	end
	sts=@variables begin
		ω(t)=ωhat; pg(t); qg(t); p̂g(t); q̂g(t); vr(t)
		nd(t); nq(t); ipar(t); iperp(t); ihpar(t); ihperp(t)
		χpar(t)=real(ia_m); χperp(t)=imag(ia_m)
		eobsd(t); eobsq(t); edampd(t); edampq(t); etotd(t); etotq(t)
		eecd(t); eecq(t); em(t); ecd(t); ecq(t)
		vrd(t); vrq(t)=abs(uc_m)
		vd(t); vq(t); vds(t); vqs(t); vmr(t)
		i1d(t)=real(i1_m); i1q(t)=imag(i1_m)
		ucd(t)=real(uc_m); ucq(t)=imag(uc_m)
		id(t)=id_m; iq(t)=iq_m; ud(t); uq(t)
	end
	eqs=[
		pg~ucd*i1d+ucq*i1q, qg~-ucd*i1q+ucq*i1d,
		p̂g~vrd*i1d+vrq*i1q, q̂g~-vrd*i1q+vrq*i1d,
		vrd~0, vr~vrq, nd~0, nq~1,
		ipar~i1q, iperp~-i1d,
		D(χpar)/ωb~wa*(ipar-χpar), D(χperp)/ωb~wa*(iperp-χperp),
		ihpar~ipar-χpar, ihperp~iperp-χperp,
		eobsd~kv*(v_m-vr),
		eobsq~kp*(pr-p̂g)+(v_m-vr),
		edampd~ra*ihperp,
		edampq~-ra*ihpar,
		etotd~eobsd+edampd, etotq~eobsq+edampq,
		eecd~i1d+etotd/kc, eecq~i1q+etotq/kc,
		em~min(sqrt(eecd^2+eecq^2),imax),
		ecd~kc*(em*eecd/sqrt(eecd^2+eecq^2)-i1d),
		ecq~kc*(em*eecq/sqrt(eecd^2+eecq^2)-i1q),
		ω~ωhat-ko*kv*(v_m-vr)/vr,
		# As in the full COBS model, auxiliary damping is applied only to the
		# inverter-voltage command and does not modify the OBS synchronization law.
		D(vrq)/ωb~ko*eobsq,
		vd~ecd+rf*i1d-lf*i1q,
		vq~vrq+ecq+rf*i1q+lf*i1d,
		vmr~ifelse(sqrt(vd^2+vq^2)<=vmax,1.0,vmax/sqrt(vd^2+vq^2)),
		vds~vd*vmr, vqs~vq*vmr,
		D(i1d)/ωb~-rf*i1d/lf+(vds-ucd)/lf+ω*i1q,
		D(i1q)/ωb~-rf*i1q/lf+(vqs-ucq)/lf-ω*i1d,
		D(ucd)/ωb~(i1d-id)/cf+ω*ucq,
		D(ucq)/ωb~(i1q-iq)/cf-ω*ucd,
		D(id)/ωb~-rg*id/lg+(ucd-ud)/lg+ω*iq,
		D(iq)/ωb~-rg*iq/lg+(ucq-uq)/lg-ω*id,
	]
	return System(eqs,t,sts,ps;name)
end

function reference_axis_phasor(::typeof(COBS_Ref),voltage,current,params::NamedTuple)
	ratio=get(params,:ratio,COBS_REF_DEFAULTS.ratio)
	scale=current_capacity_scale(real(current),imag(current),ratio)
	rg=get(params,:rg,COBS_REF_DEFAULTS.rg)/scale
	lg=get(params,:lg,COBS_REF_DEFAULTS.lg)/scale
	return voltage+(rg+1im*lg)*current
end
