function SG6Ref end

const SG6_REF_DEFAULTS = (
		ωb = 2pi*60, lg = 0.2, rg = 0.05, H = 2.656, Da = 2.0,
		Ra = 0.0043, Xl = 0.15, Xd = 1.67, Xq = 1.6,
		Xd′ = 0.265, Xq′ = 0.46, X″ = 0.205,
		Td0′ = 0.5871, Tq0′ = 0.1351, Td0″ = 0.0248, Tq0″ = 0.0267,
		S10 = 0.091, S12 = 0.4,
		id_m = 0.0, iq_m = 0.1, ud_m = 0.0, uq_m = 1.0,
		ratio = 1.0,
)

function operating_point(::typeof(SG6Ref), voltage, current, params::NamedTuple = (;))
	p = merge(SG6_REF_DEFAULTS, params,
		(id_m=real(current), iq_m=imag(current), ud_m=real(voltage), uq_m=imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	lg, rg = p.lg/scale, p.rg/scale
	H, Da, Ra = p.H*scale, p.Da*scale, p.Ra/scale
	Xl, Xd, Xq = p.Xl/scale, p.Xd/scale, p.Xq/scale
	Xd′, Xq′, X″ = p.Xd′/scale, p.Xq′/scale, p.X″/scale
	if p.S10 == 0.0 && p.S12 == 0.0
		A, B = 1.0, 0.0
	else
		ρ = sqrt(1.2*p.S12/p.S10)
		A = (ρ-1.2)/(ρ-1.0)
		B = p.S10/(1.0-A)^2
	end
	Id, Iq = real(current), imag(current)
	Vd = Id*rg-lg*Iq+real(voltage)
	Vq = Iq*rg+lg*Id+imag(voltage)
	ψ″ = hypot(Vd-X″*Iq+Ra*Id, Vq+X″*Id+Ra*Iq)
	Sat = ψ″ > A ? B*(ψ″-A)^2/ψ″ : 0.0
	ψd, ψq = Vq+Ra*Iq, -Vd-Ra*Id
	ψd″, ψq″ = ψd+X″*Id, ψq+X″*Iq
	Eq′ = ψd″+(Xd′-X″)*Id
	ψd′ = Eq′-(Xd′-Xl)*Id
	Ed′ = -ψq″-(Xq′-X″)*Iq
	ψq′ = Ed′+(Xq′-Xl)*Iq
	LadIfd = Eq′+(Xd-Xd′)*(Id+(Xd′-X″)/(Xd′-Xl)^2*(Eq′-ψd′-(Xd′-Xl)*Id))+Sat/ψ″*ψd″
	Te = ψd*Iq-ψq*Id
	states = (ω=1.0, Eq′=Eq′, ψd′=ψd′, Ed′=Ed′, ψq′=ψq′, id=Id, iq=Iq)
	parameters = (ωb=p.ωb, Xl=Xl, Xd=Xd, Xq=Xq, Xd′=Xd′, Xq′=Xq′, X″=X″,
		Ra=Ra, Td0′=p.Td0′, Tq0′=p.Tq0′, Td0″=p.Td0″, Tq0″=p.Tq0″,
		H=H, Da=Da, Efd=LadIfd, Pmech=Te, lg=lg, rg=rg, A=A, B=B,
		ratio=p.ratio, scale=scale)
	return (; states, parameters)
end

function SG6Ref(; name, params::NamedTuple = (;))
	p = merge(SG6_REF_DEFAULTS, params)
	(; ωb, lg, rg, H, Da, Ra, Xl, Xd, Xq, Xd′, Xq′, X″,
		Td0′, Tq0′, Td0″, Tq0″, S10, S12, id_m, iq_m, ud_m, uq_m) = p
	scale = current_capacity_scale(id_m, iq_m, p.ratio)
	lg /= scale;
	rg /= scale;
	H *= scale;
	Da *= scale;
	Ra /= scale
	Xl /= scale;
	Xd /= scale;
	Xq /= scale;
	Xd′ /= scale;
	Xq′ /= scale;
	X″ /= scale

	if S10 == 0.0 && S12 == 0.0
		A, B = 1.0, 0.0
	else
		ρ = sqrt(1.2*S12/S10)
		A = (ρ-1.2)/(ρ-1.0)
		B = S10/(1.0-A)^2
	end
	γqd = (Xq-Xl)/(Xd-Xl)
	vd_m = id_m*rg - lg*iq_m + ud_m
	vq_m = iq_m*rg + lg*id_m + uq_m
	Id_m, Iq_m = id_m, iq_m
	Vd_m, Vq_m = vd_m, vq_m
	ψ″_m = hypot(Vd_m - X″*Iq_m + Ra*Id_m, Vq_m + X″*Id_m + Ra*Iq_m)
	Sat_m = ψ″_m > A ? B*(ψ″_m-A)^2/ψ″_m : 0.0
	ψd_m = Vq_m + Ra*Iq_m
	ψq_m = -Vd_m - Ra*Id_m
	ψd″_m = ψd_m + X″*Id_m
	ψq″_m = ψq_m + X″*Iq_m
	Eq′_m = ψd″_m + (Xd′-X″)*Id_m
	ψd′_m = Eq′_m - (Xd′-Xl)*Id_m
	Ed′_m = -ψq″_m - (Xq′-X″)*Iq_m
	ψq′_m = Ed′_m + (Xq′-Xl)*Iq_m
	LadIfd_m = Eq′_m + (Xd-Xd′)*(Id_m + (Xd′-X″)/(Xd′-Xl)^2*(Eq′_m-ψd′_m-(Xd′-Xl)*Id_m)) + Sat_m/ψ″_m*ψd″_m
	Te_m = ψd_m*Iq_m - ψq_m*Id_m

	ps = @parameters begin
		ωb=ωb
		Xl=Xl
		Xd=Xd
		Xq=Xq
		Xd′=Xd′
		Xq′=Xq′
		X″=X″
		Ra=Ra
		Td0′=Td0′
		Tq0′=Tq0′
		Td0″=Td0″
		Tq0″=Tq0″
		H=H
		Da=Da
		Efd=LadIfd_m
		Pmech=Te_m
		lg=lg
		rg=rg
		A=A
		B=B
		ratio=p.ratio
		scale=scale
	end
	sts = @variables begin
		ω(t)=1.0
		Eq′(t)=Eq′_m
		ψd′(t)=ψd′_m
		Ed′(t)=Ed′_m
		ψq′(t)=ψq′_m
		id(t)=id_m
		iq(t)=iq_m
		ud(t)
		uq(t)
		ψd″(t)
		ψq″(t)
		ψd(t)
		ψq(t)
		ψ″(t)
		Sat(t)
		LadIfd(t)
		Te(t)
		pg(t)
		qg(t)
		Id(t)
		Iq(t)
		vd(t)
		vq(t)
		Vd(t)
		Vq(t)
	end
	eqs = [
		2H*D(ω) ~ (Pmech-Da*(ω-1))/ω - Te,
		Td0′*D(Eq′) ~ Efd-LadIfd,
		Td0″*D(ψd′) ~ Eq′-ψd′-(Xd′-Xl)*Id,
		ψd″ ~ (X″-Xl)/(Xd′-Xl)*Eq′ + (Xd′-X″)/(Xd′-Xl)*ψd′,
		Tq0′*D(Ed′) ~ -Ed′+(Xq-Xq′)*(Iq-(Xq′-X″)/(Xq′-Xl)^2*(Ed′-ψq′+(Xq′-Xl)*Iq))+(Xq-Xl)/(Xd-Xl)*Sat/ψ″*ψq″,
		Tq0″*D(ψq′) ~ Ed′-ψq′+(Xq′-Xl)*Iq,
		ψq″ ~ -((X″-Xl)/(Xq′-Xl)*Ed′+(Xq′-X″)/(Xq′-Xl)*ψq′),
		ψ″ ~ sqrt(ψd″^2+ψq″^2),
		Sat ~ ifelse(ψ″>A, B*(ψ″-A)^2/ψ″, 0.0),
		LadIfd ~ Eq′+(Xd-Xd′)*(Id+(Xd′-X″)/(Xd′-Xl)^2*(Eq′-ψd′-(Xd′-Xl)*Id))+Sat/ψ″*ψd″,
		ψd ~ ψd″-X″*Id,
		ψq ~ ψq″-X″*Iq,
		Te ~ ψd*Iq-ψq*Id,
		Vd ~ -ω*ψq-Ra*Id,
		Vq ~ ω*ψd-Ra*Iq,
		pg ~ Vd*Id+Vq*Iq,
		qg ~ -Vd*Iq+Vq*Id,
		Id ~ id,
		Iq ~ iq,
		vd ~ Vd,
		vq ~ Vq,
		D(id)/ωb ~ -rg*id/lg+(vd-ud)/lg+ω*iq,
		D(iq)/ωb ~ -rg*iq/lg+(vq-uq)/lg-ω*id,
	]
	return System(eqs, t, sts, ps; name = name)
end


function reference_axis_phasor(::typeof(SG6Ref), voltage, current, params::NamedTuple)
	ratio = get(params, :ratio, 1.0)
	scale = current_capacity_scale(real(current), imag(current), ratio)
	rg = get(params, :rg, 0.05) / scale
	lg = get(params, :lg, 0.2) / scale
	Ra = get(params, :Ra, 0.0043) / scale
	Xl = get(params, :Xl, 0.15) / scale
	Xd = get(params, :Xd, 1.67) / scale
	Xq = get(params, :Xq, 1.6) / scale
	X″ = get(params, :X″, 0.205) / scale
	S10 = get(params, :S10, 0.091)
	S12 = get(params, :S12, 0.4)

	if S10 == 0.0 && S12 == 0.0
		A, B = 1.0, 0.0
	else
		ρ = sqrt(1.2*S12/S10)
		A = (ρ-1.2)/(ρ-1.0)
		B = S10/(1.0-A)^2
	end

	internal_voltage = voltage + (rg + 1im*lg) * current
	vd, vq = real(internal_voltage), imag(internal_voltage)
	id, iq = real(current), imag(current)
	ψ″ = hypot(vd - X″*iq + Ra*id, vq + X″*id + Ra*iq)
	Sat = ψ″ > A ? B*(ψ″-A)^2/ψ″ : 0.0
	γqd = (Xq-Xl)/(Xd-Xl)
	Xq_eff = X″ + (Xq-X″)/(1.0 + γqd*Sat/ψ″)
	δ = atan(Xq_eff*iq - vd - Ra*id, Xq_eff*id + vq + Ra*iq)

	return cis(δ + π/2)
end
